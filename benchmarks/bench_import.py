#!/usr/bin/env python3
"""Measure start-up cost of the gfc CLI: wall time and peak RSS per subcommand.

Each measurement runs in a fresh interpreter, so it includes interpreter
start-up and every import the subcommand triggers, but no conversion work:

* ``lazy``  -- ``python -m genome_format_converters.cli <cmd> --help``. Since
  the dispatch table imports converter modules on first use, ``--help`` for a
  subcommand loads only the CLI; running the subcommand then loads only its
  own converter module. ``import:<cmd>`` rows import that one module, which is
  what a real run of the subcommand adds on top of the CLI.
* ``all``   -- import every converter module in one interpreter. This is the
  start-up cost of the pre-0.2 CLI, which imported all converters eagerly.

Usage::

    python benchmarks/bench_import.py --trials 5 --out import_times.tsv

Columns: mode, target, trial, wall_s, peak_rss_mb, modules (heavy libraries
present in ``sys.modules`` afterwards).
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time

HEAVY = ("Bio", "BCBio", "pandas", "numpy", "cyvcf2", "pysam")

CONVERTERS = [
    "gff3_to_gtf", "gff3_to_bed", "gff3_to_bed12", "genbank_to_gff3",
    "gff3_to_table", "gff3_to_protein", "gtf_to_gff3", "fasta_to_fastq",
    "fastq_to_fasta", "fasta_qual_to_fastq", "fastq_to_fasta_qual",
    "convert_alignment", "fasta_to_table", "stockholm_to_fasta",
    "fasta_to_stockholm", "bam_to_bed", "blast_tab_to_links", "delta_to_tab",
    "maf_to_xmfa", "hmmer_tblout_to_tsv", "orthogroups_to_fasta", "vcf_to_bed",
    "vcf_to_table", "vcf_to_consensus", "vcf_to_eigenstrat",
    "vcf_to_pseudohaploid", "vcf_to_plink", "vcf_to_snp_alignment",
    "eigenstrat_to_vcf", "tree_convert", "annotate_tree",
    "convert_all_gff_fasta_to_gbk",
]

# Benchmark tasks T1-T8 and the converter module each one runs.
TASK_MODULES = {
    "T1": "vcf_to_eigenstrat", "T2": "vcf_to_plink",
    "T3": "convert_all_gff_fasta_to_gbk", "T4": "gff3_to_gtf",
    "T5": "gff3_to_bed12", "T6": "hmmer_tblout_to_tsv",
    "T7": "orthogroups_to_fasta", "T8": "vcf_to_pseudohaploid",
}

_REPORT = ("import sys; print(','.join(m for m in %r if m in sys.modules))"
           % (HEAVY,))


def _run(code: str):
    """Run `code` in a fresh interpreter; return (wall_s, rss_mb, stdout)."""
    t0 = time.perf_counter()
    proc = subprocess.Popen([sys.executable, "-c", code],
                            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    out = proc.stdout.read().decode().strip()
    _, status, ru = os.wait4(proc.pid, 0)
    wall = time.perf_counter() - t0
    proc.returncode = os.waitstatus_to_exitcode(status)
    # ru_maxrss is KiB on Linux, bytes on macOS.
    rss = ru.ru_maxrss / (1024 * 1024 if sys.platform == "darwin" else 1024)
    if proc.returncode != 0:
        raise RuntimeError(f"failed: {code}")
    return wall, rss, out.splitlines()[-1] if out else ""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--trials", type=int, default=5)
    ap.add_argument("--out", default="-", help="TSV path, or - for stdout")
    ap.add_argument("--tasks-only", action="store_true",
                    help="Only the converters used by benchmark tasks T1-T8")
    args = ap.parse_args()

    targets = list(TASK_MODULES.values()) if args.tasks_only else CONVERTERS
    rows = []
    for trial in range(1, args.trials + 1):
        cli = ("import genome_format_converters.cli; " + _REPORT)
        rows.append(("lazy", "cli", trial, *_run(cli)))
        for mod in targets:
            code = ("import genome_format_converters.cli, importlib; "
                    f"importlib.import_module('genome_format_converters.converters.{mod}'); "
                    + _REPORT)
            try:
                rows.append(("lazy", f"import:{mod}", trial, *_run(code)))
            except RuntimeError:  # e.g. pysam/cyvcf2 missing on Windows
                rows.append(("lazy", f"import:{mod}", trial, float("nan"), float("nan"), "unavailable"))
        everything = ("import importlib\n"
                      f"for m in {targets!r}:\n"
                      "    try: importlib.import_module('genome_format_converters.converters.' + m)\n"
                      "    except ImportError: pass\n" + _REPORT)
        rows.append(("all", "all-converters", trial, *_run(everything)))

    fh = sys.stdout if args.out == "-" else open(args.out, "w")
    fh.write("mode\ttarget\ttrial\twall_s\tpeak_rss_mb\tmodules\n")
    for mode, target, trial, wall, rss, mods in rows:
        fh.write(f"{mode}\t{target}\t{trial}\t{wall:.4f}\t{rss:.2f}\t{mods}\n")
    if fh is not sys.stdout:
        fh.close()


if __name__ == "__main__":
    main()
