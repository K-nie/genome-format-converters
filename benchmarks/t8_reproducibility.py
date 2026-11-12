#!/usr/bin/env python3
"""Reproducibility and correctness checks for T8 (VCF -> pseudohaploid EIGENSTRAT).

No published tool converts VCF to pseudohaploid EIGENSTRAT directly, so T8
is reported as a correctness and reproducibility task rather than a speed
comparison. This script runs ``gfc vcf-to-pseudohaploid`` several ways and
writes one row per check:

``same_seed``        two runs with --seed 42 give byte-identical output.
``threads``          --threads N over a batch equals the serial batch.
``file_order``       a file converted alone (--pattern) equals the same file
                     converted inside the batch: the per-file seed does not
                     depend on which other files are present.
``new_seed_differs`` --seed 43 changes the .geno draws (if any hets exist).
``reference_free``   every output passes ``four_state.check_pseudohaploid``
                     against its input VCF (calls consistent with the source
                     genotypes; het draws balanced, two-sided binomial test).

The batch is built from symlinks to the input VCF under several names, so the
check works on a single VCF; each name gets its own seed (seed XOR crc32).

Usage::

    python benchmarks/t8_reproducibility.py --vcf chr22.vcf.gz \\
        --work /tmp/t8_repro --copies 4 --threads 4 --out t8_repro.tsv

Exit status is non-zero if any check fails.
"""
from __future__ import annotations

import argparse
import filecmp
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from four_state import check_pseudohaploid  # noqa: E402

EXTS = (".geno", ".snp", ".ind")


def _gfc(*args: str) -> None:
    cmd = [sys.executable, "-m", "genome_format_converters.cli",
           "vcf-to-pseudohaploid", *args]
    res = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                         text=True)
    if res.returncode != 0:
        sys.stderr.write(res.stderr)
        raise SystemExit(f"gfc failed: {' '.join(cmd)}")


def _same_tree(a: Path, b: Path, stems) -> tuple:
    for stem in stems:
        for ext in EXTS:
            fa, fb = a / (stem + ext), b / (stem + ext)
            if not fa.is_file() or not fb.is_file():
                return False, f"missing {stem}{ext}"
            if not filecmp.cmp(fa, fb, shallow=False):
                return False, f"{stem}{ext} differs"
    return True, f"{len(stems)} file sets identical"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--vcf", required=True, help="Input VCF (.vcf or .vcf.gz)")
    ap.add_argument("--work", required=True, help="Scratch directory (recreated)")
    ap.add_argument("--copies", type=int, default=4,
                    help="Number of differently named inputs in the batch")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out", default="-", help="TSV path, or - for stdout")
    args = ap.parse_args()

    vcf = Path(args.vcf).resolve()
    suffix = ".vcf.gz" if vcf.name.endswith(".vcf.gz") else ".vcf"
    work = Path(args.work)
    if work.exists():
        shutil.rmtree(work)
    inp = work / "input"
    inp.mkdir(parents=True)
    stems = [f"sample_set_{i}" for i in range(1, args.copies + 1)]
    for stem in stems:
        (inp / (stem + suffix)).symlink_to(vcf)
    pattern = "*" + suffix
    out = {name: work / name for name in
           ("serial_a", "serial_b", "threads", "alone", "seed43")}

    _gfc("--input-dir", str(inp), "--output-dir", str(out["serial_a"]),
         "--pattern", pattern, "--seed", "42", "--force")
    _gfc("--input-dir", str(inp), "--output-dir", str(out["serial_b"]),
         "--pattern", pattern, "--seed", "42", "--force")
    _gfc("--input-dir", str(inp), "--output-dir", str(out["threads"]),
         "--pattern", pattern, "--seed", "42", "--threads", str(args.threads),
         "--force")
    _gfc("--input-dir", str(inp), "--output-dir", str(out["alone"]),
         "--pattern", stems[-1] + suffix, "--seed", "42", "--force")
    _gfc("--input-dir", str(inp), "--output-dir", str(out["seed43"]),
         "--pattern", pattern, "--seed", "43", "--force")

    rows = []
    ok, why = _same_tree(out["serial_a"], out["serial_b"], stems)
    rows.append(("same_seed", ok, why))
    ok, why = _same_tree(out["serial_a"], out["threads"], stems)
    rows.append(("threads", ok, f"--threads {args.threads}: {why}"))
    ok, why = _same_tree(out["serial_a"], out["alone"], stems[-1:])
    rows.append(("file_order", ok, why))

    geno_a = (out["serial_a"] / (stems[0] + ".geno")).read_bytes()
    geno_b = (out["seed43"] / (stems[0] + ".geno")).read_bytes()
    # check_pseudohaploid reads one EIGENSTRAT set per directory, so each
    # output set is copied into its own directory first.
    v42 = []
    for d in (out["serial_a"], out["seed43"]):
        for stem in stems:
            single = work / "check" / d.name / stem
            single.mkdir(parents=True, exist_ok=True)
            for ext in EXTS:
                shutil.copy(d / (stem + ext), single / (stem + ext))
            v42.append((d.name, stem, check_pseudohaploid(single, vcf)))
    # A file with no heterozygous calls has nothing to draw, so identical
    # output under a new seed is then expected. The verdict reason reads
    # "stochastic draw; <ref>/<total> het draws REF".
    def _n_hets(reason: str) -> int:
        try:
            return int(reason.split(";")[1].split()[0].split("/")[1])
        except (IndexError, ValueError):
            return 0
    has_hets = any(_n_hets(r.reason) > 0 for _, _, r in v42)
    if geno_a != geno_b:
        rows.append(("new_seed_differs", True, "--seed 43 changes the draws"))
    elif not has_hets:
        rows.append(("new_seed_differs", True, "no heterozygous calls to draw"))
    else:
        rows.append(("new_seed_differs", False, "--seed 43 gave identical draws"))
    bad = [(d, s, r) for d, s, r in v42 if r.state != "structural"]
    if bad:
        d, s, r = bad[0]
        rows.append(("reference_free", False, f"{d}/{s}: {r.state}: {r.reason}"))
    else:
        rows.append(("reference_free", True,
                     f"{len(v42)} outputs consistent; e.g. {v42[0][2].reason}"))

    fh = sys.stdout if args.out == "-" else open(args.out, "w")
    fh.write("check\tpass\tdetail\n")
    for name, ok, why in rows:
        fh.write(f"{name}\t{'yes' if ok else 'no'}\t{why}\n")
    if fh is not sys.stdout:
        fh.close()
    sys.exit(0 if all(ok for _, ok, _ in rows) else 1)


if __name__ == "__main__":
    main()
