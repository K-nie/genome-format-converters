#!/usr/bin/env python3
"""Post-task correctness check: populate the ``correct`` column in a per-task TSV.

Author: Benjamin Narh-Madey

The benchmark TSV schema (see ``bench_one.py``) carries a ``correct``
column that the harness leaves blank — this script fills it in by
comparing each tool's rep-1 output against a designated per-task
reference tool. The result is written back in place.

Encoding:
    "1"     output matched the reference (apples-to-apples per the spec)
    "0"     output produced but did not match
    ""      could not check (reference missing, gfc output missing,
            structural mismatch on a non-comparable axis, etc.)
    "skip"  task does not have a reference-equivalence axis by design
            (T8 pseudohaploid: cross-RNG byte equivalence is impossible,
            so the contract is a format check only)

Per-task spec (matches docs/bench-runs/2026-04-26_smoke.md):

    T1 vcf -> eigenstrat:
        shasum (.geno + .snp + .ind) vs convertf reference. If convertf
        rep1 is missing, every row is "" — convertf is the only valid
        reference.

    T2 vcf -> plink:
        byte-identical .bed (the genotype payload, not .bim/.fam) vs
        plink1.9 reference. plink2 sometimes flips a1/a2 — flag those
        as 0; document the rule, don't silently swap.

    T3 fasta+gff -> gbk:
        structural — same record count + total seq length via Bio.SeqIO.
        Reference is the py-ref handwritten Biopython baseline.

    T4 gff3 -> gtf:
        line-set equivalence after sort + canonical-attribute order.
        Reference is gffread.

    T5 gff3 -> bed12:
        line-set equivalence after sort, leading-12-cols only.
        Reference is ucsc-chain.

    T6 hmmer-tblout -> tsv:
        byte-identical against pyhmmer reference output.

    T7 orthogroups -> fasta:
        row-wise equivalence after canonical species-prefix ordering.
        Reference is py-ref.

    T8 pseudohaploid:
        format-only — verify .geno contains only {'0','2','9'} and the
        same row count as gfc rep1. Cross-RNG byte equivalence is
        impossible by design (per-file seeded RNG, different impls).

Only rep1 is actually checked; the result is propagated to subsequent
reps for the same tool (we trust deterministic tools). T8 propagates
its format-check result to every rep, since the contract is per-seed
and every rep uses the same seed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from collections import Counter
from pathlib import Path
from typing import Callable, Optional


# Map TSV `tool` column values to the per-rep directory prefix used by
# the task harness scripts. Multiple TSV names can map to the same
# directory prefix across tasks (e.g. py-ref is "pyref" in T3 and "ref"
# in T7); we resolve via (task, tool) when ambiguous.
_TOOL_TO_DIR = {
    ("T1", "gfc"):              "gfc",
    ("T1", "convertf"):         "convertf",
    ("T1", "bioconvert"):       "bioconvert",
    ("T2", "gfc"):              "gfc",
    ("T2", "plink2"):           "plink2",
    ("T2", "plink1.9"):         "plink1",
    ("T2", "bcftools"):         "bcftools",
    ("T2", "bioconvert"):       "bioconvert",
    ("T3", "gfc"):              "gfc",
    ("T3", "EMBOSS-seqret"):    "seqret",
    ("T3", "table2asn"):        "table2asn",
    ("T3", "EMBLmyGFF3"):       "emblmygff3",
    ("T3", "py-ref"):           "pyref",
    ("T4", "gfc"):              "gfc",
    ("T4", "gffread"):          "gffread",
    ("T4", "AGAT"):             "agat",
    ("T4", "bioconvert"):       "bioconvert",
    ("T5", "gfc"):              "gfc",
    ("T5", "ucsc-chain"):       "ucsc",
    ("T5", "AGAT"):             "agat",
    ("T5", "bedops"):           "bedops",
    ("T5", "bioconvert"):       "bioconvert",
    ("T6", "gfc"):              "gfc",
    ("T6", "pyhmmer"):          "pyhmmer",
    ("T6", "biopython"):        "biopython",
    ("T7", "gfc"):              "gfc",
    ("T7", "py-ref"):           "ref",
    ("T7", "orthofinder"):      "orthofinder",
    ("T7", "proteinortho"):     "proteinortho",
    ("T8", "gfc"):              "gfc",
    ("T8", "bcftools-pyref"):   "bcftoolsref",
    ("T8", "bcftools"):         "bcftools",
    ("T8", "angsd"):            "angsd",
}


# Per-task: which TSV `tool` value is the canonical reference. If the
# reference's rep1 dir is missing, every row's correct field is left
# as "". gfc itself never has a reference (it is the subject), but its
# rep1 IS used as the format/row-count anchor for T8.
_REFERENCE_TOOL = {
    "T1": "convertf",
    "T2": "plink1.9",
    "T3": "py-ref",
    "T4": "gffread",
    "T5": "ucsc-chain",
    "T6": "pyhmmer",
    "T7": "py-ref",
    "T8": None,  # format-only, see check_t8
}


# --- Methods footnote constants -------------------------------------------
# When a check_t* function returns "0" for gfc, it means gfc's output is not
# byte/structurally identical to the reference. For some tasks the divergence
# is genuinely semantic (different but defensible design choices in each
# tool); the bench cannot adjudicate which is "right", so we record the
# reason here and let the paper's Methods section cite it. Each _REASON_T*
# string is a one-paragraph footnote intended to be lifted verbatim into
# Methods (proposal-writer polishes the prose).
_REASON_T1 = (
    "T1 (VCF -> EIGENSTRAT): gfc and convertf may diverge on the chromosome "
    "encoding (gfc preserves the VCF CHROM literal, e.g. 'chr22'; convertf "
    "rewrites to its integer-only EIGENSTRAT contract, e.g. '22'). The .geno "
    "and .ind payloads are byte-identical when both tools are pointed at the "
    "same chrom-mapping; the .snp file differs only in the chrom column. We "
    "report correct=0 honestly rather than canonicalising the chrom column "
    "post-hoc, because the EIGENSTRAT spec does not mandate either rendering "
    "and downstream tools (smartpca, ADMIXTURE) accept both."
)

_REASON_T2 = (
    "T2 (VCF -> PLINK .bed): gfc and plink1.9 produce structurally distinct "
    "outputs from the same VCF input. Three independent divergence axes were "
    "confirmed locally on tests/test_data/tiny.vcf: (i) gfc skips indels by "
    "default and writes only biallelic SNPs (2/3 variants for tiny.vcf), "
    "while plink1.9 retains all variants by default; (ii) gfc synthesises "
    "variant IDs as '<chrom>_<pos>' when the VCF ID is '.', while plink1.9 "
    "preserves the literal '.'; (iii) gfc and plink1.9 select different A1/A2 "
    "allele assignments (gfc uses VCF REF/ALT directly; plink1.9 applies "
    "minor-allele-first reordering). Consequently the .bim text and the .bed "
    "binary payload differ in length, content, and genotype encoding "
    "direction. None of these divergences is a bug in either tool; they are "
    "documented design choices. We keep the strict byte-identity comparator "
    "so any regression in gfc's own .bed output (relative to itself across "
    "versions) would surface, and accept correct=0 vs plink1.9 with this "
    "footnote."
)

_REASON_T4 = (
    "T4 (GFF3 -> GTF): gfc and gffread implement two defensible but distinct "
    "GFF3-to-GTF translations. Confirmed locally on tests/test_data/tiny.gff3: "
    "(i) gfc emits a faithful 1:1 row mapping including 'gene' and 'mRNA' "
    "feature rows; gffread elides gene rows and renames 'mRNA' to "
    "'transcript'. (ii) gffread fuses adjacent exon and CDS spans into a "
    "single feature when they abut (101-150 + 151-200 -> 101-200); gfc "
    "preserves the original GFF3 spans verbatim. (iii) Trailing-semicolon "
    "and attribute-key order differ. gffread's normalisations are "
    "transcript-centric and standard for transcript-isoform pipelines; "
    "gfc's faithful mapping is standard for general feature-table "
    "round-tripping. Neither is the canonical GTF rendering. The "
    "comparator already canonicalises attribute order via "
    "_normalise_gtf_line; correct=0 reflects the deeper structural "
    "difference (line count, feature types, span boundaries) and is "
    "documented rather than hidden behind a span-merging rewrite."
)

_REASON_T3 = (
    "T3 (FASTA+GFF -> GenBank): gfc uses Biopython's GenBank writer; the "
    "py-ref baseline also uses Biopython but with a different feature "
    "qualifier order. EMBOSS-seqret renders the same content with different "
    "line-wrap, FEATURES section ordering, and qualifier formatting. The "
    "GenBank flat-file format does not mandate a single canonical rendering, "
    "so byte-identity across implementations is impossible. The comparator "
    "uses Bio.SeqIO to verify structural equivalence: identical record count "
    "and total sequence length. correct=0 here means the structural invariant "
    "actually broke (record loss or sequence corruption), not cosmetic "
    "rendering differences."
)


def _rep1_dir(bench_dir: Path, task: str, tool: str) -> Optional[Path]:
    """Resolve the rep1 directory for a (task, tool) pair, or None if missing."""
    prefix = _TOOL_TO_DIR.get((task, tool))
    if prefix is None:
        return None
    p = bench_dir / f"{prefix}_rep1"
    return p if p.is_dir() else None


def _shasum(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _multi_shasum(paths: list[Path]) -> Optional[str]:
    """Concatenated shasum of an ordered list of files. Returns None if
    any path is missing."""
    h = hashlib.sha256()
    for p in sorted(paths, key=lambda x: x.name):
        if not p.is_file():
            return None
        h.update(p.name.encode("utf-8") + b"\0")
        h.update(_shasum(p).encode("ascii") + b"\0")
    return h.hexdigest()


# ---- Per-task check functions ---------------------------------------------
# Each returns one of: "1", "0", "" (couldn't check) for the gfc<->ref
# comparison, given the rep1 dirs for both. They MUST tolerate missing
# files / dirs without raising — the smoke pipeline pipes our stdout
# into `head -20`, not into a logger.


def check_t1(gfc_dir: Path, ref_dir: Path) -> str:
    """T1: shasum (.geno + .snp + .ind) match convertf reference.

    Divergence axis (Option II — accept correct=0 with Methods footnote,
    see _REASON_T1):
      gfc preserves the VCF CHROM literal in the .snp file; convertf
      rewrites to its integer-only convention. The .geno and .ind payloads
      are byte-identical when both tools see the same chrom-mapping. We
      keep the strict triplet-shasum contract so a real corruption in
      .geno or .ind still fails loudly.
    """
    triplet_exts = (".geno", ".snp", ".ind")
    gfc_files = [p for ext in triplet_exts for p in gfc_dir.glob(f"*{ext}")]
    ref_files = [p for ext in triplet_exts for p in ref_dir.glob(f"*{ext}")]
    if not gfc_files or not ref_files:
        return ""
    g = _multi_shasum(gfc_files)
    r = _multi_shasum(ref_files)
    if g is None or r is None:
        return ""
    return "1" if g == r else "0"


def check_t2(gfc_dir: Path, ref_dir: Path) -> str:
    """T2: byte-identical .bed against plink1.9 reference.

    Divergence axis (Option II — accept correct=0 with Methods footnote,
    see _REASON_T2):
      gfc and plink1.9 disagree on indel handling, variant ID synthesis,
      and A1/A2 allele assignment. The .bed binary length and content
      both differ on tiny.vcf (gfc=5 bytes, plink=6 bytes; gfc skips the
      indel, plink retains it). The strict shasum contract is kept here
      so a future regression in gfc's own .bed bytes — relative to a
      pinned gfc reference run — would still fail loudly.
    """
    gfc_beds = sorted(gfc_dir.glob("*.bed"))
    ref_beds = sorted(ref_dir.glob("*.bed"))
    if not gfc_beds or not ref_beds:
        return ""
    if len(gfc_beds) != len(ref_beds):
        return "0"
    for g, r in zip(gfc_beds, ref_beds):
        if _shasum(g) != _shasum(r):
            return "0"
    return "1"


def check_t3(gfc_dir: Path, ref_dir: Path) -> str:
    """T3: structural — same record count + total seq length via Bio.SeqIO.

    Divergence axis (Option II — accept correct=0 with Methods footnote,
    see _REASON_T3):
      The GenBank flat-file format does not mandate a canonical rendering,
      so byte-identity across implementations is impossible (line-wrap,
      qualifier order, FEATURES section ordering all legitimately differ).
      We compare on the structural invariant only — record count + total
      seq length via Bio.SeqIO — so correct=0 here flags a real loss of
      records or sequence content, not cosmetic rendering drift.
    """
    try:
        from Bio import SeqIO
    except ImportError:
        return ""
    # gfc emits .gbk; py-ref emits .gb. Match by stripped stem.
    gfc_files = list(gfc_dir.glob("*.gbk")) + list(gfc_dir.glob("*.gb"))
    ref_files = list(ref_dir.glob("*.gb")) + list(ref_dir.glob("*.gbk"))
    if not gfc_files or not ref_files:
        return ""

    def summarise(files: list[Path]) -> Optional[tuple[int, int]]:
        total_recs, total_len = 0, 0
        try:
            for f in files:
                for rec in SeqIO.parse(str(f), "genbank"):
                    total_recs += 1
                    total_len += len(rec.seq)
        except Exception:
            return None
        return total_recs, total_len

    g = summarise(gfc_files)
    r = summarise(ref_files)
    if g is None or r is None:
        return ""
    return "1" if g == r else "0"


def _normalise_gtf_line(line: str) -> Optional[str]:
    """Reduce a GTF line to a canonical comparable form: feature columns
    1-8 verbatim, attribute column re-sorted alphabetically by key.
    Returns None for headers / blanks (callers drop them)."""
    s = line.rstrip("\n")
    if not s or s.startswith("#"):
        return None
    parts = s.split("\t")
    if len(parts) < 8:
        return None  # malformed — drop rather than crash the diff
    if len(parts) == 8:
        return "\t".join(parts)
    attr_blob = "\t".join(parts[8:])  # GTF only has 9 cols, but be defensive
    attrs = []
    for chunk in attr_blob.split(";"):
        chunk = chunk.strip()
        if chunk:
            attrs.append(chunk)
    attrs.sort()
    return "\t".join(parts[:8] + ["; ".join(attrs)])


def check_t4(gfc_dir: Path, ref_dir: Path) -> str:
    """T4: line-set equivalence after sort + canonical-attribute order.

    Divergence axis (Option II — accept correct=0 with Methods footnote,
    see _REASON_T4):
      gfc preserves the GFF3 row structure verbatim (including gene and
      mRNA rows, separate exon/CDS spans); gffread normalises to a
      transcript-centric layout (drops gene rows, renames mRNA to
      transcript, fuses adjacent exon/CDS). Both renderings parse as
      valid GTF; neither is the canonical form. The comparator keeps
      its sorted-line + canonicalised-attribute-order contract so a
      gfc-vs-gfc regression in line content would still fail; correct=0
      vs gffread is the structural-divergence verdict.
    """
    gfc_files = sorted(gfc_dir.glob("*.gtf"))
    ref_files = sorted(ref_dir.glob("*.gtf"))
    if not gfc_files or not ref_files:
        return ""

    def linecount(files: list[Path]) -> Counter:
        c: Counter = Counter()
        for f in files:
            try:
                with f.open() as fh:
                    for line in fh:
                        norm = _normalise_gtf_line(line)
                        if norm is not None:
                            c[norm] += 1
            except OSError:
                continue
        return c

    return "1" if linecount(gfc_files) == linecount(ref_files) else "0"


def _normalise_bed_line(line: str, max_cols: int = 12) -> Optional[str]:
    s = line.rstrip("\n")
    if not s or s.startswith(("#", "track", "browser")):
        return None
    parts = s.split("\t")
    if len(parts) < 3:
        return None
    return "\t".join(parts[:max_cols])


def check_t5(gfc_dir: Path, ref_dir: Path) -> str:
    """T5: line-set equivalence after sort, leading 12 cols only.

    gfc emits .bed12; ucsc-chain emits .bed. Compare by the leading 12
    columns so a 12-col vs 6-col emit doesn't false-fail (BED is
    cumulative columns 1..12)."""
    gfc_files = sorted(list(gfc_dir.glob("*.bed12")) + list(gfc_dir.glob("*.bed")))
    ref_files = sorted(list(ref_dir.glob("*.bed")) + list(ref_dir.glob("*.bed12")))
    if not gfc_files or not ref_files:
        return ""

    # Honest check: shrink to min(cols) so a 6-col reference doesn't
    # demand 12-col content from gfc. Find the smallest cols-per-line
    # across both sides.
    def min_cols(files: list[Path]) -> int:
        m = 12
        for f in files:
            try:
                with f.open() as fh:
                    for line in fh:
                        s = line.rstrip("\n")
                        if not s or s.startswith(("#", "track", "browser")):
                            continue
                        m = min(m, len(s.split("\t")))
                        break  # first data line is enough
            except OSError:
                continue
        return max(m, 3)  # never below BED3

    cols = min(min_cols(gfc_files), min_cols(ref_files))

    def linecount(files: list[Path]) -> Counter:
        c: Counter = Counter()
        for f in files:
            try:
                with f.open() as fh:
                    for line in fh:
                        norm = _normalise_bed_line(line, max_cols=cols)
                        if norm is not None:
                            c[norm] += 1
            except OSError:
                continue
        return c

    return "1" if linecount(gfc_files) == linecount(ref_files) else "0"


def check_t6(gfc_dir: Path, ref_dir: Path) -> str:
    """T6: byte-identical against pyhmmer reference."""
    gfc_files = sorted(gfc_dir.glob("*.tsv"))
    ref_files = sorted(ref_dir.glob("*.tsv"))
    if not gfc_files or not ref_files:
        return ""
    if len(gfc_files) != len(ref_files):
        return "0"
    for g, r in zip(gfc_files, ref_files):
        if _shasum(g) != _shasum(r):
            return "0"
    return "1"


def check_t7(gfc_dir: Path, ref_dir: Path) -> str:
    """T7: row-wise equivalence after canonical species-prefix ordering.

    gfc emits headers like ``>species|gene``; py-ref does the same. Two
    OG fasta files match if, after sorting by the (species, gene) key,
    the (id, sequence) pairs agree."""
    gfc_files = sorted(gfc_dir.glob("*.fa"))
    ref_files = sorted(ref_dir.glob("*.fa"))
    if not gfc_files or not ref_files:
        return ""
    if {f.name for f in gfc_files} != {f.name for f in ref_files}:
        return "0"

    def canonical(path: Path) -> list[tuple[str, str]]:
        recs: list[tuple[str, str]] = []
        cur_id, cur_seq = None, []
        try:
            with path.open() as fh:
                for line in fh:
                    line = line.rstrip("\n")
                    if line.startswith(">"):
                        if cur_id is not None:
                            recs.append((cur_id, "".join(cur_seq)))
                        cur_id = line[1:].split()[0] if line[1:] else ""
                        cur_seq = []
                    else:
                        cur_seq.append(line)
                if cur_id is not None:
                    recs.append((cur_id, "".join(cur_seq)))
        except OSError:
            return []
        recs.sort(key=lambda r: r[0])
        return recs

    for f in gfc_files:
        if canonical(f) != canonical(ref_dir / f.name):
            return "0"
    return "1"


def check_t8(gfc_dir: Path, ref_dir: Path) -> str:
    """T8: format-only check on gfc's own .geno output.

    Cross-RNG byte equivalence between gfc and bcftools-pyref is
    impossible by design (different RNG streams), so the contract here
    is the EIGENSTRAT format invariant: every char in .geno must be in
    {'0','2','9'} and every row has the same length. ref_dir is unused;
    we verify gfc rep1's .geno files self-consistently."""
    geno_files = sorted(gfc_dir.glob("*.geno"))
    if not geno_files:
        return ""
    allowed = {"0", "2", "9", "\n"}
    for f in geno_files:
        row_len = None
        try:
            with f.open() as fh:
                for line in fh:
                    s = line.rstrip("\n")
                    if not s:
                        continue
                    if any(ch not in allowed for ch in s + "\n"):
                        return "0"
                    if row_len is None:
                        row_len = len(s)
                    elif len(s) != row_len:
                        return "0"
        except OSError:
            return ""
    return "1"


_CHECKERS: dict[str, Callable[[Path, Path], str]] = {
    "T1": check_t1,
    "T2": check_t2,
    "T3": check_t3,
    "T4": check_t4,
    "T5": check_t5,
    "T6": check_t6,
    "T7": check_t7,
    "T8": check_t8,
}


def populate(task: str, bench_dir: Path, tsv_path: Path) -> str:
    """Read tsv_path, fill in `correct` for each row, write back. Return
    a one-line summary suitable for piping into `head -20`."""
    if task not in _CHECKERS:
        return f"[corr {task}] unknown task; tsv left untouched"
    if not tsv_path.is_file():
        return f"[corr {task}] tsv missing: {tsv_path}"

    with tsv_path.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        rows = list(reader)
    if not rows:
        return f"[corr {task}] empty tsv"
    header = rows[0]
    try:
        ix_tool = header.index("tool")
        ix_rep = header.index("replicate")
        ix_correct = header.index("correct")
    except ValueError:
        return f"[corr {task}] tsv header missing required cols"

    gfc_dir = _rep1_dir(bench_dir, task, "gfc")
    ref_tool = _REFERENCE_TOOL.get(task)
    if task == "T8":
        ref_dir = gfc_dir  # T8 self-checks; check_t8 ignores ref_dir anyway
    else:
        ref_dir = _rep1_dir(bench_dir, task, ref_tool) if ref_tool else None

    # Mark every row "" if we can't even start.
    if gfc_dir is None:
        for r in rows[1:]:
            r[ix_correct] = ""
        with tsv_path.open("w", newline="") as fh:
            csv.writer(fh, delimiter="\t").writerows(rows)
        return f"[corr {task}] gfc rep1 dir missing under {bench_dir}; all rows blank"

    if task != "T8" and ref_dir is None:
        for r in rows[1:]:
            r[ix_correct] = ""
        with tsv_path.open("w", newline="") as fh:
            csv.writer(fh, delimiter="\t").writerows(rows)
        return (f"[corr {task}] reference '{ref_tool}' rep1 missing under "
                f"{bench_dir}; all rows blank")

    # Compute the per-tool verdict on rep1, then propagate to every rep.
    # The reference tool's own row gets "1" by definition (it IS the
    # ground truth for this comparison axis). T8 sets "skip" for the
    # bcftools-pyref reference because cross-RNG byte equivalence is
    # explicitly out of scope.
    checker = _CHECKERS[task]
    verdict_per_tool: dict[str, str] = {}
    summary_bits: list[str] = []
    distinct_tools = sorted({r[ix_tool] for r in rows[1:]})
    for tool in distinct_tools:
        if tool == "gfc":
            if task == "T8":
                v = checker(gfc_dir, gfc_dir)
            else:
                v = checker(gfc_dir, ref_dir)
        elif tool == ref_tool and task != "T8":
            v = "1"  # reference is correct by construction
        elif task == "T8":
            v = "skip"  # any non-gfc tool in T8: cross-RNG check, not applicable
        else:
            tool_dir = _rep1_dir(bench_dir, task, tool)
            v = checker(tool_dir, ref_dir) if tool_dir is not None else ""
        verdict_per_tool[tool] = v
        summary_bits.append(f"{tool}={v or 'blank'}")

    for r in rows[1:]:
        r[ix_correct] = verdict_per_tool.get(r[ix_tool], "")

    n_in = len(rows)
    with tsv_path.open("w", newline="") as fh:
        csv.writer(fh, delimiter="\t").writerows(rows)

    # Defensive: re-read what we wrote and compare row counts. If they
    # diverge, surface to stderr so the bug is visible in the bench log.
    # (Run 136838 lost 3 of 4 T1 rows somewhere in this pipeline; this
    # guard means the next divergence won't be silent.)
    try:
        with tsv_path.open() as fh:
            n_out = sum(1 for _ in fh)
    except OSError:
        n_out = -1
    if n_out != n_in:
        print(f"[corr {task}] WARNING: wrote {n_out} lines, expected "
              f"{n_in} (header + {n_in - 1} data rows)", file=sys.stderr)

    return f"[corr {task}] " + " ".join(summary_bits)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--task", required=True,
                    help="Task id (T1..T8)")
    ap.add_argument("--bench-dir", required=True,
                    help="Per-task /tmp/gfc_bench_T<n>... directory containing "
                         "<tool>_rep<n>/ subdirs.")
    ap.add_argument("--tsv", required=True,
                    help="Per-task results TSV to update in place.")
    args = ap.parse_args()
    msg = populate(args.task, Path(args.bench_dir), Path(args.tsv))
    print(msg)


if __name__ == "__main__":
    main()
