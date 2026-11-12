#!/usr/bin/env python3
"""Four-state correctness comparators for the gfc benchmark.

Each ``compare_tN(candidate_dir, reference_dir, **ctx)`` returns a
``Verdict(state, reason)`` where ``state`` is one of:

``byte``
    After line-ending normalisation, every output file is byte-identical
    to the reference.
``structural``
    Not byte-identical, but the task-specific comparator finds the same
    biological content once a *named* convention difference is allowed
    for. ``reason`` lists exactly which conventions were needed (for
    example ``chrom-label``, ``allele-order``, ``exon-fusion``).
``limitation``
    The two outputs cannot be compared on this axis (for example a BED6
    reference against BED12 output, or a missing reference file).
    ``reason`` says why.
``disagreement``
    The outputs differ in content that no documented convention explains.
    ``reason`` names the first difference found.

The rules are deliberately executable rather than written as footnotes:
a "structural" verdict is only returned when the code has checked the
content and recorded which tolerance it used. ``mutation_test.py`` checks
that the comparators return ``disagreement`` on seeded defects.
"""
from __future__ import annotations

import csv
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

BYTE = "byte"
STRUCTURAL = "structural"
LIMITATION = "limitation"
DISAGREEMENT = "disagreement"
STATES = (BYTE, STRUCTURAL, LIMITATION, DISAGREEMENT)


@dataclass
class Verdict:
    state: str
    reason: str = ""
    conventions: List[str] = field(default_factory=list)

    def __post_init__(self):
        assert self.state in STATES, self.state

    @property
    def correct_flag(self) -> str:
        """Back-compatible value for the TSV ``correct`` column."""
        return {BYTE: "1", STRUCTURAL: "1", DISAGREEMENT: "0"}.get(self.state, "")


# ---------------------------------------------------------------- helpers


def _files(d: Optional[Path], *patterns: str) -> List[Path]:
    if d is None or not Path(d).is_dir():
        return []
    out: List[Path] = []
    for p in patterns:
        out.extend(sorted(Path(d).glob(p)))
    return out


def _norm_bytes(p: Path) -> bytes:
    return p.read_bytes().replace(b"\r\n", b"\n")


def _byte_equal(cand: Sequence[Path], ref: Sequence[Path]) -> bool:
    if len(cand) != len(ref) or not cand:
        return False
    return all(_norm_bytes(a) == _norm_bytes(b) for a, b in zip(cand, ref))


def _data_lines(p: Path, skip=("#",)) -> List[str]:
    out = []
    with open(p) as fh:
        for line in fh:
            s = line.rstrip("\r\n")
            if s and not s.startswith(skip):
                out.append(s)
    return out


_CHR_PREFIX = re.compile(r"^(chr|Chr|CHR)")


def _norm_chrom(c: str) -> str:
    c = _CHR_PREFIX.sub("", c)
    return {"X": "23", "Y": "24", "MT": "26", "M": "26"}.get(c, c)


def _missing(what: str) -> Verdict:
    return Verdict(LIMITATION, f"missing {what}")


# ------------------------------------------------- T1 VCF -> EIGENSTRAT


def _read_eigenstrat(d: Path):
    geno = _files(d, "*.geno")
    snp = _files(d, "*.snp")
    ind = _files(d, "*.ind")
    if not (geno and snp and ind):
        return None
    snps = [l.split() for f in snp for l in _data_lines(f)]
    genos = [l.strip() for f in geno for l in _data_lines(f)]
    inds = [l.split()[0] for f in ind for l in _data_lines(f)]
    return snps, genos, inds, geno + snp + ind


def _flip_geno(row: str) -> str:
    return row.translate(str.maketrans("02", "20"))


def compare_t1(cand: Path, ref: Path, **_) -> Verdict:
    c = _read_eigenstrat(cand)
    r = _read_eigenstrat(ref)
    if c is None or r is None:
        return _missing("EIGENSTRAT triplet")
    c_snps, c_geno, c_ind, c_files = c
    r_snps, r_geno, r_ind, r_files = r
    if _byte_equal(c_files, r_files):
        return Verdict(BYTE)
    conv = set()
    if c_ind != r_ind:
        if sorted(c_ind) != sorted(r_ind):
            return Verdict(DISAGREEMENT, "sample sets differ")
        return Verdict(DISAGREEMENT, "sample order differs")
    # Index reference sites by (normalised chrom, position).
    r_by_pos: Dict[Tuple[str, int], Tuple[List[str], str]] = {}
    for s, g in zip(r_snps, r_geno):
        r_by_pos[(_norm_chrom(s[1]), int(s[3]))] = (s, g)
    if len(c_snps) != len(c_geno) or len(r_snps) != len(r_geno):
        return Verdict(DISAGREEMENT, ".snp and .geno row counts differ")
    if len(c_snps) != len(r_snps):
        return Verdict(DISAGREEMENT,
                       f"site count {len(c_snps)} vs reference {len(r_snps)}")
    for s, g in zip(c_snps, c_geno):
        key = (_norm_chrom(s[1]), int(s[3]))
        if key not in r_by_pos:
            return Verdict(DISAGREEMENT, f"site {s[1]}:{s[3]} absent in reference")
        rs, rg = r_by_pos[key]
        if s[1] != rs[1]:
            conv.add("chrom-label")
        if s[0] != rs[0]:
            conv.add("snp-id")
        if not math.isclose(float(s[2]), float(rs[2]), abs_tol=1e-6):
            return Verdict(DISAGREEMENT, f"genetic position differs at {s[1]}:{s[3]}")
        ca, ra = (s[4].upper(), s[5].upper()), (rs[4].upper(), rs[5].upper())
        if ca == ra:
            if g != rg:
                return Verdict(DISAGREEMENT, f"genotypes differ at {s[1]}:{s[3]}")
        elif ca == (ra[1], ra[0]):
            conv.add("allele-order")
            if _flip_geno(g) != rg:
                return Verdict(DISAGREEMENT, f"genotypes differ at {s[1]}:{s[3]} (after allele swap)")
        else:
            return Verdict(DISAGREEMENT, f"alleles differ at {s[1]}:{s[3]}")
    if [(_norm_chrom(s[1]), int(s[3])) for s in c_snps] != \
            [(_norm_chrom(s[1]), int(s[3])) for s in r_snps]:
        conv.add("site-order")
    return Verdict(STRUCTURAL, ", ".join(sorted(conv)) or "whitespace/format",
                   sorted(conv))


# ------------------------------------------------- T2 VCF -> PLINK .bed


def _read_plink(d: Path):
    bed, bim, fam = _files(d, "*.bed"), _files(d, "*.bim"), _files(d, "*.fam")
    if not (bed and bim and fam):
        return None
    bims = [l.split() for l in _data_lines(bim[0])]
    fams = [l.split() for l in _data_lines(fam[0])]
    raw = bed[0].read_bytes()
    if raw[:3] != b"\x6c\x1b\x01":
        return ("bad-magic",)
    n = len(fams)
    per = (n + 3) // 4
    body = raw[3:]
    if len(body) != per * len(bims):
        return ("bad-length",)
    # Decode to per-sample A1 dosage: 0b00 hom A1 -> 2, 0b10 het -> 1,
    # 0b11 hom A2 -> 0, 0b01 missing -> None.
    code = {0: 2, 2: 1, 3: 0, 1: None}
    genos = []
    for i in range(len(bims)):
        chunk = body[i * per:(i + 1) * per]
        row = []
        for j in range(n):
            row.append(code[(chunk[j // 4] >> (2 * (j % 4))) & 3])
        genos.append(row)
    return bims, fams, genos, bed + bim + fam


def compare_t2(cand: Path, ref: Path, **_) -> Verdict:
    c, r = _read_plink(cand), _read_plink(ref)
    if c is None or r is None:
        return _missing("PLINK .bed/.bim/.fam")
    if len(c) == 1 or len(r) == 1:
        return Verdict(DISAGREEMENT, f"malformed .bed ({c[0] if len(c) == 1 else r[0]})")
    c_bim, c_fam, c_geno, c_files = c
    r_bim, r_fam, r_geno, r_files = r
    if _byte_equal(c_files, r_files):
        return Verdict(BYTE)
    conv = set()
    if [f[1] for f in c_fam] != [f[1] for f in r_fam]:
        return Verdict(DISAGREEMENT, "sample IDs or order differ in .fam")
    if [f[:2] for f in c_fam] != [f[:2] for f in r_fam]:
        conv.add("family-id")

    r_index: Dict[Tuple[str, int], int] = {}
    for i, b in enumerate(r_bim):
        r_index.setdefault((_norm_chrom(b[0]), int(b[3])), i)
    matched = set()
    for i, b in enumerate(c_bim):
        key = (_norm_chrom(b[0]), int(b[3]))
        if key not in r_index:
            return Verdict(DISAGREEMENT, f"variant {b[0]}:{b[3]} absent in reference")
        j = r_index[key]
        matched.add(j)
        rb = r_bim[j]
        if b[0] != rb[0]:
            conv.add("chrom-label")
        if b[1] != rb[1]:
            conv.add("variant-id")
        ca, ra = (b[4], b[5]), (rb[4], rb[5])
        cg, rg = c_geno[i], r_geno[j]
        if ca == ra:
            same = cg == rg
        elif ca == (ra[1], ra[0]):
            conv.add("allele-order")
            same = [None if x is None else 2 - x for x in cg] == rg
        elif "0" in ra and (ca[0] in ra or ca[1] in ra):
            # Reference wrote a monomorphic site with A1 = '0'.
            conv.add("monomorphic-allele-code")
            if ca[1] == ra[1]:
                same = cg == rg
            else:
                same = [None if x is None else 2 - x for x in cg] == rg
        else:
            return Verdict(DISAGREEMENT, f"alleles differ at {b[0]}:{b[3]}")
        if not same:
            return Verdict(DISAGREEMENT, f"genotypes differ at {b[0]}:{b[3]}")
    extra = [r_bim[j] for j in range(len(r_bim)) if j not in matched]
    if extra:
        non_snp = [b for b in extra if len(b[4]) != 1 or len(b[5]) != 1]
        if len(non_snp) != len(extra):
            return Verdict(DISAGREEMENT,
                           f"{len(extra) - len(non_snp)} reference SNPs missing from candidate")
        conv.add("indels-excluded")
    order = [r_index[(_norm_chrom(b[0]), int(b[3]))] for b in c_bim]
    if order != sorted(order):
        conv.add("site-order")
    return Verdict(STRUCTURAL, ", ".join(sorted(conv)) or "format", sorted(conv))


# ------------------------------------------------- T3 FASTA+GFF -> GenBank


def compare_t3(cand: Path, ref: Path, **_) -> Verdict:
    try:
        from Bio import SeqIO
    except ImportError:
        return Verdict(LIMITATION, "Biopython not available")
    c_files = _files(cand, "*.gbk", "*.gb")
    r_files = _files(ref, "*.gbk", "*.gb")
    if not c_files or not r_files:
        return _missing("GenBank files")
    if _byte_equal(c_files, r_files):
        return Verdict(BYTE)

    def load(files):
        recs = {}
        for f in files:
            for rec in SeqIO.parse(str(f), "genbank"):
                feats = Counter()
                for ft in rec.features:
                    if ft.type == "source":
                        continue
                    parts = [(int(p.start), int(p.end), p.strand) for p in ft.location.parts]
                    feats[(ft.type, tuple(parts))] += 1
                recs[rec.id] = (str(rec.seq).upper(), feats)
        return recs

    try:
        c, r = load(c_files), load(r_files)
    except Exception as exc:  # unparsable GenBank is a defect, not a convention
        return Verdict(DISAGREEMENT, f"GenBank parse error: {exc}")
    if set(c) != set(r):
        return Verdict(DISAGREEMENT, "record IDs differ")
    for rid in c:
        if c[rid][0] != r[rid][0]:
            return Verdict(DISAGREEMENT, f"sequence differs in {rid}")
        if c[rid][1] != r[rid][1]:
            return Verdict(DISAGREEMENT, f"feature set differs in {rid}")
    return Verdict(STRUCTURAL, "rendering (line wrap, qualifier order)",
                   ["rendering"])


# ------------------------------------------------- T4 GFF3 -> GTF


def _gtf_attr(attrs: str, key: str) -> Optional[str]:
    m = re.search(rf'(?:^|;)\s*{key}\s+"([^"]*)"', attrs)
    return m.group(1) if m else None


def _merge(iv: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    out: List[List[int]] = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1] + 1:  # overlapping or abutting
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def _gtf_model(files: Sequence[Path]):
    tx = defaultdict(lambda: {"exon": [], "CDS": [], "chrom": None, "strand": None})
    for f in files:
        for line in _data_lines(f):
            p = line.split("\t")
            if len(p) < 9:
                continue
            tid = _gtf_attr(p[8], "transcript_id")
            if tid is None or p[2] not in ("exon", "CDS"):
                continue
            t = tx[tid]
            t["chrom"], t["strand"] = p[0], p[6]
            t[p[2]].append((int(p[3]), int(p[4])))
    return tx


def compare_t4(cand: Path, ref: Path, **_) -> Verdict:
    c_files, r_files = _files(cand, "*.gtf"), _files(ref, "*.gtf")
    if not c_files or not r_files:
        return _missing("GTF files")
    if _byte_equal(c_files, r_files):
        return Verdict(BYTE)
    conv = set()
    c_lines = [l for f in c_files for l in _data_lines(f)]
    r_lines = [l for f in r_files for l in _data_lines(f)]
    if any(l.split("\t")[2] == "gene" for l in c_lines) != \
            any(l.split("\t")[2] == "gene" for l in r_lines):
        conv.add("gene-rows")
    c, r = _gtf_model(c_files), _gtf_model(r_files)
    if set(c) != set(r):
        missing = sorted(set(r) ^ set(c))[:3]
        return Verdict(DISAGREEMENT, f"transcript sets differ (e.g. {missing})")
    for tid in c:
        a, b = c[tid], r[tid]
        if (a["chrom"], a["strand"]) != (b["chrom"], b["strand"]):
            return Verdict(DISAGREEMENT, f"{tid}: chrom/strand differ")
        for kind in ("exon", "CDS"):
            if sorted(a[kind]) == sorted(b[kind]):
                continue
            if _merge(a[kind]) == _merge(b[kind]):
                conv.add("exon-fusion")
                continue
            return Verdict(DISAGREEMENT, f"{tid}: {kind} coordinates differ")
    norm = lambda ls: sorted(re.sub(r";\s*$", "", l) for l in ls)  # noqa: E731
    if norm(c_lines) != norm(r_lines):
        conv.add("attribute-format")
    return Verdict(STRUCTURAL, ", ".join(sorted(conv)) or "ordering", sorted(conv))


# ------------------------------------------------- T5 GFF3 -> BED12


def _bed_rows(files: Sequence[Path]) -> List[List[str]]:
    return [l.split("\t") for f in files
            for l in _data_lines(f, skip=("#", "track", "browser"))]


def compare_t5(cand: Path, ref: Path, **_) -> Verdict:
    c_files = _files(cand, "*.bed12", "*.bed")
    r_files = _files(ref, "*.bed12", "*.bed")
    if not c_files or not r_files:
        return _missing("BED files")
    if _byte_equal(c_files, r_files):
        return Verdict(BYTE)
    c, r = _bed_rows(c_files), _bed_rows(r_files)
    if not c or not r:
        return Verdict(DISAGREEMENT, "one side has no BED rows")
    if len(r[0]) < 12 <= len(c[0]):
        return Verdict(LIMITATION, f"reference is BED{len(r[0])}, candidate BED12")
    conv = set()

    def key(row):
        blocks = tuple(int(x) for x in row[10].rstrip(",").split(",")) if len(row) >= 12 else ()
        starts = tuple(int(x) for x in row[11].rstrip(",").split(",")) if len(row) >= 12 else ()
        return (row[0], int(row[1]), int(row[2]), row[5] if len(row) > 5 else ".",
                blocks, starts)

    c_by = {}
    for row in c:
        c_by.setdefault(row[3] if len(row) > 3 else key(row), row)
    r_by = {}
    for row in r:
        r_by.setdefault(row[3] if len(row) > 3 else key(row), row)
    if set(c_by) != set(r_by):
        return Verdict(DISAGREEMENT, "feature names differ")
    for name, row in c_by.items():
        rr = r_by[name]
        if key(row) != key(rr):
            return Verdict(DISAGREEMENT, f"{name}: coordinates or blocks differ")
        if len(row) >= 8 and len(rr) >= 8 and (row[6], row[7]) != (rr[6], rr[7]):
            if int(rr[6]) == int(rr[7]) or int(row[6]) == int(row[7]):
                conv.add("thick-span-noncoding")
            else:
                return Verdict(DISAGREEMENT, f"{name}: thickStart/thickEnd differ")
        if len(row) > 4 and row[4] != rr[4]:
            conv.add("score")
        if len(row) > 8 and row[8] != rr[8]:
            conv.add("itemRgb")
    if [r_[3] for r_ in c] != [r_[3] for r_ in r]:
        conv.add("row-order")
    return Verdict(STRUCTURAL, ", ".join(sorted(conv)) or "format", sorted(conv))


# ------------------------------------------------- T6 HMMER tblout -> TSV


def compare_t6(cand: Path, ref: Path, **_) -> Verdict:
    c_files, r_files = _files(cand, "*.tsv"), _files(ref, "*.tsv")
    if not c_files or not r_files:
        return _missing("TSV files")
    if _byte_equal(c_files, r_files):
        return Verdict(BYTE)

    def rows(files):
        out = []
        for f in files:
            with open(f, newline="") as fh:
                out.extend(list(csv.reader(fh, delimiter="\t")))
        return out

    c, r = rows(c_files), rows(r_files)
    if len(c) != len(r):
        return Verdict(DISAGREEMENT, f"row count {len(c)} vs {len(r)}")
    conv = set()
    for i, (a, b) in enumerate(zip(c, r)):
        if len(a) != len(b):
            return Verdict(DISAGREEMENT, f"row {i}: column count differs")
        for x, y in zip(a, b):
            if x == y:
                continue
            try:
                if math.isclose(float(x), float(y), rel_tol=1e-6):
                    conv.add("number-format")
                    continue
            except ValueError:
                pass
            if x.strip() == y.strip():
                conv.add("whitespace")
                continue
            return Verdict(DISAGREEMENT, f"row {i}: {x!r} vs {y!r}")
    return Verdict(STRUCTURAL, ", ".join(sorted(conv)), sorted(conv))


# ------------------------------------------------- T7 orthogroups -> FASTA


def _fasta(path: Path) -> List[Tuple[str, str]]:
    recs, cur, seq = [], None, []
    for line in open(path):
        line = line.rstrip("\r\n")
        if line.startswith(">"):
            if cur is not None:
                recs.append((cur, "".join(seq).upper()))
            cur, seq = line[1:].split()[0] if line[1:] else "", []
        elif line:
            seq.append(line)
    if cur is not None:
        recs.append((cur, "".join(seq).upper()))
    return recs


def compare_t7(cand: Path, ref: Path, **_) -> Verdict:
    c_files = _files(cand, "*.fa", "*.fasta")
    r_files = _files(ref, "*.fa", "*.fasta")
    if not c_files or not r_files:
        return _missing("FASTA files")
    if _byte_equal(c_files, r_files):
        return Verdict(BYTE)
    if {f.name for f in c_files} != {f.name for f in r_files}:
        return Verdict(DISAGREEMENT, "orthogroup file sets differ")
    conv = set()
    r_map = {f.name: f for f in r_files}
    for f in c_files:
        a, b = _fasta(f), _fasta(r_map[f.name])
        if a == b:
            conv.add("line-wrap")
            continue
        if sorted(a) == sorted(b):
            conv.add("record-order")
            continue
        return Verdict(DISAGREEMENT, f"{f.name}: records differ")
    return Verdict(STRUCTURAL, ", ".join(sorted(conv)), sorted(conv))


# ------------------------------------------------- T8 pseudohaploid


def _vcf_gt_classes(vcf: Path):
    """Map (norm chrom, pos) -> list of gt classes for biallelic SNPs."""
    import gzip
    opener = gzip.open if open(vcf, "rb").read(2) == b"\x1f\x8b" else open
    sites = {}
    with opener(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if "," in p[4] or len(p[3]) != 1 or len(p[4]) != 1:
                continue
            cls = []
            for s in p[9:]:
                g = re.split(r"[/|]", s.split(":")[0])
                if any(a == "." for a in g):
                    cls.append("miss")
                elif all(a == "0" for a in g):
                    cls.append("ref")
                elif all(a == "1" for a in g):
                    cls.append("alt")
                else:
                    cls.append("het")
            sites[(_norm_chrom(p[0]), int(p[1]))] = (p[3].upper(), p[4].upper(), cls)
    return sites


def check_pseudohaploid(out_dir: Path, input_vcf: Path, alpha: float = 1e-6) -> Verdict:
    """Reference-free check of a pseudohaploid EIGENSTRAT against its input.

    Every call must be consistent with the source genotype (hom-ref -> the
    REF allele, hom-alt -> ALT, missing -> 9, het -> either), and het draws
    must be balanced: a two-sided binomial test of #REF picks among hets
    at ``alpha`` flags a biased sampler.
    """
    e = _read_eigenstrat(out_dir)
    if e is None or not Path(input_vcf).is_file():
        return _missing("pseudohaploid output or input VCF")
    snps, genos, _inds, _ = e
    src = _vcf_gt_classes(Path(input_vcf))
    het_ref = het_total = 0
    for s, g in zip(snps, genos):
        key = (_norm_chrom(s[1]), int(s[3]))
        if key not in src:
            return Verdict(DISAGREEMENT, f"site {s[1]}:{s[3]} not in input")
        ref_a, alt_a, cls = src[key]
        major = s[4].upper()
        if {s[4].upper(), s[5].upper()} != {ref_a, alt_a}:
            return Verdict(DISAGREEMENT, f"alleles differ at {s[1]}:{s[3]}")
        if len(g) != len(cls):
            return Verdict(DISAGREEMENT, "sample count differs from input")
        ref_code = "2" if major == ref_a else "0"
        alt_code = "0" if ref_code == "2" else "2"
        for call, c in zip(g, cls):
            expected = {"ref": ref_code, "alt": alt_code, "miss": "9"}.get(c)
            if expected is not None and call != expected:
                return Verdict(DISAGREEMENT, f"call {call} inconsistent with {c} at {s[1]}:{s[3]}")
            if c == "het":
                if call not in "02":
                    return Verdict(DISAGREEMENT, f"het drawn as {call} at {s[1]}:{s[3]}")
                het_total += 1
                het_ref += call == ref_code
    if het_total:
        p = _binom_two_sided(het_ref, het_total)
        if p < alpha:
            return Verdict(DISAGREEMENT,
                           f"het draws biased: {het_ref}/{het_total} REF (p={p:.2g})")
    return Verdict(STRUCTURAL, f"stochastic draw; {het_ref}/{het_total} het draws REF",
                   ["stochastic-draw"])


def _binom_two_sided(k: int, n: int) -> float:
    # Normal approximation with continuity correction; adequate for the
    # large n of real runs and conservative for small n.
    if n == 0:
        return 1.0
    mu, sd = n / 2.0, math.sqrt(n) / 2.0
    z = (abs(k - mu) - 0.5) / sd
    return math.erfc(max(z, 0.0) / math.sqrt(2))


def compare_t8(cand: Path, ref: Path, input_vcf: Optional[Path] = None, **_) -> Verdict:
    if input_vcf is None:
        return Verdict(LIMITATION, "input VCF needed for the pseudohaploid check")
    return check_pseudohaploid(cand, input_vcf)


COMPARATORS = {
    "T1": compare_t1, "T2": compare_t2, "T3": compare_t3, "T4": compare_t4,
    "T5": compare_t5, "T6": compare_t6, "T7": compare_t7, "T8": compare_t8,
}
