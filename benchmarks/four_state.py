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


class _Mismatch(Exception):
    """A structural defect found while streaming (row counts, .bed length)."""


def _same_bytes(cand: Sequence[Path], ref: Sequence[Path]) -> bool:
    """Byte equality (CRLF = LF for text) without reading whole files."""
    import filecmp
    from itertools import zip_longest
    if len(cand) != len(ref) or not cand:
        return False
    for a, b in zip(cand, ref):
        if a.suffix == ".bed":
            if not filecmp.cmp(a, b, shallow=False):
                return False
            continue
        with open(a, "rb") as fa, open(b, "rb") as fb:
            for la, lb in zip_longest(fa, fb):
                if la is None or lb is None or la.rstrip(b"\r\n") != lb.rstrip(b"\r\n"):
                    return False
    return True


def _pair_stream(c_items, r_items, key, loose):
    """Pair two record streams, holding only out-of-step records in memory.

    Items are ``(index, fields, payload)``. Records are paired on ``key``
    (chrom, pos and allele set), so the two biallelic records 1000 Genomes
    writes for a split multi-allelic site are told apart. When both files
    list sites in the same order, which converters of one VCF do, nothing
    is buffered and memory stays flat on whole chromosomes. Records left
    over are paired on ``loose`` (chrom, pos) so allele mismatches and
    monomorphic '0' alleles are still reported; any remaining ones come out
    as ``(c, None)`` or ``(None, r)``.
    """
    from collections import deque
    from itertools import zip_longest
    c_pend: Dict = {}
    r_pend: Dict = {}

    def pop(d, k):
        q = d[k]
        x = q.popleft()
        if not q:
            del d[k]
        return x

    for c, r in zip_longest(c_items, r_items):
        if c is not None and r is not None:
            kc, kr = key(c), key(r)
            if kc == kr and kc not in c_pend and kr not in r_pend:
                yield c, r
                continue
        else:
            kc = key(c) if c is not None else None
            kr = key(r) if r is not None else None
        if c is not None:
            if kc in r_pend:
                yield c, pop(r_pend, kc)
            else:
                c_pend.setdefault(kc, deque()).append(c)
        if r is not None:
            if kr in c_pend:
                yield pop(c_pend, kr), r
            else:
                r_pend.setdefault(kr, deque()).append(r)
    r_loose: Dict = {}
    for q in r_pend.values():
        for r in q:
            r_loose.setdefault(loose(r), deque()).append(r)
    for q in c_pend.values():
        for c in q:
            lk = loose(c)
            yield c, (pop(r_loose, lk) if lk in r_loose else None)
    for q in r_loose.values():
        for r in q:
            yield None, r


def _site_key(ci: int, pi: int, a1: int, a2: int):
    def key(item):
        f = item[1]
        return (_norm_chrom(f[ci]), int(f[pi]), frozenset((f[a1].upper(), f[a2].upper())))

    def loose(item):
        f = item[1]
        return (_norm_chrom(f[ci]), int(f[pi]))
    return key, loose


def _in_order(pairs_idx) -> bool:
    """True if reference indices, taken in candidate order, only increase."""
    last = -1
    for _ci, ri in sorted(pairs_idx):
        if ri < last:
            return False
        last = ri
    return True


def _count_data_lines(files: Sequence[Path]) -> int:
    n = 0
    for f in files:
        with open(f, "rb") as fh:
            for line in fh:
                if line.strip() and not line.startswith(b"#"):
                    n += 1
    return n


# ------------------------------------------------- T1 VCF -> EIGENSTRAT


def _read_eigenstrat(d: Path):
    """Whole-file reader, kept for small inputs (tests, mutation helpers)."""
    geno = _files(d, "*.geno")
    snp = _files(d, "*.snp")
    ind = _files(d, "*.ind")
    if not (geno and snp and ind):
        return None
    snps = [l.split() for f in snp for l in _data_lines(f)]
    genos = [l.strip() for f in geno for l in _data_lines(f)]
    inds = [l.split()[0] for f in ind for l in _data_lines(f)]
    return snps, genos, inds, geno + snp + ind


def _eig_items(snp_files: Sequence[Path], geno_files: Sequence[Path]):
    """Yield (index, .snp fields, .geno row bytes), one site at a time."""
    idx = 0
    if len(snp_files) != len(geno_files):
        raise _Mismatch(".snp and .geno file counts differ")
    for sf, gf in zip(snp_files, geno_files):
        with open(sf) as fs_, open(gf, "rb") as fg:
            for line in fs_:
                if not line.strip() or line.startswith("#"):
                    continue
                g = fg.readline()
                if not g:
                    raise _Mismatch(".snp and .geno row counts differ")
                yield idx, line.split(), g.rstrip(b"\r\n")
                idx += 1
            if fg.readline().strip():
                raise _Mismatch(".snp and .geno row counts differ")


_FLIP02 = bytes.maketrans(b"02", b"20")


def _flip_geno(row):
    if isinstance(row, bytes):
        return row.translate(_FLIP02)
    return row.translate(str.maketrans("02", "20"))


def compare_t1(cand: Path, ref: Path, **_) -> Verdict:
    """T1 comparator. Streams both EIGENSTRAT sets; memory does not grow
    with the number of sites when both list them in the same order."""
    c_snp, c_geno, c_ind = (_files(cand, p) for p in ("*.snp", "*.geno", "*.ind"))
    r_snp, r_geno, r_ind = (_files(ref, p) for p in ("*.snp", "*.geno", "*.ind"))
    if not (c_snp and c_geno and c_ind and r_snp and r_geno and r_ind):
        return _missing("EIGENSTRAT triplet")
    if _same_bytes(c_geno + c_snp + c_ind, r_geno + r_snp + r_ind):
        return Verdict(BYTE)
    ci = [l.split()[0] for f in c_ind for l in _data_lines(f)]
    ri = [l.split()[0] for f in r_ind for l in _data_lines(f)]
    if ci != ri:
        if sorted(ci) != sorted(ri):
            return Verdict(DISAGREEMENT, "sample sets differ")
        return Verdict(DISAGREEMENT, "sample order differs")
    nc, nr = _count_data_lines(c_snp), _count_data_lines(r_snp)
    conv = set()
    key, loose = _site_key(1, 3, 4, 5)
    idx_pairs = []
    # Sites present on one side only are collected (not a stop at the first)
    # so the verdict names them: "site count" alone does not say which tool
    # dropped what.
    only_c: List[str] = []
    only_r: List[str] = []
    try:
        for c, r in _pair_stream(_eig_items(c_snp, c_geno), _eig_items(r_snp, r_geno),
                                 key, loose):
            if r is None or c is None:
                f = (c or r)[1]
                (only_r if c is None else only_c).append(f"{f[1]}:{f[3]} {f[4]}/{f[5]}")
                continue
            (i, s, g), (j, rs, rg) = c, r
            idx_pairs.append((i, j))
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
                if g.translate(_FLIP02) != rg:
                    return Verdict(DISAGREEMENT,
                                   f"genotypes differ at {s[1]}:{s[3]} (after allele swap)")
            else:
                return Verdict(DISAGREEMENT, f"alleles differ at {s[1]}:{s[3]}")
    except _Mismatch as e:
        return Verdict(DISAGREEMENT, str(e))
    if only_c and not only_r:
        # EIGENSOFT convertf keeps one record per position: of two biallelic
        # records at one site (a split multi-allelic, e.g. 1000 Genomes
        # chr22:43236421 G>A and G>T) it silently drops the second. A
        # candidate site counts as such a record when the reference has a
        # different record at the same position; that is a reference
        # limitation, not a candidate error.
        positions = {(_norm_chrom(x.split(":")[0]), int(x.split(":")[1].split()[0]))
                     for x in only_c}
        ref_pos = set()
        for f in r_snp:
            with open(f) as fh:
                for line in fh:
                    t = line.split()
                    if len(t) > 3 and (_norm_chrom(t[1]), int(t[3])) in positions:
                        ref_pos.add((_norm_chrom(t[1]), int(t[3])))
        if ref_pos == positions:
            conv.add("reference-drops-colocated")
            if not _in_order(idx_pairs):
                conv.add("site-order")
            return Verdict(STRUCTURAL,
                           f"reference keeps one record per position; {len(only_c)} "
                           f"co-located record(s) only in candidate (e.g. "
                           f"{', '.join(only_c[:3])})", sorted(conv))
    if only_c or only_r:
        bits = [f"site count {nc} vs reference {nr}"]
        if only_c:
            bits.append(f"{len(only_c)} only in candidate (e.g. {', '.join(only_c[:3])})")
        if only_r:
            bits.append(f"{len(only_r)} only in reference (e.g. {', '.join(only_r[:3])})")
        return Verdict(DISAGREEMENT, "; ".join(bits))
    if not _in_order(idx_pairs):
        conv.add("site-order")
    return Verdict(STRUCTURAL, ", ".join(sorted(conv)) or "whitespace/format",
                   sorted(conv))


# ------------------------------------------------- T2 VCF -> PLINK .bed


def _read_plink(d: Path):
    """Whole-file reader, kept for small inputs (tests, mutation helpers)."""
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
    code = {0: 2, 2: 1, 3: 0, 1: None}
    genos = []
    for i in range(len(bims)):
        chunk = body[i * per:(i + 1) * per]
        genos.append([code[(chunk[j // 4] >> (2 * (j % 4))) & 3] for j in range(n)])
    return bims, fams, genos, bed + bim + fam


def _plink_items(bim: Path, bed: Path, per: int):
    """Yield (index, .bim fields, raw .bed bytes for the variant)."""
    with open(bim) as fb, open(bed, "rb") as fd:
        fd.read(3)
        idx = 0
        for line in fb:
            if not line.strip():
                continue
            chunk = fd.read(per)
            if len(chunk) != per:
                raise _Mismatch("malformed .bed (bad-length)")
            yield idx, line.split(), chunk
            idx += 1


def _dosage(chunk: bytes, n: int):
    """Per-sample A1 dosage: 00 hom A1 -> 2, 10 het -> 1, 11 hom A2 -> 0,
    01 missing -> -1."""
    import numpy as np
    a = np.frombuffer(chunk, dtype=np.uint8)
    codes = np.stack([(a >> s) & 3 for s in (0, 2, 4, 6)], axis=1).reshape(-1)[:n]
    return np.array([2, -1, 1, 0], dtype=np.int8)[codes]


def _swap_dosage(x):
    import numpy as np
    return np.where(x < 0, x, 2 - x)


def compare_t2(cand: Path, ref: Path, **_) -> Verdict:
    """T2 comparator. Streams the .bed files variant by variant."""
    import numpy as np
    cb, cm, cf = _files(cand, "*.bed"), _files(cand, "*.bim"), _files(cand, "*.fam")
    rb, rm, rf = _files(ref, "*.bed"), _files(ref, "*.bim"), _files(ref, "*.fam")
    if not (cb and cm and cf and rb and rm and rf):
        return _missing("PLINK .bed/.bim/.fam")
    c_fam = [l.split() for l in _data_lines(cf[0])]
    r_fam = [l.split() for l in _data_lines(rf[0])]
    n = len(c_fam)
    per = (n + 3) // 4
    for bed, bim, fam in ((cb[0], cm[0], c_fam), (rb[0], rm[0], r_fam)):
        with open(bed, "rb") as fh:
            if fh.read(3) != b"\x6c\x1b\x01":
                return Verdict(DISAGREEMENT, "malformed .bed (bad-magic)")
        p = (len(fam) + 3) // 4
        if bed.stat().st_size - 3 != p * _count_data_lines([bim]):
            return Verdict(DISAGREEMENT, "malformed .bed (bad-length)")
    if _same_bytes(cb + cm + cf, rb + rm + rf):
        return Verdict(BYTE)
    conv = set()
    if [f[1] for f in c_fam] != [f[1] for f in r_fam]:
        return Verdict(DISAGREEMENT, "sample IDs or order differ in .fam")
    if [f[:2] for f in c_fam] != [f[:2] for f in r_fam]:
        conv.add("family-id")
    key, loose = _site_key(0, 3, 4, 5)
    idx_pairs = []
    extra_snp = extra_other = 0
    try:
        for c, r in _pair_stream(_plink_items(cm[0], cb[0], per),
                                 _plink_items(rm[0], rb[0], per), key, loose):
            if r is None:
                b = c[1]
                return Verdict(DISAGREEMENT, f"variant {b[0]}:{b[3]} absent in reference")
            if c is None:
                rb_ = r[1]
                if len(rb_[4]) == 1 and len(rb_[5]) == 1:
                    extra_snp += 1
                else:
                    extra_other += 1
                continue
            (i, b, cg_raw), (j, rb_, rg_raw) = c, r
            idx_pairs.append((i, j))
            if b[0] != rb_[0]:
                conv.add("chrom-label")
            if b[1] != rb_[1]:
                conv.add("variant-id")
            ca, ra = (b[4], b[5]), (rb_[4], rb_[5])
            cg, rg = _dosage(cg_raw, n), _dosage(rg_raw, n)
            if ca == ra:
                same = np.array_equal(cg, rg)
            elif ca == (ra[1], ra[0]):
                conv.add("allele-order")
                same = np.array_equal(_swap_dosage(cg), rg)
            elif "0" in ra and (ca[0] in ra or ca[1] in ra):
                # Reference wrote a monomorphic site with A1 = '0'.
                conv.add("monomorphic-allele-code")
                same = np.array_equal(cg if ca[1] == ra[1] else _swap_dosage(cg), rg)
            else:
                return Verdict(DISAGREEMENT, f"alleles differ at {b[0]}:{b[3]}")
            if not same:
                return Verdict(DISAGREEMENT, f"genotypes differ at {b[0]}:{b[3]}")
    except _Mismatch as e:
        return Verdict(DISAGREEMENT, str(e))
    if extra_snp:
        return Verdict(DISAGREEMENT, f"{extra_snp} reference SNPs missing from candidate")
    if extra_other:
        conv.add("indels-excluded")
    if not _in_order(idx_pairs):
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

    # Load each side on its own so the verdict says whose file is broken.
    try:
        r = load(r_files)
    except Exception as exc:
        return Verdict(LIMITATION, f"reference GenBank unreadable: {exc}")
    try:
        c = load(c_files)
    except Exception as exc:  # unparsable GenBank is a defect, not a convention
        return Verdict(DISAGREEMENT, f"GenBank parse error: {exc}")
    if set(c) != set(r):
        return Verdict(DISAGREEMENT, "record IDs differ")
    # Biopython's GenBank writer drops sub_features, so a BCBio-based
    # reference keeps only top-level features (tRNA, gene) while gfc and
    # seqret also write the child exon/mRNA/CDS rows. Feature types the
    # reference never writes are reported, not judged; every feature type
    # the reference does write must match exactly.
    r_types = {t for rid in r for (t, _loc) in r[rid][1]}
    c_types = {t for rid in c for (t, _loc) in c[rid][1]}
    extra_types = sorted(c_types - r_types)
    for rid in c:
        if c[rid][0] != r[rid][0]:
            return Verdict(DISAGREEMENT, f"sequence differs in {rid}")
        cf = Counter({k: v for k, v in c[rid][1].items() if k[0] in r_types})
        if cf != r[rid][1]:
            return Verdict(DISAGREEMENT, f"feature set differs in {rid}")
    conv = ["rendering"]
    reason = "rendering (line wrap, qualifier order)"
    if extra_types:
        conv.append("reference-omits-children")
        reason += f"; reference writes no {', '.join(extra_types)} features (not checked)"
    return Verdict(STRUCTURAL, reason, sorted(conv))


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


_ID_PREFIX = re.compile(r"^(?:transcript|gene|mRNA|rna|CDS|exon):")


def _set_diff_reason(what: str, c: set, r: set) -> str:
    only_c, only_r = sorted(c - r), sorted(r - c)
    bits = [f"{what} differ"]
    if only_c:
        bits.append(f"{len(only_c)} only in candidate (e.g. {only_c[:3]})")
    if only_r:
        bits.append(f"{len(only_r)} only in reference (e.g. {only_r[:3]})")
    return "; ".join(bits)


def _strip_prefix_keys(d: dict) -> Optional[dict]:
    """Re-key by ID without an Ensembl type prefix ('transcript:ENST...' ->
    'ENST...'); None if that would merge two keys."""
    out = {}
    for k, v in d.items():
        k2 = _ID_PREFIX.sub("", k)
        if k2 in out:
            return None
        out[k2] = v
    return out


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
        c2, r2 = _strip_prefix_keys(c), _strip_prefix_keys(r)
        if c2 is not None and r2 is not None and set(c2) == set(r2):
            conv.add("id-prefix")
            c, r = c2, r2
        else:
            return Verdict(DISAGREEMENT, _set_diff_reason("transcript sets", set(c), set(r)))
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
        c2, r2 = _strip_prefix_keys(c_by), _strip_prefix_keys(r_by)
        if c2 is not None and r2 is not None and set(c2) == set(r2):
            conv.add("id-prefix")
            c_by, r_by = c2, r2
        else:
            # Tools differ in which attribute names a row (ID vs Name);
            # match rows on coordinates and blocks instead.
            ck = Counter(key(row) for row in c)
            rk = Counter(key(row) for row in r)
            # Only when the naming differs throughout: one renamed row in an
            # otherwise matching file is a defect, not a convention.
            unnamed = len(set(c_by) - set(r_by)) / max(len(c_by), 1)
            if ck == rk and max(ck.values()) == 1 and unnamed >= 0.9:
                conv.add("name-attribute")
                c_by = {key(row): row for row in c}
                r_by = {key(row): row for row in r}
            else:
                return Verdict(DISAGREEMENT,
                               _set_diff_reason("feature names", set(c_by), set(r_by)))
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
    if "name-attribute" in conv:
        if [key(x) for x in c] != [key(x) for x in r]:
            conv.add("row-order")
    elif [_ID_PREFIX.sub("", r_[3]) for r_ in c] != [_ID_PREFIX.sub("", r_[3]) for r_ in r]:
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


def _iter_vcf_sites(vcf: Path):
    """Yield (chrom, pos, REF, ALT, gt) for biallelic SNPs, streaming.

    ``gt`` is an int array per sample: 0 hom-ref, 1 het, 2 missing, 3
    hom-alt (cyvcf2's gt_types coding). Uses cyvcf2 when installed, else a
    pure-Python parser; neither holds more than one record in memory.
    """
    import numpy as np
    try:
        import cyvcf2
    except ImportError:  # pragma: no cover - exercised only without cyvcf2
        cyvcf2 = None
    if cyvcf2 is not None:
        with cyvcf2.VCF(str(vcf)) as fh:
            for v in fh:
                if len(v.ALT) != 1 or len(v.REF) != 1 or len(v.ALT[0]) != 1:
                    continue
                yield v.CHROM, v.POS, v.REF.upper(), v.ALT[0].upper(), v.gt_types
        return
    import gzip
    opener = gzip.open if open(vcf, "rb").read(2) == b"\x1f\x8b" else open
    with opener(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if "," in p[4] or len(p[3]) != 1 or len(p[4]) != 1:
                continue
            gt = []
            for smp in p[9:]:
                a = re.split(r"[/|]", smp.split(":")[0])
                if any(x == "." for x in a):
                    gt.append(2)
                elif all(x == "0" for x in a):
                    gt.append(0)
                elif all(x == "1" for x in a):
                    gt.append(3)
                else:
                    gt.append(1)
            yield p[0], int(p[1]), p[3].upper(), p[4].upper(), np.array(gt)


def check_pseudohaploid(out_dir: Path, input_vcf: Path, alpha: float = 1e-6) -> Verdict:
    """Reference-free check of a pseudohaploid EIGENSTRAT against its input.

    Every call must be consistent with the source genotype (hom-ref -> the
    REF allele, hom-alt -> ALT, missing -> 9, het -> either), and het draws
    must be balanced: a two-sided binomial test of #REF picks among hets
    at ``alpha`` flags a biased sampler.

    The output and the VCF are read in lockstep, one site at a time, with
    numpy per site, so memory stays flat on whole-chromosome inputs (the
    earlier version held every genotype of the VCF as Python strings,
    about 20 GB for 1000 Genomes chr22). Output sites must appear in input
    order; input sites the converter filtered out are skipped. Sites are
    matched on chrom, position and allele set, so split multi-allelic
    records at one position are told apart.
    """
    import numpy as np
    out_dir = Path(out_dir)
    geno = _files(out_dir, "*.geno")
    snp = _files(out_dir, "*.snp")
    ind = _files(out_dir, "*.ind")
    if not (geno and snp and ind) or not Path(input_vcf).is_file():
        return _missing("pseudohaploid output or input VCF")
    n_ind = len(_data_lines(ind[0]))
    src = _iter_vcf_sites(Path(input_vcf))
    het_ref = het_total = 0
    n_sites = 0
    with open(snp[0]) as sf, open(geno[0], "rb") as gf:
        for sline in sf:
            if not sline.strip() or sline.startswith("#"):
                continue
            s = sline.split()
            g = np.frombuffer(gf.readline().rstrip(b"\r\n"), dtype=np.uint8)
            chrom, pos = _norm_chrom(s[1]), int(s[3])
            major, minor = s[4].upper(), s[5].upper()
            for vchrom, vpos, ref_a, alt_a, gt in src:
                if (_norm_chrom(vchrom) == chrom and vpos == pos
                        and {ref_a, alt_a} == {major, minor}):
                    break
            else:
                return Verdict(DISAGREEMENT,
                               f"site {s[1]}:{s[3]} not in input (or out of order)")
            if len(g) != len(gt) or len(g) != n_ind:
                return Verdict(DISAGREEMENT,
                               f"sample count differs from input at {s[1]}:{s[3]}")
            ref_code = ord("2") if major == ref_a else ord("0")
            alt_code = ord("0") if ref_code == ord("2") else ord("2")
            for cls, want in ((0, ref_code), (3, alt_code), (2, ord("9"))):
                bad = (gt == cls) & (g != want)
                if bad.any():
                    call = chr(g[bad.argmax()])
                    name = {0: "ref", 3: "alt", 2: "miss"}[cls]
                    return Verdict(DISAGREEMENT,
                                   f"call {call} inconsistent with {name} at {s[1]}:{s[3]}")
            het = gt == 1
            hg = g[het]
            if ((hg != ord("0")) & (hg != ord("2"))).any():
                return Verdict(DISAGREEMENT, f"het drawn as non-0/2 at {s[1]}:{s[3]}")
            het_total += int(het.sum())
            het_ref += int((hg == ref_code).sum())
            n_sites += 1
        if gf.readline().strip():
            return Verdict(DISAGREEMENT, ".geno has more rows than .snp")
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
