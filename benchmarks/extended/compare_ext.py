"""Four-state comparators for the extended benchmark (E01-E24).

Same contract as ``benchmarks/four_state.py``: each comparator takes the
candidate's output directory and the reference (another tool's output
directory, or the task's input directory for a reference-free check) and
returns ``Verdict(state, reason, conventions)``. A structural verdict lists
every convention it needed; anything else that differs is a disagreement,
and the reason names the first difference.

Files are paired by stem (``reads.fasta`` with ``reads.fasta``). Records are
paired by identifier, so record order is a convention, not an error.
"""
from __future__ import annotations

import gzip
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from four_state import (BYTE, DISAGREEMENT, LIMITATION, STRUCTURAL, Verdict,
                        _merge, _same_bytes)

# ------------------------------------------------------------------ helpers


def _files(d: Path, exts: Sequence[str]) -> Dict[str, Path]:
    """Map stem -> file for files in ``d`` with one of ``exts``."""
    out: Dict[str, Path] = {}
    if not Path(d).is_dir():
        return out
    for p in sorted(Path(d).iterdir()):
        if p.name.startswith("."):
            continue
        for e in sorted(exts, key=len, reverse=True):
            if p.name.endswith(e):
                out[p.name[: -len(e)]] = p
                break
    return out


def _open(p: Path):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p, encoding="latin-1")


def _pairs(cand: Path, ref: Path, cext, rext=None):
    """Pair candidate and reference files by stem. Returns (pairs, verdict)
    where verdict is set when the file sets cannot be paired."""
    c, r = _files(cand, cext), _files(ref, rext or cext)
    if not r:
        return None, Verdict(LIMITATION, f"no reference files ({', '.join(rext or cext)})")
    if not c:
        return None, Verdict(DISAGREEMENT, f"no candidate output ({', '.join(cext)})")
    if set(c) != set(r):
        miss = sorted(set(r) - set(c))[:3]
        extra = sorted(set(c) - set(r))[:3]
        return None, Verdict(DISAGREEMENT, f"file sets differ: missing {miss}, extra {extra}")
    return [(c[k], r[k]) for k in sorted(r)], None


class Diff(Exception):
    """A content difference no convention explains."""


class Conv:
    """Collects the conventions a comparison needed."""

    def __init__(self):
        self.used: List[str] = []

    def add(self, name: str) -> None:
        if name not in self.used:
            self.used.append(name)


def _finish(conv: Conv, byte_equal: bool, what: str = "") -> Verdict:
    if byte_equal:
        return Verdict(BYTE, "")
    if not conv.used:
        conv.add("rendering")
    return Verdict(STRUCTURAL, what, list(conv.used))


def _run(fn, cand, ref, cext, rext=None, **kw) -> Verdict:
    pairs, v = _pairs(cand, ref, cext, rext)
    if v:
        return v
    conv = Conv()
    try:
        notes = [fn(c, r, conv, **kw) for c, r in pairs]
    except Diff as exc:
        return Verdict(DISAGREEMENT, str(exc))
    byte = all(_same_bytes([c], [r]) for c, r in pairs)
    return _finish(conv, byte, "; ".join(n for n in notes if n))


def _compare_maps(c: dict, r: dict, what: str, conv: Conv, eq=None, label="record") -> None:
    """Compare two id -> value maps; ``eq(cv, rv, conv)`` may record
    conventions and returns False on a real difference."""
    if set(c) != set(r):
        miss = sorted(set(r) - set(c))
        extra = sorted(set(c) - set(r))
        raise Diff(f"{what}: {len(miss)} {label}s only in reference (e.g. {miss[:3]}), "
                   f"{len(extra)} only in candidate (e.g. {extra[:3]})")
    n_bad, first = 0, None
    for k in r:
        ok = eq(c[k], r[k], conv) if eq else c[k] == r[k]
        if not ok:
            n_bad += 1
            first = first or k
    if n_bad:
        raise Diff(f"{what}: {n_bad} of {len(r)} {label}s differ (first: {first}: "
                   f"candidate {str(c[first])[:80]!r} vs reference {str(r[first])[:80]!r})")


# ------------------------------------------------------------- sequences


def _read_fasta(p: Path) -> List[Tuple[str, str, str]]:
    out, head, buf = [], None, []
    with _open(p) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if head is not None:
                    out.append((head.split()[0] if head.split() else "", head, "".join(buf)))
                head, buf = line[1:].strip(), []
            elif head is not None:
                buf.append(line.strip())
    if head is not None:
        out.append((head.split()[0] if head.split() else "", head, "".join(buf)))
    return out


def _id_map(recs, what: str, conv: Conv, order_ref=None) -> dict:
    m = {}
    for rid, head, val in recs:
        if rid in m:
            raise Diff(f"{what}: duplicate record id {rid!r}")
        m[rid] = (head, val)
    return m


def _cmp_records(c_recs, r_recs, what: str, conv: Conv, seq_eq) -> None:
    c, r = _id_map(c_recs, what, conv), _id_map(r_recs, what, conv)
    if [x[0] for x in c_recs] != [x[0] for x in r_recs] and set(c) == set(r):
        conv.add("record-order")
    if any(c[k][0] != r[k][0] for k in r if k in c):
        conv.add("description")
    _compare_maps({k: v[1] for k, v in c.items()}, {k: v[1] for k, v in r.items()},
                  what, conv, seq_eq)


def _seq_eq(a: str, b: str, conv: Conv) -> bool:
    if a == b:
        return True
    if a.upper() == b.upper():
        conv.add("case")
        return True
    return False


def _fasta_pair(c: Path, r: Path, conv: Conv, **_) -> str:
    _cmp_records(_read_fasta(c), _read_fasta(r), c.name, conv, _seq_eq)
    return ""


def compare_fasta(cand, ref, **_):
    return _run(_fasta_pair, cand, ref, [".fasta", ".fa", ".fas"])


def _read_fastq(p: Path):
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    with _open(p) as fh:
        try:
            return [(t.split()[0] if t.split() else "", t, (s, q))
                    for t, s, q in FastqGeneralIterator(fh)]
        except ValueError as exc:
            raise Diff(f"{p.name}: not valid FASTQ ({exc})")


def _fastq_pair(c: Path, r: Path, conv: Conv, fixed_quality_ok: str = "", **_) -> str:
    def eq(a, b, conv):
        if not _seq_eq(a[0], b[0], conv):
            return False
        if a[1] == b[1]:
            return True
        # FASTA -> FASTQ has no input quality, so each tool invents one.
        if fixed_quality_ok and len(set(a[1])) <= 1 and len(set(b[1])) <= 1 \
                and len(a[1]) == len(b[1]):
            conv.add("fixed-quality-value")
            return True
        return False
    _cmp_records(_read_fastq(c), _read_fastq(r), c.name, conv, eq)
    return ""


def compare_fastq(cand, ref, **kw):
    return _run(_fastq_pair, cand, ref, [".fastq", ".fq"], **kw)


def _qual_values(p: Path, head: str, buf) -> tuple:
    text = " ".join(buf).split()
    try:
        return tuple(int(x) for x in text)
    except ValueError:
        raise Diff(f"{p.name}: record {head.split()[0]} has non-numeric quality "
                   f"{' '.join(text)[:20]!r} (QUAL holds Phred scores as integers)")


def _read_qual(p: Path):
    out, head, buf = [], None, []
    with _open(p) as fh:
        for line in fh:
            if line.startswith(">"):
                if head is not None:
                    out.append((head.split()[0], head, _qual_values(p, head, buf)))
                head, buf = line[1:].strip(), []
            elif head is not None:
                buf.append(line)
    if head is not None:
        out.append((head.split()[0], head, _qual_values(p, head, buf)))
    return out


def compare_fasta_qual(cand, ref, **_):
    v1 = _run(_fasta_pair, cand, ref, [".fasta", ".fa"])
    if v1.state == DISAGREEMENT or v1.state == LIMITATION:
        return Verdict(v1.state, "FASTA part: " + v1.reason)

    def qual_pair(c, r, conv, **_):
        _cmp_records(_read_qual(c), _read_qual(r), c.name, conv, None)
        return ""
    v2 = _run(qual_pair, cand, ref, [".qual"])
    if v2.state in (DISAGREEMENT, LIMITATION):
        return Verdict(v2.state, "QUAL part: " + v2.reason)
    if v1.state == BYTE and v2.state == BYTE:
        return v1
    convs = list(dict.fromkeys(v1.conventions + v2.conventions))
    return Verdict(STRUCTURAL, "", convs or ["rendering"])


def _table_pair(c: Path, r: Path, conv: Conv, **_) -> str:
    def rows(p):
        out = []
        with _open(p) as fh:
            for line in fh:
                f = line.rstrip("\r\n").split("\t")
                if len(f) >= 2:
                    out.append(f)
        return out
    cr, rr = rows(c), rows(r)
    for side in (cr, rr):
        if side and side[0][0].lower() in ("id", "name", "#id") and side[0][1].lower().startswith("seq"):
            del side[0]
            conv.add("header-row")
    def recs(rows_, side):
        out = []
        for f in rows_:
            name = f[0]
            if len(name.split()) > 1:
                conv.add("full-header-as-name")
            out.append((name.split()[0] if name.split() else "", name, f[1]))
        return out
    if any(len(f) > 2 and all(x == "" for x in f[2:]) for f in rr + cr):
        conv.add("trailing-empty-column")
    _cmp_records(recs(cr, "c"), recs(rr, "r"), c.name, conv, _seq_eq)
    return ""


def compare_seq_table(cand, ref, **_):
    return _run(_table_pair, cand, ref, [".tsv", ".txt", ".tab"])


def _bed6_pair(c: Path, r: Path, conv: Conv, **_) -> str:
    def rows(p):
        out = Counter()
        with _open(p) as fh:
            for line in fh:
                if line.startswith(("#", "track", "browser")) or not line.strip():
                    continue
                f = line.rstrip("\r\n").split("\t")
                out[tuple(f[:6])] += 1
        return out
    cr, rr = rows(c), rows(r)
    if cr == rr:
        return ""
    # Read names: some tools append /1 or /2.
    strip = lambda cnt: Counter({(a, b, cc, re.sub(r"/[12]$", "", d), e, g): n
                                 for (a, b, cc, d, e, g), n in cnt.items()})
    if strip(cr) == strip(rr):
        conv.add("read-name-suffix")
        return ""
    only_r = rr - cr
    only_c = cr - rr
    raise Diff(f"{c.name}: {sum(only_r.values())} rows only in reference (e.g. {list(only_r)[:2]}), "
               f"{sum(only_c.values())} only in candidate (e.g. {list(only_c)[:2]})")


def compare_bed6(cand, ref, **_):
    return _run(_bed6_pair, cand, ref, [".bed"])


# ------------------------------------------------------------- annotation


def _attrs(s: str) -> Dict[str, str]:
    out = {}
    for part in s.strip().strip(";").split(";"):
        part = part.strip()
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def _gff_lines(p: Path):
    with _open(p) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                if line.startswith("##FASTA"):
                    break
                continue
            f = line.rstrip("\r\n").split("\t")
            if len(f) == 9:
                yield f


_PREFIX = re.compile(r"^(?:transcript|gene|mRNA|rna|CDS|cds|exon|agat-[a-z]+-\d+):")


def _norm_id(x: str, conv: Conv) -> str:
    y = _PREFIX.sub("", x)
    if y != x:
        conv.add("id-prefix")
    return y


def _gff3_models(p: Path, conv: Conv):
    """Transcript id -> (seqid, strand, merged exons, merged CDS)."""
    exons, cds, info = defaultdict(list), defaultdict(list), {}
    for f in _gff_lines(p):
        a = _attrs(f[8])
        typ = f[2]
        if typ in ("exon", "CDS"):
            for par in a.get("Parent", "").split(","):
                if par:
                    (exons if typ == "exon" else cds)[par].append((int(f[3]), int(f[4])))
                    info.setdefault(par, (f[0], f[6]))
    out = {}
    for tid in set(exons) | set(cds):
        ex, cd = exons.get(tid, []), cds.get(tid, [])
        m_ex, m_cd = _merge(ex), _merge(cd)
        if len(m_ex) < len(ex) or len(m_cd) < len(cd):
            conv.add("abutting-exons-merged")
        out[_norm_id(tid, conv)] = (info[tid][0], info[tid][1], tuple(m_ex or m_cd), tuple(m_cd))
    return out


def _gff3_models_pair(c: Path, r: Path, conv: Conv, **_) -> str:
    cm, rm = _gff3_models(c, conv), _gff3_models(r, conv)
    _compare_maps(cm, rm, c.name, conv, label="transcript")
    return f"{len(rm)} transcripts"


def compare_gff3_models(cand, ref, **_):
    return _run(_gff3_models_pair, cand, ref, [".gff3", ".gff"])


def _input_features(gff: Path):
    """Every feature line of an input GFF3 as (seqid, type, start, end, strand, ID)."""
    return [(f[0], f[2], int(f[3]), int(f[4]), f[6] if f[6] in "+-" else ".",
             _attrs(f[8]).get("ID", "")) for f in _gff_lines(gff)]


def _merged_note(exp: Counter, got: Counter, key_id) -> str:
    """Explain missing rows when features that share an ID (multi-line CDS)
    came out as one row spanning all their segments."""
    by_id = defaultdict(list)
    for k, n in exp.items():
        if key_id(k):
            by_id[key_id(k)].extend([k] * n)
    multi = {i: ks for i, ks in by_id.items() if len(ks) > 1}
    return f"; {len(multi)} input IDs are shared by several lines (e.g. {list(multi)[:2]})" if multi else ""


def compare_gff3_bed6(cand, ref, input_dir=None, **_):
    """gfc BED6 against the input GFF3: one row per feature line, 0-based
    start, ID (or the type when there is no ID) as name, strand."""
    pairs, v = _pairs(cand, ref, [".bed"], [".gff3", ".gff"])
    if v:
        return v
    conv = Conv()
    for c, g in pairs:
        exp = Counter((s, a - 1, b, i or t, st) for s, t, a, b, st, i in _input_features(g))
        got = Counter()
        with _open(c) as fh:
            for line in fh:
                f = line.rstrip("\r\n").split("\t")
                if len(f) >= 6:
                    got[(f[0], int(f[1]), int(f[2]), f[3], f[5])] += 1
        if got != exp:
            miss, extra = exp - got, got - exp
            return Verdict(DISAGREEMENT,
                           f"{c.name}: {sum(miss.values())} of {sum(exp.values())} input features have no "
                           f"matching row (e.g. {list(miss)[:2]}), {sum(extra.values())} rows match no "
                           f"feature (e.g. {list(extra)[:2]})" + _merged_note(exp, got, lambda k: k[3]))
    conv.add("one row per input feature")
    return Verdict(STRUCTURAL, "checked against the input GFF3", conv.used)


def compare_gff3_table(cand, ref, input_dir=None, **_):
    """gfc feature table against the input GFF3: one row per feature line
    with seqid, type, 1-based start and end, strand and ID."""
    pairs, v = _pairs(cand, ref, [".tsv"], [".gff3", ".gff"])
    if v:
        return v
    for c, g in pairs:
        exp = Counter(_input_features(g))
        got = Counter()
        with _open(c) as fh:
            head = fh.readline().rstrip("\r\n").split("\t")
            ix = {h: i for i, h in enumerate(head)}
            need = ("seqid", "type", "start", "end", "strand", "ID")
            if not all(h in ix for h in need):
                return Verdict(DISAGREEMENT, f"{c.name}: header lacks {need}")
            for line in fh:
                f = line.rstrip("\r\n").split("\t")
                got[(f[ix["seqid"]], f[ix["type"]], int(f[ix["start"]]), int(f[ix["end"]]),
                     f[ix["strand"]], f[ix["ID"]])] += 1
        if got != exp:
            miss, extra = exp - got, got - exp
            return Verdict(DISAGREEMENT,
                           f"{c.name}: {sum(miss.values())} of {sum(exp.values())} input features have no "
                           f"matching row (e.g. {list(miss)[:2]}), {sum(extra.values())} rows match no "
                           f"feature (e.g. {list(extra)[:2]})" + _merged_note(exp, got, lambda k: k[5]))
    return Verdict(STRUCTURAL, "checked against the input GFF3", ["one row per input feature"])


def _protein_pair(c: Path, r: Path, conv: Conv, cds_parent=None, **_) -> str:
    def recs(p):
        out = []
        for rid, head, seq in _read_fasta(p):
            k = _norm_id(rid, conv)
            k = (cds_parent or {}).get(k, k)
            s = seq.upper()
            if s.endswith(("*", ".")):
                s = s[:-1]
                conv.add("terminal-stop-symbol")
            if "." in s:
                s = s.replace(".", "*")
                conv.add("stop-symbol")
            out.append((k, head, s))
        return out
    cr, rr = recs(c), recs(r)
    if len({x[0] for x in cr}) != len(cr) or len({x[0] for x in rr}) != len(rr):
        raise Diff(f"{c.name}: duplicate protein ids after normalization")
    _cmp_records(cr, rr, c.name, conv, None)
    return f"{len(rr)} proteins"


def compare_protein(cand, ref, input_dir=None, **_):
    cds_parent = {}
    for g in _files(input_dir, [".gff3", ".gff"]).values() if input_dir else []:
        for f in _gff_lines(g):
            if f[2] == "CDS":
                a = _attrs(f[8])
                if "ID" in a and "Parent" in a:
                    cds_parent[_PREFIX.sub("", a["ID"])] = _PREFIX.sub("", a["Parent"].split(",")[0])
    return _run(_protein_pair, cand, ref, [".faa", ".fa", ".fasta"], cds_parent=cds_parent)


def _genbank_gff3_pair(c: Path, r: Path, conv: Conv, **_) -> str:
    def feats(p):
        cnt = Counter()
        for f in _gff_lines(p):
            cnt[(f[0], f[2], int(f[3]), int(f[4]), f[6])] += 1
        return cnt
    cf, rf = feats(c), feats(r)
    if cf == rf:
        return ""
    def by_type(cnt):
        out = Counter()
        for k, n in cnt.items():
            out[k[1]] += n
        return out
    miss, extra = rf - cf, cf - rf
    raise Diff(f"{c.name}: {sum(miss.values())} reference features missing "
               f"(by type {dict(by_type(miss).most_common(4))}), {sum(extra.values())} extra "
               f"(by type {dict(by_type(extra).most_common(4))})")


def compare_genbank_gff3(cand, ref, **_):
    return _run(_genbank_gff3_pair, cand, ref, [".gff3", ".gff"])


_RECORD_TYPES = {"region", "source", "remark", "chromosome", "databank_entry"}


def _genbank_input_pair(c: Path, g: Path, conv: Conv, **_) -> str:
    """GFF3 against the GenBank input it came from: every feature location
    part, as (seqid, type, start, end, strand), once."""
    from Bio import SeqIO
    parts, spans, lengths, n_joined = Counter(), Counter(), {}, 0
    st = {1: "+", -1: "-"}
    with _open(g) as fh:
        for rec in SeqIO.parse(fh, "genbank"):
            lengths[rec.id] = len(rec.seq)
            for f in rec.features:
                if f.type == "source" or f.location is None:
                    continue
                ps = f.location.parts
                n_joined += len(ps) > 1
                for q in ps:
                    parts[(rec.id, f.type, int(q.start) + 1, int(q.end), st.get(q.strand, "."))] += 1
                spans[(rec.id, f.type, int(f.location.start) + 1, int(f.location.end),
                       st.get(f.location.strand, "."))] += 1
    got = Counter()
    # Accession without its version (NC_001133 for NC_001133.9) is a
    # naming convention, not a content difference.
    bare = {k.rsplit(".", 1)[0]: k for k in lengths if "." in k}
    for f in _gff_lines(c):
        if f[0] not in lengths and f[0] in bare:
            conv.add("accession without version")
            f[0] = bare[f[0]]
        key = (f[0], f[2], int(f[3]), int(f[4]), f[6])
        if f[2] in _RECORD_TYPES and key[2] == 1 and key[3] == lengths.get(f[0]):
            conv.add("record-level line")
            continue
        got[key] += 1
    if got == parts:
        if n_joined:
            conv.add("joined locations as one line per part")
        return f"{sum(parts.values())} feature parts"
    gs, ps = {k[0] for k in got}, {k[0] for k in parts}
    if got and not gs & ps:
        raise Diff(f"{c.name}: sequence names differ from the GenBank records "
                   f"(e.g. {sorted(gs)[:2]} vs {sorted(ps)[:2]})")
    if got == spans and n_joined:
        raise Diff(f"{c.name}: {n_joined} joined (spliced) locations written as one line "
                   "spanning their introns")
    miss, extra = parts - got, got - parts
    def by_type(cnt):
        out = Counter()
        for k, n in cnt.items():
            out[k[1]] += n
        return dict(out.most_common(4))
    raise Diff(f"{c.name}: {sum(miss.values())} input feature parts missing (by type "
               f"{by_type(miss)}), {sum(extra.values())} extra (by type {by_type(extra)})")


def compare_genbank_vs_input(cand, ref, **_):
    return _run(_genbank_input_pair, cand, ref, [".gff3", ".gff"], [".gbk", ".gb", ".gbff"])


# ------------------------------------------------------------- alignments and trees


def _read_alignment(p: Path, fmt: str):
    from Bio import AlignIO
    fmts = [fmt] + (["phylip", "phylip-sequential"] if fmt.startswith("phylip") else [])
    last = None
    for f in fmts:
        try:
            with _open(p) as fh:
                aln = AlignIO.read(fh, f)
            return [(rec.id, rec.id, str(rec.seq)) for rec in aln]
        except Exception as exc:  # try the next flavour
            last = exc
    raise Diff(f"{p.name}: cannot parse as {fmt}: {last}")


def _name_prefix(a: str, b: str) -> bool:
    """One name is a prefix of the other, allowing punctuation to have been
    replaced by '_' (EMBOSS writes 'NAME/1-90' as 'NAME_1-90')."""
    a, b = re.sub(r"[^A-Za-z0-9-]", "_", a), re.sub(r"[^A-Za-z0-9-]", "_", b)
    return a == b[: len(a)] or b == a[: len(b)]


def _alignment_pair(c: Path, r: Path, conv: Conv, cand_fmt="fasta", ref_fmt="fasta", **_) -> str:
    cr, rr = _read_alignment(c, cand_fmt), _read_alignment(r, ref_fmt)
    if len(cr) != len(rr):
        raise Diff(f"{c.name}: {len(cr)} rows against {len(rr)}")
    def norm(s):
        t = s.replace(".", "-")
        if t != s:
            conv.add("gap-character")
        if t != t.upper():
            conv.add("case")
        return t.upper()
    cids, rids = [x[0] for x in cr], [x[0] for x in rr]
    dup = len(set(cids)) < len(cids) or len(set(rids)) < len(rids)
    if all(_name_prefix(x, y) for x, y in zip(cids, rids)):
        # Rows line up by position. Compare by position, so names that
        # strict PHYLIP's 10-character cut made identical cannot hide a row.
        if cids != rids:
            conv.add("name-truncation")
            if not all(x == y[: len(x)] or y == x[: len(y)] for x, y in zip(cids, rids)):
                conv.add("name-character substitution")
        if len(set(cids)) < len(cids):
            conv.add("duplicate names after truncation")
        elif dup:
            conv.add("reference has duplicate names after truncation")
        bad = [i for i, (x, y) in enumerate(zip(cr, rr)) if norm(x[2]) != norm(y[2])]
        if bad:
            i = bad[0]
            raise Diff(f"{c.name}: {len(bad)} of {len(rr)} rows differ (first: row {i + 1}, "
                       f"{cids[i]!r}: {norm(cr[i][2])[:60]!r} vs {norm(rr[i][2])[:60]!r})")
        return ""
    if set(cids) == set(rids) and not dup:
        conv.add("record-order")
        _compare_maps({k: norm(x) for k, _, x in cr}, {k: norm(x) for k, _, x in rr},
                      c.name, conv, label="row")
        return ""
    i = next(i for i, (x, y) in enumerate(zip(cids, rids)) if not _name_prefix(x, y))
    raise Diff(f"{c.name}: names differ (first at row {i + 1}: {cids[i]!r} vs {rids[i]!r})")


def compare_alignment(cand, ref, cand_fmt="fasta", ref_fmt="fasta", **_):
    ext = {"fasta": [".fasta", ".fa"], "stockholm": [".sto", ".stockholm", ".stk"],
           "phylip": [".phylip", ".phy"], "phylip-relaxed": [".phylip", ".phy"]}
    return _run(_alignment_pair, cand, ref, ext[cand_fmt], ext[ref_fmt],
                cand_fmt=cand_fmt, ref_fmt=ref_fmt)


def _splits(tree) -> dict:
    """Each clade's leaf set -> sorted (branch length, support) pairs."""
    out: dict = {}
    for cl in tree.find_clades():
        if cl is tree.root:
            continue
        # Bio.Phylo keeps NEXUS quotes on translated names; drop them.
        leaves = frozenset((t.name or "").strip("'") for t in cl.get_terminals())
        out.setdefault(leaves, []).append((cl.branch_length, cl.confidence))
    return {k: sorted(v, key=lambda x: (x[0] is None, x[0] or 0)) for k, v in out.items()}


def _num_close(a, b, rel=1e-6):
    if a is None or b is None:
        return (a or 0) == (b or 0)
    return abs(float(a) - float(b)) <= rel * max(abs(float(a)), abs(float(b)), 1e-9)


def _tree_diff(a, b, conv: Conv) -> str:
    """Empty when trees ``a`` and ``b`` agree on leaves, clades and branch
    lengths (relative tolerance 1e-6); support values that differ only by
    rounding to two decimals are recorded as a convention."""
    sa, sb = _splits(a), _splits(b)
    la = {t.name.strip("'") for t in a.get_terminals() if t.name}
    lb = {t.name.strip("'") for t in b.get_terminals() if t.name}
    if la != lb:
        return f"leaf sets differ ({sorted(la ^ lb)[:3]})"
    if set(sa) != set(sb):
        return f"topology differs ({len(set(sa) ^ set(sb))} clades not shared)"
    for k in sb:
        for (ba, ca), (bb, cb) in zip(sa[k], sb[k]):
            if not _num_close(ba, bb):
                return f"branch length {ba} vs {bb}"
            if not _num_close(ca, cb):
                if ca is not None and cb is not None and abs(float(ca) - float(cb)) <= 0.005 + 1e-9:
                    conv.add("support rounding")
                else:
                    return f"support value {ca} vs {cb}"
    return ""


def _read_tree(p: Path):
    from Bio import Phylo
    for fmt in ("nexus", "newick", "phyloxml"):
        try:
            with _open(p) as fh:
                trees = list(Phylo.parse(fh, fmt))
            if trees:
                return trees
        except Exception:
            continue
    raise Diff(f"{p.name}: no tree could be parsed")


def _tree_pair(c: Path, r: Path, conv: Conv, **_) -> str:
    ct, rt = _read_tree(c), _read_tree(r)
    if len(ct) != len(rt):
        raise Diff(f"{c.name}: {len(ct)} trees against {len(rt)}")
    for i, (a, b) in enumerate(zip(ct, rt)):
        why = _tree_diff(a, b, conv)
        if why:
            raise Diff(f"{c.name} tree {i + 1}: {why}")
    return ""


def compare_tree(cand, ref, **_):
    return _run(_tree_pair, cand, ref, [".nexus", ".nex", ".tre"])


def compare_tree_vs_input(cand, ref, **_):
    """Converted tree against the input Newick: same leaves, clades and
    branch lengths."""
    return _run(_tree_pair, cand, ref, [".nexus", ".nex", ".tre"], [".nwk"])


def compare_annotated_tree(cand, ref, input_dir=None, **_):
    """gfc NEXUS against the input tree and alignment: same taxa, every
    matrix row equal to the alignment row, same tree."""
    from Bio.Nexus import Nexus
    gap_conv = sec_conv = False
    conv = Conv()
    pairs, v = _pairs(cand, ref, [".nex"], [".nwk"])
    if v:
        return v
    for c, nwk in pairs:
        aln = {rid: s for rid, _, s in _read_fasta(nwk.with_suffix(".fasta"))}
        try:
            nx = Nexus.Nexus(str(c))
        except Exception as exc:
            return Verdict(DISAGREEMENT, f"{c.name}: not readable as NEXUS: {exc}")
        mat = {k: str(s) for k, s in (nx.matrix or {}).items()}
        if set(mat) != set(aln):
            return Verdict(DISAGREEMENT, f"{c.name}: taxa differ from the alignment "
                           f"(e.g. {sorted(set(mat) ^ set(aln))[:3]})")
        def norm(x):
            return x.upper().replace(".", "-")
        if any("." in aln[k] for k in aln):
            gap_conv = True
        bad = [k for k in aln if norm(mat[k]) != norm(aln[k])]
        if bad and all(norm(mat[k]) == norm(aln[k]).replace("U", "X") for k in aln):
            sec_conv, bad = True, []
        if bad:
            return Verdict(DISAGREEMENT, f"{c.name}: {len(bad)} matrix rows differ from the "
                           f"alignment (first {bad[0]})")
        from Bio import Phylo
        try:
            out_tree = _read_tree(c)[0]
        except Diff as exc:
            return Verdict(DISAGREEMENT, str(exc))
        in_tree = Phylo.read(str(nwk), "newick")
        why = _tree_diff(out_tree, in_tree, conv)
        if why:
            return Verdict(DISAGREEMENT, f"{c.name}: tree differs from the input tree: {why}")
    return Verdict(STRUCTURAL, "checked against the input tree and alignment",
                   ["taxon translate table"] + (["gap-character"] if gap_conv else [])
                   + (["selenocysteine as X"] if sec_conv else []) + conv.used)


# ------------------------------------------------------------- comparative genomics


def compare_blast_links(cand, ref, input_dir=None, **_):
    """One link row per BLAST hit, with the same names, coordinates and identity."""
    pairs, v = _pairs(cand, ref, [".links.tsv"], [".tab", ".tsv", ".tabular"])
    if v:
        return v
    conv = Conv()
    for c, b in pairs:
        exp = Counter()
        with _open(b) as fh:
            for line in fh:
                f = line.rstrip("\r\n").split("\t")
                if len(f) >= 12 and not line.startswith("#"):
                    exp[(f[0], int(f[6]), int(f[7]), f[1], int(f[8]), int(f[9]), float(f[2]))] += 1
        got = Counter()
        with _open(c) as fh:
            first = fh.readline()
            if not first.split("\t")[1].strip().lstrip("-").isdigit():
                conv.add("header-row")
            else:
                fh.seek(0)
            for line in fh:
                f = line.rstrip("\r\n").split("\t")
                got[(f[0], int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5]), float(f[6]))] += 1
        if got != exp:
            miss, extra = exp - got, got - exp
            return Verdict(DISAGREEMENT, f"{c.name}: {sum(miss.values())} hits without a link "
                           f"(e.g. {list(miss)[:2]}), {sum(extra.values())} links without a hit "
                           f"(e.g. {list(extra)[:2]})")
    return Verdict(STRUCTURAL, "checked against the input BLAST table", conv.used or ["header-row"])


def _coords(p: Path, cand: bool):
    """Alignments as (ref, rstart, rend, qry, qstart, qend, strand)."""
    out = Counter()
    with _open(p) as fh:
        lines = [l.rstrip("\r\n").split("\t") for l in fh if l.strip()]
    if cand:
        head = lines[0]
        ix = {h.strip().lower(): i for i, h in enumerate(head)}
        need = ("ref_name", "ref_start", "ref_end", "query_name", "query_start", "query_end")
        if not all(k in ix for k in need):
            raise Diff(f"{p.name}: header lacks {need}")
        has_strand = "strand" in ix
        for f in lines[1:]:
            qs, qe = int(f[ix["query_start"]]), int(f[ix["query_end"]])
            st = f[ix["strand"]] if has_strand else None
            out[(f[ix["ref_name"]], int(f[ix["ref_start"]]), int(f[ix["ref_end"]]),
                 f[ix["query_name"]], min(qs, qe), max(qs, qe), st)] += 1
        return out, has_strand
    for f in lines:  # show-coords -T -H: S1 E1 S2 E2 LEN1 LEN2 %IDY REF QRY
        s1, e1, s2, e2 = (int(x) for x in f[:4])
        st = "+" if (s2 <= e2) == (s1 <= e1) else "-"
        out[(f[7], min(s1, e1), max(s1, e1), f[8], min(s2, e2), max(s2, e2), st)] += 1
    return out, True


def _delta_pair(c: Path, r: Path, conv: Conv, **_) -> str:
    cc, has_strand = _coords(c, True)
    rc, _ = _coords(r, False)
    if not has_strand:
        n_rev = sum(n for k, n in rc.items() if k[6] == "-")
        rc = Counter({k[:6] + (None,): n for k, n in rc.items()})
        if cc == rc and n_rev:
            raise Diff(f"{c.name}: coordinates match, but the orientation of {n_rev} reverse "
                       "alignments is not recorded")
    if cc != rc:
        miss, extra = rc - cc, cc - rc
        raise Diff(f"{c.name}: {sum(miss.values())} of {sum(rc.values())} alignments only in "
                   f"show-coords (e.g. {list(miss)[:1]}), {sum(extra.values())} only in candidate "
                   f"(e.g. {list(extra)[:1]})")
    conv.add("reverse alignments as a strand column")
    return f"{sum(rc.values())} alignments"


def compare_delta(cand, ref, **_):
    return _run(_delta_pair, cand, ref, [".tsv"], [".coords"])


_XMFA_HEAD = re.compile(r"^>\s*(\d+):(\d+)-(\d+)\s+([+-])(?:\s+(.*))?$")


def compare_xmfa(cand, ref, input_dir=None, **_):
    """gfc XMFA against the input MAF: one block per MAF block, Mauve headers
    ("> index:start-end strand name"), forward-strand 1-based coordinates
    and the same aligned sequences."""
    from Bio import AlignIO
    pairs, v = _pairs(cand, ref, [".xmfa"], [".maf"])
    if v:
        return v
    for c, m in pairs:
        exp = []
        with _open(m) as fh:
            for aln in AlignIO.parse(fh, "maf"):
                blk = []
                for rec in aln:
                    a = rec.annotations
                    size, start, src = a["size"], a["start"], a["srcSize"]
                    if a.get("strand") in (-1, "-"):
                        lo, hi, st = src - start - size + 1, src - start, "-"
                    else:
                        lo, hi, st = start + 1, start + size, "+"
                    blk.append((rec.id, lo, hi, st, str(rec.seq).upper()))
                exp.append(blk)
        got, blk, names = [], [], {}
        head, buf = None, []
        with _open(c) as fh:
            for line in fh:
                line = line.rstrip("\r\n")
                if line.startswith("#") or not line:
                    continue
                if line.startswith(">") or line.startswith("="):
                    if head is not None:
                        blk.append(head + ("".join(buf).upper(),))
                        head, buf = None, []
                    if line.startswith("="):
                        got.append(blk)
                        blk = []
                        continue
                    mh = _XMFA_HEAD.match(line)
                    if not mh:
                        return Verdict(DISAGREEMENT, f"{c.name}: header {line[:60]!r} is not "
                                       "'> index:start-end strand name' (XMFA)")
                    idx, lo, hi, st, nm = mh.groups()
                    nm = (nm or "").split()[0] if nm else ""
                    if names.setdefault(idx, nm) != nm:
                        return Verdict(DISAGREEMENT, f"{c.name}: sequence index {idx} names "
                                       f"both {names[idx]!r} and {nm!r}")
                    head = (nm or idx, int(lo), int(hi), st)
                else:
                    buf.append(line.strip())
        if len(got) != len(exp):
            return Verdict(DISAGREEMENT, f"{c.name}: {len(got)} blocks against {len(exp)} MAF blocks")
        for i, (g, e) in enumerate(zip(got, exp)):
            if g != e:
                return Verdict(DISAGREEMENT, f"{c.name}: block {i + 1} differs "
                               f"(XMFA {[x[:4] for x in g][:2]} vs MAF {[x[:4] for x in e][:2]})")
    return Verdict(STRUCTURAL, "checked against the input MAF", ["XMFA rendering"])


# ------------------------------------------------------------- variants


def _vcf_bed_pair(c: Path, r: Path, conv: Conv, **_) -> str:
    def rows(p):
        out, ncol = Counter(), set()
        with _open(p) as fh:
            for line in fh:
                f = line.rstrip("\r\n").split("\t")
                if len(f) >= 3 and not line.startswith(("#", "track")):
                    try:
                        out[(f[0], int(f[1]), int(f[2]))] += 1
                    except ValueError:
                        raise Diff(f"{p.name}: non-integer coordinates in {line[:60]!r}")
                    ncol.add(len(f))
        return out, ncol
    (cr, cn), (rr, _) = rows(c), rows(r)
    if cr != rr:
        miss, extra = rr - cr, cr - rr
        raise Diff(f"{c.name}: {sum(miss.values())} intervals only in reference (e.g. {list(miss)[:2]}), "
                   f"{sum(extra.values())} only in candidate (e.g. {list(extra)[:2]})")
    if cn == {3}:
        conv.add("BED3 (no name column)")
    else:
        conv.add("name column (ID_ALT against ID)")
    return f"{sum(rr.values())} intervals"


def compare_vcf_bed(cand, ref, **_):
    return _run(_vcf_bed_pair, cand, ref, [".bed"])


def _vcf_table_pair(c: Path, r: Path, conv: Conv, **_) -> str:
    def rows(p):
        with _open(p) as fh:
            lines = [l.rstrip("\r\n").split("\t") for l in fh if l.strip()]
        head, body = lines[0], lines[1:]
        names = [re.sub(r"^#?\s*\[\d+\]", "", h).strip().lstrip("#") for h in head]
        return names, body
    cn, cb = rows(c)
    rn, rb = rows(r)
    if [x.upper() for x in cn[: len(rn)]] != [x.upper() for x in rn]:
        raise Diff(f"{c.name}: columns {cn} against {rn}")
    if cn != rn:
        conv.add("header labels")
    if len(cb) != len(rb):
        raise Diff(f"{c.name}: {len(cb)} rows against {len(rb)}")
    for i, (a, b) in enumerate(zip(cb, rb)):
        a = a[: len(b)]
        if a == b:
            continue
        for j, (x, y) in enumerate(zip(a, b)):
            if x == y:
                continue
            try:
                if float(x) == float(y):
                    conv.add("number format")
                    continue
            except ValueError:
                pass
            raise Diff(f"{c.name}: row {i + 1} column {rn[j]}: {x!r} vs {y!r}")
    return f"{len(rb)} rows"


def compare_vcf_table(cand, ref, **_):
    return _run(_vcf_table_pair, cand, ref, [".tsv"])


def _consensus_pair(c: Path, r: Path, conv: Conv, **_) -> str:
    def recs(p):
        out = []
        for rid, head, seq in _read_fasta(p):
            out.append((rid, head, seq))
        return out
    _cmp_records(recs(c), recs(r), c.name, conv, _seq_eq)
    return ""


def compare_consensus(cand, ref, **_):
    """gfc writes <stem>_consensus.fa; the bcftools loop writes <stem>.fa."""
    c = _files(cand, ["_consensus.fa", ".fa", ".fasta"])
    r = _files(ref, [".fa"])
    if not r:
        return Verdict(LIMITATION, "no reference consensus")
    if set(c) != set(r):
        return Verdict(DISAGREEMENT, f"file sets differ: {sorted(c)} vs {sorted(r)}")
    conv = Conv()
    try:
        for k in r:
            _consensus_pair(c[k], r[k], conv)
    except Diff as exc:
        return Verdict(DISAGREEMENT, str(exc))
    return _finish(conv, all(_same_bytes([c[k]], [r[k]]) for k in r))


def _phylip_rows(p: Path) -> Dict[str, str]:
    with _open(p) as fh:
        lines = [l.rstrip("\r\n") for l in fh if l.strip()]
    n, m = (int(x) for x in lines[0].split()[:2])
    rows = {}
    order = []
    body = lines[1:]
    for line in body[:n]:
        name, seq = line.split(None, 1) if " " in line or "\t" in line else (line[:10], line[10:])
        rows[name] = seq.replace(" ", "")
        order.append(name)
    # interleaved continuation blocks
    for i, line in enumerate(body[n:]):
        rows[order[i % n]] += line.replace(" ", "")
    for k, s in rows.items():
        if len(s) != m:
            raise Diff(f"{p.name}: row {k} has {len(s)} of {m} columns")
    return rows


def _snp_aln_pair(c: Path, r: Path, conv: Conv, **_) -> str:
    cr, rr = _phylip_rows(c), _phylip_rows(r)
    if set(cr) != set(rr):
        strip = {re.sub(r"^(.+)_\1$", r"\1", k): v for k, v in rr.items()}
        if set(strip) == set(cr):
            conv.add("sample-id (FID_IID)")
            rr = strip
    _compare_maps(cr, rr, c.name, conv, eq=lambda a, b, conv: _seq_eq(a, b, conv), label="sample")
    return f"{len(rr)} samples"


def compare_snp_alignment(cand, ref, **_):
    return _run(_snp_aln_pair, cand, ref, [".phylip", ".phy"])


def _vcf_calls(p: Path):
    """(chrom, pos) -> (allele set, {sample: sorted base genotype})."""
    out, samples = {}, []
    with _open(p) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\r\n").split("\t")
            if line.startswith("#"):
                samples = f[9:]
                continue
            alleles = [f[3]] + f[4].split(",")
            fmt = f[8].split(":")
            gi = fmt.index("GT")
            gts = []
            for s in f[9:]:
                g = s.split(":")[gi]
                gts.append(tuple(sorted(alleles[int(x)] if x not in (".", "") else "."
                                        for x in re.split(r"[/|]", g))))
            out[(f[0], int(f[1]))] = (f[3], frozenset(alleles), tuple(gts))
    return samples, out


def _vcf_gt_pair(c: Path, r: Path, conv: Conv, **_) -> str:
    cs, cc = _vcf_calls(c)
    rs, rc = _vcf_calls(r)
    if cs != rs:
        strip = [re.sub(r"^(.+)_\1$", r"\1", x) for x in rs]
        if strip == cs:
            conv.add("sample-id (FID_IID)")
        else:
            raise Diff(f"{c.name}: samples differ (e.g. {cs[:2]} vs {rs[:2]})")
    def eq(a, b, conv):
        if a[1] != b[1] or a[2] != b[2]:
            return False
        if a[0] != b[0]:
            conv.add("REF allele choice")
        return True
    _compare_maps(cc, rc, c.name, conv, eq, label="site")
    return f"{len(rc)} sites"


def compare_vcf_genotypes(cand, ref, **_):
    return _run(_vcf_gt_pair, cand, ref, [".vcf", ".vcf.gz"])


COMPARATORS = {
    "fasta": compare_fasta, "fastq": compare_fastq, "fasta_qual": compare_fasta_qual,
    "seq_table": compare_seq_table, "bed6": compare_bed6,
    "gff3_models": compare_gff3_models, "gff3_bed6": compare_gff3_bed6,
    "gff3_table": compare_gff3_table, "protein": compare_protein,
    "genbank_gff3": compare_genbank_gff3, "genbank_vs_input": compare_genbank_vs_input, "alignment": compare_alignment,
    "tree": compare_tree, "tree_vs_input": compare_tree_vs_input, "annotated_tree": compare_annotated_tree,
    "blast_links": compare_blast_links, "delta": compare_delta, "xmfa": compare_xmfa,
    "vcf_bed": compare_vcf_bed, "vcf_table": compare_vcf_table,
    "consensus": compare_consensus, "snp_alignment": compare_snp_alignment,
    "vcf_genotypes": compare_vcf_genotypes,
}
