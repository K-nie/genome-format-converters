#!/usr/bin/env python3
"""Mutation test of the four-state comparators on real benchmark data.

``mutation_test.py`` plants defects in gfc's output on small synthetic
inputs. This script does the same on the real inputs of the benchmark, so
the comparators are tested on files with real-world structure (thousands of
samples, multi-exon transcripts, transposable elements, many records):

1. run gfc on a real input;
2. compare the output with an exact copy (expects ``byte``);
3. compare it with copies rewritten to another valid convention
   (expects ``structural``);
4. plant one defect at a randomly chosen site, sample, record or feature,
   ``--reps`` times per defect type with a seeded RNG, and compare with the
   original (expects ``disagreement``).

Every case records where the defect was planted, so a reader can check it.
Defects that cannot be planted in a given input (for example, a missing
call in a VCF that has none) are reported as ``not applicable``, not as
passes.

Inputs (all optional; each enables the tasks it feeds):

  --vcf FILE             T1, T2 and T8 (the first --sites biallelic SNPs)
  --annotation-dir DIR   T4 and T5 (every *.gff3 in DIR); repeatable
  --genome-dir DIR       T3 (paired <stem>.fasta/.fa + <stem>.gff3 in DIR)

Example (cluster)::

    python benchmarks/mutation_test_real.py \\
        --vcf data/1kg/chr22.biallelic.vcf.gz --sites 20000 \\
        --annotation-dir data/annotations/hs_chr22 \\
        --annotation-dir data/annotations/sc_r64 \\
        --annotation-dir data/y1000plus \\
        --genome-dir data/annotations/sc_r64_genome \\
        --reps 5 --seed 1 --out results/mutation_real.tsv
"""
from __future__ import annotations

import argparse
import gzip
import random
import re
import shutil
import sys
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, List, Optional

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import four_state as fs  # noqa: E402
from mutation_test import (  # noqa: E402
    copy_dir, eig_allele_swap, gfc, plink_allele_swap)

NA = "not applicable"


@dataclass
class Case:
    dataset: str
    task: str
    kind: str
    name: str
    expected: str
    got: str = ""
    where: str = ""
    reason: str = ""
    seconds: float = 0.0

    @property
    def ok(self) -> bool:
        return self.got == self.expected


class NotApplicable(Exception):
    pass


# ------------------------------------------------------------- helpers


def _read(f: Path) -> List[str]:
    return f.read_text().splitlines()


def _write(f: Path, lines: List[str]) -> None:
    f.write_text("\n".join(lines) + "\n")


def _one(d: Path, pat: str, rng: random.Random) -> Path:
    files = sorted(d.glob(pat))
    if not files:
        raise NotApplicable(f"no {pat} file")
    return rng.choice(files)


def _toggle_chr(c: str) -> str:
    return c[3:] if c.lower().startswith("chr") else "chr" + c


# ------------------------------------------------------------- EIGENSTRAT (T1)


def eig_chrom_label(d, rng):
    f = next(d.glob("*.snp"))
    out = []
    for l in _read(f):
        p = l.split()
        p[1] = _toggle_chr(p[1])
        out.append("\t".join(p))
    _write(f, out)
    return "every chromosome label"


def eig_allele_order(d, rng):
    eig_allele_swap(d)
    return "every site"


def eig_flip_one_call(d, rng):
    geno = next(d.glob("*.geno"))
    g = _read(geno)
    for _ in range(1000):
        i = rng.randrange(len(g))
        cols = [k for k, c in enumerate(g[i]) if c in "02"]
        if cols:
            j = rng.choice(cols)
            g[i] = g[i][:j] + ("2" if g[i][j] == "0" else "0") + g[i][j + 1:]
            _write(geno, g)
            return f"site {i + 1}, sample {j + 1}"
    raise NotApplicable("no homozygous call found")


def eig_pos_off_by_one(d, rng):
    snp = next(d.glob("*.snp"))
    s = _read(snp)
    i = rng.randrange(len(s))
    p = s[i].split()
    p[3] = str(int(p[3]) + 1)
    s[i] = "\t".join(p)
    _write(snp, s)
    return f"site {i + 1}"


def eig_swap_alleles_only(d, rng):
    snp = next(d.glob("*.snp"))
    s = _read(snp)
    i = rng.randrange(len(s))
    p = s[i].split()
    p[4], p[5] = p[5], p[4]
    s[i] = "\t".join(p)
    _write(snp, s)
    return f"site {i + 1}"


def eig_drop_site(d, rng):
    n = len(_read(next(d.glob("*.snp"))))
    i = rng.randrange(n)
    for pat in ("*.snp", "*.geno"):
        f = next(d.glob(pat))
        l = _read(f)
        del l[i]
        _write(f, l)
    return f"site {i + 1}"


def eig_swap_samples(d, rng):
    ind = next(d.glob("*.ind"))
    geno = _read(next(d.glob("*.geno")))
    l = _read(ind)
    for _ in range(1000):
        a, b = rng.sample(range(len(l)), 2)
        if any(row[a] != row[b] for row in geno):
            l[a], l[b] = l[b], l[a]
            _write(ind, l)
            return f"samples {a + 1} and {b + 1}"
    raise NotApplicable("no two samples with different genotypes")


# ------------------------------------------------------------- PLINK (T2)


def _plink(d):
    bed, bim, fam = next(d.glob("*.bed")), next(d.glob("*.bim")), next(d.glob("*.fam"))
    n = len(_read(fam))
    return bed, bim, n, (n + 3) // 4


def plink_chrom_label(d, rng):
    f = next(d.glob("*.bim"))
    out = []
    for l in _read(f):
        p = l.split("\t")
        p[0] = _toggle_chr(p[0])
        out.append("\t".join(p))
    _write(f, out)
    return "every chromosome label"


def plink_allele_order(d, rng):
    plink_allele_swap(d)
    return "every variant"


def plink_flip_genotype(d, rng):
    bed, bim, n, per = _plink(d)
    n_var = len(_read(bim))
    raw = bytearray(bed.read_bytes())
    for _ in range(10000):
        v, s = rng.randrange(n_var), rng.randrange(n)
        i, shift = 3 + v * per + s // 4, 2 * (s % 4)
        c = (raw[i] >> shift) & 3
        if c in (0, 3):
            raw[i] ^= (3 << shift)
            bed.write_bytes(bytes(raw))
            return f"variant {v + 1}, sample {s + 1}"
    raise NotApplicable("no homozygous call found")


def plink_pos_shift(d, rng):
    bim = next(d.glob("*.bim"))
    l = _read(bim)
    i = rng.randrange(len(l))
    p = l[i].split("\t")
    p[3] = str(int(p[3]) + 1)
    l[i] = "\t".join(p)
    _write(bim, l)
    return f"variant {i + 1}"


def plink_swap_bim_only(d, rng):
    bim = next(d.glob("*.bim"))
    l = _read(bim)
    i = rng.randrange(len(l))
    p = l[i].split("\t")
    p[4], p[5] = p[5], p[4]
    l[i] = "\t".join(p)
    _write(bim, l)
    return f"variant {i + 1}"


# ------------------------------------------------------------- pseudohaploid (T8)


class VcfClasses:
    """Genotype class per (site, sample) of the input VCF, loaded once."""

    def __init__(self, vcf: Path):
        self.sites = {}
        for c, pos, r, a, gt in fs._iter_vcf_sites(vcf):
            self.sites[(fs._norm_chrom(c), int(pos), frozenset((r, a)))] = gt.copy()  # cyvcf2 may reuse the buffer

    def get(self, p):
        # Key on alleles too: a VCF can hold two records at one position.
        return self.sites[(fs._norm_chrom(p[1]), int(p[3]), frozenset((p[4].upper(), p[5].upper())))]

    def rows(self, snp_lines):
        for i, s in enumerate(snp_lines):
            p = s.split()
            yield i, self.get(p)


def _ph_set(d, cls: VcfClasses, rng, want: int, new: Callable[[str], str], label):
    snp = _read(next(d.glob("*.snp")))
    geno_f = next(d.glob("*.geno"))
    g = _read(geno_f)
    order = list(range(len(snp)))
    rng.shuffle(order)
    for i in order:
        p = snp[i].split()
        gt = cls.get(p)
        cols = (gt == want).nonzero()[0].tolist()
        if cols:
            j = rng.choice(cols)
            g[i] = g[i][:j] + new(g[i][j]) + g[i][j + 1:]
            _write(geno_f, g)
            return f"site {i + 1}, sample {j + 1}"
    raise NotApplicable(f"input has no {label} call")


def ph_homref_flip(d, cls, rng):
    return _ph_set(d, cls, rng, 0, lambda c: "0" if c == "2" else "2", "hom-ref")


def ph_homalt_flip(d, cls, rng):
    return _ph_set(d, cls, rng, 3, lambda c: "2" if c == "0" else "0", "hom-alt")


def ph_missing_filled(d, cls, rng):
    return _ph_set(d, cls, rng, 2, lambda c: "2", "missing")


def ph_het_as_one(d, cls, rng):
    return _ph_set(d, cls, rng, 1, lambda c: "1", "heterozygous")


def ph_biased_sampler(d, cls, rng):
    snp = _read(next(d.glob("*.snp")))
    geno_f = next(d.glob("*.geno"))
    g = _read(geno_f)
    for i, gt in cls.rows(snp):
        row = list(g[i])
        for k in (gt == 1).nonzero()[0].tolist():
            row[k] = "2"
        g[i] = "".join(row)
    _write(geno_f, g)
    return "every heterozygous call drawn as REF"


# ------------------------------------------------------------- GTF (T4)


def _gtf_rows(f):
    return [l.split("\t") for l in _read(f)]


def _tid(p):
    m = re.search(r'transcript_id "([^"]+)"', p[8]) if len(p) >= 9 else None
    return m.group(1) if m else None


def gtf_gffread_style(d, rng):
    for f in d.glob("*.gtf"):
        out = []
        for p in _gtf_rows(f):
            if len(p) < 9:
                out.append("\t".join(p))
                continue
            if p[2] == "gene":
                continue
            if p[2] == "mRNA":
                p[2] = "transcript"
            p[8] = p[8].rstrip().rstrip(";")
            out.append("\t".join(p))
        _write(f, out)
    return "all files"


def gtf_split_exon(d, rng):
    f = _one(d, "*.gtf", rng)
    rows = _gtf_rows(f)
    cand = [i for i, p in enumerate(rows)
            if len(p) >= 9 and p[2] == "exon" and int(p[4]) - int(p[3]) > 10]
    if not cand:
        raise NotApplicable("no exon longer than 10 bp")
    i = rng.choice(cand)
    p = rows[i]
    mid = (int(p[3]) + int(p[4])) // 2
    rows[i:i + 1] = [p[:4] + [str(mid)] + p[5:], p[:3] + [str(mid + 1)] + p[4:]]
    _write(f, ["\t".join(p) for p in rows])
    return f"{f.name}: exon of {_tid(p)}"


def _gtf_pick(d, rng, pred, edit, what):
    f = _one(d, "*.gtf", rng)
    rows = _gtf_rows(f)
    cand = [i for i, p in enumerate(rows) if len(p) >= 9 and pred(p)]
    if not cand:
        raise NotApplicable(f"no {what} in {f.name}")
    i = rng.choice(cand)
    tid = _tid(rows[i])
    rows[i] = edit(rows[i])
    _write(f, ["\t".join(p) for p in rows])
    return f"{f.name}: {what} of {tid}"


def gtf_exon_shift(d, rng):
    return _gtf_pick(d, rng, lambda p: p[2] == "exon",
                     lambda p: p[:3] + [str(int(p[3]) + 1)] + p[4:], "exon")


def gtf_cds_trim(d, rng):
    return _gtf_pick(d, rng, lambda p: p[2] == "CDS" and int(p[4]) - int(p[3]) > 3,
                     lambda p: p[:4] + [str(int(p[4]) - 3)] + p[5:], "CDS")


def _gtf_transcripts(f):
    tids = []
    for p in _gtf_rows(f):
        t = _tid(p)
        if t and t not in tids:
            tids.append(t)
    return tids


def gtf_strand_flip(d, rng):
    f = _one(d, "*.gtf", rng)
    tids = _gtf_transcripts(f)
    if not tids:
        raise NotApplicable(f"no transcript in {f.name}")
    t = rng.choice(tids)
    rows = _gtf_rows(f)
    for p in rows:
        if _tid(p) == t:
            p[6] = "-" if p[6] == "+" else "+"
    _write(f, ["\t".join(p) for p in rows])
    return f"{f.name}: transcript {t}"


def gtf_drop_transcript(d, rng):
    f = _one(d, "*.gtf", rng)
    tids = _gtf_transcripts(f)
    if not tids:
        raise NotApplicable(f"no transcript in {f.name}")
    t = rng.choice(tids)
    _write(f, ["\t".join(p) for p in _gtf_rows(f) if _tid(p) != t])
    return f"{f.name}: transcript {t}"


# ------------------------------------------------------------- BED12 (T5)


def bed_score_and_order(d, rng):
    for f in d.glob("*.bed12"):
        rows = [l.split("\t") for l in _read(f) if l]
        for r in rows:
            r[4] = "1000"
            r[8] = "255,0,0"
        rows.reverse()
        _write(f, ["\t".join(r) for r in rows])
    return "all files"


def _bed_pick(d, rng, pred, edit, what):
    f = _one(d, "*.bed12", rng)
    rows = [l.split("\t") for l in _read(f) if l]
    cand = [i for i, r in enumerate(rows) if pred(r)]
    if not cand:
        raise NotApplicable(f"no {what} row in {f.name}")
    i = rng.choice(cand)
    edit(rows[i])
    _write(f, ["\t".join(r) for r in rows])
    return f"{f.name}: {rows[i][3]}"


def bed_start_off_by_one(d, rng):
    return _bed_pick(d, rng, lambda r: True,
                     lambda r: r.__setitem__(1, str(int(r[1]) + 1)), "any")


def bed_block_size(d, rng):
    def edit(r):
        s = r[10].rstrip(",").split(",")
        k = max(range(len(s)), key=lambda x: int(s[x]))
        s[k] = str(int(s[k]) - 1)
        r[10] = ",".join(s) + ","
    return _bed_pick(d, rng, lambda r: max(map(int, r[10].rstrip(",").split(","))) > 1,
                     edit, "block")


def bed_strand(d, rng):
    return _bed_pick(d, rng, lambda r: r[5] in "+-",
                     lambda r: r.__setitem__(5, "-" if r[5] == "+" else "+"), "stranded")


def bed_thick(d, rng):
    return _bed_pick(d, rng, lambda r: int(r[7]) - int(r[6]) > 3,
                     lambda r: r.__setitem__(6, str(int(r[6]) + 3)), "coding")


def bed_rename_one(d, rng):
    return _bed_pick(d, rng, lambda r: True,
                     lambda r: r.__setitem__(3, r[3] + "_renamed"), "any")


# ------------------------------------------------------------- GenBank (T3)

_FEAT = re.compile(r"^     (\S+)\s+((?:complement\(|join\(|order\()*)<?(\d+)")


def gbk_point_mutation(d, rng):
    f = _one(d, "*.gbk", rng)
    lines = _read(f)
    seq_rows = []
    in_seq = False
    for i, l in enumerate(lines):
        if l.startswith("ORIGIN"):
            in_seq = True
        elif l.startswith("//"):
            in_seq = False
        elif in_seq and re.match(r"^\s*\d+ [a-z]", l):
            seq_rows.append(i)
    if not seq_rows:
        raise NotApplicable(f"no sequence in {f.name}")
    i = rng.choice(seq_rows)
    l = lines[i]
    pos = [k for k, c in enumerate(l) if c in "acgt"]
    if not pos:
        raise NotApplicable("no a/c/g/t base")
    k = rng.choice(pos)
    lines[i] = l[:k] + ("a" if l[k] != "a" else "c") + l[k + 1:]
    _write(f, lines)
    return f"{f.name}: line {i + 1}, base {l.split()[0]}+"


def gbk_shift_feature(d, rng):
    f = _one(d, "*.gbk", rng)
    lines = _read(f)
    cand = [i for i, l in enumerate(lines) if _FEAT.match(l) and not l.startswith("     source")]
    if not cand:
        raise NotApplicable(f"no feature in {f.name}")
    i = rng.choice(cand)
    m = _FEAT.match(lines[i])
    a, b = m.start(3), m.end(3)
    lines[i] = lines[i][:a] + str(int(m.group(3)) + 1) + lines[i][b:]
    _write(f, lines)
    return f"{f.name}: {m.group(1)} at {m.group(3)}"


def gbk_drop_feature(d, rng):
    f = _one(d, "*.gbk", rng)
    lines = _read(f)
    cand = [i for i, l in enumerate(lines) if _FEAT.match(l) and not l.startswith("     source")]
    if not cand:
        raise NotApplicable(f"no feature in {f.name}")
    i = rng.choice(cand)
    m = _FEAT.match(lines[i])
    j = i + 1
    while j < len(lines) and lines[j].startswith(" " * 21):
        j += 1
    del lines[i:j]
    _write(f, lines)
    return f"{f.name}: {m.group(1)} at {m.group(3)}"


# ------------------------------------------------------------- driver


@dataclass
class Runner:
    work: Path
    reps: int
    rng: random.Random
    cases: List[Case] = field(default_factory=list)
    n: int = 0

    def _dir(self, tag):
        self.n += 1
        return self.work / f"c{self.n:04d}_{tag}"

    def check(self, ds, task, kind, name, expected, cand, ref, where="", **ctx):
        t0 = time.perf_counter()
        v = fs.COMPARATORS[task](cand, ref, **ctx)
        c = Case(ds, task, kind, name, expected, v.state, where, v.reason,
                 time.perf_counter() - t0)
        self.cases.append(c)
        print(f"[{'ok' if c.ok else 'UNEXPECTED'}] {ds} {task} {kind} {name} "
              f"({where}): {c.got}", file=sys.stderr, flush=True)

    def na(self, ds, task, kind, name, expected, why):
        self.cases.append(Case(ds, task, kind, name, expected, NA, "", why))
        print(f"[n/a] {ds} {task} {name}: {why}", file=sys.stderr, flush=True)

    def suite(self, ds, task, base: Path, conventions, mutations, **ctx):
        self.check(ds, task, "identity", "exact copy", fs.BYTE, base,
                   copy_dir(base, self._dir(f"{task}_id")), **ctx)
        for fn in conventions:
            ref = copy_dir(base, self._dir(f"{task}_{fn.__name__}"))
            try:
                where = fn(ref, self.rng)
            except NotApplicable as e:
                self.na(ds, task, "convention", fn.__name__, fs.STRUCTURAL, str(e))
                continue
            self.check(ds, task, "convention", fn.__name__, fs.STRUCTURAL, base, ref,
                       where, **ctx)
            shutil.rmtree(ref, ignore_errors=True)
        for fn in mutations:
            for _ in range(self.reps):
                cand = copy_dir(base, self._dir(f"{task}_{fn.__name__}"))
                try:
                    where = fn(cand, self.rng)
                except NotApplicable as e:
                    self.na(ds, task, "mutation", fn.__name__, fs.DISAGREEMENT, str(e))
                    shutil.rmtree(cand, ignore_errors=True)
                    break
                self.check(ds, task, "mutation", fn.__name__, fs.DISAGREEMENT, cand, base,
                           where, **ctx)
                shutil.rmtree(cand, ignore_errors=True)


def slice_vcf(src: Path, dst: Path, n_sites: int) -> int:
    """First n biallelic SNP records of src (header kept), plain text."""
    op = gzip.open if open(src, "rb").read(2) == b"\x1f\x8b" else open
    kept = 0
    with op(src, "rt") as fi, open(dst, "w") as fo:
        for line in fi:
            if line.startswith("#"):
                fo.write(line)
                continue
            p = line.split("\t", 6)
            if len(p[3]) == 1 and len(p[4]) == 1 and p[4] in "ACGT":
                fo.write(line)
                kept += 1
                if kept >= n_sites:
                    break
    return kept


def run(args, work: Path) -> List[Case]:
    R = Runner(work, args.reps, random.Random(args.seed))
    inp = work / "inputs"
    inp.mkdir(parents=True)

    if args.vcf:
        vdir = inp / "vcf"
        vdir.mkdir()
        vcf = vdir / "slice.vcf"
        n = slice_vcf(Path(args.vcf), vcf, args.sites)
        ds = f"{Path(args.vcf).name} ({n} SNPs)"
        print(f"[vcf] {ds}", file=sys.stderr)
        gfc("vcf-to-eigenstrat", "--input-dir", str(vdir), "--output-dir", str(work / "t1"), "--force")
        R.suite(ds, "T1", work / "t1", [eig_chrom_label, eig_allele_order],
                [eig_flip_one_call, eig_pos_off_by_one, eig_swap_alleles_only,
                 eig_drop_site, eig_swap_samples])
        gfc("vcf-to-plink", "--input-dir", str(vdir), "--output-dir", str(work / "t2"), "--force")
        R.suite(ds, "T2", work / "t2",
                [plink_allele_order, plink_chrom_label],
                [plink_flip_genotype, plink_pos_shift, plink_swap_bim_only])
        gfc("vcf-to-pseudohaploid", "--input-dir", str(vdir), "--output-dir", str(work / "t8"),
            "--seed", "42", "--force")
        cls = VcfClasses(vcf)
        t8 = work / "t8"
        R.check(ds, "T8", "identity", "gfc output vs input", fs.STRUCTURAL, t8, None,
                input_vcf=vcf)
        for fn in (ph_homref_flip, ph_homalt_flip, ph_missing_filled, ph_het_as_one,
                   ph_biased_sampler):
            reps = 1 if fn is ph_biased_sampler else args.reps
            for _ in range(reps):
                cand = copy_dir(t8, R._dir(f"T8_{fn.__name__}"))
                try:
                    where = fn(cand, cls, R.rng)
                except NotApplicable as e:
                    R.na(ds, "T8", "mutation", fn.__name__, fs.DISAGREEMENT, str(e))
                    shutil.rmtree(cand, ignore_errors=True)
                    break
                R.check(ds, "T8", "mutation", fn.__name__, fs.DISAGREEMENT, cand, None,
                        where, input_vcf=vcf)
                shutil.rmtree(cand, ignore_errors=True)

    for k, adir in enumerate(args.annotation_dir or []):
        adir = Path(adir)
        ds = f"{adir.name} ({len(list(adir.glob('*.gff3')))} GFF3)"
        t4, t5 = work / f"t4_{k}", work / f"t5_{k}"
        gfc("gff3-to-gtf", "--input-dir", str(adir), "--output-dir", str(t4),
            "--pattern", "*.gff3", "--force")
        R.suite(ds, "T4", t4, [gtf_gffread_style, gtf_split_exon],
                [gtf_exon_shift, gtf_strand_flip, gtf_drop_transcript, gtf_cds_trim])
        gfc("gff3-to-bed12", "--input-dir", str(adir), "--output-dir", str(t5),
            "--pattern", "*.gff3", "--force")
        R.suite(ds, "T5", t5, [bed_score_and_order],
                [bed_start_off_by_one, bed_block_size, bed_strand, bed_thick, bed_rename_one])
        bed6 = copy_dir(t5, R._dir("T5_bed6"))
        for f in bed6.glob("*.bed12"):
            _write(f, ["\t".join(l.split("\t")[:6]) for l in _read(f) if l])
        R.check(ds, "T5", "limitation", "BED6 reference", fs.LIMITATION, t5, bed6)

    if args.genome_dir:
        gdir = Path(args.genome_dir)
        ds = f"{gdir.name} (FASTA + GFF3)"
        t3 = work / "t3"
        gfc("fasta-gff-to-gbk", "--input-dir", str(gdir), "--output-dir", str(t3), "--force")
        R.suite(ds, "T3", t3, [], [gbk_point_mutation, gbk_shift_feature, gbk_drop_feature])
    return R.cases


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf")
    ap.add_argument("--sites", type=int, default=20000)
    ap.add_argument("--annotation-dir", action="append")
    ap.add_argument("--genome-dir")
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", help="write the case table as TSV here")
    ap.add_argument("--keep", help="working directory to keep (default: temporary)")
    args = ap.parse_args()
    if not (args.vcf or args.annotation_dir or args.genome_dir):
        ap.error("give at least one of --vcf, --annotation-dir, --genome-dir")
    work = Path(args.keep) if args.keep else Path(tempfile.mkdtemp(prefix="gfc_mutreal_"))
    if args.keep:
        shutil.rmtree(work, ignore_errors=True)
        work.mkdir(parents=True)
    cases = run(args, work)
    head = ("dataset", "task", "kind", "case", "expected", "got", "ok", "where", "seconds", "reason")
    rows = [head] + [(c.dataset, c.task, c.kind, c.name, c.expected, c.got,
                      "n/a" if c.got == NA else ("yes" if c.ok else "NO"), c.where,
                      f"{c.seconds:.2f}", c.reason.replace("\t", " ")) for c in cases]
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text("\n".join("\t".join(map(str, r)) for r in rows) + "\n")
    run_cases = [c for c in cases if c.got != NA]
    bad = [c for c in run_cases if not c.ok]
    muts = [c for c in run_cases if c.kind == "mutation"]
    print(f"\n{len(run_cases)} cases run ({len(cases) - len(run_cases)} not applicable), "
          f"{len(bad)} unexpected; defects caught {sum(c.ok for c in muts)}/{len(muts)}")
    for c in bad:
        print(f"  UNEXPECTED {c.dataset} {c.task} {c.name} ({c.where}): {c.got}: {c.reason}")
    if not args.keep:
        shutil.rmtree(work, ignore_errors=True)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
