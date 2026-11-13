#!/usr/bin/env python3
"""Mutation test for the four-state correctness comparators.

For each benchmark task this script

1. builds a small synthetic input and runs gfc on it;
2. compares gfc's output with an exact copy (expects ``byte``);
3. compares it with copies rewritten to a *different but valid*
   convention, the kind other tools use (expects ``structural``);
4. compares *mutated* copies carrying a seeded defect against the
   original (expects ``disagreement``).

A comparator that calls a seeded defect "structural" would be hiding real
bugs behind convention tolerance; this test is the evidence that it does
not. Run it with ``python benchmarks/mutation_test.py`` (prints a table and
writes ``--out`` TSV) or through ``tests/test_four_state.py``.
"""
from __future__ import annotations

import argparse
import random
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, List

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import four_state as fs  # noqa: E402

GFC = [sys.executable, "-c", "from genome_format_converters.cli import main; main()"]


def gfc(*args: str) -> None:
    subprocess.run([*GFC, *args], check=True, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)


@dataclass
class Case:
    task: str
    kind: str        # identity | convention | mutation
    name: str
    expected: str
    got: str = ""
    reason: str = ""

    @property
    def ok(self) -> bool:
        return self.got == self.expected


# ------------------------------------------------------------- fixtures


def make_vcf(path: Path, n_sites=120, n_samples=14, seed=3) -> None:
    rnd = random.Random(seed)
    gts = ["0|0"] * 4 + ["0|1", "1|0", "1|1", "./."]
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n##contig=<ID=chr22>\n"
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                + "\t".join(f"S{i}" for i in range(n_samples)) + "\n")
        for i in range(n_sites):
            r, a = rnd.sample("ACGT", 2)
            if i % 25 == 7:
                a += "G"  # an indel that gfc skips
            f.write(f"chr22\t{1000 + 37 * i}\t.\t{r}\t{a}\t.\tPASS\t.\tGT\t"
                    + "\t".join(rnd.choice(gts) for _ in range(n_samples)) + "\n")


def make_gff3(path: Path, n_genes=8, seed=5) -> None:
    rnd = random.Random(seed)
    lines = ["##gff-version 3"]
    pos = 100
    for g in range(n_genes):
        strand = "+" if g % 2 == 0 else "-"
        n_ex = 1 + g % 3
        exons, s = [], pos
        for _ in range(n_ex):
            ln = rnd.randint(60, 200)
            exons.append((s, s + ln - 1))
            s = s + ln + rnd.randint(30, 120)
        start, end = exons[0][0], exons[-1][1]
        gid, tid = f"gene{g}", f"mRNA{g}"
        lines.append(f"ctg1\t.\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gid};Name={gid}")
        lines.append(f"ctg1\t.\tmRNA\t{start}\t{end}\t.\t{strand}\t.\tID={tid};Parent={gid}")
        for k, (a, b) in enumerate(exons):
            lines.append(f"ctg1\t.\texon\t{a}\t{b}\t.\t{strand}\t.\tID={tid}.e{k};Parent={tid}")
        for k, (a, b) in enumerate(exons):
            lines.append(f"ctg1\t.\tCDS\t{a}\t{b}\t.\t{strand}\t0\tID={tid}.c{k};Parent={tid}")
        pos = end + 300
    path.write_text("\n".join(lines) + "\n")
    return pos


def make_fasta(path: Path, length: int, seed=9, name="ctg1") -> None:
    rnd = random.Random(seed)
    seq = "".join(rnd.choice("ACGT") for _ in range(length))
    path.write_text(f">{name}\n" + "\n".join(seq[i:i + 60] for i in range(0, length, 60)) + "\n")


def make_tblout(path: Path, n=30, seed=11) -> None:
    rnd = random.Random(seed)
    out = ["# target name  accession  query name  accession  E-value  score  bias  "
           "E-value  score  bias  exp reg clu  ov env dom rep inc description of target",
           "#------------"]
    for i in range(n):
        e = 10 ** -rnd.randint(3, 40)
        sc = rnd.uniform(20, 500)
        out.append(f"prot{i:03d}  -  PF{rnd.randint(1, 9999):05d}.{rnd.randint(1, 20)}  -  "
                   f"{e:.1e}  {sc:.1f}  {rnd.uniform(0, 5):.1f}  {e:.1e}  {sc:.1f}  "
                   f"{rnd.uniform(0, 5):.1f}  1.0  1  0  0  1  1  1  1  domain of protein {i}")
    out.append("#\n# Program: hmmscan")
    path.write_text("\n".join(out) + "\n")


def make_orthogroups(d: Path, seed=13) -> Path:
    rnd = random.Random(seed)
    fasta_dir = d / "fastas"
    fasta_dir.mkdir()
    species = ["sp_a", "sp_b", "sp_c"]
    genes = {s: [f"{s}_g{i}" for i in range(12)] for s in species}
    for s in species:
        with open(fasta_dir / f"{s}.fa", "w") as f:
            for g in genes[s]:
                seq = "".join(rnd.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(rnd.randint(80, 160)))
                f.write(f">{g}\n{seq}\n")
    og = d / "Orthogroups.tsv"
    with open(og, "w") as f:
        f.write("Orthogroup\t" + "\t".join(species) + "\n")
        for k in range(6):
            cols = [", ".join(genes[s][2 * k:2 * k + 2]) for s in species]
            f.write(f"OG{k:07d}\t" + "\t".join(cols) + "\n")
    return og


# ------------------------------------------------------------- rewriting


def copy_dir(src: Path, dst: Path) -> Path:
    shutil.copytree(src, dst)
    return dst


def rewrite(d: Path, glob: str, fn: Callable[[str], str]) -> None:
    for f in d.glob(glob):
        f.write_text(fn(f.read_text()))


def rewrite_lines(d: Path, glob: str, fn: Callable[[int, List[str]], List[str]]) -> None:
    for f in d.glob(glob):
        lines = f.read_text().splitlines()
        f.write_text("\n".join(fn(0, lines)) + "\n")


def _first_data(lines: List[str]) -> int:
    for i, l in enumerate(lines):
        if l and not l.startswith("#"):
            return i
    raise ValueError("no data line")


# -- EIGENSTRAT


def eig_chrom_label(d):
    rewrite(d, "*.snp", lambda t: re.sub(r"\tchr22\t", "\t22\t", t))


def eig_allele_swap(d):
    snp = next(d.glob("*.snp"))
    geno = next(d.glob("*.geno"))
    s = snp.read_text().splitlines()
    g = geno.read_text().splitlines()
    s = ["\t".join(p[:4] + [p[5], p[4]]) for p in (l.split("\t") for l in s)]
    g = [fs._flip_geno(x) for x in g]
    snp.write_text("\n".join(s) + "\n")
    geno.write_text("\n".join(g) + "\n")


def eig_flip_one_call(d):
    geno = next(d.glob("*.geno"))
    g = geno.read_text().splitlines()
    row = g[5]
    i = next(k for k, c in enumerate(row) if c in "02")
    g[5] = row[:i] + ("2" if row[i] == "0" else "0") + row[i + 1:]
    geno.write_text("\n".join(g) + "\n")


def eig_pos_off_by_one(d):
    snp = next(d.glob("*.snp"))
    s = snp.read_text().splitlines()
    p = s[10].split("\t")
    p[3] = str(int(p[3]) + 1)
    s[10] = "\t".join(p)
    snp.write_text("\n".join(s) + "\n")


def eig_swap_alleles_only(d):
    snp = next(d.glob("*.snp"))
    s = snp.read_text().splitlines()
    p = s[3].split("\t")
    p[4], p[5] = p[5], p[4]
    s[3] = "\t".join(p)
    snp.write_text("\n".join(s) + "\n")


def eig_drop_site(d):
    for pat in ("*.snp", "*.geno"):
        f = next(d.glob(pat))
        l = f.read_text().splitlines()
        del l[20]
        f.write_text("\n".join(l) + "\n")


def eig_swap_samples(d):
    ind = next(d.glob("*.ind"))
    l = ind.read_text().splitlines()
    l[0], l[1] = l[1], l[0]
    ind.write_text("\n".join(l) + "\n")


# -- PLINK


def _bed_parts(d):
    bed, bim, fam = next(d.glob("*.bed")), next(d.glob("*.bim")), next(d.glob("*.fam"))
    n = len(fam.read_text().splitlines())
    return bed, bim, (n + 3) // 4


def plink_allele_swap(d):
    bed, bim, per = _bed_parts(d)
    raw = bytearray(bed.read_bytes())
    for i in range(3, len(raw)):
        b = raw[i]
        out = 0
        for k in range(4):
            c = (b >> (2 * k)) & 3
            c = {0: 3, 3: 0}.get(c, c)
            out |= c << (2 * k)
        raw[i] = out
    # Padding bits in the last byte of each row were 00 and are now 11;
    # PLINK ignores them, and so does the comparator.
    bed.write_bytes(bytes(raw))
    lines = [l.split("\t") for l in bim.read_text().splitlines()]
    bim.write_text("\n".join("\t".join(p[:4] + [p[5], p[4]]) for p in lines) + "\n")


def plink_chrom_label(d):
    rewrite(d, "*.bim", lambda t: re.sub(r"(?m)^chr22\t", "22\t", t))


def plink_flip_genotype(d):
    bed, _bim, per = _bed_parts(d)
    raw = bytearray(bed.read_bytes())
    i = 3 + per * 4
    for shift in range(0, 8, 2):
        c = (raw[i] >> shift) & 3
        if c in (0, 3):
            raw[i] ^= (3 << shift)
            break
    bed.write_bytes(bytes(raw))


def plink_pos_shift(d):
    bim = next(d.glob("*.bim"))
    l = bim.read_text().splitlines()
    p = l[2].split("\t")
    p[3] = str(int(p[3]) + 1)
    l[2] = "\t".join(p)
    bim.write_text("\n".join(l) + "\n")


def plink_swap_bim_only(d):
    bim = next(d.glob("*.bim"))
    l = bim.read_text().splitlines()
    p = l[4].split("\t")
    p[4], p[5] = p[5], p[4]
    l[4] = "\t".join(p)
    bim.write_text("\n".join(l) + "\n")


# -- GTF


def gtf_gffread_style(d):
    def fn(t):
        out = []
        for line in t.splitlines():
            p = line.split("\t")
            if len(p) < 9:
                out.append(line)
                continue
            if p[2] == "gene":
                continue  # gffread elides gene rows
            if p[2] == "mRNA":
                p[2] = "transcript"
            p[8] = p[8].rstrip().rstrip(";")
            out.append("\t".join(p))
        return "\n".join(out) + "\n"
    rewrite(d, "*.gtf", fn)


def gtf_split_exon(d):
    """Split one exon into two abutting halves (the inverse of fusion)."""
    def fn(t):
        out, done = [], False
        for line in t.splitlines():
            p = line.split("\t")
            if not done and len(p) >= 9 and p[2] == "exon" and int(p[4]) - int(p[3]) > 10:
                mid = (int(p[3]) + int(p[4])) // 2
                out.append("\t".join(p[:4] + [str(mid)] + p[5:]))
                out.append("\t".join(p[:3] + [str(mid + 1)] + p[4:]))
                done = True
            else:
                out.append(line)
        return "\n".join(out) + "\n"
    rewrite(d, "*.gtf", fn)


def _gtf_edit(d, pred, edit):
    def fn(t):
        out, done = [], False
        for line in t.splitlines():
            p = line.split("\t")
            if not done and len(p) >= 9 and pred(p):
                r = edit(p)
                done = True
                if r is None:
                    continue
                p = r
            out.append("\t".join(p) if isinstance(p, list) else line)
        return "\n".join(out) + "\n"
    rewrite(d, "*.gtf", fn)


def gtf_exon_shift(d):
    _gtf_edit(d, lambda p: p[2] == "exon",
              lambda p: p[:3] + [str(int(p[3]) + 1)] + p[4:])


def gtf_strand_flip(d):
    def fn(t):
        return re.sub(r'(\t)\+(\t[^\n]*transcript_id "mRNA2")', r"\1-\2", t)
    rewrite(d, "*.gtf", fn)


def gtf_drop_transcript(d):
    rewrite(d, "*.gtf", lambda t: "\n".join(
        l for l in t.splitlines() if 'transcript_id "mRNA3"' not in l) + "\n")


def gtf_cds_trim(d):
    _gtf_edit(d, lambda p: p[2] == "CDS",
              lambda p: p[:4] + [str(int(p[4]) - 3)] + p[5:])


# -- BED12


def bed_score_and_order(d):
    def fn(t):
        rows = [l.split("\t") for l in t.splitlines() if l]
        for r in rows:
            r[4] = "1000"
            r[8] = "255,0,0"
        rows.reverse()
        return "\n".join("\t".join(r) for r in rows) + "\n"
    rewrite(d, "*.bed12", fn)


def _bed_edit(d, fn_row):
    def fn(t):
        rows = [l.split("\t") for l in t.splitlines() if l]
        fn_row(rows[2])
        return "\n".join("\t".join(r) for r in rows) + "\n"
    rewrite(d, "*.bed12", fn)


def bed_start_off_by_one(d):
    _bed_edit(d, lambda r: r.__setitem__(1, str(int(r[1]) + 1)))


def bed_block_size(d):
    def f(r):
        sizes = r[10].rstrip(",").split(",")
        sizes[0] = str(int(sizes[0]) - 1)
        r[10] = ",".join(sizes) + ","
    _bed_edit(d, f)


def bed_strand(d):
    _bed_edit(d, lambda r: r.__setitem__(5, "-" if r[5] == "+" else "+"))


def bed_thick(d):
    _bed_edit(d, lambda r: r.__setitem__(6, str(int(r[6]) + 3)))


# -- TSV (HMMER)


def tsv_number_format(d):
    def fn(t):
        lines = t.splitlines()
        out = [lines[0]]
        for l in lines[1:]:
            p = l.split("\t")
            p = [re.sub(r"^(\d+)\.(\d)e-(\d+)$", r"\1.\g<2>0e-\3", x) for x in p]
            out.append("\t".join(p))
        return "\n".join(out) + "\n"
    rewrite(d, "*.tsv", fn)


def tsv_score_change(d):
    def fn(t):
        lines = t.splitlines()
        p = lines[3].split("\t")
        p[5] = str(float(p[5]) + 1.0)
        lines[3] = "\t".join(p)
        return "\n".join(lines) + "\n"
    rewrite(d, "*.tsv", fn)


def tsv_drop_row(d):
    rewrite(d, "*.tsv", lambda t: "\n".join(
        l for i, l in enumerate(t.splitlines()) if i != 4) + "\n")


def tsv_swap_columns(d):
    def fn(t):
        out = []
        for l in t.splitlines():
            p = l.split("\t")
            p[0], p[2] = p[2], p[0]
            out.append("\t".join(p))
        return "\n".join(out) + "\n"
    rewrite(d, "*.tsv", fn)


# -- FASTA (orthogroups)


def fa_rewrap(d):
    for f in d.glob("*.fa"):
        recs = fs._fasta(f)
        f.write_text("".join(f">{n}\n" + "\n".join(s[i:i + 50] for i in range(0, len(s), 50)) + "\n"
                             for n, s in reversed(recs)))


def fa_point_mutation(d):
    f = sorted(d.glob("*.fa"))[0]
    lines = f.read_text().splitlines()
    s = lines[1]
    lines[1] = ("A" if s[0] != "A" else "C") + s[1:]
    f.write_text("\n".join(lines) + "\n")


def fa_drop_record(d):
    f = sorted(d.glob("*.fa"))[1]
    recs = fs._fasta(f)[1:]
    f.write_text("".join(f">{n}\n{s}\n" for n, s in recs))


def fa_mislabel(d):
    f = sorted(d.glob("*.fa"))[2]
    f.write_text(re.sub(r">sp_a\|", ">sp_b|", f.read_text(), count=1))


# -- GenBank


def gbk_point_mutation(d):
    f = next(d.glob("*.gbk"))
    t = f.read_text()
    i = t.index("ORIGIN")
    m = re.search(r"\d+ ([acgt])", t[i:])
    j = i + m.start(1)
    t = t[:j] + ("a" if t[j] != "a" else "c") + t[j + 1:]
    f.write_text(t)


def gbk_shift_feature(d):
    f = next(d.glob("*.gbk"))
    t = f.read_text()
    t = re.sub(r"(     CDS             (?:complement\()?(?:join\()?)(\d+)",
               lambda m: m.group(1) + str(int(m.group(2)) + 1), t, count=1)
    f.write_text(t)


def gbk_drop_feature(d):
    f = next(d.glob("*.gbk"))
    lines = f.read_text().splitlines()
    out, skipping = [], False
    dropped = False
    for l in lines:
        if not dropped and l.startswith("     exon"):
            skipping, dropped = True, True
            continue
        if skipping and l.startswith("                     "):
            continue
        skipping = False
        out.append(l)
    f.write_text("\n".join(out) + "\n")


# -- pseudohaploid


def _gt_classes(vcf):
    """(norm chrom, pos) -> (REF, ALT, class per sample) for a small VCF."""
    names = {0: "ref", 1: "het", 2: "miss", 3: "alt"}
    return {(fs._norm_chrom(c), int(pos)): (r, a, [names[int(x)] for x in gt])
            for c, pos, r, a, gt in fs._iter_vcf_sites(vcf)}


def ph_homref_flip(d, vcf):
    src = _gt_classes(vcf)
    snp = next(d.glob("*.snp")).read_text().splitlines()
    geno_f = next(d.glob("*.geno"))
    g = geno_f.read_text().splitlines()
    for i, s in enumerate(snp):
        p = s.split("\t")
        cls = src[(fs._norm_chrom(p[1]), int(p[3]))][2]
        if "ref" in cls:
            j = cls.index("ref")
            g[i] = g[i][:j] + ("0" if g[i][j] == "2" else "2") + g[i][j + 1:]
            break
    geno_f.write_text("\n".join(g) + "\n")


def ph_missing_filled(d, vcf):
    src = _gt_classes(vcf)
    snp = next(d.glob("*.snp")).read_text().splitlines()
    geno_f = next(d.glob("*.geno"))
    g = geno_f.read_text().splitlines()
    for i, s in enumerate(snp):
        p = s.split("\t")
        cls = src[(fs._norm_chrom(p[1]), int(p[3]))][2]
        if "miss" in cls:
            j = cls.index("miss")
            g[i] = g[i][:j] + "2" + g[i][j + 1:]
            break
    geno_f.write_text("\n".join(g) + "\n")


def ph_biased_sampler(d, vcf):
    """Every het resolved to REF: a sampler that is not random."""
    src = _gt_classes(vcf)
    snp = next(d.glob("*.snp")).read_text().splitlines()
    geno_f = next(d.glob("*.geno"))
    g = geno_f.read_text().splitlines()
    for i, s in enumerate(snp):
        p = s.split("\t")
        cls = src[(fs._norm_chrom(p[1]), int(p[3]))][2]
        g[i] = "".join("2" if c == "het" else x for x, c in zip(g[i], cls))
    geno_f.write_text("\n".join(g) + "\n")


def ph_het_as_one(d, vcf):
    src = _gt_classes(vcf)
    snp = next(d.glob("*.snp")).read_text().splitlines()
    geno_f = next(d.glob("*.geno"))
    g = geno_f.read_text().splitlines()
    for i, s in enumerate(snp):
        p = s.split("\t")
        cls = src[(fs._norm_chrom(p[1]), int(p[3]))][2]
        if "het" in cls:
            j = cls.index("het")
            g[i] = g[i][:j] + "1" + g[i][j + 1:]
            break
    geno_f.write_text("\n".join(g) + "\n")


# ------------------------------------------------------------- driver


def run(work: Path) -> List[Case]:
    cases: List[Case] = []
    inp = work / "in"
    inp.mkdir()

    def check(task, kind, name, expected, cand, ref, **ctx):
        v = fs.COMPARATORS[task](cand, ref, **ctx)
        cases.append(Case(task, kind, name, expected, v.state, v.reason))

    def suite(task, base: Path, conventions, mutations, **ctx):
        check(task, "identity", "exact copy", fs.BYTE,
              base, copy_dir(base, work / f"{task}_identity"), **ctx)
        for fn in conventions:
            ref = copy_dir(base, work / f"{task}_{fn.__name__}")
            fn(ref)
            check(task, "convention", fn.__name__, fs.STRUCTURAL, base, ref, **ctx)
        for fn in mutations:
            cand = copy_dir(base, work / f"{task}_{fn.__name__}")
            fn(cand)
            check(task, "mutation", fn.__name__, fs.DISAGREEMENT, cand, base, **ctx)

    # Variant tasks
    vcf_dir = inp / "vcf"
    vcf_dir.mkdir()
    vcf = vcf_dir / "toy.vcf"
    make_vcf(vcf)
    gfc("vcf-to-eigenstrat", "--input-dir", str(vcf_dir), "--output-dir", str(work / "t1"), "--force")
    suite("T1", work / "t1", [eig_chrom_label, eig_allele_swap],
          [eig_flip_one_call, eig_pos_off_by_one, eig_swap_alleles_only,
           eig_drop_site, eig_swap_samples])

    gfc("vcf-to-plink", "--input-dir", str(vcf_dir), "--output-dir", str(work / "t2"), "--force")
    suite("T2", work / "t2", [plink_allele_swap, plink_chrom_label],
          [plink_flip_genotype, plink_pos_shift, plink_swap_bim_only])

    vcf_big_dir = inp / "vcf_big"
    vcf_big_dir.mkdir()
    vcf_big = vcf_big_dir / "toy_big.vcf"
    make_vcf(vcf_big, n_sites=400, n_samples=40, seed=21)
    gfc("vcf-to-pseudohaploid", "--input-dir", str(vcf_big_dir), "--output-dir",
        str(work / "t8"), "--seed", "42", "--force")
    t8 = work / "t8"
    check("T8", "identity", "gfc output vs input", fs.STRUCTURAL, t8, None, input_vcf=vcf_big)
    for fn in (ph_homref_flip, ph_missing_filled, ph_biased_sampler, ph_het_as_one):
        cand = copy_dir(t8, work / f"T8_{fn.__name__}")
        fn(cand, vcf_big)
        check("T8", "mutation", fn.__name__, fs.DISAGREEMENT, cand, None, input_vcf=vcf_big)

    # Annotation tasks
    ann = inp / "ann"
    ann.mkdir()
    end = make_gff3(ann / "toy.gff3")
    gfc("gff3-to-gtf", "--input-dir", str(ann), "--output-dir", str(work / "t4"), "--force")
    suite("T4", work / "t4", [gtf_gffread_style, gtf_split_exon],
          [gtf_exon_shift, gtf_strand_flip, gtf_drop_transcript, gtf_cds_trim])

    gfc("gff3-to-bed12", "--input-dir", str(ann), "--output-dir", str(work / "t5"), "--force")
    suite("T5", work / "t5", [bed_score_and_order],
          [bed_start_off_by_one, bed_block_size, bed_strand, bed_thick])

    gb_in = inp / "gb"
    gb_in.mkdir()
    shutil.copy(ann / "toy.gff3", gb_in / "toy.gff3")
    make_fasta(gb_in / "toy.fasta", end + 100)
    gfc("fasta-gff-to-gbk", "--input-dir", str(gb_in), "--output-dir", str(work / "t3"), "--force")
    suite("T3", work / "t3", [],
          [gbk_point_mutation, gbk_shift_feature, gbk_drop_feature])

    hm = inp / "hmm"
    hm.mkdir()
    make_tblout(hm / "toy.tblout")
    gfc("hmmer-tblout-to-tsv", "--hmmer-format", "tblout", "--input-dir", str(hm),
        "--output-dir", str(work / "t6"), "--force")
    suite("T6", work / "t6", [tsv_number_format],
          [tsv_score_change, tsv_drop_row, tsv_swap_columns])

    og_dir = inp / "og"
    og_dir.mkdir()
    og = make_orthogroups(og_dir)
    gfc("orthogroups-to-fasta", "--orthogroups", str(og), "--fasta-dir",
        str(og_dir / "fastas"), "--output-dir", str(work / "t7"), "--force")
    suite("T7", work / "t7", [fa_rewrap],
          [fa_point_mutation, fa_drop_record, fa_mislabel])

    # Limitation: a BED6 reference against BED12 output.
    bed6 = copy_dir(work / "t5", work / "T5_bed6")
    rewrite(bed6, "*.bed12", lambda t: "\n".join(
        "\t".join(l.split("\t")[:6]) for l in t.splitlines() if l) + "\n")
    check("T5", "limitation", "BED6 reference", fs.LIMITATION, work / "t5", bed6)
    return cases


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", help="write the case table as TSV here")
    ap.add_argument("--keep", help="keep the working directory here")
    args = ap.parse_args()
    work = Path(args.keep) if args.keep else Path(tempfile.mkdtemp(prefix="gfc_mut_"))
    if args.keep:
        shutil.rmtree(work, ignore_errors=True)
        work.mkdir(parents=True)
    cases = run(work)
    rows = [("task", "kind", "case", "expected", "got", "ok", "reason")]
    rows += [(c.task, c.kind, c.name, c.expected, c.got, "yes" if c.ok else "NO", c.reason)
             for c in cases]
    w = [max(len(str(r[i])) for r in rows) for i in range(6)]
    for r in rows:
        print("  ".join(str(x).ljust(w[i]) for i, x in enumerate(r[:6])) + "  " + r[6])
    if args.out:
        Path(args.out).write_text("\n".join("\t".join(map(str, r)) for r in rows) + "\n")
    n_bad = sum(not c.ok for c in cases)
    n_mut = sum(c.kind == "mutation" for c in cases)
    caught = sum(c.kind == "mutation" and c.ok for c in cases)
    print(f"\n{len(cases)} cases, {n_bad} unexpected; mutations caught {caught}/{n_mut}")
    if not args.keep:
        shutil.rmtree(work, ignore_errors=True)
    return 1 if n_bad else 0


if __name__ == "__main__":
    sys.exit(main())
