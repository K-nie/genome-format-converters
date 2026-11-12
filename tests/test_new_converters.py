"""eigenstrat-to-vcf and vcf-to-snp-alignment."""

import subprocess
import sys
from pathlib import Path

import pytest

cyvcf2 = pytest.importorskip("cyvcf2")

BENCH = Path(__file__).resolve().parents[1] / "benchmarks"
sys.path.insert(0, str(BENCH))
import mutation_test  # noqa: E402


def gfc(*args):
    subprocess.run([sys.executable, "-m", "genome_format_converters.cli", *args],
                   check=True, capture_output=True)


@pytest.fixture
def toy_vcf(tmp_path):
    d = tmp_path / "vcf"
    d.mkdir()
    mutation_test.make_vcf(d / "toy.vcf", n_sites=60, n_samples=9)
    return d


def test_eigenstrat_vcf_round_trip(tmp_path, toy_vcf):
    gfc("vcf-to-eigenstrat", "--input-dir", str(toy_vcf), "--output-dir", str(tmp_path / "e1"))
    gfc("eigenstrat-to-vcf", "--input-dir", str(tmp_path / "e1"), "--output-dir", str(tmp_path / "v"))
    out_vcf = tmp_path / "v" / "toy.vcf"
    recs = list(cyvcf2.VCF(str(out_vcf)))  # parses as valid VCF
    assert len(recs) == len((tmp_path / "e1" / "toy.snp").read_text().splitlines())
    gfc("vcf-to-eigenstrat", "--input-dir", str(tmp_path / "v"), "--output-dir", str(tmp_path / "e2"))
    for ext in (".geno", ".ind"):
        assert (tmp_path / "e1" / f"toy{ext}").read_bytes() == (tmp_path / "e2" / f"toy{ext}").read_bytes()
    s1 = [l.split("\t") for l in (tmp_path / "e1" / "toy.snp").read_text().splitlines()]
    s2 = [l.split("\t") for l in (tmp_path / "e2" / "toy.snp").read_text().splitlines()]
    assert [r[1:] for r in s1] == [r[1:] for r in s2]


def test_eigenstrat_to_vcf_haploid(tmp_path, toy_vcf):
    gfc("vcf-to-pseudohaploid", "--input-dir", str(toy_vcf), "--output-dir", str(tmp_path / "p"),
        "--seed", "1")
    gfc("eigenstrat-to-vcf", "--input-dir", str(tmp_path / "p"), "--output-dir", str(tmp_path / "v"),
        "--haploid")
    geno = (tmp_path / "p" / "toy.geno").read_text().split()
    for rec, row in zip(cyvcf2.VCF(str(tmp_path / "v" / "toy.vcf")), geno):
        expect = [{"2": 0, "0": 3, "9": 2}[c] for c in row]
        assert list(rec.gt_types) == expect


def test_rejects_packed_geno(tmp_path):
    d = tmp_path / "eig"
    d.mkdir()
    (d / "x.geno").write_bytes(b"GENO 2 3 abc\x00\x01")
    (d / "x.snp").write_text("s1\t1\t0.0\t10\tA\tC\n")
    (d / "x.ind").write_text("a\tU\tp\nb\tU\tp\n")
    r = subprocess.run([sys.executable, "-m", "genome_format_converters.cli", "eigenstrat-to-vcf",
                        "--input-dir", str(d), "--output-dir", str(tmp_path / "o")],
                       capture_output=True, text=True)
    assert r.returncode != 0 and "packed" in r.stderr


def _write_vcf(path, rows, samples=("a", "b", "c")):
    head = ("##fileformat=VCFv4.2\n##contig=<ID=1>\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
    path.write_text(head + "".join(f"1\t{p}\t.\t{r}\t{a}\t.\t.\t.\tGT\t" + "\t".join(g) + "\n"
                                   for p, r, a, g in rows))


ROWS = [
    (10, "A", "G", ["0/0", "0/1", "1/1"]),
    (20, "C", "T", ["./.", "0|0", "1|1"]),
    (30, "A", "AT", ["0/0", "0/1", "1/1"]),   # indel: skipped
    (40, "G", "C", ["0/0", "0/0", "0/0"]),    # invariant
]


def test_snp_alignment_phylip_iupac(tmp_path):
    d = tmp_path / "in"
    d.mkdir()
    _write_vcf(d / "s.vcf", ROWS)
    gfc("vcf-to-snp-alignment", "--input-dir", str(d), "--output-dir", str(tmp_path / "o"))
    lines = (tmp_path / "o" / "s.phy").read_text().splitlines()
    assert lines[0] == "3 3"
    assert [l.split() for l in lines[1:]] == [["a", "ANG"], ["b", "RCG"], ["c", "GTG"]]


def test_snp_alignment_options(tmp_path):
    d = tmp_path / "in"
    d.mkdir()
    _write_vcf(d / "s.vcf", ROWS)
    gfc("vcf-to-snp-alignment", "--input-dir", str(d), "--output-dir", str(tmp_path / "o"),
        "--format", "fasta", "--het", "missing", "--min-samples", "3", "--drop-invariant")
    text = (tmp_path / "o" / "s.fasta").read_text()
    assert text == ">a\nA\n>b\nN\n>c\nG\n"
    gfc("vcf-to-snp-alignment", "--input-dir", str(d), "--output-dir", str(tmp_path / "n"),
        "--format", "nexus")
    nex = (tmp_path / "n" / "s.nex").read_text()
    assert "NTAX=3 NCHAR=3" in nex and nex.startswith("#NEXUS")
