"""The four-state comparators: conventions -> structural, defects -> disagreement."""

import csv
import shutil
import sys
from pathlib import Path

import pytest

pytest.importorskip("cyvcf2")

BENCH = Path(__file__).resolve().parents[1] / "benchmarks"
sys.path.insert(0, str(BENCH))

import four_state as fs  # noqa: E402
import mutation_test  # noqa: E402
import check_correctness  # noqa: E402


@pytest.fixture(scope="module")
def cases(tmp_path_factory):
    return mutation_test.run(tmp_path_factory.mktemp("mut"))


def test_every_case_gets_expected_state(cases):
    bad = [(c.task, c.name, c.expected, c.got, c.reason) for c in cases if not c.ok]
    assert not bad


def test_every_task_has_mutations(cases):
    tasks = {c.task for c in cases if c.kind == "mutation"}
    assert tasks == {f"T{i}" for i in range(1, 9)}


def test_populate_writes_four_state_sidecar(tmp_path):
    work = tmp_path / "w"
    work.mkdir()
    vcf_dir = work / "in"
    vcf_dir.mkdir()
    mutation_test.make_vcf(vcf_dir / "toy.vcf")
    bench = tmp_path / "bench"
    mutation_test.gfc("vcf-to-eigenstrat", "--input-dir", str(vcf_dir),
                      "--output-dir", str(bench / "gfc_rep1"), "--force")
    shutil.copytree(bench / "gfc_rep1", bench / "convertf_rep1")
    mutation_test.eig_chrom_label(bench / "convertf_rep1")
    tsv = tmp_path / "T1.tsv"
    with tsv.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["task", "tool", "version", "replicate", "wall_s",
                    "peak_rss_mb", "exit_code", "correct", "notes"])
        for tool in ("gfc", "convertf"):
            w.writerow(["T1", tool, "x", "1", "1", "1", "0", "", ""])
    check_correctness.populate("T1", bench, tsv)
    rows = list(csv.DictReader(tsv.open(), delimiter="\t"))
    assert {r["tool"]: r["correct"] for r in rows} == {"gfc": "1", "convertf": "1"}
    side = list(csv.DictReader(tsv.with_suffix(".verdicts.tsv").open(), delimiter="\t"))
    gfc_row = next(r for r in side if r["tool"] == "gfc")
    assert gfc_row["state"] == fs.STRUCTURAL
    assert gfc_row["reason"] == "chrom-label"


def test_populate_sidecar_names_conventions_and_reference(tmp_path):
    test_populate_writes_four_state_sidecar(tmp_path)
    side = list(csv.DictReader((tmp_path / "T1.verdicts.tsv").open(), delimiter="\t"))
    by_tool = {r["tool"]: r for r in side}
    assert by_tool["gfc"]["conventions"] == "chrom-label"
    assert by_tool["convertf"]["state"] == "reference"


# ---- split multi-allelic sites: two biallelic records at one position ----

_DUP_VCF = """##fileformat=VCFv4.2
##contig=<ID=22>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc\td\te\tf
22\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1\t0|1\t0|0\t1|0
22\t200\t.\tG\tT\t.\tPASS\t.\tGT\t0|1\t0|0\t0|0\t1|1\t0|1\t0|0
22\t200\t.\tG\tA\t.\tPASS\t.\tGT\t0|0\t1|1\t0|1\t0|0\t0|0\t0|1
22\t300\t.\tC\tT\t.\tPASS\t.\tGT\t1|1\t0|1\t0|0\t0|0\t1|0\t0|1
"""


@pytest.fixture()
def dup_vcf(tmp_path):
    d = tmp_path / "in"
    d.mkdir()
    (d / "d.vcf").write_text(_DUP_VCF)
    return d


def test_t1_split_multiallelic_is_not_a_disagreement(dup_vcf, tmp_path):
    out = tmp_path / "eig"
    mutation_test.gfc("vcf-to-eigenstrat", "--input-dir", str(dup_vcf),
                      "--output-dir", str(out), "--force")
    ref = tmp_path / "ref"
    shutil.copytree(out, ref)
    mutation_test.eig_allele_swap(ref)
    v = fs.compare_t1(out, ref)
    assert v.state == fs.STRUCTURAL, v.reason
    assert "allele-order" in v.conventions


def test_t2_split_multiallelic_is_not_a_disagreement(dup_vcf, tmp_path):
    out = tmp_path / "pl"
    mutation_test.gfc("vcf-to-plink", "--input-dir", str(dup_vcf),
                      "--output-dir", str(out), "--force")
    ref = tmp_path / "ref"
    shutil.copytree(out, ref)
    mutation_test.plink_allele_swap(ref)
    v = fs.compare_t2(out, ref)
    assert v.state == fs.STRUCTURAL, v.reason


def test_t2_still_catches_a_flipped_genotype_at_a_split_site(dup_vcf, tmp_path):
    out = tmp_path / "pl"
    mutation_test.gfc("vcf-to-plink", "--input-dir", str(dup_vcf),
                      "--output-dir", str(out), "--force")
    ref = tmp_path / "ref"
    shutil.copytree(out, ref)
    # Flip sample a at the third record (22:200 G>A): 00 (hom A1) <-> 11.
    bed = next(ref.glob("*.bed"))
    raw = bytearray(bed.read_bytes())
    i = 3 + 2 * 2  # 6 samples -> 2 bytes per variant
    raw[i] ^= 0b11
    bed.write_bytes(bytes(raw))
    v = fs.compare_t2(out, ref)
    assert v.state == fs.DISAGREEMENT and "22:200" in v.reason, v.reason


def test_pseudohaploid_check_handles_split_sites(dup_vcf, tmp_path):
    out = tmp_path / "ph"
    mutation_test.gfc("vcf-to-pseudohaploid", "--input-dir", str(dup_vcf),
                      "--output-dir", str(out), "--seed", "1", "--force")
    v = fs.check_pseudohaploid(out, dup_vcf / "d.vcf")
    assert v.state == fs.STRUCTURAL, v.reason
    # A hom-ref call turned into the ALT allele must be caught.
    geno = next(out.glob("*.geno"))
    rows = geno.read_text().splitlines()
    rows[2] = ("0" if rows[2][0] == "2" else "2") + rows[2][1:]  # sample a, 22:200 G>A, hom-ref
    geno.write_text("\n".join(rows) + "\n")
    assert fs.check_pseudohaploid(out, dup_vcf / "d.vcf").state == fs.DISAGREEMENT


def test_t8_reproducibility_with_threads_over_symlinks(dup_vcf, tmp_path):
    import subprocess
    res = subprocess.run(
        [sys.executable, str(BENCH / "t8_reproducibility.py"), "--vcf",
         str(dup_vcf / "d.vcf"), "--work", str(tmp_path / "w"), "--copies", "4",
         "--threads", "4", "--out", str(tmp_path / "r.tsv")],
        capture_output=True, text=True, timeout=300)
    assert res.returncode == 0, res.stderr
    rows = list(csv.DictReader((tmp_path / "r.tsv").open(), delimiter="\t"))
    assert all(r["pass"] == "yes" for r in rows), rows
