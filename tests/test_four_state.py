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
