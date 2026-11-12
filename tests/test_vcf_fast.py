"""The fast VCF reader must agree exactly with cyvcf2."""

import gzip
import shutil
from pathlib import Path

import numpy as np
import pytest

cyvcf2 = pytest.importorskip("cyvcf2")

from genome_format_converters.converters import _vcf_fast  # noqa: E402
from genome_format_converters.converters.vcf_to_pseudohaploid import _PickBuffer  # noqa: E402

DATA = Path(__file__).parent / "test_data"
FILES = ["tiny.vcf", "gt_edge_cases.vcf"]


def _records(path):
    vcf = cyvcf2.VCF(str(path))
    return list(_vcf_fast.iter_records(path, vcf))


def _assert_same(path):
    fast = _records(path)
    ref = list(cyvcf2.VCF(str(path)))
    assert len(fast) == len(ref)
    for f, r in zip(fast, ref):
        assert (f.CHROM, f.POS, f.REF, f.ALT) == (r.CHROM, r.POS, r.REF, r.ALT)
        assert f.ID == r.ID
        assert f.INFO.get("AA") == r.INFO.get("AA")
        np.testing.assert_array_equal(f.gt_types, r.gt_types)


@pytest.mark.parametrize("name", FILES)
def test_fast_reader_matches_cyvcf2(name):
    _assert_same(DATA / name)


@pytest.mark.parametrize("name", FILES)
def test_fast_reader_matches_cyvcf2_gzip(tmp_path, name):
    gz = tmp_path / (name + ".gz")
    with open(DATA / name, "rb") as src, gzip.open(gz, "wb") as dst:
        shutil.copyfileobj(src, dst)
    _assert_same(gz)


def test_fast_path_used_for_regular_gt_block():
    block = b"0|0\t0|1\t1/0\t1|1\t./."
    np.testing.assert_array_equal(
        _vcf_fast._decode_gt_block(block, 5), [0, 1, 1, 3, 2]
    )


@pytest.mark.parametrize("block", [b"0/.\t0|1", b"0\t0|1", b"0|2\t0|1", b"0|0:3\t0|1"])
def test_irregular_blocks_defer_to_cyvcf2(block):
    assert _vcf_fast._decode_gt_block(block, 2) is None


@pytest.mark.parametrize("n", [1, 7, 2504, 2505])
def test_pick_buffer_matches_per_site_draws(n):
    a = np.random.default_rng(42)
    expected = np.concatenate(
        [a.integers(0, 2, size=n, dtype=np.int32) for _ in range(600)]
    )
    buf = _PickBuffer(np.random.default_rng(42), n)
    got = np.concatenate([buf.take(n).copy() for _ in range(600)])
    np.testing.assert_array_equal(got, expected)


def test_threads_output_matches_sequential(tmp_path):
    import subprocess
    import sys

    in_dir = tmp_path / "in"
    in_dir.mkdir()
    for i in range(3):
        shutil.copy(DATA / "gt_edge_cases.vcf", in_dir / f"part{i}.vcf")
    outs = {}
    for t in (1, 3):
        out = tmp_path / f"out{t}"
        subprocess.run([sys.executable, "-m", "genome_format_converters.cli",
                        "vcf-to-pseudohaploid", "--input-dir", str(in_dir),
                        "--output-dir", str(out), "--seed", "5", "--threads", str(t)],
                       check=True, capture_output=True)
        outs[t] = {p.name: p.read_bytes() for p in out.iterdir()}
    assert outs[1] == outs[3]
    assert len(outs[3]) == 9
