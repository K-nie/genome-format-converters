"""Pytest fixtures shared across converter tests."""

from pathlib import Path
import shutil

import pytest

FIXTURE_DIR = Path(__file__).parent / "test_data"


@pytest.fixture
def fixtures() -> Path:
    """Absolute path to the bundled test-data directory."""
    return FIXTURE_DIR


@pytest.fixture
def tmp_input_dir(tmp_path: Path, request) -> Path:
    """Stage a copy of a requested fixture file inside a per-test input
    directory. Parametrise via `@pytest.mark.parametrize('tmp_input_dir',
    ['tiny.vcf', 'tiny.gff3'], indirect=True)` or set an attribute via
    `request.param`.
    """
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    target = getattr(request, "param", None)
    if target:
        src = FIXTURE_DIR / target
        shutil.copy(src, in_dir / src.name)
    return in_dir


@pytest.fixture
def tmp_output_dir(tmp_path: Path) -> Path:
    out = tmp_path / "out"
    return out  # intentionally NOT created — let the converter do it
