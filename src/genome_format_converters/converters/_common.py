"""Shared helpers for converter modules.

Centralises:
    * stderr logging (so stdout stays clean for pipeable output)
    * compound-suffix stripping (`sample.vcf.gz` → `sample`)
    * two-column TSV map loading (used by pop-map, sex-map, chrom-map)
    * output-directory setup with optional overwrite guard
    * fatal-error helpers so converters fail loudly instead of silently
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, Iterable, Optional


# ---------- logging ---------------------------------------------------------

# Controlled via the CLI (--quiet, --verbose) by pointing these at no-op
# callables or swapping in tee'd variants.
_VERBOSE = False
_QUIET = False


def set_verbosity(verbose: bool = False, quiet: bool = False) -> None:
    global _VERBOSE, _QUIET
    _VERBOSE = verbose
    _QUIET = quiet


def log_info(msg: str) -> None:
    """Progress line. Suppressed under --quiet. Goes to stderr."""
    if not _QUIET:
        print(msg, file=sys.stderr)


def log_warn(msg: str) -> None:
    """Warning. Always shown unless --quiet. Goes to stderr."""
    if not _QUIET:
        print(f"[warn] {msg}", file=sys.stderr)


def log_debug(msg: str) -> None:
    """Debug-only line. Shown only under --verbose."""
    if _VERBOSE and not _QUIET:
        print(f"[debug] {msg}", file=sys.stderr)


def die(msg: str, exit_code: int = 2) -> None:
    """Fatal error — print to stderr and exit non-zero."""
    print(f"[error] {msg}", file=sys.stderr)
    sys.exit(exit_code)


def require_pysam():
    """Lazy-import pysam.

    htslib does not build cleanly on Windows, so upstream `pysam` does not
    publish Windows wheels on PyPI. We therefore leave pysam as a
    platform-conditional dependency (installed on Linux/macOS, omitted on
    Windows) and defer its import into the specific subcommands that need
    it. Windows users invoking a VCF/BAM/BCF subcommand get a readable
    error pointing them at WSL or conda-forge rather than a cryptic
    ImportError during module load.
    """
    try:
        import pysam  # noqa: F401 — imported for side-effects of availability check
        return pysam
    except ImportError as exc:
        die(
            "This subcommand needs the `pysam` package, which isn't available "
            "on Windows via pip. Options:\n"
            "  1. Use Windows Subsystem for Linux (WSL) and `pip install "
            "genome-format-converters` inside the WSL shell.\n"
            "  2. Install pysam from conda-forge: "
            "`conda install -c conda-forge pysam`.\n"
            "  3. If you're on Linux/macOS and seeing this, reinstall "
            "genome-format-converters to pull in pysam:\n"
            "     pip install --force-reinstall genome-format-converters\n"
            f"Underlying error: {exc}"
        )


# ---------- filename / extension handling -----------------------------------

def strip_compound_suffix(path: Path, exts: Iterable[str]) -> str:
    """Return the path's stem with any recognised compound suffix stripped.

    `Path.with_suffix` only handles the final suffix, which silently produces
    wrong output names for `.vcf.gz` / `.tar.gz` / etc. This function checks
    the full name against an explicit list of extensions and strips the first
    match.
    """
    name = path.name
    # Sort longest first so `.vcf.gz` matches before `.gz`.
    for ext in sorted(exts, key=len, reverse=True):
        if name.endswith(ext):
            return name[: -len(ext)]
    return path.stem


# ---------- two-column map files --------------------------------------------

def load_two_col_map(path: Optional[str]) -> Dict[str, str]:
    """Load a whitespace-delimited two-column TSV into a dict.

    Blank lines and `#`-prefixed comment lines are ignored. Rows with fewer
    than two columns are skipped silently. Returns an empty dict when `path`
    is None, which lets callers treat the map as optional.
    """
    if path is None:
        return {}
    mapping: Dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            mapping[parts[0]] = parts[1]
    return mapping


# ---------- output directory + overwrite guard ------------------------------

def prepare_output_dir(output_dir: str, force: bool = False) -> Path:
    """Create the output directory if needed; optionally refuse to clobber.

    With `force=False` (the default when called from the CLI without
    `--force`), an existing non-empty directory raises SystemExit(2). This
    prevents silent overwrites of previous results.
    """
    out_path = Path(output_dir)
    if out_path.exists() and any(out_path.iterdir()) and not force:
        die(
            f"Output directory {out_path} exists and is non-empty. "
            "Pass --force to overwrite, or choose a fresh --output-dir."
        )
    out_path.mkdir(parents=True, exist_ok=True)
    return out_path


# ---------- input globbing --------------------------------------------------

def iter_input_files(input_dir: Path, exts: Iterable[str],
                     pattern: Optional[str] = None) -> Iterable[Path]:
    """Yield files in `input_dir` matching any of `exts`, optionally
    restricted by a user-provided glob `pattern`. Case-insensitive on the
    extension match."""
    lowered = {e.lower() for e in exts}
    seen = set()
    if pattern:
        candidates = sorted(input_dir.glob(pattern))
    else:
        candidates = sorted(input_dir.iterdir())
    for path in candidates:
        if not path.is_file():
            continue
        name_lower = path.name.lower()
        if any(name_lower.endswith(e) for e in lowered):
            if path not in seen:
                seen.add(path)
                yield path
