"""Genome Format Converters package.

A uniform command-line toolkit for converting common bioinformatics file
formats. Exposes the `gfc` CLI (see genome_format_converters.cli:main).
"""

from importlib.metadata import PackageNotFoundError, version as _pkg_version

try:
    __version__ = _pkg_version("genome-format-converters")
except PackageNotFoundError:
    # Package isn't installed (e.g. running from a source checkout without
    # `pip install -e .`). Fall back to a local sentinel.
    __version__ = "0.0.0+unknown"
