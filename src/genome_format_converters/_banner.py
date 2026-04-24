"""ASCII banner for the ``gfc`` CLI.

Shown on bare ``gfc`` invocation and at the top of ``--help``. Suppressed
when stdout is not a terminal (so redirection, pipes, and CI logs stay
clean), and suppressed unconditionally under ``--quiet``.
"""

from __future__ import annotations

import sys

from genome_format_converters import __version__


# Monospace block-letter G F C. 6 lines, 36 characters wide.
_BLOCK = r"""   ____    _____    ____
  / ___|  |  ___|  / ___|
 | |  _   | |_    | |
 | |_| |  |  _|   | |___
  \____|  |_|      \____|"""


def get_banner(version: str = __version__) -> str:
    """Return the full banner string, including tagline."""
    return (
        f"{_BLOCK}\n\n"
        f"   genome-format-converters v{version}\n"
        f"   Uniform CLI for 30 bioinformatics file-format conversions.\n"
        f"   Benjamin Narh-Madey · Hittinger Lab · UW-Madison\n"
        f"   https://github.com/K-nie/genome-format-converters\n"
    )


def print_banner(force: bool = False) -> None:
    """Print the banner to stderr if stderr is a terminal (or `force=True`)."""
    if force or sys.stderr.isatty():
        print(get_banner(), file=sys.stderr)
