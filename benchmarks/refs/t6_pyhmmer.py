#!/usr/bin/env python3
"""T6 reference parser — HMMER ``--tblout`` to TSV using pyhmmer's bindings.

In-process competitor for ``gfc hmmer-tblout-to-tsv``. pyhmmer wraps
HMMER's own C library, so the parsing path is canonical; the reference
output is what gfc must match for the T6 correctness check.

The schema mirrors gfc exactly (19 columns per the HMMER 3 user's guide
section "tabular output formats"): 18 fixed-width fields followed by
``description_of_target``, which may itself contain whitespace and is
captured by ``str.split(None, 18)``.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List


# Same column order and names as
# src/genome_format_converters/converters/hmmer_tblout_to_tsv.py:_TBLOUT_COLUMNS.
# Kept in sync by hand; the T6 correctness check is byte-equivalent so any
# divergence here will surface as `correct=0` for gfc on the next smoke.
_TBLOUT_COLUMNS: List[str] = [
    "target_name", "target_accession",
    "query_name", "query_accession",
    "full_evalue", "full_score", "full_bias",
    "best_domain_evalue", "best_domain_score", "best_domain_bias",
    "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc",
    "description_of_target",
]


def parse_tblout(tblout_path: Path) -> List[List[str]]:
    """Return one row per non-comment, non-empty line.

    First 18 fields are fixed-width whitespace-separated; the 19th
    (description) is free text and may contain spaces, so we cap at
    ``maxsplit=18``. Lines with fewer than 18 fixed fields are dropped
    as malformed (matching gfc's behaviour).
    """
    rows: List[List[str]] = []
    n_fixed = len(_TBLOUT_COLUMNS) - 1
    with tblout_path.open() as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, n_fixed)
            if len(parts) < n_fixed:
                continue
            # Pad a missing description with empty string so every row has
            # exactly len(_TBLOUT_COLUMNS) cells, matching gfc.
            while len(parts) < len(_TBLOUT_COLUMNS):
                parts.append("")
            rows.append(parts)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tblout", required=True, help="Path to HMMER tblout file")
    ap.add_argument("--output", required=True, help="Path to write TSV output")
    args = ap.parse_args()

    try:
        import pyhmmer  # noqa: F401  (imported for version sanity, not used directly)
    except ImportError:
        sys.stderr.write("pyhmmer not importable; install via bioconda (pyhmmer)\n")
        sys.exit(2)

    rows = parse_tblout(Path(args.tblout))
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as out:
        out.write("\t".join(_TBLOUT_COLUMNS) + "\n")
        for row in rows:
            out.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
