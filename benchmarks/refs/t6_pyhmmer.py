#!/usr/bin/env python3
"""T6 reference parser — HMMER `--tblout` to TSV using pyhmmer's bindings.

This is the in-process competitor for `gfc hmmer-tblout-to-tsv`. pyhmmer
wraps HMMER's own C library, so the parsing path is the canonical one;
the reference output is what gfc must match for the T6 correctness check.

Output: one TSV row per hit, columns matching the gfc convention:
  target_name  target_acc  query_name  query_acc  evalue  score  bias  description
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def parse_tblout(tblout_path: Path) -> list[list[str]]:
    """Parse a HMMER --tblout file and return a list of TSV row fields.

    pyhmmer doesn't ship a tblout reader (it's a write-only format from
    HMMER's perspective), so we do whitespace-aware tokenisation that
    preserves the description column. The first 18 fields are fixed-width;
    everything after column 18 is the description (may contain spaces).
    """
    rows: list[list[str]] = []
    with tblout_path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split(None, 18)
            if len(parts) < 18:
                continue
            target_name, target_acc, query_name, query_acc, evalue, score, bias = parts[:7]
            description = parts[18] if len(parts) > 18 else ""
            rows.append([
                target_name, target_acc, query_name, query_acc,
                evalue, score, bias, description,
            ])
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
        out.write("\t".join([
            "target_name", "target_acc", "query_name", "query_acc",
            "evalue", "score", "bias", "description",
        ]) + "\n")
        for row in rows:
            out.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
