#!/usr/bin/env python3
"""Concatenate every `results/raw/T*.tsv` into `results/raw/all.tsv`.

Run after ``run_bench.sh``. Figures (`plot_wall_time.py`,
`plot_memory.py`) read from `all.tsv`.
"""
from pathlib import Path
import sys

RAW = Path(__file__).resolve().parent.parent / "results" / "raw"

def main() -> int:
    files = sorted(RAW.glob("T*.tsv"))
    if not files:
        print(f"[error] no T*.tsv under {RAW}. Did you run run_bench.sh?",
              file=sys.stderr)
        return 1
    out = RAW / "all.tsv"
    header = None
    rows = []
    for f in files:
        lines = f.read_text().splitlines()
        if not lines:
            continue
        if header is None:
            header = lines[0]
        rows.extend(l for l in lines[1:] if l.strip())
    if header is None:
        print(f"[warn] every T*.tsv under {RAW} was empty; nothing to merge.",
              file=sys.stderr)
        return 0
    with open(out, "w") as fh:
        fh.write(header + "\n")
        for r in rows:
            fh.write(r + "\n")
    print(f"[done] merged {len(files)} files → {out} "
          f"({len(rows)} rows).", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
