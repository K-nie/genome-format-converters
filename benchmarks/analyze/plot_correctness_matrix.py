#!/usr/bin/env python3
"""Figure 5 — correctness heatmap, task x tool.

One cell per (task, tool); colour encodes whether check_correctness.py
flagged that tool's output as matching the per-task reference.

Cell encoding:
  1.0 = correct (matches reference)
  0.0 = incorrect (output diverged from reference)
  0.5 = couldn't check (no reference output, or comparator skipped)
  NaN = tool didn't produce a row for this task (cell stays blank)

Paper utility: a single heatmap is the cleanest way to claim
"gfc agrees with the canonical reference on tasks T1..T8" without
forcing the reader to read a 40-row TSV.
"""
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

RAW = Path(__file__).resolve().parent.parent / "results" / "raw" / "all.tsv"
FIG = Path(__file__).resolve().parent.parent / "results" / "figures"


def _correct_to_score(s: str) -> float:
    """Map the correct column's string to a numeric score for the heatmap."""
    if s == "1":
        return 1.0
    if s == "0":
        return 0.0
    # "" / "skip" / anything else
    return 0.5


def main() -> int:
    if not RAW.exists():
        print(f"[error] {RAW} missing.", file=sys.stderr)
        return 1
    df = pd.read_csv(RAW, sep="\t")
    if "correct" not in df.columns:
        print("[error] no `correct` column in all.tsv; "
              "did check_correctness.py run?", file=sys.stderr)
        return 1
    df["correct"] = df["correct"].fillna("").astype(str)
    # Use rep1 as the canonical correctness verdict (deterministic tools
    # have the same flag for every rep; check_correctness.py only writes
    # rep1 then propagates).
    df1 = df[df["replicate"] == 1].copy()
    df1["score"] = df1["correct"].apply(_correct_to_score)

    matrix = df1.pivot_table(index="task", columns="tool", values="score",
                             aggfunc="first")

    FIG.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="paper")

    fig, ax = plt.subplots(figsize=(max(8, 0.7 * matrix.shape[1] + 4),
                                    max(4, 0.5 * matrix.shape[0] + 2)))
    cmap = sns.color_palette(["#d62728", "#bdbdbd", "#2ca02c"], n_colors=3)
    sns.heatmap(matrix, ax=ax, cmap=cmap, vmin=0.0, vmax=1.0,
                cbar_kws={"ticks": [0.0, 0.5, 1.0],
                          "label": "0 = mismatch | 0.5 = unchecked | 1 = match"},
                linewidths=0.5, linecolor="white",
                annot=True, fmt=".0%", annot_kws={"size": 9})
    ax.set_xlabel("tool")
    ax.set_ylabel("benchmark task")
    ax.set_title("Output correctness — gfc and competitors vs. per-task reference")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"correctness_matrix.{ext}", dpi=300,
                    bbox_inches="tight")
    print(f"[done] wrote {FIG / 'correctness_matrix.png'} + .pdf",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
