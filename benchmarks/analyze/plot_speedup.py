#!/usr/bin/env python3
"""Figure 7 — speedup ratio bars (gfc vs each competitor, per task).

For every task and every competitor, plots the ratio
``competitor_wall / gfc_wall``:

  > 1   means gfc is FASTER than the competitor (bar above 1 line)
  = 1   means tied
  < 1   means gfc is SLOWER than the competitor (bar below 1 line)

A horizontal reference line at y=1 makes the win/loss boundary visible
at a glance. Y-axis is log so that 14x slower (T2 vs plink2) and 13x
faster (T4 vs AGAT) read symmetrically.

Paper utility: this is the single chart that says "where does gfc
win, where does it lose, by how much" without requiring the reader to
do mental arithmetic against the wall_time bar chart.

Failed competitor reps are excluded; tasks with no competitor produce
no bars.
"""
from pathlib import Path
import sys

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

RAW = Path(__file__).resolve().parent.parent / "results" / "raw" / "all.tsv"
FIG = Path(__file__).resolve().parent.parent / "results" / "figures"


def main() -> int:
    if not RAW.exists():
        print(f"[error] {RAW} missing.", file=sys.stderr)
        return 1
    df = pd.read_csv(RAW, sep="\t")
    df = df[df["exit_code"] == 0].copy()
    if df.empty:
        print("[error] no successful reps in all.tsv", file=sys.stderr)
        return 1

    # Mean per (task, tool).
    means = df.groupby(["task", "tool"], as_index=False)["wall_s"].mean()

    # Pivot so we can divide each competitor row by gfc.
    wide = means.pivot(index="task", columns="tool", values="wall_s")
    if "gfc" not in wide.columns:
        print("[error] no gfc rows in all.tsv", file=sys.stderr)
        return 1

    gfc_wall = wide["gfc"]
    competitors = wide.drop(columns=["gfc"])

    # speedup = competitor_wall / gfc_wall.  > 1 = gfc faster, < 1 = competitor faster.
    speedup = competitors.divide(gfc_wall, axis=0)

    # Long form for seaborn.
    long = (speedup.reset_index()
                   .melt(id_vars="task", var_name="tool", value_name="speedup_x")
                   .dropna(subset=["speedup_x"]))
    if long.empty:
        print("[error] no competitor data after pivoting", file=sys.stderr)
        return 1

    FIG.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="paper")

    fig, ax = plt.subplots(figsize=(10, 5))
    sns.barplot(data=long, x="task", y="speedup_x", hue="tool", ax=ax)
    ax.axhline(1.0, color="black", linestyle="--", linewidth=0.8)
    ax.set_yscale("log")
    ax.set_ylabel("speedup ratio = competitor_wall / gfc_wall (log scale)")
    ax.set_xlabel("benchmark task")
    ax.set_title("gfc speedup vs. competitors — bars above the dashed line "
                 "mean gfc is faster")
    ax.legend(title="competitor", loc="best", frameon=False, fontsize=8)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"speedup.{ext}", dpi=300, bbox_inches="tight")
    print(f"[done] wrote {FIG / 'speedup.png'} + .pdf", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
