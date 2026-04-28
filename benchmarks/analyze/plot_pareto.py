#!/usr/bin/env python3
"""Figure 6 — speed vs. memory Pareto scatter.

One point per (task, tool) on a 2D plane: x = wall_s (log scale),
y = peak_rss_mb (log scale). Lower-left = better tradeoff. Points are
coloured by tool and labelled by task to make per-task tools easy to
trace.

Paper utility: where the bar charts say "gfc is competitive on speed
and dominant on memory," the Pareto plot SHOWS it — readers see
plink2's chr22 datapoint sit way up the y-axis (1300 MB) while gfc
sits in the lower left, even though plink2 wins on the x-axis. That
visual is the headline argument for "memory-efficient alternative."

Failed reps are excluded; only successful (exit_code == 0) replicates
contribute.
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

    # Mean per (task, tool) so each tool plots one point per task, not
    # `n_replicates` overlapping dots.
    g = df.groupby(["task", "tool"], as_index=False).agg(
        wall_s=("wall_s", "mean"),
        peak_rss_mb=("peak_rss_mb", "mean"),
    )

    FIG.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="paper")

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.scatterplot(data=g, x="wall_s", y="peak_rss_mb",
                    hue="tool", style="tool", s=140, ax=ax,
                    alpha=0.85, edgecolor="black", linewidth=0.5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("wall time (s, log scale)")
    ax.set_ylabel("peak RSS (MB, log scale)")
    ax.set_title("Speed vs. memory tradeoff — one point per (task, tool); "
                 "lower-left is better")

    # Annotate each point with its task label so the reader can trace
    # which dots are which task without cross-referencing the legend.
    for _, row in g.iterrows():
        ax.annotate(row["task"], (row["wall_s"], row["peak_rss_mb"]),
                    xytext=(4, 4), textcoords="offset points",
                    fontsize=7, alpha=0.7)

    ax.legend(title="tool", loc="best", frameon=False, fontsize=8)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"pareto.{ext}", dpi=300, bbox_inches="tight")
    print(f"[done] wrote {FIG / 'pareto.png'} + .pdf", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
