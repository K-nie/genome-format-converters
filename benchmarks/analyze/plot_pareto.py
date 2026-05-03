#!/usr/bin/env python3
# Author: Benjamin Narh-Madey
"""Figure 6 — speed vs. memory Pareto scatter.

One point per (task, tool) on a 2D log-log plane:
  x = mean wall_s across replicates (lower is faster)
  y = mean peak_rss_mb across replicates (lower uses less memory)

Lower-left = better trade-off. Each point is coloured by tool via the
shared `_style.TOOL_COLORS` mapping (same colour for the same tool in
every figure) and labelled with its task code so dots can be traced
without cross-referencing the legend.

Paper utility: where the bar charts say "gfc is competitive on speed
and dominant on memory," the Pareto plot SHOWS it — the chr22 plink2
+ convertf points sit ~22x up the y-axis at ~1300 MB while every gfc
point clusters at ~60 MB across all eight tasks. That horizontal
band IS the memory-efficiency story.

Failed replicates are excluded (exit_code == 0 only).
"""
from pathlib import Path
import sys

import pandas as pd
import matplotlib.pyplot as plt

from _style import (
    color_for_tool,
    double_col_size,
    save_figure,
    set_bioinformatics_style,
)

RAW = Path(__file__).resolve().parent.parent / "results" / "raw" / "all.tsv"


def main() -> int:
    if not RAW.exists():
        print(f"[error] {RAW} missing.", file=sys.stderr)
        return 1
    df = pd.read_csv(RAW, sep="\t")
    df = df[df["exit_code"] == 0].copy()
    if df.empty:
        print("[error] no successful reps in all.tsv", file=sys.stderr)
        return 1

    g = df.groupby(["task", "tool"], as_index=False).agg(
        wall_s=("wall_s", "mean"),
        peak_rss_mb=("peak_rss_mb", "mean"),
    )

    set_bioinformatics_style()

    fig, ax = plt.subplots(figsize=double_col_size(height_in=5.0))

    tools = ["gfc"] + sorted(t for t in g["tool"].unique() if t != "gfc")
    for tool in tools:
        sub = g[g["tool"] == tool]
        ax.scatter(sub["wall_s"], sub["peak_rss_mb"],
                   c=color_for_tool(tool), s=80, label=tool,
                   alpha=0.85, edgecolor="black", linewidth=0.5)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("wall time (s, log scale) — lower is faster")
    ax.set_ylabel("peak resident set size (MB, log scale) — lower uses less memory")
    ax.set_title("Speed vs. memory trade-off across all (task, tool) pairs; "
                 "lower-left is better")

    # Per-point task label — small offset so dot + label do not overlap.
    for _, row in g.iterrows():
        ax.annotate(row["task"], (row["wall_s"], row["peak_rss_mb"]),
                    xytext=(4, 4), textcoords="offset points",
                    fontsize=6, alpha=0.7)

    ax.legend(title="tool", loc="upper left",
              bbox_to_anchor=(1.02, 1.0), frameon=False)
    fig.tight_layout()
    save_figure(fig, "pareto")
    print(f"[done] wrote pareto.png + .pdf", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
