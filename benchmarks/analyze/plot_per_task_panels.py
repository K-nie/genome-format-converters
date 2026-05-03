#!/usr/bin/env python3
# Author: Benjamin Narh-Madey
"""Supplementary figure — small multiples, one panel per task (2x4 grid).

Where the grouped wall-time bar chart compresses 8 tasks into one
busy x-axis, this lays out a 2x4 grid of mini-panels — one per task,
each showing every tool's mean wall_s as bars (left axis, log scale)
and mean peak_rss_mb as dots on a twinned right axis (also log).

Tool colours are pinned to the shared `_style.TOOL_COLORS` mapping,
so "gfc" is the same colour in every panel — readers can scan
across panels without re-decoding the legend each time.

Successful replicates only.
"""
from pathlib import Path
import sys

import pandas as pd
import matplotlib.pyplot as plt

from _style import (
    color_for_tool,
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

    means = df.groupby(["task", "tool"], as_index=False).agg(
        wall_s=("wall_s", "mean"),
        peak_rss_mb=("peak_rss_mb", "mean"),
    )

    tasks = sorted(means["task"].unique())
    n_tasks = len(tasks)
    n_cols = 4
    n_rows = (n_tasks + n_cols - 1) // n_cols

    set_bioinformatics_style()

    # 180 mm double-column wide; tall enough to fit 2 rows of square-ish panels.
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(180.0 / 25.4, 5.5))
    axes = axes.flatten() if n_rows * n_cols > 1 else [axes]

    for i, task in enumerate(tasks):
        ax = axes[i]
        sub = means[means["task"] == task].sort_values("wall_s")
        x = list(range(len(sub)))
        bar_colors = [color_for_tool(t) for t in sub["tool"]]
        ax.bar(x, sub["wall_s"], color=bar_colors, alpha=0.85,
               edgecolor="black", linewidth=0.4)
        ax.set_xticks(x)
        ax.set_xticklabels(sub["tool"], rotation=45, ha="right", fontsize=6)
        ax.set_ylabel("wall (s)")
        ax.set_yscale("log")
        ax.set_title(task)
        ax.tick_params(axis="y", labelsize=6)

        ax2 = ax.twinx()
        ax2.plot(x, sub["peak_rss_mb"], "ko-",
                 markersize=3, linewidth=0.7, alpha=0.75)
        ax2.set_yscale("log")
        ax2.set_ylabel("RSS (MB)", fontsize=7)
        ax2.tick_params(axis="y", labelsize=6)
        ax2.grid(False)

    for j in range(len(tasks), len(axes)):
        axes[j].axis("off")

    fig.suptitle("Per-task detail: bars = mean wall time (left axis); "
                 "dots = mean peak RSS (right axis); both log scale")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save_figure(fig, "per_task_panels")
    print(f"[done] wrote per_task_panels.png + .pdf", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
