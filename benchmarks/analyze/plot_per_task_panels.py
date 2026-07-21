#!/usr/bin/env python3
"""Supplementary figure — small multiples, one panel per task.

Where the grouped wall-time bar chart compresses 8 tasks into one
busy x-axis, this lays out a 2x4 grid of mini-panels — one per task,
each showing every tool's wall_s and peak_rss_mb side by side as
twinned bars (left axis: wall, right axis: RSS). Easier to read for
single-task comparisons.

Paper utility: supplementary "S2 — per-task detail" — when a reviewer
asks "what's actually going on with T2?" this is the figure to point
at without forcing them to filter the headline bar chart visually.

Successful reps only.
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

    means = df.groupby(["task", "tool"], as_index=False).agg(
        wall_s=("wall_s", "mean"),
        peak_rss_mb=("peak_rss_mb", "mean"),
    )

    tasks = sorted(means["task"].unique())
    n_tasks = len(tasks)
    n_cols = 4
    n_rows = (n_tasks + n_cols - 1) // n_cols

    FIG.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="paper")

    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(3.2 * n_cols, 3.2 * n_rows))
    axes = axes.flatten() if n_rows * n_cols > 1 else [axes]

    palette = sns.color_palette("husl", n_colors=12)
    tool_color = {t: palette[i % len(palette)]
                  for i, t in enumerate(sorted(means["tool"].unique()))}

    for i, task in enumerate(tasks):
        ax = axes[i]
        sub = means[means["task"] == task].sort_values("wall_s")
        x = range(len(sub))
        ax.bar(x, sub["wall_s"], color=[tool_color[t] for t in sub["tool"]],
               alpha=0.85)
        ax.set_xticks(list(x))
        ax.set_xticklabels(sub["tool"], rotation=45, ha="right", fontsize=7)
        ax.set_ylabel("wall (s)", fontsize=8)
        ax.set_yscale("log")
        ax.set_title(task, fontsize=10)
        ax.tick_params(axis="y", labelsize=7)

        # Twinned axis for RSS.
        ax2 = ax.twinx()
        ax2.plot(x, sub["peak_rss_mb"], "ko-",
                 markersize=4, linewidth=0.8, alpha=0.7,
                 label="peak RSS (MB)")
        ax2.set_yscale("log")
        ax2.set_ylabel("RSS (MB)", fontsize=8)
        ax2.tick_params(axis="y", labelsize=7)
        ax2.grid(False)

    # Hide any unused panels.
    for j in range(len(tasks), len(axes)):
        axes[j].axis("off")

    fig.suptitle("Per-task detail — bars: wall time (left axis); "
                 "dots: peak RSS (right axis)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"per_task_panels.{ext}", dpi=300,
                    bbox_inches="tight")
    print(f"[done] wrote {FIG / 'per_task_panels.png'} + .pdf",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
