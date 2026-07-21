#!/usr/bin/env python3
"""Supplementary figure — cumulative compute time per tool, stacked by task.

For each tool, the bar height is the sum of mean wall_s across every
task that tool ran. gfc's bar covers all 8 segments (one per task);
single-purpose competitors cover only 1-2. Visualises the
"one tool, no overhead, handles every task" story.

Paper utility: supplementary "S3 — total cost of running the bench" or
"workflow comparison" — the implicit argument is that even where a
single competitor wins on one task, the user pays the cost of
installing, configuring, and learning N tools to cover what gfc does
alone.

Failed reps excluded; tasks where a tool didn't participate contribute
zero to that tool's stack.
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

    # Mean wall per (task, tool); pivot so tasks are stacked on top of
    # each tool's column.
    means = df.groupby(["task", "tool"], as_index=False)["wall_s"].mean()
    wide = means.pivot(index="tool", columns="task", values="wall_s").fillna(0.0)

    # Sort tools by total cumulative wall so the figure reads
    # "biggest stack on the right."
    wide = wide.assign(_total=wide.sum(axis=1)).sort_values("_total")
    totals = wide.pop("_total")

    FIG.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="paper")
    palette = sns.color_palette("tab10", n_colors=len(wide.columns))

    fig, ax = plt.subplots(figsize=(10, 5))
    # Use integer x-positions explicitly so the per-bar text annotation
    # below has a numeric x to anchor against. seaborn / matplotlib's
    # categorical bar handling is fine for the bars themselves but
    # `ax.text(tool, total, ...)` with a string x crashes inside
    # `convert_xunits` on recent matplotlib.
    x_positions = list(range(len(wide.index)))
    bottom = pd.Series(0.0, index=wide.index)
    for i, task in enumerate(wide.columns):
        heights = wide[task]
        ax.bar(x_positions, heights.values, bottom=bottom.values,
               label=task, color=palette[i], edgecolor="white",
               linewidth=0.4)
        bottom = bottom + heights
    ax.set_xticks(x_positions)
    ax.set_xticklabels(wide.index)

    # Annotate each bar's total at the top.
    for x, (tool, total) in zip(x_positions, totals.items()):
        ax.text(x, float(total), f"{float(total):.0f}s",
                ha="center", va="bottom", fontsize=8, alpha=0.8)

    ax.set_ylabel("cumulative wall time across tasks (s)")
    ax.set_xlabel("tool")
    ax.set_title("Total compute cost per tool, stacked by task — "
                 "gfc's stack covers every task; single-purpose tools cover one")
    ax.legend(title="task", loc="best", frameon=False, fontsize=8,
              ncol=2)
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"compute_cost.{ext}", dpi=300, bbox_inches="tight")
    print(f"[done] wrote {FIG / 'compute_cost.png'} + .pdf", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
