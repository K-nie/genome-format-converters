#!/usr/bin/env python3
# Author: Benjamin Narh-Madey
"""Supplementary figure — cumulative compute time per tool, stacked by task.

For each tool, the bar height is the sum of mean wall_s across every
task that tool ran. gfc's bar covers all 8 segments (one per task);
single-purpose competitors cover only 1-2 segments.

The framing in the title is explicit about the methodology: each tool's
bar is the sum of wall time across the tasks it can perform, not a
like-for-like single-task comparison. Reading the figure as "gfc is
slower than plink2" misses the point — plink2's bar covers only T2;
gfc's bar covers T1..T8. The implicit user-cost argument is "to cover
what gfc does in one CLI you would install N tools whose summed
runtime is the sum of their bars."

Failed reps excluded; tasks where a tool did not participate
contribute zero to that tool's stack. Task-segment colours come from
a 8-element pull from the shared Okabe-Ito palette so they match
the per-task colours used elsewhere in the figure suite.
"""
from pathlib import Path
import sys

import pandas as pd
import matplotlib.pyplot as plt

from _style import (
    PALETTE,
    double_col_size,
    save_figure,
    set_bioinformatics_style,
)

RAW = Path(__file__).resolve().parent.parent / "results" / "raw" / "all.tsv"

# Stable T1..T8 -> Okabe-Ito colour. Same task = same colour wherever
# tasks are colour-coded (only this figure today, but pinned for
# forward consistency).
TASK_COLORS = {
    "T1": PALETTE["orange"],
    "T2": PALETTE["sky_blue"],
    "T3": PALETTE["bluish_green"],
    "T4": PALETTE["yellow"],
    "T5": PALETTE["blue"],
    "T6": PALETTE["vermillion"],
    "T7": PALETTE["reddish_purple"],
    "T8": "#999999",   # 8th slot beyond Okabe-Ito's 7-non-black colours.
}


def main() -> int:
    if not RAW.exists():
        print(f"[error] {RAW} missing.", file=sys.stderr)
        return 1
    df = pd.read_csv(RAW, sep="\t")
    df = df[df["exit_code"] == 0].copy()
    if df.empty:
        print("[error] no successful reps in all.tsv", file=sys.stderr)
        return 1

    means = df.groupby(["task", "tool"], as_index=False)["wall_s"].mean()
    wide = means.pivot(index="tool", columns="task", values="wall_s").fillna(0.0)
    wide = wide.assign(_total=wide.sum(axis=1)).sort_values("_total")
    totals = wide.pop("_total")

    set_bioinformatics_style()

    fig, ax = plt.subplots(figsize=double_col_size(height_in=4.5))
    x_positions = list(range(len(wide.index)))
    bottom = pd.Series(0.0, index=wide.index)
    for task in wide.columns:
        heights = wide[task]
        ax.bar(x_positions, heights.values, bottom=bottom.values,
               label=task, color=TASK_COLORS.get(task, "#444444"),
               edgecolor="white", linewidth=0.4)
        bottom = bottom + heights
    ax.set_xticks(x_positions)
    ax.set_xticklabels(wide.index, rotation=45, ha="right")

    for x, (tool, total) in zip(x_positions, totals.items()):
        ax.text(x, float(total), f"{float(total):.0f}s",
                ha="center", va="bottom", fontsize=6, alpha=0.85)

    ax.set_ylabel("cumulative wall time across tasks (s)")
    ax.set_xlabel("tool")
    ax.set_title("Total compute cost per tool, stacked by task — each tool's bar "
                 "is the sum of wall time across the tasks it can perform; "
                 "gfc covers all 8")
    ax.legend(title="task", loc="upper left",
              bbox_to_anchor=(1.02, 1.0), frameon=False, ncol=1)
    fig.tight_layout()
    save_figure(fig, "compute_cost")
    print(f"[done] wrote compute_cost.png + .pdf", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
