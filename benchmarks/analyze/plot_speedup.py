#!/usr/bin/env python3
# Author: Benjamin Narh-Madey
"""Figure 7 — wall-time ratio (gfc vs each competitor, per task).

For every task and every competitor, plots the ratio
  wall_time_ratio = competitor_wall / gfc_wall

  > 1  : gfc is faster than the competitor (bar above 1 line)
  = 1  : tied
  < 1  : gfc is slower than the competitor (bar below 1 line)

A horizontal reference line at y = 1 makes the parity boundary
visible. Y-axis is log so symmetric speed/slowdown deltas read
symmetrically; e.g. T2 plink2 (0.07x = gfc 14x slower) and T5
AGAT (13.3x = gfc 13x faster) are equidistant from the parity line.

The earlier "speedup" framing in the title overstated the result:
of the 12 competitor pairs benchmarked here, only 3 are clean wins
(audit 2026-05-03). The honest framing is "wall-time ratio" and
"which side is faster."

Statistical annotation (Phase 3, 2026-05-03): bars where the paired
Wilcoxon test on log10(wall_s) survives Benjamini-Hochberg FDR
correction at q = 0.05 (across all task x competitor x metric
tests) get a trailing asterisk. p_adj_bh values are pulled from
``benchmarks/results/figures/stats_pairs.tsv``; if the stats file
is absent the plot still renders without annotations and warns.

Correctness annotation: tasks where any tool's output diverged
structurally from the per-task reference (`correct == "0"`) get a
trailing asterisk on the x-axis label; the title points to the
correctness matrix figure for the per-task footnote.

Failed competitor reps excluded; tasks with no competitor produce
no bars.
"""
from pathlib import Path
import sys

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from _style import (
    color_for_tool,
    double_col_size,
    save_figure,
    set_bioinformatics_style,
)

RAW = Path(__file__).resolve().parent.parent / "results" / "raw" / "all.tsv"
STATS = Path(__file__).resolve().parent.parent / "results" / "figures" / "stats_pairs.tsv"


def main() -> int:
    if not RAW.exists():
        print(f"[error] {RAW} missing.", file=sys.stderr)
        return 1
    df = pd.read_csv(RAW, sep="\t")
    df = df[df["exit_code"] == 0].copy()
    if df.empty:
        print("[error] no successful reps in all.tsv", file=sys.stderr)
        return 1

    if "correct" in df.columns:
        df["correct"] = df["correct"].fillna("").astype(str)
        diverged_tasks = sorted(set(df.loc[df["correct"] == "0", "task"]))
    else:
        diverged_tasks = []

    means = df.groupby(["task", "tool"], as_index=False)["wall_s"].mean()
    wide = means.pivot(index="task", columns="tool", values="wall_s")
    if "gfc" not in wide.columns:
        print("[error] no gfc rows in all.tsv", file=sys.stderr)
        return 1

    gfc_wall = wide["gfc"]
    competitors = wide.drop(columns=["gfc"])
    ratio = competitors.divide(gfc_wall, axis=0)

    long = (ratio.reset_index()
                 .melt(id_vars="task", var_name="tool",
                       value_name="ratio_x")
                 .dropna(subset=["ratio_x"]))
    if long.empty:
        print("[error] no competitor data after pivoting", file=sys.stderr)
        return 1

    # Pull per-(task, competitor) BH-adjusted p-values for the wall_s
    # metric; build a lookup keyed on (task, tool) -> bool significant.
    sig_lookup: dict[tuple[str, str], bool] = {}
    n_tests = 0
    if STATS.exists():
        stats_df = pd.read_csv(STATS, sep="\t")
        wall_stats = stats_df[stats_df["metric"] == "wall_s"].copy()
        n_tests = len(stats_df)
        for _, row in wall_stats.iterrows():
            sig_lookup[(row["task"], row["competitor"])] = bool(
                pd.notna(row["p_adj_bh"]) and row["p_adj_bh"] < 0.05)
    else:
        print(f"[warn] {STATS} missing — no significance asterisks rendered.",
              file=sys.stderr)

    set_bioinformatics_style()

    competitor_order = sorted(long["tool"].unique())
    palette = {t: color_for_tool(t) for t in competitor_order}

    fig, ax = plt.subplots(figsize=double_col_size(height_in=4.5))
    sns.barplot(data=long, x="task", y="ratio_x", hue="tool",
                hue_order=competitor_order, palette=palette, ax=ax)
    ax.axhline(1.0, color="black", linestyle="--", linewidth=0.8)
    ax.set_yscale("log")
    ax.set_ylabel("wall-time ratio (competitor / gfc, log scale)\n"
                  "> 1 : gfc faster   |   < 1 : competitor faster")
    ax.set_xlabel("benchmark task")

    xticks = list(range(len(ax.get_xticklabels())))
    ax.set_xticks(xticks)
    new_xticklabels = []
    for txt in ax.get_xticklabels():
        task = txt.get_text()
        new_xticklabels.append(f"{task}*" if task in diverged_tasks else task)
    ax.set_xticklabels(new_xticklabels)

    # Annotate significant bars with an asterisk above the bar.
    # seaborn's modern barplot returns dense containers (only the
    # bars that have data for that hue, not one bar per x position),
    # so we can't index by task position. Map each bar back to its
    # task via the x-tick at the bar's x-centre instead.
    if sig_lookup:
        x_tasks = [t.get_text().rstrip("*") for t in ax.get_xticklabels()]
        x_centres = ax.get_xticks()
        # Tolerance: half a bar group's worth of x.
        bar_w = 0.0
        for cont in ax.containers:
            for bar in cont:
                bar_w = max(bar_w, bar.get_width())
        tol = max(bar_w, 0.5)
        for h, hue in enumerate(competitor_order):
            if h >= len(ax.containers):
                continue
            for bar in ax.containers[h]:
                xc = bar.get_x() + bar.get_width() / 2
                # Find nearest x-tick.
                idx = int(min(range(len(x_centres)),
                              key=lambda j: abs(x_centres[j] - xc)))
                if idx >= len(x_tasks):
                    continue
                if abs(x_centres[idx] - xc) > tol:
                    continue
                task = x_tasks[idx]
                if not sig_lookup.get((task, hue), False):
                    continue
                height = bar.get_height()
                if not (height and height > 0):
                    continue
                ax.text(xc, height * 1.06, "*",
                        ha="center", va="bottom",
                        fontsize=11, color="black",
                        clip_on=False)

    title = ("Wall-time ratio (competitor / gfc) per benchmark task; "
             "dashed line at 1.0 is parity")
    if sig_lookup:
        title += (f"\n* p_BH < 0.05 (paired Wilcoxon on log10 wall_s; "
                  f"BH-FDR across {n_tests} tests)")
    if diverged_tasks:
        title += (f"\n** = at least one tool's output flagged as "
                  f"structural-divergence (see correctness_matrix figure)")
    ax.set_title(title)

    ax.legend(title="competitor", loc="upper left",
              bbox_to_anchor=(1.02, 1.0), frameon=False)
    fig.tight_layout()
    save_figure(fig, "speedup")
    print(f"[done] wrote speedup.png + .pdf", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
