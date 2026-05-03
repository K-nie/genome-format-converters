#!/usr/bin/env python3
# Author: Benjamin Narh-Madey
"""Figure 4 — peak resident memory per benchmark task, grouped by tool.

Bar height is mean peak RSS (MB) across N replicates; error bars
encode the min/max range so reviewers see the spread directly rather
than inferring it from sd. Y-axis is log because the chr22 plink2 +
convertf points sit ~22x above gfc.

Phase 3 stats annotation (2026-05-03): the headline RAM ratio (the
single largest competitor / gfc median ratio across all (task,
competitor) pairs) is annotated on the figure with its bootstrap
95 % CI from ``stats_pairs.tsv``. This stops the bare "22x" claim
from dropping into the paper without uncertainty: the actual
number is 22.09x [95% CI 21.70-22.46] (or whatever shifts in a
re-run). If ``stats_pairs.tsv`` is absent the annotation is
skipped silently.

Correctness annotation matches `plot_wall_time.py`: tasks where any
tool's output diverged structurally from the per-task reference get
a trailing asterisk and the title points to the correctness matrix
for context (no auto-generated "X disagreed" inline text — that
framing read as an indictment of tools whose divergences are
documented design choices).
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


def _headline_ram_pair(stats_path: Path) -> dict | None:
    """Return the (task, competitor, median, lo, hi) tuple for the
    single largest peak_rss_mb median ratio in stats_pairs.tsv.
    None if the file is absent or empty."""
    if not stats_path.exists():
        return None
    sdf = pd.read_csv(stats_path, sep="\t")
    rss = sdf[sdf["metric"] == "peak_rss_mb"].copy()
    rss = rss.dropna(subset=["median_ratio"])
    if rss.empty:
        return None
    top = rss.loc[rss["median_ratio"].idxmax()]
    return {
        "task": str(top["task"]),
        "competitor": str(top["competitor"]),
        "median_ratio": float(top["median_ratio"]),
        "ci_lo": float(top["ci_lo_95"]),
        "ci_hi": float(top["ci_hi_95"]),
    }


def main() -> int:
    if not RAW.exists():
        print(f"[error] {RAW} missing.", file=sys.stderr)
        return 1
    df = pd.read_csv(RAW, sep="\t")
    n_total = len(df)
    failed = df[df["exit_code"] != 0]
    df_ok = df[df["exit_code"] == 0].copy()
    if df_ok.empty:
        print("[error] no successful replicates in all.tsv", file=sys.stderr)
        return 1

    if "correct" in df_ok.columns:
        df_ok["correct"] = df_ok["correct"].fillna("").astype(str)
        diverged_tasks = sorted(set(df_ok.loc[df_ok["correct"] == "0", "task"]))
    else:
        diverged_tasks = []

    set_bioinformatics_style()

    tool_order = ["gfc"] + sorted(t for t in df_ok["tool"].unique() if t != "gfc")
    palette = {t: color_for_tool(t) for t in tool_order}

    fig, ax = plt.subplots(figsize=double_col_size(height_in=4.5))
    sns.barplot(data=df_ok, x="task", y="peak_rss_mb", hue="tool",
                hue_order=tool_order, palette=palette,
                errorbar=lambda v: (v.min(), v.max()),
                ax=ax)
    ax.set_yscale("log")
    ax.set_ylabel("peak resident set size (MB, log scale)\n"
                  "max across N replicates per (task, tool)")
    ax.set_xlabel("benchmark task")

    xticks = list(range(len(ax.get_xticklabels())))
    ax.set_xticks(xticks)
    new_xticklabels = []
    for txt in ax.get_xticklabels():
        task = txt.get_text()
        new_xticklabels.append(f"{task}*" if task in diverged_tasks else task)
    ax.set_xticklabels(new_xticklabels)

    # Annotate the headline RAM ratio with bootstrap 95 % CI (Phase 3).
    headline = _headline_ram_pair(STATS)
    if headline is not None:
        x_tasks = [t.get_text().rstrip("*") for t in ax.get_xticklabels()]
        # Find the bar for (headline['task'], headline['competitor']).
        if headline["competitor"] in tool_order:
            hue_idx = tool_order.index(headline["competitor"])
            if (hue_idx < len(ax.containers)
                    and headline["task"] in x_tasks):
                container = ax.containers[hue_idx]
                bar_idx = x_tasks.index(headline["task"])
                if bar_idx < len(container):
                    bar = container[bar_idx]
                    height = bar.get_height()
                    if height and height > 0:
                        annot = (f"{headline['median_ratio']:.1f}x vs gfc\n"
                                 f"[95% CI {headline['ci_lo']:.1f}-"
                                 f"{headline['ci_hi']:.1f}]")
                        ax.annotate(
                            annot,
                            xy=(bar.get_x() + bar.get_width() / 2, height),
                            xytext=(0, 14),
                            textcoords="offset points",
                            ha="center", va="bottom",
                            fontsize=8, color="black",
                            arrowprops=dict(arrowstyle="-",
                                            color="black", lw=0.6),
                        )

    n_reps = df_ok.groupby(["task", "tool"]).size().min()
    title = (f"Peak resident memory per benchmark task "
             f"(bar = mean, whiskers = min..max, n={n_reps} replicates)")
    if headline is not None:
        title += (f"\nHeadline RAM ratio: {headline['competitor']} on "
                  f"{headline['task']} = "
                  f"{headline['median_ratio']:.2f}x gfc "
                  f"[95% CI {headline['ci_lo']:.2f}-{headline['ci_hi']:.2f}, "
                  f"10000-resample bootstrap]")
    if diverged_tasks:
        title += (f"\n* = at least one tool's output flagged as "
                  f"structural-divergence (see correctness_matrix figure)")
    if not failed.empty:
        title += f"\n{len(failed)}/{n_total} replicates excluded (non-zero exit)"
    ax.set_title(title)

    ax.legend(title="tool", loc="upper left",
              bbox_to_anchor=(1.02, 1.0), frameon=False)
    fig.tight_layout()
    save_figure(fig, "peak_memory")
    print(f"[done] wrote peak_memory.png + .pdf", file=sys.stderr)
    if headline is not None:
        print(f"[info] headline RAM ratio: {headline['competitor']} on "
              f"{headline['task']} = {headline['median_ratio']:.2f}x "
              f"[95% CI {headline['ci_lo']:.2f}-{headline['ci_hi']:.2f}]",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
