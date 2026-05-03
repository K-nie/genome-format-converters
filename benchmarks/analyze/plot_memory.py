#!/usr/bin/env python3
# Author: Benjamin Narh-Madey
"""Figure 4 — peak resident memory per benchmark task, grouped by tool.

Bar height is mean peak RSS (MB) across N replicates; error bars
encode the min/max range so reviewers see the spread directly rather
than inferring it from sd. Y-axis is log because the chr22 plink2 +
convertf points sit ~22x above gfc.

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

    n_reps = df_ok.groupby(["task", "tool"]).size().min()
    title = (f"Peak resident memory per benchmark task "
             f"(bar = mean, whiskers = min..max, n={n_reps} replicates)")
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
    return 0


if __name__ == "__main__":
    sys.exit(main())
