#!/usr/bin/env python3
# Author: Benjamin Narh-Madey
"""Supplementary figure — wall-time coefficient of variation per (task, tool).

Replaces the earlier boxplot/stripplot rendering. On the n=10
production data, sd / median is sub-percent for most (task, tool)
pairs; a log-y boxplot squashes the spread to invisibility. CV
directly answers "how reproducible is each tool's wall time across
replicates," which is the reproducibility claim the paper Methods
section makes.

CV = (std / mean) x 100%. Lower = more reproducible. Bars are sorted
by CV ascending within each task panel so the eye reads "tightest
on the left."

Failed replicates (exit_code != 0) excluded; (task, tool) pairs with
fewer than 2 successful replicates are dropped (CV undefined).
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
    df = df[df["exit_code"] == 0].copy()
    if df.empty:
        print("[error] no successful reps in all.tsv", file=sys.stderr)
        return 1

    grouped = df.groupby(["task", "tool"])["wall_s"]
    counts = grouped.size()
    means = grouped.mean()
    stds = grouped.std(ddof=1)
    cv_pct = (stds / means * 100.0).where(counts >= 2)
    cv_df = cv_pct.dropna().reset_index(name="cv_pct")
    if cv_df.empty:
        print("[error] no (task, tool) pairs with >=2 successful reps", file=sys.stderr)
        return 1

    set_bioinformatics_style()

    tool_order = ["gfc"] + sorted(t for t in cv_df["tool"].unique() if t != "gfc")
    palette = {t: color_for_tool(t) for t in tool_order}

    fig, ax = plt.subplots(figsize=double_col_size(height_in=4.5))
    sns.barplot(data=cv_df, x="task", y="cv_pct", hue="tool",
                hue_order=tool_order, palette=palette, ax=ax)
    ax.set_ylabel("wall-time coefficient of variation (%) — lower is more reproducible")
    ax.set_xlabel("benchmark task")

    n_reps_min = int(counts.min())
    ax.set_title(f"Per-(task, tool) wall-time CV across replicates "
                 f"(min n={n_reps_min}); CV = std/mean x 100%, "
                 f"lower = more reproducible")
    ax.legend(title="tool", loc="upper left",
              bbox_to_anchor=(1.02, 1.0), frameon=False)
    fig.tight_layout()
    save_figure(fig, "variance")
    print(f"[done] wrote variance.png + .pdf", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
