#!/usr/bin/env python3
"""Supplementary figure — cross-replicate variance per (task, tool).

Boxplot of wall_s across replicates, grouped by task on x and tool by
hue. With reps >= 2 the box collapses to a line and the median dot is
the only visible mark; with reps >= 5 (production) the boxes carry
real information about run-to-run jitter.

Paper utility: paper claims "tight wall-time variance, reproducible
across reps." This figure is the data behind that claim, suitable for
a supplementary "S1 — measurement reproducibility" panel.

Failed reps (exit_code != 0) excluded; tasks/tools with fewer than 2
successful reps render as a single point (still informative — shows
where reps are missing).
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

    FIG.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="paper")

    fig, ax = plt.subplots(figsize=(11, 5))
    sns.boxplot(data=df, x="task", y="wall_s", hue="tool", ax=ax,
                showfliers=True, fliersize=3, linewidth=0.8)
    # Overlay individual reps as small dots so n-rep is visible.
    sns.stripplot(data=df, x="task", y="wall_s", hue="tool", ax=ax,
                  dodge=True, size=2.2, color="black", alpha=0.6,
                  legend=False)
    ax.set_yscale("log")
    ax.set_ylabel("wall time (s, log scale)")
    ax.set_xlabel("benchmark task")
    n_reps = df.groupby(["task", "tool"]).size().min()
    ax.set_title(f"Per-(task, tool) wall-time distribution — boxes show "
                 f"replicate spread (min n={n_reps})")
    ax.legend(title="tool", loc="best", frameon=False, fontsize=8)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"variance.{ext}", dpi=300, bbox_inches="tight")
    print(f"[done] wrote {FIG / 'variance.png'} + .pdf", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
