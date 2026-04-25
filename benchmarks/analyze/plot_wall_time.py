#!/usr/bin/env python3
"""Figure 3 — wall-time per benchmark task, grouped by tool.

Reads ``benchmarks/results/raw/all.tsv`` and writes
``benchmarks/results/figures/wall_time.{png,pdf}``.
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
        print(f"[error] {RAW} missing. Run run_bench.sh + collect_results.py first.",
              file=sys.stderr)
        return 1
    df = pd.read_csv(RAW, sep="\t")
    # Successful replicates only.
    df = df[df["exit_code"] == 0].copy()
    if df.empty:
        print("[error] no successful replicates in all.tsv", file=sys.stderr)
        return 1

    FIG.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="paper")

    fig, ax = plt.subplots(figsize=(9, 5))
    sns.barplot(data=df, x="task", y="wall_s", hue="tool",
                errorbar="sd", ax=ax)
    ax.set_yscale("log")
    ax.set_ylabel("wall time (s, log scale)")
    ax.set_xlabel("benchmark task")
    ax.set_title("gfc vs. competitors — wall-clock time (mean ± sd, n replicates)")
    ax.legend(title="tool", loc="best", frameon=False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"wall_time.{ext}", dpi=300, bbox_inches="tight")
    print(f"[done] wrote {FIG / 'wall_time.png'} + .pdf", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
