#!/usr/bin/env python3
"""Figure 4 — peak resident memory per benchmark task, grouped by tool."""
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
    FIG.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="paper")

    fig, ax = plt.subplots(figsize=(9, 5))
    sns.barplot(data=df, x="task", y="peak_rss_mb", hue="tool",
                errorbar="sd", ax=ax)
    ax.set_ylabel("peak RSS (MB)")
    ax.set_xlabel("benchmark task")
    ax.set_title("Peak resident memory — gfc vs. competitors")
    ax.legend(title="tool", loc="best", frameon=False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"peak_memory.{ext}", dpi=300, bbox_inches="tight")
    print(f"[done] wrote {FIG / 'peak_memory.png'} + .pdf", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
