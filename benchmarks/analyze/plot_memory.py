#!/usr/bin/env python3
"""Figure 4 — peak resident memory per benchmark task, grouped by tool.

Same correctness annotations as plot_wall_time.py: tasks where any
tool's output diverged from the reference get a `*` next to the task
label, and the title appends a one-line note naming which tool
diverged. Y-axis is log-scale because the chr22 plink2 datapoint sits
~25× above gfc on T2.
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
    n_total = len(df)
    failed = df[df["exit_code"] != 0]
    df_ok = df[df["exit_code"] == 0].copy()
    if df_ok.empty:
        print("[error] no successful replicates in all.tsv", file=sys.stderr)
        return 1

    if "correct" in df_ok.columns:
        df_ok["correct"] = df_ok["correct"].fillna("").astype(str)
        bad_tools = (df_ok[df_ok["correct"] == "0"]
                     .groupby("task")["tool"].apply(lambda s: ",".join(sorted(set(s)))))
    else:
        bad_tools = pd.Series(dtype=str)

    FIG.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="paper")

    fig, ax = plt.subplots(figsize=(10, 5))
    sns.barplot(data=df_ok, x="task", y="peak_rss_mb", hue="tool",
                errorbar="sd", ax=ax)
    ax.set_yscale("log")
    ax.set_ylabel("peak RSS (MB, log scale)")
    ax.set_xlabel("benchmark task")

    new_xticklabels = []
    for txt in ax.get_xticklabels():
        task = txt.get_text()
        new_xticklabels.append(f"{task}*" if task in bad_tools.index else task)
    ax.set_xticklabels(new_xticklabels)

    title = "Peak resident memory — gfc vs. competitors"
    if not bad_tools.empty:
        notes = "; ".join(f"{t}: {tools} disagreed with reference"
                          for t, tools in bad_tools.items())
        title += f"\n* = output-mismatch flagged ({notes})"
    if not failed.empty:
        title += f"\n{len(failed)}/{n_total} replicates excluded (non-zero exit)"
    ax.set_title(title, fontsize=10)
    ax.legend(title="tool", loc="best", frameon=False, fontsize=8)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"peak_memory.{ext}", dpi=300, bbox_inches="tight")
    print(f"[done] wrote {FIG / 'peak_memory.png'} + .pdf", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
