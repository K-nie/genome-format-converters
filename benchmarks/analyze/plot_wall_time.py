#!/usr/bin/env python3
"""Figure 3 — wall-time per benchmark task, grouped by tool.

Reads ``benchmarks/results/raw/all.tsv`` and writes
``benchmarks/results/figures/wall_time.{png,pdf}``.

Correctness-aware: when a tool's row for a task carries ``correct=0``
the task label gets a trailing ``*`` and the figure title appends a
note naming which tool diverged. Failed reps (``exit_code != 0``) are
dropped from aggregation but the count appears in the title so the
reader can see how many didn't run.
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
    n_total = len(df)
    failed = df[df["exit_code"] != 0]
    df_ok = df[df["exit_code"] == 0].copy()
    if df_ok.empty:
        print("[error] no successful replicates in all.tsv", file=sys.stderr)
        return 1

    # Tasks with at least one tool whose output diverged from the per-task
    # reference. Empty / "skip" / "1" all count as fine.
    if "correct" in df_ok.columns:
        df_ok["correct"] = df_ok["correct"].fillna("").astype(str)
        bad_tools = (df_ok[df_ok["correct"] == "0"]
                     .groupby("task")["tool"].apply(lambda s: ",".join(sorted(set(s)))))
    else:
        bad_tools = pd.Series(dtype=str)

    FIG.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="paper")

    fig, ax = plt.subplots(figsize=(10, 5))
    sns.barplot(data=df_ok, x="task", y="wall_s", hue="tool",
                errorbar="sd", ax=ax)
    ax.set_yscale("log")
    ax.set_ylabel("wall time (s, log scale)")
    ax.set_xlabel("benchmark task")

    # Annotate task labels with `*` when at least one tool's output
    # diverged from the per-task reference (so the speed bars don't
    # silently hide an output mismatch).
    new_xticklabels = []
    for txt in ax.get_xticklabels():
        task = txt.get_text()
        new_xticklabels.append(f"{task}*" if task in bad_tools.index else task)
    ax.set_xticklabels(new_xticklabels)

    title = "gfc vs. competitors — wall-clock time (mean ± sd)"
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
        fig.savefig(FIG / f"wall_time.{ext}", dpi=300, bbox_inches="tight")
    print(f"[done] wrote {FIG / 'wall_time.png'} + .pdf", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
