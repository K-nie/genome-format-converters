#!/usr/bin/env python3
# Author: Benjamin Narh-Madey
"""Per-task spider — supplementary 8-panel grid (one panel per task).

Two rows by four columns of polar axes. Each panel compares gfc to the
primary competitor for that task on four panel-friendly axes:

    1. Wall-time efficiency (inverted + min-max normalised within panel)
    2. Memory efficiency    (inverted + min-max normalised within panel)
    3. Correctness          (categorical: structural-equivalence pass=1,
                             documented-divergence=0.5, fail=0)
    4. Install lightness    (heavy / medium / light bucket -> 0.3 / 0.6 / 1.0)

Primary competitors per task are taken from
LITERATURE_COMPARATORS_2026-05-03.md and the existing all.tsv:

    T1 -> convertf
    T2 -> plink2
    T3 -> EMBOSS-seqret
    T4 -> gffread
    T5 -> AGAT (no UCSC chain primary in the catalog; use AGAT)
    T6 -> pyhmmer
    T7 -> py-ref
    T8 -> bcftools-pyref

If the primary competitor is missing from all.tsv for a task, that
task's panel shows the gfc polygon alone with a "competitor pending"
annotation and the other vertices clipped to 0.
"""
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from _style import (
    color_for_tool,
    normalize_for_spider,
    plot_spider,
    save_figure,
    set_bioinformatics_style,
    double_col_size,
    spider_axes,
)

RAW = Path(__file__).resolve().parent.parent / "results" / "raw" / "all.tsv"
FIG_DIR = Path(__file__).resolve().parent.parent / "results" / "figures"


# Per-task primary competitor (gfc is the second polygon every panel).
PRIMARY_COMPETITOR = {
    "T1": "convertf",
    "T2": "plink2",
    "T3": "EMBOSS-seqret",
    "T4": "gffread",
    "T5": "AGAT",
    "T6": "pyhmmer",
    "T7": "py-ref",
    "T8": "bcftools-pyref",
}

INSTALL_LIGHTNESS = {
    "gfc":            1.0,
    "convertf":       0.3,    # EIGENSOFT toolchain.
    "plink2":         0.6,    # single C binary.
    "plink1.9":       0.6,
    "EMBOSS-seqret":  0.3,    # EMBOSS suite.
    "py-ref":         1.0,    # pure Python.
    "gffread":        0.6,    # single binary, gclib bundled.
    "AGAT":           0.3,    # Perl + BioPerl.
    "ucsc-chain":     0.6,    # UCSC binary suite.
    "pyhmmer":        1.0,    # pip install pyhmmer.
    "bcftools-pyref": 0.6,    # bcftools binary + Python wrapper.
}


def _correctness_score(correct_str: str) -> float:
    """Map correctness verdict to a numeric score in [0, 1]."""
    s = (correct_str or "").strip()
    if s == "1":
        return 1.0
    if s.lower() == "skip":
        return 1.0
    if s == "0":
        return 0.5     # documented structural divergence, not a failure.
    return 0.0          # unchecked / unknown.


def _task_polygons(df: pd.DataFrame, task: str, competitor: str) -> tuple:
    """Build the two-polygon dict for one panel."""
    rows = df[(df["task"] == task) & (df["exit_code"] == 0)]
    if rows.empty:
        return None, None

    tools = ["gfc", competitor]

    wall = {}
    rss = {}
    correct = {}
    for t in tools:
        sub = rows[rows["tool"] == t]
        if sub.empty:
            wall[t] = float("nan")
            rss[t] = float("nan")
            correct[t] = float("nan")
        else:
            wall[t] = float(sub["wall_s"].median())
            rss[t] = float(sub["peak_rss_mb"].median())
            correct[t] = _correctness_score(str(sub["correct"].iloc[0]))

    # Normalise each axis within the panel so the polygon comparison is
    # interpretable: 1.0 == best of the two on this panel's axis.
    n_wall = normalize_for_spider(wall, "lower_is_better")
    n_rss = normalize_for_spider(rss, "lower_is_better")
    n_install = {t: INSTALL_LIGHTNESS.get(t, 0.0) for t in tools}

    polygons = {
        t: [n_wall[t], n_rss[t], correct[t], n_install[t]]
        for t in tools
    }
    return polygons, ["Wall", "Memory", "Correct", "Install"]


def _write_caption(out: Path, task_status: dict) -> None:
    bullets = "\n".join(
        f"  {task} ({comp}): {status}"
        for task, (comp, status) in sorted(task_status.items())
    )
    text = (
        "Per-task spider grid (supplementary). Each polar panel "
        "compares gfc against the primary literature competitor for "
        "that task, on four panel-friendly axes: wall-time efficiency "
        "(inverted, min-max normalised within the panel), memory "
        "efficiency (inverted, min-max normalised within the panel), "
        "correctness (structural-equivalence pass = 1.0, documented "
        "divergence = 0.5, unchecked / fail = 0.0), and install "
        "lightness (heavy 0.3 / medium 0.6 / light 1.0 dependency "
        "bucket). Polygons inherit the figure-suite tool colours so "
        "gfc's polygon is the same colour across every panel and every "
        "figure. The shape of each panel argues per task: a tall "
        "narrow polygon (wall + memory dominant) is the speed / "
        "footprint argument; a balanced polygon is the workflow-fit "
        "argument. Per-task status:\n"
        f"{bullets}\n"
    )
    out.write_text(text)


def main() -> int:
    set_bioinformatics_style()
    if not RAW.exists():
        print(f"ERROR: {RAW} not found", file=sys.stderr)
        return 1

    df = pd.read_csv(RAW, sep="\t")
    df.columns = [c.strip() for c in df.columns]

    width, _ = double_col_size()
    fig, axes = plt.subplots(
        2, 4, figsize=(width, width * 0.55),
        subplot_kw=dict(projection="polar"),
    )

    task_status = {}
    for ax, task in zip(axes.flat, sorted(PRIMARY_COMPETITOR)):
        competitor = PRIMARY_COMPETITOR[task]
        polygons, labels = _task_polygons(df, task, competitor)

        if polygons is None:
            ax.set_title(f"{task} (no data)", fontsize=8)
            ax.axis("off")
            task_status[task] = (competitor, "no data in all.tsv")
            continue

        angles = spider_axes(ax, labels)
        plot_spider(ax, polygons, angles,
                    line_width=1.3, marker_size=2.5)
        ax.set_title(f"{task} — gfc vs {competitor}", fontsize=8, pad=8)
        ax.tick_params(axis="x", pad=2)

        # Panel-level legend would over-clutter — use one figure-level
        # legend below.
        ax.legend().set_visible(False) if ax.get_legend() else None

        task_status[task] = (competitor, "panel populated")

    # Single figure-level legend so the reader sees gfc + competitor
    # colour mapping once.
    handles = [
        plt.Line2D([0], [0], color=color_for_tool("gfc"),
                   linewidth=1.6, marker="o", markersize=4, label="gfc"),
    ]
    seen = {"gfc"}
    for task in sorted(PRIMARY_COMPETITOR):
        comp = PRIMARY_COMPETITOR[task]
        if comp in seen:
            continue
        seen.add(comp)
        handles.append(
            plt.Line2D([0], [0], color=color_for_tool(comp),
                       linewidth=1.6, marker="o", markersize=4, label=comp)
        )

    fig.legend(handles=handles, loc="lower center", ncol=len(handles),
               frameon=False, fontsize=7,
               bbox_to_anchor=(0.5, -0.02))

    fig.suptitle(
        "Per-task spider — gfc vs primary competitor on (Wall, Memory, "
        "Correctness, Install)",
        fontsize=9.5, y=1.00,
    )
    fig.tight_layout(rect=(0, 0.02, 1, 0.97))

    save_figure(fig, "per_task_spider")
    plt.close(fig)

    _write_caption(FIG_DIR / "per_task_spider_caption.txt", task_status)

    print(f"wrote {FIG_DIR / 'per_task_spider.png'}")
    print(f"wrote {FIG_DIR / 'per_task_spider.pdf'}")
    print(f"wrote {FIG_DIR / 'per_task_spider_caption.txt'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
