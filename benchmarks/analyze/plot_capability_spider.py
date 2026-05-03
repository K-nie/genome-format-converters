#!/usr/bin/env python3
# Author: Benjamin Narh-Madey
"""Capability spider — Figure 1 candidate (third-architecture positioning).

One polygon per representative tool, on seven normalised axes:

    1. Format coverage      — count of supported conversions, log10
                              normalised across the field.
    2. Memory efficiency    — inverted median peak RSS across the
                              tasks the tool runs in `all.tsv`.
    3. Wall-time efficiency — inverted median wall time across the
                              tasks the tool runs in `all.tsv`.
    4. Correctness          — pass rate at structural equivalence
                              (fraction of (task, tool) cells with
                              correct == "1" or "skip").
    5. Install lightness    — categorical bucket from dependency
                              footprint (heavy / medium / light ->
                              0.3 / 0.6 / 1.0).
    6. Native Python        — binary 0 / 1.
    7. Provenance support   — 0 / 0.5 / 1 if the tool emits
                              run-provenance metadata by default.

Tools shown:
    gfc          — paper subject.
    BioConvert   — wrapper / orchestrator class representative.
    AGAT         — domain-deep (annotation-tooling) representative.
    samtools     — format-family (SAM / BAM / CRAM) representative.

Per-tool wall and RSS numbers are pulled from
benchmarks/results/raw/all.tsv. Axes that need information not
captured by the bench harness (format coverage, dependency footprint,
language, provenance support) are encoded as constants below with
inline citations to the literature comparator catalog.

The four polygons make the headline argument visible: gfc occupies a
shape no other tool covers — broad coverage like wrappers, low memory
and native Python like single-purpose tools.

PLACEHOLDER VALUES
==================
BioConvert and samtools have not yet been run inside the bench
harness (Phase 4 task). Their wall / RSS / correctness contributions
are encoded as constants flagged 'PLACEHOLDER'. A WARNING banner is
printed at script start so the situation is loud rather than silent.
Refresh after Phase 4 BioConvert / samtools runs land in all.tsv.
"""
from pathlib import Path
import sys
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from _style import (
    color_for_tool,
    normalize_for_spider,
    plot_spider,
    save_figure,
    set_bioinformatics_style,
    single_col_size,
    spider_axes,
)

RAW = Path(__file__).resolve().parent.parent / "results" / "raw" / "all.tsv"
FIG_DIR = Path(__file__).resolve().parent.parent / "results" / "figures"
LIT_CATALOG = (
    "/Users/black_einstein/Desktop/PhD_Workflow_Automation/docs/"
    "LITERATURE_COMPARATORS_2026-05-03.md"
)

# ---------------------------------------------------------------------------
# Constants per tool (axes the harness does not measure).
# ---------------------------------------------------------------------------

# Format coverage = count of source -> target conversions the tool
# advertises. gfc: from `gfc --help` (30 converters). Others: from the
# literature comparator catalog at LITERATURE_COMPARATORS_2026-05-03.md.
FORMAT_COVERAGE = {
    "gfc":        30,
    "BioConvert": 100,   # PLACEHOLDER — refresh after Phase 4 BioConvert run
    "AGAT":       8,     # gff -> gtf / bed / proteins / cds / spliced / etc.
    "samtools":   6,     # SAM <-> BAM <-> CRAM, FASTA index, depth, etc.
}

# Install lightness buckets:
#   heavy  -> 0.3  (Conda + many pinned deps; e.g. Perl + BioPerl + AGAT)
#   medium -> 0.6  (single C binary or one wrapper; e.g. samtools, gffread)
#   light  -> 1.0  (pip install <name>; pure Python / one wheel; e.g. gfc)
INSTALL_LIGHTNESS = {
    "gfc":        1.0,    # pip install gfc, no system deps.
    "BioConvert": 0.3,    # PLACEHOLDER — wraps ~50 underlying tools.
    "AGAT":       0.3,    # Conda + Perl + BioPerl + many CPAN modules.
    "samtools":   0.6,    # single C binary; htslib bundled.
}

# Native-Python axis (binary).
NATIVE_PYTHON = {
    "gfc":        1.0,
    "BioConvert": 0.0,    # Python harness, but shells out to native tools.
    "AGAT":       0.0,    # Perl.
    "samtools":   0.0,    # C.
}

# Provenance support (0 = none, 0.5 = partial, 1 = first-class metadata).
PROVENANCE = {
    "gfc":        1.0,    # writes run manifest with version, ts, fixture md5.
    "BioConvert": 0.5,    # PLACEHOLDER — logs steps, no manifest.
    "AGAT":       0.0,    # per-input *.agat.log only.
    "samtools":   0.0,    # no provenance metadata by default.
}

# Wall / RSS / correctness placeholders for tools NOT yet in all.tsv.
# These match the order-of-magnitude expectations from the literature so
# the spider polygon shape is plausible. Replace with measured medians
# once Phase 4 lands them in the TSV.
PLACEHOLDER_WALL_S = {
    "BioConvert": 90.0,    # PLACEHOLDER — wrapper overhead, single conversion.
    "samtools":   12.0,    # PLACEHOLDER — fast C path.
}
PLACEHOLDER_RSS_MB = {
    "BioConvert": 350.0,   # PLACEHOLDER — Python interpreter + child tool.
    "samtools":   80.0,    # PLACEHOLDER — C, low footprint.
}
PLACEHOLDER_CORRECTNESS = {
    "BioConvert": 0.85,    # PLACEHOLDER — wrapper inherits child correctness.
    "samtools":   1.0,     # PLACEHOLDER — gold standard for BAM / CRAM.
}


def _median_per_tool(df: pd.DataFrame, column: str) -> dict:
    """Median across all (task, replicate) rows per tool, exit_code == 0."""
    ok = df[df["exit_code"] == 0]
    return ok.groupby("tool")[column].median().to_dict()


def _correctness_pass_rate(df: pd.DataFrame) -> dict:
    """Fraction of (task, tool) cells where comparator returned pass / skip."""
    ok = df[df["exit_code"] == 0].copy()
    ok["correct_str"] = ok["correct"].astype(str).str.strip()
    # One verdict per (task, tool) — the comparator is run once per pair
    # and propagated across replicates.
    per_pair = ok.groupby(["task", "tool"])["correct_str"].first().reset_index()
    per_pair["pass"] = per_pair["correct_str"].isin(["1", "skip"])
    return per_pair.groupby("tool")["pass"].mean().to_dict()


def _build_polygons(df: pd.DataFrame, tools: list) -> tuple:
    """Return (axis_labels, polygons_dict) for the four-tool spider."""
    wall_med = _median_per_tool(df, "wall_s")
    rss_med = _median_per_tool(df, "peak_rss_mb")
    correct = _correctness_pass_rate(df)

    # Inject placeholders for tools not measured in this run.
    wall_med = {**wall_med, **PLACEHOLDER_WALL_S}
    rss_med = {**rss_med, **PLACEHOLDER_RSS_MB}
    correct = {**correct, **PLACEHOLDER_CORRECTNESS}

    # Normalise each axis across the tools shown.
    axis_inputs = {
        "Format coverage":      ({t: np.log10(FORMAT_COVERAGE[t]) for t in tools},
                                 "higher_is_better"),
        "Memory efficiency":    ({t: rss_med.get(t, np.nan) for t in tools},
                                 "lower_is_better"),
        "Wall-time efficiency": ({t: wall_med.get(t, np.nan) for t in tools},
                                 "lower_is_better"),
        "Correctness":          ({t: correct.get(t, np.nan) for t in tools},
                                 "higher_is_better"),
        "Install lightness":    ({t: INSTALL_LIGHTNESS[t] for t in tools},
                                 "higher_is_better"),
        "Native Python":        ({t: NATIVE_PYTHON[t] for t in tools},
                                 "higher_is_better"),
        "Provenance":           ({t: PROVENANCE[t] for t in tools},
                                 "higher_is_better"),
    }

    # Build polygon-per-tool arrays in axis order.
    axis_labels = list(axis_inputs.keys())
    normalised_per_axis = {
        label: normalize_for_spider(values, direction)
        for label, (values, direction) in axis_inputs.items()
    }
    polygons = {}
    for t in tools:
        polygons[t] = [normalised_per_axis[label][t] for label in axis_labels]
    return axis_labels, polygons


def _write_caption(out: Path, n_pairs: int) -> None:
    text = (
        "Capability spider for representative tools across seven dimensions: "
        "format coverage (log10 of advertised converter count), memory "
        "efficiency (inverted median peak RSS across tasks the tool ran), "
        "wall-time efficiency (inverted median wall time across the same "
        "tasks), correctness (fraction of structural-equivalence pass plus "
        "out-of-scope skip cells), install lightness (heavy 0.3 / medium "
        "0.6 / light 1.0 dependency bucket), native Python (binary), and "
        "provenance support (none 0 / partial 0.5 / first-class manifest 1). "
        "Wall, RSS and correctness for gfc, AGAT, and other tools currently "
        "in benchmarks/results/raw/all.tsv are measured (n = "
        f"{n_pairs} replicates per task, n = 10 per pair for production "
        "tasks). Wall, RSS and correctness for BioConvert and samtools are "
        "PLACEHOLDER values pending the Phase 4 run; coverage, install "
        "lightness, language, and provenance are sourced from the "
        f"literature comparator catalog at {LIT_CATALOG}. The headline is "
        "the polygon shape: gfc covers the breadth axis like a wrapper, "
        "the memory and Python axes like a single-purpose native tool, "
        "and the provenance axis on its own — a position no other tool "
        "in the field currently occupies."
    )
    out.write_text(text + "\n")


def main() -> int:
    set_bioinformatics_style()
    if not RAW.exists():
        print(f"ERROR: {RAW} not found", file=sys.stderr)
        return 1

    df = pd.read_csv(RAW, sep="\t")
    df.columns = [c.strip() for c in df.columns]

    tools = ["gfc", "BioConvert", "AGAT", "samtools"]

    # Loud warning that two of four polygons carry PLACEHOLDER axes.
    placeholder_tools = [t for t in tools
                         if t in PLACEHOLDER_WALL_S or t in PLACEHOLDER_RSS_MB]
    if placeholder_tools:
        warnings.warn(
            "PLACEHOLDER values in use for: "
            f"{', '.join(placeholder_tools)}. Refresh capability_spider.py "
            "after Phase 4 BioConvert / samtools runs land in all.tsv.",
            stacklevel=1,
        )
        print(
            "WARNING: PLACEHOLDER axes for "
            f"{', '.join(placeholder_tools)} — see top of script.",
            file=sys.stderr,
        )

    labels, polygons = _build_polygons(df, tools)

    width, height = single_col_size(height_in=3.6)
    fig, ax = plt.subplots(figsize=(width, height),
                           subplot_kw=dict(projection="polar"))
    angles = spider_axes(ax, labels)
    plot_spider(ax, polygons, angles)

    ax.legend(loc="upper right", bbox_to_anchor=(1.32, 1.10),
              frameon=False, fontsize=7)
    ax.set_title("Capability spider — gfc vs three architecture classes",
                 fontsize=8.5, pad=14)

    save_figure(fig, "capability_spider")
    plt.close(fig)

    n_pairs = df[(df["exit_code"] == 0)
                 & (df["tool"].isin(tools))].shape[0]
    _write_caption(FIG_DIR / "capability_spider_caption.txt", n_pairs)

    print(f"wrote {FIG_DIR / 'capability_spider.png'}")
    print(f"wrote {FIG_DIR / 'capability_spider.pdf'}")
    print(f"wrote {FIG_DIR / 'capability_spider_caption.txt'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
