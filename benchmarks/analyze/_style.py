#!/usr/bin/env python3
# Author: Benjamin Narh-Madey
"""Shared figure style for the gfc-bench manuscript.

Centralises:
  - Bioinformatics journal-conformant rcParams (sans-serif, point sizes,
    300 dpi defaults).
  - Okabe-Ito colour-blind-safe 8-colour palette.
  - Stable tool->colour mapping so "gfc" is the same colour in every
    figure (anchored to black; competitors get the remaining slots).
  - Single- and double-column figure sizes in inches converted from the
    Bioinformatics standard 88 mm / 180 mm widths.
  - save_figure() helper that writes both PNG (300 dpi raster) and PDF
    (vector for typesetter) under benchmarks/results/figures/.

Every plot_*.py script in this directory imports from this module so a
palette / size / dpi change here propagates to all eight figures.
"""
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt


# Okabe-Ito 8-colour palette — colour-blind-safe, distinguishable in
# greyscale, recommended by Wong (Nature Methods 2011) and adopted by
# Bioinformatics editorial style guidance for accessible figures.
PALETTE = {
    "black":          "#000000",
    "orange":         "#E69F00",
    "sky_blue":       "#56B4E9",
    "bluish_green":   "#009E73",
    "yellow":         "#F0E442",
    "blue":           "#0072B2",
    "vermillion":     "#D55E00",
    "reddish_purple": "#CC79A7",
}

# Stable tool -> colour mapping. gfc is anchored to black (paper subject
# is foregrounded visually); competitors are assigned in a fixed order so
# the same tool appears in the same colour across every figure. Tools
# that are not yet in the benchmark TSVs but are scheduled for Phase 4
# (e.g. BioConvert) keep their slot reserved here.
TOOL_COLORS = {
    # Paper subject — anchored to black so it never moves.
    "gfc":              PALETTE["black"],

    # Tools currently in benchmarks/results/raw/all.tsv.
    "convertf":         PALETTE["orange"],
    "plink2":           PALETTE["sky_blue"],
    "plink1.9":         PALETTE["bluish_green"],
    "EMBOSS-seqret":    PALETTE["blue"],
    "py-ref":           PALETTE["reddish_purple"],
    "gffread":          PALETTE["vermillion"],
    "AGAT":             PALETTE["yellow"],
    "ucsc-chain":       PALETTE["sky_blue"],
    "pyhmmer":          PALETTE["bluish_green"],
    "bcftools-pyref":   PALETTE["orange"],

    # Phase-4 reservations (consistent with feature_matrix.py columns).
    # Assigned from the same palette so adding them later does not break
    # any existing figure's colour scheme.
    "samtools":         PALETTE["bluish_green"],
    "bcftools":         PALETTE["orange"],
    "seqkit":           PALETTE["blue"],
    "bedtools":         PALETTE["vermillion"],
    "gff3toembl":       PALETTE["reddish_purple"],
    "OrthoFinder":      PALETTE["sky_blue"],
    "BioConvert":       PALETTE["yellow"],
    "biopython.convert": PALETTE["bluish_green"],
}

# Fallback colour for any unmapped tool name. Picked outside the
# Okabe-Ito set deliberately so a missing entry is visually obvious.
_UNKNOWN_TOOL_COLOR = "#777777"


def color_for_tool(tool: str) -> str:
    """Return the Okabe-Ito colour for `tool`, or the fallback grey."""
    return TOOL_COLORS.get(tool, _UNKNOWN_TOOL_COLOR)


def palette_for_tools(tools: Iterable[str]) -> list:
    """Return a colour list aligned to the order of `tools`."""
    return [color_for_tool(t) for t in tools]


def set_bioinformatics_style() -> None:
    """Apply Bioinformatics journal rcParams to the global matplotlib state.

    Call once at the top of `main()` before any figure / axis is created.
    Sets sans-serif font (Helvetica preferred, Arial / DejaVu fallback),
    Bioinformatics-conformant point sizes, and 300 dpi save default.
    """
    mpl.rcParams.update({
        "figure.dpi":          150,
        "savefig.dpi":         300,
        "savefig.bbox":        "tight",
        "pdf.fonttype":        42,        # TrueType, editable in Illustrator.
        "ps.fonttype":         42,
        "font.family":         "sans-serif",
        "font.sans-serif":     ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size":            8,        # Bioinformatics body text point.
        "axes.titlesize":       9,
        "axes.labelsize":       8,
        "xtick.labelsize":      7,
        "ytick.labelsize":      7,
        "legend.fontsize":      7,
        "legend.title_fontsize": 7,
        "axes.spines.top":      False,
        "axes.spines.right":    False,
        "axes.grid":            True,
        "grid.alpha":           0.3,
        "grid.linewidth":       0.4,
        "axes.axisbelow":       True,
    })


# Bioinformatics column widths: 88 mm single, 180 mm double.
# 1 in = 25.4 mm.
def single_col_size(height_in: float = 3.0) -> tuple:
    """Return (width_in, height_in) for a Bioinformatics single-column figure."""
    return (88.0 / 25.4, height_in)


def double_col_size(height_in: float = 4.5) -> tuple:
    """Return (width_in, height_in) for a Bioinformatics double-column figure."""
    return (180.0 / 25.4, height_in)


def save_figure(fig: plt.Figure, name: str,
                fmts: Sequence[str] = ("png", "pdf")) -> None:
    """Write `fig` to benchmarks/results/figures/<name>.<fmt> for each fmt.

    Resolves the output directory relative to this module so callers do
    not need to compute paths.
    """
    out_dir = Path(__file__).resolve().parent.parent / "results" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in fmts:
        fig.savefig(out_dir / f"{name}.{ext}", dpi=300, bbox_inches="tight")
