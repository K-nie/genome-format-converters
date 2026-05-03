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


# ---------------------------------------------------------------------------
# Spider / radar plot helpers (Phase 2 follow-on)
# ---------------------------------------------------------------------------
#
# Spider plots map a small set of normalised metrics onto a polar grid so
# that each tool / class becomes a polygon. The shape of that polygon is
# the headline visual. Three rules govern the helpers below:
#
#   1. All axes must be normalised to [0, 1] BEFORE plotting. The caller
#      passes raw metric values to normalize_for_spider() and gets back a
#      dict of polygons ready for plot_spider().
#   2. Direction matters. Wall time and peak RSS are 'lower_is_better';
#      correctness and format coverage are 'higher_is_better'. The
#      normaliser flips the lower-is-better axes so '1.0 = best on this
#      axis' across the whole plot.
#   3. Tool colours flow from TOOL_COLORS so a tool's polygon is the
#      same colour in the spider as in every other figure.
import math

import numpy as np


def normalize_for_spider(values: dict, direction: str = "higher_is_better") -> dict:
    """Normalise a dict-of-tool->raw-value to dict-of-tool->[0,1].

    Parameters
    ----------
    values : dict
        Mapping from tool / class name to raw metric value (float). NaN
        and None are coerced to 0 in the output.
    direction : str
        'higher_is_better' -> normalise by dividing each value by the
            maximum across tools (so the best gets 1.0, others scale
            linearly). Use for correctness, format coverage, etc.
        'lower_is_better' -> invert first, then min-max normalise so the
            smallest raw value (best) maps to 1.0 and the largest (worst)
            to 0.0. Use for wall time, peak RSS.

    Returns
    -------
    dict
        Same keys as `values`; values are floats in [0, 1].
    """
    if not values:
        return {}

    # Strip NaN / None so they do not pollute the max/min.
    clean = {}
    for k, v in values.items():
        if v is None:
            continue
        try:
            fv = float(v)
        except (TypeError, ValueError):
            continue
        if math.isnan(fv):
            continue
        clean[k] = fv

    if not clean:
        return {k: 0.0 for k in values}

    if direction == "higher_is_better":
        denom = max(clean.values())
        if denom <= 0:
            return {k: 0.0 for k in values}
        out = {}
        for k in values:
            v = clean.get(k)
            out[k] = 0.0 if v is None else max(0.0, min(1.0, v / denom))
        return out

    if direction == "lower_is_better":
        lo = min(clean.values())
        hi = max(clean.values())
        rng = hi - lo
        out = {}
        for k in values:
            v = clean.get(k)
            if v is None:
                out[k] = 0.0
            elif rng <= 0:
                # All tools tied — give every present tool 1.0.
                out[k] = 1.0
            else:
                out[k] = max(0.0, min(1.0, 1.0 - (v - lo) / rng))
        return out

    raise ValueError(f"unknown direction: {direction!r}")


def spider_axes(ax, labels, ylim: tuple = (0.0, 1.0)) -> "np.ndarray":
    """Configure a polar Axes for a spider / radar plot.

    Parameters
    ----------
    ax : matplotlib polar Axes
        Created via ``plt.subplots(subplot_kw=dict(projection='polar'))``.
    labels : sequence of str
        N angular axis labels. They are placed evenly around the circle
        starting at the top (12 o'clock) and going clockwise.
    ylim : (float, float)
        Radial range. Default (0, 1) matches the normaliser output.

    Returns
    -------
    np.ndarray
        Length-N array of angles (radians) for the angular axes; pass
        this to plot_spider() as the ``angles`` argument.
    """
    n = len(labels)
    # Angles for each axis, evenly spaced; start at the top (pi/2) and
    # go clockwise so axis 1 sits at 12 o'clock and the rest read
    # left-to-right naturally.
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    angles = (np.pi / 2) - angles  # rotate so first axis is at top.
    angles = angles % (2 * np.pi)

    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)  # clockwise, matching axis enumeration.
    ax.set_xticks(np.linspace(0, 2 * np.pi, n, endpoint=False))
    ax.set_xticklabels(labels, fontsize=7)

    ax.set_ylim(*ylim)
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(["", "", "", "", ""])  # avoid radial clutter.
    ax.grid(True, alpha=0.4, linewidth=0.4)
    ax.spines["polar"].set_visible(True)
    ax.spines["polar"].set_linewidth(0.5)
    ax.spines["polar"].set_color("#999999")

    return angles


def plot_spider(ax, polygons: dict, angles: "np.ndarray",
                fill_alpha: float = 0.15, line_width: float = 1.6,
                marker: str = "o", marker_size: float = 3.0,
                color_map: dict = None) -> None:
    """Render one or more polygons onto a configured polar Axes.

    Parameters
    ----------
    ax : matplotlib polar Axes
        Must have been initialised with spider_axes() first so the
        angular ticks / radial limits are correct.
    polygons : dict
        Mapping from polygon label (tool name, class name, etc.) to a
        sequence of N normalised values in [0, 1].
    angles : np.ndarray
        Length-N angle array returned by spider_axes(). The polygon is
        closed by appending the first vertex to the end of the path so
        the spider line connects all the way around.
    fill_alpha : float
        Alpha for the polygon interior fill. Low (default 0.15) so
        overlapping polygons remain visible.
    line_width : float
        Width of the polygon outline in points.
    marker : str
        Marker style at each vertex. Set to "" to suppress markers.
    marker_size : float
        Marker size in points.
    color_map : dict | None
        Optional override mapping from polygon label to colour. If None,
        falls back to TOOL_COLORS via color_for_tool() so polygons match
        the rest of the figure suite.
    """
    # Close the polygon path so the line returns to its starting vertex.
    closed_angles = np.concatenate([angles, angles[:1]])

    for name, values in polygons.items():
        if color_map is not None and name in color_map:
            c = color_map[name]
        else:
            c = color_for_tool(name)

        v = np.asarray(list(values), dtype=float)
        v_closed = np.concatenate([v, v[:1]])

        ax.plot(closed_angles, v_closed,
                color=c, linewidth=line_width,
                marker=marker, markersize=marker_size,
                label=name)
        ax.fill(closed_angles, v_closed,
                color=c, alpha=fill_alpha)
