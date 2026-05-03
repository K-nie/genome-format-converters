#!/usr/bin/env python3
# Author: Benjamin Narh-Madey
"""Architecture spider — discussion-section figure.

Three polygons, one per architecture class, on the same seven axes
used by plot_capability_spider.py:

    1. Format coverage      (log10 normalised)
    2. Memory efficiency    (inverted median peak RSS)
    3. Wall-time efficiency (inverted median wall time)
    4. Correctness          (structural-equivalence pass rate)
    5. Install lightness    (heavy / medium / light bucket)
    6. Native Python        (binary)
    7. Provenance support   (none / partial / first-class)

The class polygons are computed as the unweighted mean of the
constituent tools' raw values, then normalised across classes. This
keeps the spider readable: instead of N tool polygons, the reader sees
three competing design philosophies and where gfc lands relative to
both.

CLASSES
=======
    Domain-deep — samtools, bcftools, AGAT, gffread.
                  One-format-family experts. High wall-time efficiency,
                  high memory efficiency, low coverage breadth.
    Wrappers    — BioConvert, EMBOSS-seqret, biopython.convert.
                  Coverage-first orchestrators. Broad reach, but pay
                  Python-startup + child-tool overhead per call.
    gfc         — single-tool, native-Python, manifest-emitting design.

The third polygon is a single tool (gfc) so the comparison stays
honest: gfc's polygon is its measured shape from all.tsv, not a
class-aggregate average.

PLACEHOLDER VALUES
==================
BioConvert, EMBOSS-seqret, biopython.convert, samtools, bcftools,
gffread are not all yet in all.tsv. Wall, RSS, and correctness for
those tools are encoded as constants below flagged 'PLACEHOLDER'. A
WARNING banner prints at script start. Refresh after Phase 4 lands.
"""
from pathlib import Path
import sys
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from _style import (
    PALETTE,
    normalize_for_spider,
    plot_spider,
    save_figure,
    set_bioinformatics_style,
    double_col_size,
    spider_axes,
)

RAW = Path(__file__).resolve().parent.parent / "results" / "raw" / "all.tsv"
FIG_DIR = Path(__file__).resolve().parent.parent / "results" / "figures"


# ---------------------------------------------------------------------------
# Class membership.
# ---------------------------------------------------------------------------
DOMAIN_DEEP = ["samtools", "bcftools", "AGAT", "gffread"]
WRAPPERS = ["BioConvert", "EMBOSS-seqret", "biopython.convert"]


# ---------------------------------------------------------------------------
# Per-tool constants (axes the harness does not measure).
# ---------------------------------------------------------------------------
FORMAT_COVERAGE = {
    "gfc":               30,
    "samtools":          6,
    "bcftools":          5,
    "AGAT":              8,
    "gffread":           4,
    "BioConvert":        100,    # PLACEHOLDER
    "EMBOSS-seqret":     45,     # PLACEHOLDER
    "biopython.convert": 25,     # PLACEHOLDER
}

INSTALL_LIGHTNESS = {
    "gfc":               1.0,
    "samtools":          0.6,
    "bcftools":          0.6,
    "AGAT":              0.3,
    "gffread":           0.6,
    "BioConvert":        0.3,    # PLACEHOLDER
    "EMBOSS-seqret":     0.3,    # PLACEHOLDER (EMBOSS suite is heavy)
    "biopython.convert": 1.0,    # pip install biopython
}

NATIVE_PYTHON = {
    "gfc":               1.0,
    "samtools":          0.0,
    "bcftools":          0.0,
    "AGAT":              0.0,
    "gffread":           0.0,
    "BioConvert":        0.0,    # shells out to native tools
    "EMBOSS-seqret":     0.0,    # C
    "biopython.convert": 1.0,
}

PROVENANCE = {
    "gfc":               1.0,
    "samtools":          0.0,
    "bcftools":          0.0,
    "AGAT":              0.0,
    "gffread":           0.0,
    "BioConvert":        0.5,    # PLACEHOLDER — logs steps
    "EMBOSS-seqret":     0.0,
    "biopython.convert": 0.0,
}

# Wall / RSS / correctness placeholders for tools not in all.tsv.
PLACEHOLDER_WALL_S = {
    "samtools":          12.0,
    "bcftools":          15.0,
    "gffread":           4.0,
    "BioConvert":        90.0,
    "EMBOSS-seqret":     60.0,
    "biopython.convert": 70.0,
}
PLACEHOLDER_RSS_MB = {
    "samtools":          80.0,
    "bcftools":          90.0,
    "gffread":           60.0,
    "BioConvert":        350.0,
    "EMBOSS-seqret":     200.0,
    "biopython.convert": 220.0,
}
PLACEHOLDER_CORRECTNESS = {
    "samtools":          1.0,
    "bcftools":          1.0,
    "gffread":           1.0,
    "BioConvert":        0.85,
    "EMBOSS-seqret":     0.9,
    "biopython.convert": 0.9,
}


def _median_per_tool(df: pd.DataFrame, column: str) -> dict:
    ok = df[df["exit_code"] == 0]
    return ok.groupby("tool")[column].median().to_dict()


def _correctness_pass_rate(df: pd.DataFrame) -> dict:
    ok = df[df["exit_code"] == 0].copy()
    ok["correct_str"] = ok["correct"].astype(str).str.strip()
    per_pair = ok.groupby(["task", "tool"])["correct_str"].first().reset_index()
    per_pair["pass"] = per_pair["correct_str"].isin(["1", "skip"])
    return per_pair.groupby("tool")["pass"].mean().to_dict()


def _class_mean(values: dict, members: list) -> float:
    """Unweighted mean across class members; skip absent / NaN entries."""
    xs = [values[m] for m in members
          if m in values and values[m] is not None and not _isnan(values[m])]
    if not xs:
        return float("nan")
    return float(np.mean(xs))


def _isnan(x) -> bool:
    try:
        return bool(np.isnan(float(x)))
    except (TypeError, ValueError):
        return False


def _build_polygons(df: pd.DataFrame) -> tuple:
    wall_med = {**_median_per_tool(df, "wall_s"), **PLACEHOLDER_WALL_S}
    rss_med = {**_median_per_tool(df, "peak_rss_mb"), **PLACEHOLDER_RSS_MB}
    correct = {**_correctness_pass_rate(df), **PLACEHOLDER_CORRECTNESS}

    # Per-axis raw values for each polygon (gfc + 2 class aggregates).
    polygons_raw = {
        "gfc": {
            "Format coverage":      np.log10(FORMAT_COVERAGE["gfc"]),
            "Memory efficiency":    rss_med.get("gfc", float("nan")),
            "Wall-time efficiency": wall_med.get("gfc", float("nan")),
            "Correctness":          correct.get("gfc", float("nan")),
            "Install lightness":    INSTALL_LIGHTNESS["gfc"],
            "Native Python":        NATIVE_PYTHON["gfc"],
            "Provenance":           PROVENANCE["gfc"],
        },
        "Domain-deep": {
            "Format coverage":      _class_mean(
                {k: np.log10(v) for k, v in FORMAT_COVERAGE.items()},
                DOMAIN_DEEP),
            "Memory efficiency":    _class_mean(rss_med, DOMAIN_DEEP),
            "Wall-time efficiency": _class_mean(wall_med, DOMAIN_DEEP),
            "Correctness":          _class_mean(correct, DOMAIN_DEEP),
            "Install lightness":    _class_mean(INSTALL_LIGHTNESS, DOMAIN_DEEP),
            "Native Python":        _class_mean(NATIVE_PYTHON, DOMAIN_DEEP),
            "Provenance":           _class_mean(PROVENANCE, DOMAIN_DEEP),
        },
        "Wrappers": {
            "Format coverage":      _class_mean(
                {k: np.log10(v) for k, v in FORMAT_COVERAGE.items()},
                WRAPPERS),
            "Memory efficiency":    _class_mean(rss_med, WRAPPERS),
            "Wall-time efficiency": _class_mean(wall_med, WRAPPERS),
            "Correctness":          _class_mean(correct, WRAPPERS),
            "Install lightness":    _class_mean(INSTALL_LIGHTNESS, WRAPPERS),
            "Native Python":        _class_mean(NATIVE_PYTHON, WRAPPERS),
            "Provenance":           _class_mean(PROVENANCE, WRAPPERS),
        },
    }

    axis_labels = list(polygons_raw["gfc"].keys())
    directions = {
        "Format coverage":      "higher_is_better",
        "Memory efficiency":    "lower_is_better",
        "Wall-time efficiency": "lower_is_better",
        "Correctness":          "higher_is_better",
        "Install lightness":    "higher_is_better",
        "Native Python":        "higher_is_better",
        "Provenance":           "higher_is_better",
    }

    # Normalise per axis across the three polygons.
    polygons_norm = {name: [] for name in polygons_raw}
    for label in axis_labels:
        per_axis = {name: polygons_raw[name][label] for name in polygons_raw}
        n = normalize_for_spider(per_axis, directions[label])
        for name in polygons_raw:
            polygons_norm[name].append(n[name])

    return axis_labels, polygons_norm


def _write_caption(out: Path) -> None:
    text = (
        "Architecture spider for the third-architecture argument. Three "
        "polygons compare gfc to two existing design philosophies: "
        "domain-deep tools (samtools, bcftools, AGAT, gffread), each "
        "expert in one format family, and wrappers (BioConvert, "
        "EMBOSS-seqret, biopython.convert), which trade per-call "
        "overhead for breadth of coverage. Class polygons are unweighted "
        "means across the constituent tools on each of the seven axes "
        "(format coverage in log10, inverted median wall and RSS, "
        "structural-equivalence pass rate, install-lightness bucket, "
        "native-Python binary, and provenance support). gfc is plotted "
        "as a single-tool polygon — its values are not aggregated. The "
        "wall, RSS and correctness contributions for tools not yet in "
        "the harness (BioConvert, samtools, bcftools, gffread, "
        "EMBOSS-seqret, biopython.convert) are PLACEHOLDER constants "
        "and will be replaced with measured medians after Phase 4. The "
        "figure supports the manuscript Discussion claim that gfc "
        "represents a third architecture: native Python, broad coverage, "
        "low memory, and first-class provenance metadata in one polygon."
    )
    out.write_text(text + "\n")


def main() -> int:
    set_bioinformatics_style()
    if not RAW.exists():
        print(f"ERROR: {RAW} not found", file=sys.stderr)
        return 1

    df = pd.read_csv(RAW, sep="\t")
    df.columns = [c.strip() for c in df.columns]

    placeholder_tools = sorted(set(PLACEHOLDER_WALL_S) | set(PLACEHOLDER_RSS_MB))
    warnings.warn(
        "PLACEHOLDER axes in use for: "
        f"{', '.join(placeholder_tools)}. Refresh after Phase 4 runs land.",
        stacklevel=1,
    )
    print(
        "WARNING: PLACEHOLDER axes for "
        f"{', '.join(placeholder_tools)} — see top of script.",
        file=sys.stderr,
    )

    labels, polygons = _build_polygons(df)

    width, height = double_col_size(height_in=4.5)
    fig, ax = plt.subplots(figsize=(width, height),
                           subplot_kw=dict(projection="polar"))
    angles = spider_axes(ax, labels)

    # Class polygons get distinct colours (NOT TOOL_COLORS — these are
    # aggregate philosophies, not single tools). gfc keeps its anchored
    # black so the paper-subject signal stays consistent across figures.
    color_map = {
        "gfc":         PALETTE["black"],
        "Domain-deep": PALETTE["sky_blue"],
        "Wrappers":    PALETTE["vermillion"],
    }
    plot_spider(ax, polygons, angles, color_map=color_map)

    ax.legend(loc="upper right", bbox_to_anchor=(1.30, 1.10),
              frameon=False, fontsize=8)
    ax.set_title(
        "Architecture spider — gfc vs domain-deep vs wrapper aggregates",
        fontsize=9, pad=14)

    save_figure(fig, "architecture_spider")
    plt.close(fig)

    _write_caption(FIG_DIR / "architecture_spider_caption.txt")
    print(f"wrote {FIG_DIR / 'architecture_spider.png'}")
    print(f"wrote {FIG_DIR / 'architecture_spider.pdf'}")
    print(f"wrote {FIG_DIR / 'architecture_spider_caption.txt'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
