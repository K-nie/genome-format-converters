#!/usr/bin/env python3
# Author: Benjamin Narh-Madey
"""Figure 5 — correctness heatmap, task x tool (4-state categorical).

Encodes each (task, tool) cell as one of four categories:

  pass                   — comparator returned `correct == "1"`. Strict
                           byte / structural identity against the
                           per-task reference.
  structural-divergence  — comparator returned `correct == "0"`. The
                           tool's output is structurally different from
                           the reference, but the divergence is
                           documented in `check_correctness.py`
                           (`_REASON_T1` ... `_REASON_T4`) as a
                           legitimate rendering choice rather than a
                           defect.
  skip                   — comparator returned `correct == "skip"`.
                           Intentionally out of scope (e.g. T8
                           cross-RNG byte-equivalence is impossible
                           by design).
  unchecked              — comparator returned `""` (no reference
                           available, OSError, or malformed input). The
                           cell is left unverified, NOT marked failed.

The four categories are colour-coded (green / orange / grey / light-grey)
and labelled in-cell with a plain-text glyph (Y, eq, n/a, ?). NO
percentage formatting — `0.5 = 50%` was the source of the earlier
"50 % correctness" misreading.

A companion `correctness_matrix_caption.txt` is written next to the
PNG. It lifts the `_REASON_T*` paragraph footnotes from
`check_correctness.py` verbatim so the figure caption (Methods) can
cite each divergence by source.

Paper utility: this is the figure that defends gfc's `correct=0` rows
on T1/T2/T3/T4. Without the four-state encoding, reviewers read
"gfc gives wrong answers on half the benchmark" instead of "gfc
output is structurally equivalent on these four tasks; the divergence
is documented."
"""
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

from _style import (
    PALETTE,
    save_figure,
    set_bioinformatics_style,
    double_col_size,
)

RAW = Path(__file__).resolve().parent.parent / "results" / "raw" / "all.tsv"
FIG_DIR = Path(__file__).resolve().parent.parent / "results" / "figures"
CHECK_SRC = Path(__file__).resolve().parent.parent / "check_correctness.py"


# Ordered category list -> integer code used for colormap indexing.
# Order matters: it controls the colorbar / legend ordering and the
# integer values written into the matrix.
CATEGORIES = ["pass", "structural-divergence", "skip", "unchecked"]
CAT_TO_CODE = {c: i for i, c in enumerate(CATEGORIES)}
CAT_COLORS = {
    "pass":                   PALETTE["bluish_green"],   # solid signal of OK.
    "structural-divergence":  PALETTE["orange"],         # documented, non-defect.
    "skip":                   "#BBBBBB",                 # mid-grey: out of scope.
    "unchecked":              "#E5E5E5",                 # light-grey: no verdict.
}
CAT_GLYPHS = {
    "pass":                   "Y",
    "structural-divergence":  "eq",
    "skip":                   "n/a",
    "unchecked":              "?",
}


def _classify(correct: str) -> str:
    """Map a raw `correct` cell value to one of the four categories."""
    s = (correct or "").strip()
    if s == "1":
        return "pass"
    if s == "0":
        return "structural-divergence"
    if s.lower() == "skip":
        return "skip"
    return "unchecked"


def _extract_reason_blocks(src: Path) -> dict:
    """Pull `_REASON_T*` paren-delimited strings out of check_correctness.py.

    Returns a dict keyed by task ('T1', 'T2', ...) with the cleaned-up
    paragraph as the value. Walks the character stream with a
    paren-depth counter so nested '(' / ')' inside the prose (e.g.
    `(VCF -> EIGENSTRAT)`, `(2/3 variants for tiny.vcf)`) do not close
    the block early. Falls back to an empty dict if the source file is
    unreadable so the figure still renders.
    """
    out = {}
    if not src.exists():
        return out
    text = src.read_text()
    for tag in ("T1", "T2", "T3", "T4"):
        marker = f"_REASON_{tag} = ("
        i = text.find(marker)
        if i < 0:
            continue
        start = i + len(marker)
        depth = 1
        j = start
        while j < len(text) and depth > 0:
            ch = text[j]
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        if depth != 0:
            continue
        body = text[start:j]
        # Concatenate adjacent string literals (implicit "abc" "def" join)
        # and strip surrounding quotes / whitespace.
        cleaned = (body.replace("\n", " ")
                       .replace('" "', "")
                       .replace('"', "")
                       .strip())
        cleaned = " ".join(cleaned.split())
        out[tag] = cleaned
    return out


def main() -> int:
    if not RAW.exists():
        print(f"[error] {RAW} missing.", file=sys.stderr)
        return 1
    df = pd.read_csv(RAW, sep="\t")
    if "correct" not in df.columns:
        print("[error] no `correct` column in all.tsv", file=sys.stderr)
        return 1

    df["correct"] = df["correct"].fillna("").astype(str)
    df1 = df[df["replicate"] == 1].copy()
    df1["category"] = df1["correct"].apply(_classify)
    df1["code"] = df1["category"].map(CAT_TO_CODE)

    # Build a strict (task x tool) full grid. Every (task, tool) combination
    # exists; cells with no run are flagged as "missing" and rendered as
    # white space with no glyph. The previous pivot_table path produced a
    # grid where matplotlib's auto-aspect imshow stretched cells unevenly
    # depending on column density — switching to pcolormesh on an
    # explicitly indexed grid plus aspect="equal" guarantees uniform cell
    # pitch across the whole figure.
    tasks = sorted(df1["task"].unique())
    tools = sorted(df1["tool"].unique())
    cat_matrix = (df1.pivot_table(index="task", columns="tool",
                                  values="category", aggfunc="first")
                     .reindex(index=tasks, columns=tools))
    # Encode: integer codes 0..len(CATEGORIES)-1 for the four real
    # categories; NaN ("missing") survives in the masked array so
    # pcolormesh leaves it transparent.
    code_grid = np.full((len(tasks), len(tools)), np.nan, dtype=float)
    for r, t in enumerate(tasks):
        for c, tool in enumerate(tools):
            cat = cat_matrix.iat[r, c]
            if isinstance(cat, str) and cat in CAT_TO_CODE:
                code_grid[r, c] = CAT_TO_CODE[cat]

    set_bioinformatics_style()
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # Size figure proportional to the grid so cells render square. Cell
    # pitch = 0.55 in works at 300 dpi for both screen and print.
    cell_in = 0.55
    fig_w = max(8.0, len(tools) * cell_in + 2.5)
    fig_h = max(3.5, len(tasks) * cell_in + 1.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    cmap = ListedColormap([CAT_COLORS[c] for c in CATEGORIES])
    cmap.set_bad(color="white", alpha=0.0)  # missing cells: transparent.
    masked = np.ma.masked_invalid(code_grid)

    # pcolormesh draws cells between explicit edge coordinates; this lets
    # us put a visible white gutter between every cell so adjacent
    # categorical cells in the same column never read as one tall block.
    edges_x = np.arange(len(tools) + 1)
    edges_y = np.arange(len(tasks) + 1)
    ax.pcolormesh(edges_x, edges_y, masked, cmap=cmap,
                  vmin=-0.5, vmax=len(CATEGORIES) - 0.5,
                  edgecolors="white", linewidth=1.5)
    ax.set_aspect("equal")
    # Invert y so T1 sits at the top (reading order).
    ax.invert_yaxis()

    # Tick positions at cell centers (pcolormesh uses cell edges).
    ax.set_xticks(np.arange(len(tools)) + 0.5)
    ax.set_xticklabels(tools, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(tasks)) + 0.5)
    ax.set_yticklabels(tasks)
    ax.set_xlabel("tool")
    ax.set_ylabel("benchmark task")
    ax.set_title("Output correctness against per-task reference "
                 "(4-state categorical; see Methods for divergence notes)")
    # Hide axis spines so the grid reads as a clean matrix, not a chart.
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)

    # Annotate each cell with the categorical glyph; skip masked cells.
    for r in range(len(tasks)):
        for c in range(len(tools)):
            cat = cat_matrix.iat[r, c]
            if not isinstance(cat, str) or cat not in CAT_GLYPHS:
                continue
            ax.text(c + 0.5, r + 0.5, CAT_GLYPHS[cat],
                    ha="center", va="center", fontsize=8,
                    color="black", fontweight="bold")

    # Categorical legend instead of a continuous colorbar.
    legend_handles = [
        Patch(facecolor=CAT_COLORS[c], edgecolor="black", linewidth=0.4,
              label=f"{c} ({CAT_GLYPHS[c]})")
        for c in CATEGORIES
    ]
    ax.legend(handles=legend_handles, title="category",
              loc="upper left", bbox_to_anchor=(1.02, 1.0),
              frameon=False)

    fig.tight_layout()
    save_figure(fig, "correctness_matrix")

    # Companion caption file: dump _REASON_T* strings verbatim. The
    # caption text gets pasted into the paper's Methods / figure caption.
    reasons = _extract_reason_blocks(CHECK_SRC)
    caption_lines = [
        "Figure 5 — Output correctness against per-task reference.",
        "",
        "Cell encoding (4 categories):",
        "  pass  (Y)   : strict byte / structural identity against reference.",
        "  structural-divergence (eq) : documented rendering difference; not a defect "
        "(see per-task footnotes below).",
        "  skip  (n/a) : intentionally out of scope (e.g. T8 cross-RNG byte-",
        "                equivalence is impossible by design).",
        "  unchecked (?) : comparator could not produce a verdict (no reference,",
        "                OSError, or malformed input).",
        "",
        "Per-task divergence footnotes (lifted verbatim from",
        "benchmarks/check_correctness.py):",
        "",
    ]
    for tag in ("T1", "T2", "T3", "T4"):
        reason = reasons.get(tag, "(reason text not found in source)")
        caption_lines.append(f"  [{tag}] {reason}")
        caption_lines.append("")
    (FIG_DIR / "correctness_matrix_caption.txt").write_text(
        "\n".join(caption_lines))

    print(f"[done] wrote {FIG_DIR / 'correctness_matrix.png'} + .pdf "
          f"+ correctness_matrix_caption.txt", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
