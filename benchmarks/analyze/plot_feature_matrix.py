#!/usr/bin/env python3
# Author: Benjamin Narh-Madey
"""Figure 2 — feature coverage matrix (gfc vs. common bioinformatics CLIs).

Static, hand-curated data — reviewed and updated before paper
submission. Renders both a markdown table (for paper appendix /
README; plain-text glyphs Y / ~ / . in place of emoji) and a heatmap
figure for the paper itself.

Style: shared `_style` rcParams (sans-serif, Bioinformatics point
sizes, 300 dpi save). Heatmap uses a 3-step blue palette with
explicit vmin/vmax so the three categories (none, scripted,
first-class) read as three distinct cell shades at print resolution.
"""
from pathlib import Path
import sys

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from _style import (
    save_figure,
    set_bioinformatics_style,
)

FIG_DIR = Path(__file__).resolve().parent.parent / "results" / "figures"

# Rows = tasks (broadly grouped). Cols = tools. Value: 1 = first-class,
# 0.5 = possible with scripting, 0 = not supported. Edit before submission.
COVERAGE: dict = {
    #                                     gfc  seqkit  samtools  bcftools  AGAT  gffread  plink2  convertf  EMBOSS  UCSC   pileupCaller
    "FASTA <-> FASTQ":                   [ 1.0, 1.0,    0,        0,        0,    0,       0,      0,        0.5,    0,     0 ],
    "FASTA -> table":                    [ 1.0, 1.0,    0,        0,        0,    0,       0,      0,        0.5,    0,     0 ],
    "GFF3 -> GTF":                       [ 1.0, 0,      0,        0,        1.0,  1.0,     0,      0,        0,      0.5,   0 ],
    "GTF -> GFF3":                       [ 1.0, 0,      0,        0,        1.0,  0.5,     0,      0,        0,      0,     0 ],
    "GFF3 -> BED6":                      [ 1.0, 0,      0,        0,        1.0,  0,       0,      0,        0,      0.5,   0 ],
    "GFF3 -> BED12":                     [ 1.0, 0,      0,        0,        1.0,  0.5,     0,      0,        0,      1.0,   0 ],
    "GFF3 -> protein FASTA":             [ 1.0, 0,      0,        0,        1.0,  1.0,     0,      0,        1.0,    0,     0 ],
    "FASTA + GFF -> GenBank":            [ 1.0, 0,      0,        0,        0,    0,       0,      0,        1.0,    0,     0 ],
    "GenBank -> GFF3":                   [ 1.0, 0,      0,        0,        0,    0,       0,      0,        1.0,    0,     0 ],
    "BAM/SAM -> BED":                    [ 1.0, 0,      0.5,      0,        0,    0,       0,      0,        0,      0,     0 ],
    "VCF -> BED":                        [ 1.0, 0,      0,        0.5,      0,    0,       0,      0,        0,      0,     0 ],
    "VCF -> table (TSV)":                [ 1.0, 0,      0,        1.0,      0,    0,       0,      0,        0,      0,     0 ],
    "VCF -> consensus FASTA":            [ 1.0, 0,      0,        1.0,      0,    0,       0,      0,        0,      0,     0 ],
    "VCF -> EIGENSTRAT":                 [ 1.0, 0,      0,        0.5,      0,    0,       0,      1.0,      0,      0,     0.5 ],
    "VCF -> pseudohaploid EIGENSTRAT":   [ 1.0, 0,      0,        0.5,      0,    0,       0,      0,        0,      0,     1.0 ],
    "VCF -> PLINK binary":               [ 1.0, 0,      0,        0.5,      0,    0,       1.0,    0,        0,      0,     0 ],
    "MAF -> XMFA":                       [ 1.0, 0,      0,        0,        0,    0,       0,      0,        0,      0,     0 ],
    "MUMmer delta -> tabular":           [ 1.0, 0,      0,        0,        0,    0,       0,      0,        0,      0,     0 ],
    "BLAST outfmt6 -> link TSV":         [ 1.0, 0,      0,        0,        0,    0,       0,      0,        0,      0,     0 ],
    "Stockholm <-> FASTA alignment":     [ 1.0, 0,      0,        0,        0,    0,       0,      0,        0.5,    0,     0 ],
    "HMMER tblout -> TSV":               [ 1.0, 0,      0,        0,        0,    0,       0,      0,        0,      0,     0 ],
    "OrthoFinder Orthogroups -> FASTAs": [ 1.0, 0,      0,        0,        0,    0,       0,      0,        0,      0,     0 ],
    "Tree Newick <-> NEXUS <-> PhyloXML": [ 1.0, 0,      0,        0,        0,    0,       0,      0,        0,      0,     0 ],
    "Tree + alignment -> annotated NEX": [ 1.0, 0,      0,        0,        0,    0,       0,      0,        0,      0,     0 ],
}

COLUMNS = ["gfc", "seqkit", "samtools", "bcftools", "AGAT", "gffread",
           "plink2", "convertf", "EMBOSS", "UCSC", "pileupCaller"]


def main() -> int:
    df = pd.DataFrame.from_dict(COVERAGE, orient="index", columns=COLUMNS)

    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # Markdown table for paper appendix / README. Plain-text glyphs:
    #   Y = first-class support, ~ = possible with scripting, . = none.
    md_lines = ["| Task | " + " | ".join(COLUMNS) + " |",
                "|" + "|".join(["---"] * (len(COLUMNS) + 1)) + "|"]
    for task, row in df.iterrows():
        cells = []
        for v in row.tolist():
            if v == 1.0:
                cells.append("Y")
            elif v == 0.5:
                cells.append("~")
            else:
                cells.append(".")
        md_lines.append(f"| {task} | " + " | ".join(cells) + " |")
    md_path = FIG_DIR / "feature_matrix.md"
    md_path.write_text("\n".join(md_lines) + "\n")

    # Heatmap figure. Discrete 3-step palette (white -> mid blue ->
    # solid blue) with explicit vmin/vmax so the three categories
    # render as three visually distinct cells at print resolution.
    set_bioinformatics_style()
    fig, ax = plt.subplots(figsize=(180.0 / 25.4, 7.0))
    cmap = sns.color_palette(["#FFFFFF", "#9ECAE1", "#08519C"], n_colors=3)
    sns.heatmap(df, ax=ax, cmap=cmap, vmin=0, vmax=1.0,
                cbar_kws={"label": "support (0 = none, 0.5 = scripted, 1 = first-class)",
                          "ticks": [0.0, 0.5, 1.0]},
                linewidths=0.5, linecolor="lightgrey",
                square=False)
    ax.set_title("Feature coverage: gfc vs. common bioinformatics CLIs")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    plt.setp(ax.get_yticklabels(), rotation=0)
    fig.tight_layout()
    save_figure(fig, "feature_matrix")
    print(f"[done] wrote {md_path} + feature_matrix.png/.pdf", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
