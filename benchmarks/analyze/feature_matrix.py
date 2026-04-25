#!/usr/bin/env python3
"""Figure 2 — feature coverage matrix (gfc vs. common bioinformatics CLIs).

Static data — intended to be reviewed by hand and updated before paper
submission. Renders both a markdown table (for README / appendix) and a
heatmap figure (for the paper).
"""
from pathlib import Path
import sys

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

FIG = Path(__file__).resolve().parent.parent / "results" / "figures"

# Rows = tasks (broadly grouped). Cols = tools. Value: 1 = first-class,
# 0.5 = possible with scripting, 0 = not supported. Edit before submission.
COVERAGE: dict = {
    #                                     gfc  seqkit  samtools  bcftools  AGAT  gffread  plink2  convertf  EMBOSS  UCSC   pileupCaller
    "FASTA ↔ FASTQ":                    [ 1.0, 1.0,    0,        0,        0,    0,       0,      0,        0.5,    0,     0 ],
    "FASTA → table":                     [ 1.0, 1.0,    0,        0,        0,    0,       0,      0,        0.5,    0,     0 ],
    "GFF3 → GTF":                        [ 1.0, 0,      0,        0,        1.0,  1.0,     0,      0,        0,      0.5,   0 ],
    "GTF → GFF3":                        [ 1.0, 0,      0,        0,        1.0,  0.5,     0,      0,        0,      0,     0 ],
    "GFF3 → BED6":                       [ 1.0, 0,      0,        0,        1.0,  0,       0,      0,        0,      0.5,   0 ],
    "GFF3 → BED12":                      [ 1.0, 0,      0,        0,        1.0,  0.5,     0,      0,        0,      1.0,   0 ],
    "GFF3 → protein FASTA":              [ 1.0, 0,      0,        0,        1.0,  1.0,     0,      0,        1.0,    0,     0 ],
    "FASTA + GFF → GenBank":             [ 1.0, 0,      0,        0,        0,    0,       0,      0,        1.0,    0,     0 ],
    "GenBank → GFF3":                    [ 1.0, 0,      0,        0,        0,    0,       0,      0,        1.0,    0,     0 ],
    "BAM/SAM → BED":                     [ 1.0, 0,      0.5,      0,        0,    0,       0,      0,        0,      0,     0 ],
    "VCF → BED":                         [ 1.0, 0,      0,        0.5,      0,    0,       0,      0,        0,      0,     0 ],
    "VCF → table (TSV)":                 [ 1.0, 0,      0,        1.0,      0,    0,       0,      0,        0,      0,     0 ],
    "VCF → consensus FASTA":             [ 1.0, 0,      0,        1.0,      0,    0,       0,      0,        0,      0,     0 ],
    "VCF → EIGENSTRAT":                  [ 1.0, 0,      0,        0.5,      0,    0,       0,      1.0,      0,      0,     0.5 ],
    "VCF → pseudohaploid EIGENSTRAT":    [ 1.0, 0,      0,        0.5,      0,    0,       0,      0,        0,      0,     1.0 ],
    "VCF → PLINK binary":                [ 1.0, 0,      0,        0.5,      0,    0,       1.0,    0,        0,      0,     0 ],
    "MAF → XMFA":                        [ 1.0, 0,      0,        0,        0,    0,       0,      0,        0,      0,     0 ],
    "MUMmer delta → tabular":            [ 1.0, 0,      0,        0,        0,    0,       0,      0,        0,      0,     0 ],
    "BLAST outfmt6 → link TSV":          [ 1.0, 0,      0,        0,        0,    0,       0,      0,        0,      0,     0 ],
    "Stockholm ↔ FASTA alignment":       [ 1.0, 0,      0,        0,        0,    0,       0,      0,        0.5,    0,     0 ],
    "HMMER tblout → TSV":                [ 1.0, 0,      0,        0,        0,    0,       0,      0,        0,      0,     0 ],
    "OrthoFinder Orthogroups → FASTAs":  [ 1.0, 0,      0,        0,        0,    0,       0,      0,        0,      0,     0 ],
    "Tree Newick ↔ NEXUS ↔ PhyloXML":    [ 1.0, 0,      0,        0,        0,    0,       0,      0,        0,      0,     0 ],
    "Tree + alignment → annotated NEX":  [ 1.0, 0,      0,        0,        0,    0,       0,      0,        0,      0,     0 ],
}

COLUMNS = ["gfc", "seqkit", "samtools", "bcftools", "AGAT", "gffread",
           "plink2", "convertf", "EMBOSS", "UCSC", "pileupCaller"]


def main() -> int:
    df = pd.DataFrame.from_dict(COVERAGE, orient="index", columns=COLUMNS)

    FIG.mkdir(parents=True, exist_ok=True)

    # Markdown table, for paper appendix / README.
    md_lines = ["| Task | " + " | ".join(COLUMNS) + " |",
                "|" + "|".join(["---"] * (len(COLUMNS) + 1)) + "|"]
    for task, row in df.iterrows():
        cells = []
        for v in row.tolist():
            if v == 1.0:
                cells.append("✅")
            elif v == 0.5:
                cells.append("◐")
            else:
                cells.append("·")
        md_lines.append(f"| {task} | " + " | ".join(cells) + " |")
    md_path = FIG / "feature_matrix.md"
    md_path.write_text("\n".join(md_lines) + "\n")

    # Heatmap for the paper figure.
    sns.set_theme(style="whitegrid", context="paper")
    fig, ax = plt.subplots(figsize=(10, 8))
    cmap = sns.color_palette("Blues", n_colors=3)
    sns.heatmap(df, ax=ax, cmap=cmap, vmin=0, vmax=1.0,
                cbar_kws={"label": "support (0 = none, 0.5 = scripted, 1 = first-class)"},
                linewidths=0.5, linecolor="lightgrey")
    ax.set_title("Feature coverage — gfc vs. common bioinformatics CLIs")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"feature_matrix.{ext}", dpi=300, bbox_inches="tight")
    print(f"[done] wrote {md_path} + feature_matrix.png/.pdf", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
