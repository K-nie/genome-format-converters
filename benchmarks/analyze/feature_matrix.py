#!/usr/bin/env python3
"""Figure 1: what gfc converts, and which other tools cover the same ground.

Panel a draws gfc's 32 conversions as a format graph, one box per domain.
Panel b is a coverage matrix: for each gfc conversion, does a commonly used
tool do the same job?

Ratings
    1    one documented command, or a documented pipeline inside the same suite
    0.5  possible only with a second tool, or with a major restriction or loss
    0    no documented route

Every non-zero cell carries the command or the restriction that justifies it
(EVIDENCE below), and the script writes them to coverage_evidence.tsv so the
matrix can be checked cell by cell. Ratings were checked against each tool's
documentation, and the EMBOSS, easel and BioConvert cells were also run on
small inputs (EMBOSS 6.6.0, easel 0.49, BioConvert 1.2.0, plink2 2.0,
bcftools 1.21, AGAT 1.4, gffread 0.12).

Usage: feature_matrix.py [--out DIR]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Wedge, Circle

HERE = Path(__file__).resolve().parent
DEFAULT_OUT = HERE.parent / "results" / "figures"

INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#8a8984"
GRID = "#e4e3df"
GFC = "#2a78d6"
BOX = "#f6f5f2"

plt.rcParams.update({
    "font.size": 7, "axes.titlesize": 8.5, "text.color": INK,
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 600,
})

TOOLS = ["gfc", "BioConvert", "seqkit", "EMBOSS", "samtools", "bedtools",
         "bcftools", "plink2", "convertf", "AGAT", "gffread", "UCSC", "easel"]

# (domain, row label, gfc subcommand)
ROWS = [
    ("Sequences & reads", "FASTA → FASTQ", "fasta-to-fastq"),
    ("Sequences & reads", "FASTQ → FASTA", "fastq-to-fasta"),
    ("Sequences & reads", "FASTA + QUAL → FASTQ", "fasta-qual-to-fastq"),
    ("Sequences & reads", "FASTQ → FASTA + QUAL", "fastq-to-fasta-qual"),
    ("Sequences & reads", "FASTA → TSV", "fasta-to-table"),
    ("Sequences & reads", "BAM/SAM/CRAM → BED6", "bam-to-bed"),
    ("Annotation", "GFF3 → GTF", "gff3-to-gtf"),
    ("Annotation", "GTF → GFF3", "gtf-to-gff3"),
    ("Annotation", "GFF3 → BED6", "gff3-to-bed"),
    ("Annotation", "GFF3 → BED12", "gff3-to-bed12"),
    ("Annotation", "GFF3 → TSV", "gff3-to-table"),
    ("Annotation", "GFF3 + FASTA → protein", "gff3-to-protein"),
    ("Annotation", "GFF3 + FASTA → GenBank", "fasta-gff-to-gbk"),
    ("Annotation", "GenBank → GFF3", "genbank-to-gff3"),
    ("Alignments", "FASTA/PHYLIP/NEXUS/Clustal", "convert-alignment"),
    ("Alignments", "Stockholm → FASTA", "stockholm-to-fasta"),
    ("Alignments", "FASTA → Stockholm", "fasta-to-stockholm"),
    ("Comparative & search", "BLAST outfmt 6 → links", "blast-to-links"),
    ("Comparative & search", "MUMmer delta → TSV", "delta-to-tab"),
    ("Comparative & search", "MAF → XMFA", "maf-to-xmfa"),
    ("Comparative & search", "HMMER tblout → TSV", "hmmer-tblout-to-tsv"),
    ("Comparative & search", "Orthogroups → FASTA", "orthogroups-to-fasta"),
    ("Variants", "VCF → BED", "vcf-to-bed"),
    ("Variants", "VCF → TSV", "vcf-to-table"),
    ("Variants", "VCF + FASTA → consensus", "vcf-to-consensus"),
    ("Variants", "VCF → EIGENSTRAT", "vcf-to-eigenstrat"),
    ("Variants", "VCF → pseudohaploid EIG.", "vcf-to-pseudohaploid"),
    ("Variants", "VCF → PLINK", "vcf-to-plink"),
    ("Variants", "VCF → SNP alignment", "vcf-to-snp-alignment"),
    ("Variants", "EIGENSTRAT → VCF", "eigenstrat-to-vcf"),
    ("Phylogenetics", "Newick/NEXUS/PhyloXML", "tree-convert"),
    ("Phylogenetics", "Tree + alignment → NEXUS", "annotate-tree"),
]

# {(subcommand, tool): (rating, evidence)}. gfc is 1 everywhere; cells not
# listed are 0.
EVIDENCE: dict[tuple[str, str], tuple[float, str]] = {
    # sequences & reads
    ("fasta-to-fastq", "BioConvert"): (1, "bioconvert fasta2fastq"),
    ("fasta-to-fastq", "EMBOSS"): (1, "seqret -osformat fastq (tested; fixed quality)"),
    ("fastq-to-fasta", "BioConvert"): (1, "bioconvert fastq2fasta"),
    ("fastq-to-fasta", "seqkit"): (1, "seqkit fq2fa"),
    ("fastq-to-fasta", "EMBOSS"): (1, "seqret -sformat fastq -osformat fasta"),
    ("fastq-to-fasta", "samtools"): (1, "samtools import | samtools fasta"),
    ("fasta-qual-to-fastq", "BioConvert"): (1, "bioconvert fasta_qual2fastq"),
    ("fastq-to-fasta-qual", "BioConvert"): (1, "bioconvert fastq2fasta_qual"),
    ("fasta-to-table", "seqkit"): (1, "seqkit fx2tab"),
    ("bam-to-bed", "bedtools"): (1, "bedtools bamtobed (BAM input only)"),
    ("bam-to-bed", "samtools"): (0.5, "no BED writer; samtools view piped to awk"),
    # annotation
    ("gff3-to-gtf", "BioConvert"): (1, "bioconvert gff32gtf"),
    ("gff3-to-gtf", "AGAT"): (1, "agat_convert_sp_gff2gtf.pl"),
    ("gff3-to-gtf", "gffread"): (1, "gffread -T"),
    ("gff3-to-gtf", "UCSC"): (1, "gff3ToGenePred | genePredToGtf"),
    ("gtf-to-gff3", "AGAT"): (1, "agat_convert_sp_gxf2gxf.pl"),
    ("gtf-to-gff3", "gffread"): (1, "gffread (GFF3 is the default output)"),
    ("gff3-to-bed", "AGAT"): (0.5, "agat_convert_sp_gff2bed.pl writes BED12; cut to BED6"),
    ("gff3-to-bed", "gffread"): (0.5, "gffread --bed writes BED12; cut to BED6"),
    ("gff3-to-bed", "UCSC"): (0.5, "gff3ToGenePred | genePredToBed gives BED12 only"),
    ("gff3-to-bed12", "AGAT"): (1, "agat_convert_sp_gff2bed.pl"),
    ("gff3-to-bed12", "gffread"): (1, "gffread --bed"),
    ("gff3-to-bed12", "UCSC"): (1, "gff3ToGenePred | genePredToBed"),
    ("gff3-to-table", "AGAT"): (1, "agat_convert_sp_gff2tsv.pl"),
    ("gff3-to-table", "gffread"): (1, "gffread --table"),
    ("gff3-to-protein", "AGAT"): (1, "agat_sp_extract_sequences.pl -p"),
    ("gff3-to-protein", "gffread"): (1, "gffread -y"),
    ("gff3-to-protein", "EMBOSS"): (0.5, "extractfeat -join then transeq; multi-step, frame by hand"),
    ("fasta-gff-to-gbk", "EMBOSS"): (1, "seqret -feature -osformat genbank (benchmarked, T3)"),
    ("genbank-to-gff3", "BioConvert"): (1, "bioconvert genbank2gff3"),
    ("genbank-to-gff3", "EMBOSS"): (0.5, "seqret -feature -offormat gff3 documented; wrote no "
                                         "features or crashed in our tests (6.6.0)"),
    # alignments
    ("convert-alignment", "BioConvert"): (1, "bioconvert fasta2phylip, nexus2phylip, clustal2fasta, ..."),
    ("convert-alignment", "EMBOSS"): (1, "seqret -osformat phylip/nexus/clustal/fasta"),
    ("convert-alignment", "easel"): (0.5, "esl-reformat afa/phylip/clustal; no NEXUS"),
    ("stockholm-to-fasta", "BioConvert"): (1, "bioconvert stockholm2clustal then clustal2fasta"),
    ("stockholm-to-fasta", "EMBOSS"): (1, "seqret -sformat stockholm (tested)"),
    ("stockholm-to-fasta", "easel"): (1, "esl-reformat afa"),
    ("fasta-to-stockholm", "BioConvert"): (1, "bioconvert fasta2clustal then clustal2stockholm"),
    ("fasta-to-stockholm", "easel"): (1, "esl-reformat stockholm"),
    # comparative & search
    ("orthogroups-to-fasta", "seqkit"): (0.5, "seqkit grep -f per orthogroup, in a shell loop"),
    # variants
    ("vcf-to-bed", "BioConvert"): (1, "bioconvert vcf2bed"),
    ("vcf-to-bed", "bcftools"): (1, "bcftools query -f '%CHROM\\t%POS0\\t%END\\n'"),
    ("vcf-to-table", "bcftools"): (1, "bcftools query -H -f ..."),
    ("vcf-to-consensus", "bcftools"): (1, "bcftools consensus -s SAMPLE -f ref.fa"),
    ("vcf-to-eigenstrat", "plink2"): (1, "plink2 --vcf ... --export eig"),
    ("vcf-to-eigenstrat", "convertf"): (0.5, "no VCF input; needs plink to make PED/BED first (T1)"),
    ("vcf-to-plink", "BioConvert"): (1, "bioconvert vcf2bplink (wraps plink)"),
    ("vcf-to-plink", "plink2"): (1, "plink2 --vcf ... --make-bed"),
    ("vcf-to-snp-alignment", "plink2"): (0.5, "plink2 --export phylip; PHYLIP only, IUPAC hets"),
    ("eigenstrat-to-vcf", "plink2"): (1, "plink2 --eiggeno/--eigind/--eigsnp ... --export vcf"),
    ("eigenstrat-to-vcf", "convertf"): (0.5, "convertf to PED, then plink --recode vcf"),
    # phylogenetics
    ("tree-convert", "BioConvert"): (1, "bioconvert newick2nexus, nexus2phyloxml, ..."),
}

DOMAIN_ORDER = ["Sequences & reads", "Annotation", "Alignments",
                "Comparative & search", "Variants", "Phylogenetics"]


def rating(cmd: str, tool: str) -> float:
    if tool == "gfc":
        return 1.0
    return float(EVIDENCE.get((cmd, tool), (0, ""))[0])


# ------------------------------------------------------------ panel a

# One box per domain: nodes at local (x, y) in a 3 x 2.4 box, and edges
# (src, dst, both_ways, label).
GRAPH = {
    "Sequences & reads": (
        {"FASTA": (0.5, 2.0), "FASTQ": (2.5, 2.0), "FASTA+QUAL": (2.5, 1.0),
         "TSV": (0.5, 1.0), "BAM/SAM/CRAM": (0.5, 0.15), "BED6": (2.5, 0.15)},
        [("FASTA", "FASTQ", True, ""), ("FASTA+QUAL", "FASTQ", True, ""),
         ("FASTA", "TSV", False, ""), ("BAM/SAM/CRAM", "BED6", False, "")]),
    "Annotation": (
        {"GFF3": (1.4, 1.1), "GTF": (2.6, 2.05), "BED6": (2.75, 1.4),
         "BED12": (2.75, 0.8), "TSV": (2.6, 0.15), "GenBank": (0.25, 2.05),
         "protein": (0.25, 0.15)},
        [("GFF3", "GTF", True, ""), ("GFF3", "BED6", False, ""),
         ("GFF3", "BED12", False, ""), ("GFF3", "TSV", False, ""),
         ("GFF3", "GenBank", True, "+FASTA"), ("GFF3", "protein", False, "+FASTA")]),
    "Alignments": (
        {"Clustal": (0.6, 1.95), "PHYLIP": (2.4, 1.95), "FASTA": (0.6, 1.2),
         "NEXUS": (2.4, 1.2), "Stockholm": (0.6, 0.05)},
        [("Stockholm", "FASTA", True, "")]),
    "Comparative & search": (
        {"BLAST": (0.5, 2.15), "links": (2.5, 2.15), "delta": (0.5, 1.6),
         "coords TSV": (2.5, 1.6), "MAF": (0.5, 1.05), "XMFA": (2.5, 1.05),
         "tblout": (0.5, 0.5), "hits TSV": (2.5, 0.5),
         "orthogroups": (0.5, -0.05), "per-OG FASTA": (2.5, -0.05)},
        [("BLAST", "links", False, ""), ("delta", "coords TSV", False, ""),
         ("MAF", "XMFA", False, ""), ("tblout", "hits TSV", False, ""),
         ("orthogroups", "per-OG FASTA", False, "+FASTAs")]),
    "Variants": (
        {"VCF": (1.35, 1.1), "BED": (0.2, 2.05), "TSV": (1.35, 2.15),
         "consensus": (2.6, 2.05), "EIGENSTRAT": (2.75, 1.1),
         "pseudohaploid": (2.6, 0.15), "PLINK": (1.35, 0.05),
         "SNP aln": (0.2, 0.15)},
        [("VCF", "BED", False, ""), ("VCF", "TSV", False, ""),
         ("VCF", "consensus", False, "+FASTA"), ("VCF", "EIGENSTRAT", True, ""),
         ("VCF", "pseudohaploid", False, ""), ("VCF", "PLINK", False, ""),
         ("VCF", "SNP aln", False, "")]),
    "Phylogenetics": (
        {"Newick": (0.5, 2.0), "NEXUS": (2.5, 2.0), "PhyloXML": (1.5, 1.1),
         "tree": (0.5, 0.15), "NEXUS+aln": (2.5, 0.15)},
        [("Newick", "NEXUS", True, ""), ("Newick", "PhyloXML", True, ""),
         ("NEXUS", "PhyloXML", True, ""), ("tree", "NEXUS+aln", False, "+alignment")]),
}


def _node(ax, x, y, text):
    ax.text(x, y, text, ha="center", va="center", fontsize=6.3, color=INK,
            zorder=3, bbox=dict(boxstyle="round,pad=0.28,rounding_size=0.18",
                                fc="white", ec=GFC, lw=0.7))


def _edge_end(a, b, text):
    """Point where the segment a->b leaves the node box drawn around a."""
    hw, hh = 0.075 + 0.052 * len(text), 0.16
    dx, dy = b[0] - a[0], b[1] - a[1]
    t = min(hw / abs(dx) if dx else 1e9, hh / abs(dy) if dy else 1e9)
    return a[0] + dx * t, a[1] + dy * t


def draw_graph(ax) -> None:
    ax.set_axis_off()
    bw, bh, gx, gy = 3.4, 3.0, 0.35, 0.3
    for i, dom in enumerate(DOMAIN_ORDER):
        col, row = i % 3, 1 - i // 3
        ox, oy = col * (bw + gx), row * (bh + gy)
        ax.add_patch(FancyBboxPatch((ox, oy), bw, bh,
                                    boxstyle="round,pad=0,rounding_size=0.12",
                                    fc=BOX, ec=GRID, lw=0.6, zorder=0))
        n = sum(1 for r in ROWS if r[0] == dom)
        ax.text(ox + 0.12, oy + bh - 0.12, dom, ha="left", va="top",
                fontsize=7.2, fontweight="bold", color=INK)
        ax.text(ox + bw - 0.12, oy + bh - 0.12, f"{n}", ha="right", va="top",
                fontsize=7.2, color=INK_2)
        nodes, edges = GRAPH[dom]
        pos = {k: (ox + 0.2 + x, oy + 0.3 + y * 0.97) for k, (x, y) in nodes.items()}
        if dom == "Alignments":
            # convert-alignment converts any of the four formats to any other
            (x0, y0), (x1, y1) = pos["FASTA"], pos["PHYLIP"]
            ax.add_patch(FancyBboxPatch((x0 - 0.55, y0 - 0.3), x1 - x0 + 1.1,
                                        y1 - y0 + 0.6,
                                        boxstyle="round,pad=0,rounding_size=0.1",
                                        fc="none", ec=GFC, lw=0.6, ls=(0, (3, 2)),
                                        zorder=1))
            ax.text((x0 + x1) / 2, (y0 + y1) / 2, "any ↔ any", ha="center",
                    va="center", fontsize=5.6, color=INK_2, style="italic")
        for src, dst, both, label in edges:
            p1 = _edge_end(pos[src], pos[dst], src)
            p2 = _edge_end(pos[dst], pos[src], dst)
            (x1, y1), (x2, y2) = pos[src], pos[dst]
            ax.add_patch(FancyArrowPatch(
                p1, p2, arrowstyle="<|-|>" if both else "-|>",
                mutation_scale=6, lw=0.7, color=INK_2, zorder=1,
                shrinkA=1.5, shrinkB=1.5))
            if label:
                ax.text((x1 + x2) / 2, (y1 + y2) / 2 + 0.09, label, ha="center",
                        va="bottom", fontsize=5.4, color=INK_2, style="italic",
                        zorder=2)
        for k, (x, y) in pos.items():
            _node(ax, x, y, k)
    ax.set_xlim(-0.05, 3 * bw + 2 * gx + 0.05)
    ax.set_ylim(-0.05, 2 * bh + gy + 0.05)
    ax.set_aspect("equal")


# ------------------------------------------------------------ panel b

def draw_matrix(ax) -> None:
    ax.set_axis_off()
    nt = len(TOOLS)
    y = 0.0
    ys = []
    seps = []
    prev = None
    for dom, label, cmd in ROWS:
        if dom != prev:
            if prev is not None:
                y += 0.5
            seps.append((y, dom))
            y += 1.0
            prev = dom
        ys.append(y)
        y += 1.0
    total_y = y + 0.5

    # gfc column band
    ax.add_patch(FancyBboxPatch((-0.45, -0.2), 0.9, total_y + 0.6,
                                boxstyle="round,pad=0,rounding_size=0.2",
                                fc="#e8f0fb", ec="none", zorder=0))
    for j, tool in enumerate(TOOLS):
        ax.text(j, -0.55, tool, rotation=50, ha="left", va="bottom",
                rotation_mode="anchor", fontsize=6.6,
                color=GFC if tool == "gfc" else INK,
                fontweight="bold" if tool == "gfc" else "normal")
    for sy, dom in seps:
        ax.text(-11.4, sy, dom, ha="left", va="center", fontsize=6.6,
                fontweight="bold", color=INK)
        ax.plot([-11.5, nt - 0.5], [sy - 0.6, sy - 0.6], color=GRID, lw=0.6,
                clip_on=False)
    for (dom, label, cmd), yy in zip(ROWS, ys):
        ax.text(-1.0, yy, label, ha="right", va="center", fontsize=6.3, color=INK)
        others = 0
        anyroute = any(rating(cmd, t) for t in TOOLS[1:])
        for j, tool in enumerate(TOOLS):
            r = rating(cmd, tool)
            col = GFC if tool == "gfc" else INK_2
            if tool != "gfc" and r == 1:
                others += 1
            if r == 1:
                ax.add_patch(Circle((j, yy), 0.3, fc=col, ec=col, lw=0.6, zorder=2))
            elif r == 0.5:
                ax.add_patch(Circle((j, yy), 0.3, fc="white", ec=col, lw=0.6, zorder=2))
                ax.add_patch(Wedge((j, yy), 0.3, 90, 270, fc=col, ec=col, lw=0.6, zorder=3))
            else:
                ax.add_patch(Circle((j, yy), 0.05, fc=GRID, ec="none", zorder=1))
        ax.text(nt - 0.2, yy, str(others) if anyroute else "only gfc",
                ha="left", va="center", fontsize=5.8,
                color=GFC if not anyroute else INK_2,
                fontweight="bold" if not anyroute else "normal")
    ax.text(nt + 0.4, -2.4, "other tools\nwith a full route", ha="left",
            va="bottom", fontsize=5.8, color=INK_2)

    # column totals: full routes (half routes)
    ty = total_y
    ax.plot([-11.5, nt - 0.5], [ty - 0.55, ty - 0.55], color=GRID, lw=0.6,
            clip_on=False)
    ax.text(-1.0, ty, "full routes", ha="right", va="center",
            fontsize=6.3, color=INK)
    ax.text(-1.0, ty + 0.9, "partial routes", ha="right", va="center",
            fontsize=6.3, color=INK_2)
    for j, tool in enumerate(TOOLS):
        full = sum(rating(c, tool) == 1 for _, _, c in ROWS)
        half = sum(rating(c, tool) == 0.5 for _, _, c in ROWS)
        kw = dict(ha="center", va="center", fontsize=6,
                  color=GFC if tool == "gfc" else INK_2,
                  fontweight="bold" if tool == "gfc" else "normal")
        ax.text(j, ty, str(full), **kw)
        ax.text(j, ty + 0.9, str(half) if half else "", **kw)

    # legend
    ly0 = ty + 2.2
    lx = -11.0
    for k, (r, text) in enumerate([(1, "one command, or a pipeline inside the suite"),
                                   (0.5, "needs a second tool, or a major restriction"),
                                   (0, "no documented route")]):
        x, ly = lx, ly0 + k * 0.9
        if r == 1:
            ax.add_patch(Circle((x, ly), 0.3, fc=INK_2, ec=INK_2, clip_on=False))
        elif r == 0.5:
            ax.add_patch(Circle((x, ly), 0.3, fc="white", ec=INK_2, lw=0.6, clip_on=False))
            ax.add_patch(Wedge((x, ly), 0.3, 90, 270, fc=INK_2, ec=INK_2, clip_on=False))
        else:
            ax.add_patch(Circle((x, ly), 0.05, fc=GRID, ec="none", clip_on=False))
        ax.text(x + 0.55, ly, text, ha="left", va="center", fontsize=6.2, color=INK_2)

    ax.set_xlim(-11.5, nt + 2.6)
    ax.set_ylim(ly0 + 2.4, -4.2)
    ax.set_aspect("equal")


# ------------------------------------------------------------ output

def write_tables(out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "coverage_evidence.tsv", "w") as fh:
        fh.write("domain\tconversion\tgfc_subcommand\ttool\trating\tevidence\n")
        for dom, label, cmd in ROWS:
            fh.write(f"{dom}\t{label}\t{cmd}\tgfc\t1\tgfc {cmd}\n")
            for tool in TOOLS[1:]:
                if (cmd, tool) in EVIDENCE:
                    r, ev = EVIDENCE[(cmd, tool)]
                    fh.write(f"{dom}\t{label}\t{cmd}\t{tool}\t{r:g}\t{ev}\n")
    sym = {1.0: "●", 0.5: "◐", 0.0: ""}
    lines = ["| Conversion | " + " | ".join(TOOLS) + " |",
             "|" + "---|" * (len(TOOLS) + 1)]
    for dom, label, cmd in ROWS:
        lines.append(f"| {label} | " + " | ".join(sym[rating(cmd, t)] for t in TOOLS) + " |")
    (out / "feature_matrix.md").write_text("\n".join(lines) + "\n")


def _save(fig, out: Path, name: str) -> None:
    for ext in ("pdf", "png"):
        fig.savefig(out / f"{name}.{ext}", bbox_inches="tight",
                    dpi=600 if ext == "png" else None)
    plt.close(fig)
    print(f"[fig ] {out / name}.pdf", file=sys.stderr)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    a = ap.parse_args(argv)
    assert all(t in TOOLS for _, t in EVIDENCE), "unknown tool in EVIDENCE"
    assert {c for c, _ in EVIDENCE} <= {r[2] for r in ROWS}, "unknown row in EVIDENCE"
    write_tables(a.out)

    fig = plt.figure(figsize=(7.2, 9.6))
    ax_a = fig.add_axes([0.0, 0.70, 1.0, 0.29])
    ax_b = fig.add_axes([0.02, 0.0, 0.96, 0.68])
    draw_graph(ax_a)
    draw_matrix(ax_b)
    fig.text(0.0, 0.995, "a", fontsize=10, fontweight="bold", va="top")
    fig.text(0.0, 0.685, "b", fontsize=10, fontweight="bold", va="top")
    a.out.mkdir(parents=True, exist_ok=True)
    _save(fig, a.out, "fig_coverage")
    return 0


if __name__ == "__main__":
    sys.exit(main())
