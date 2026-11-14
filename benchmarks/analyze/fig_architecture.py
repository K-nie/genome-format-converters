#!/usr/bin/env python3
"""Schematic of how gfc runs a conversion, drawn from the code.

    command line  ->  dispatch (cli.py)  ->  per-file conversion  ->  outputs
                                          (one converter module,
                                           one parsing backend)

What each box says is taken from the source:
  * cli.py: --input/--output is run as a one-file batch (--pattern equal to the
    file name, then renamed to the requested output); converter modules are
    imported on first use (_LazyModule); --threads N runs one worker process per
    input file (_run_parallel), each the same CLI with --threads 1 and stdin
    closed.
  * converters/*.py: every converter except annotate-tree and
    orthogroups-to-fasta exposes batch_convert(input_dir, output_dir, ...,
    force, pattern), which lists files with _common.iter_input_files and calls
    _convert_file(in, out) on each.
  * Backends, by import: Biopython SeqIO/AlignIO/Phylo; BCBio.GFF; streaming
    line parsers (gff3_to_gtf, gff3_to_bed12, gtf_to_gff3); pysam/htslib;
    _vcf_fast (numpy GT decode, cyvcf2 fallback) + _eigenstrat_common;
    standard library only.

Usage: fig_architecture.py [--out DIR]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

HERE = Path(__file__).resolve().parent
DEFAULT_OUT = HERE.parent / "results" / "figures"

INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#d9d8d3"
GFC = "#2a78d6"
GFC_BG = "#e8f0fb"
BOX = "#f6f5f2"
MONO = "DejaVu Sans Mono"

plt.rcParams.update({"font.size": 7, "text.color": INK,
                     "pdf.fonttype": 42, "ps.fonttype": 42})

# (family, backend line, subcommands)
FAMILIES = [
    ("Sequences, alignments, trees",
     "Biopython SeqIO / AlignIO / Phylo",
     "fasta\u2194fastq, fasta+qual\u2194fastq, fasta-to-table, convert-alignment,\n"
     "stockholm\u2194fasta, maf-to-xmfa, tree-convert, annotate-tree"),
    ("Annotation, streaming",
     "line parser, two passes,\nno object model",
     "gff3-to-gtf, gff3-to-bed12, gtf-to-gff3"),
    ("Annotation, feature model",
     "BCBio.GFF, or a GFF3 reader,\n+ Biopython SeqFeature",
     "gff3-to-bed, gff3-to-table, gff3-to-protein,\ngenbank-to-gff3, fasta-gff-to-gbk"),
    ("Reads and variant tables",
     "pysam (htslib)",
     "bam-to-bed, vcf-to-bed, vcf-to-table,\nvcf-to-consensus, eigenstrat-to-vcf"),
    ("Population-genetics VCF",
     "numpy GT decode of VCF text;\ncyvcf2 for BCF and odd records",
     "vcf-to-eigenstrat, vcf-to-pseudohaploid,\nvcf-to-plink, vcf-to-snp-alignment"),
    ("Tool-output tables",
     "Python standard library only",
     "blast-to-links, delta-to-tab,\nhmmer-tblout-to-tsv, orthogroups-to-fasta"),
]
COUNTS = [11, 3, 5, 5, 4, 4]
assert sum(COUNTS) == 32


# Axes units are inches (the figure is drawn at its print size), and y is
# measured down from the top edge, so text sizes and box sizes can be checked
# against each other directly.
W, H = 7.2, 6.72


def box(ax, x, t, w, h, fc=BOX, ec=GRID, lw=0.6, r=0.06, ls="-", z=0):
    ax.add_patch(FancyBboxPatch((x, H - t - h), w, h, fc=fc, ec=ec, lw=lw,
                                ls=ls, boxstyle=f"round,pad=0,rounding_size={r}",
                                zorder=z))


def txt(ax, x, t, s, size=6.0, color=INK, mono=False, bold=False,
        ha="left", va="top", z=3):
    ax.text(x, H - t, s, ha=ha, va=va, fontsize=size, color=color,
            family=MONO if mono else None, fontweight="bold" if bold else None,
            linespacing=1.3, zorder=z)


def arrow(ax, p, q, color=INK, lw=0.9, ls="-"):
    ax.add_patch(FancyArrowPatch((p[0], H - p[1]), (q[0], H - q[1]),
                                 arrowstyle="-|>", mutation_scale=8, lw=lw,
                                 color=color, ls=ls, zorder=4,
                                 shrinkA=0, shrinkB=0))


def title(ax, x, t, step, text):
    txt(ax, x, t, step, size=7.5, color=GFC, bold=True)
    txt(ax, x + 0.15, t, text, size=7.5, bold=True)


def draw(ax) -> None:
    ax.set_axis_off()
    pad = 0.1

    # 1  command line ---------------------------------------------------
    x1, t1, w1, h1 = 0.0, 0.0, 3.6, 1.55
    box(ax, x1, t1, w1, h1)
    title(ax, x1 + pad, t1 + pad, "1", "Command line")
    txt(ax, x1 + 1.55, t1 + pad + 0.01, "gfc <subcommand> ...", size=6.2, mono=True)
    modes = [
        ("Single file", "one-file batch, renamed to O", "--input F --output O"),
        ("Batch", "every matching file in D", "--input-dir D --output-dir O\n[--pattern G]"),
        ("Parallel batch", "one worker per file (see 5)", "... --threads N"),
    ]
    t = t1 + 0.38
    for name, note, flags in modes:
        box(ax, x1 + pad, t, w1 - 2 * pad, 0.36, fc="white", ec=GRID, r=0.04, z=1)
        txt(ax, x1 + 0.18, t + 0.05, name, size=6.2, bold=True)
        txt(ax, x1 + 0.18, t + 0.19, note, size=5.6, color=INK_2)
        txt(ax, x1 + 1.75, t + 0.05, flags, size=5.6, mono=True, color=GFC)
        t += 0.39

    # 2  dispatch -------------------------------------------------------
    x2, t2, w2, h2 = 3.85, 0.0, 3.35, 1.55
    box(ax, x2, t2, w2, h2)
    title(ax, x2 + pad, t2 + pad, "2", "Dispatch (cli.py)")
    steps = [
        ("argparse", "32 subcommands with shared I/O and batch flags"),
        ("lazy import", "loads only the chosen converter and its libraries"),
        ("batch_convert", "(in_dir, out_dir, <options>, force, pattern)"),
    ]
    t = t2 + 0.38
    for i, (head, body) in enumerate(steps):
        box(ax, x2 + pad, t, w2 - 2 * pad, 0.36, fc="white", ec=GRID, r=0.04, z=1)
        txt(ax, x2 + 0.18, t + 0.05, head, size=6.0, mono=True, color=GFC)
        txt(ax, x2 + 0.18, t + 0.2, body, size=5.6, color=INK_2,
            mono=head == "batch_convert")
        t += 0.39
    arrow(ax, (x1 + w1 + 0.02, 0.8), (x2 - 0.02, 0.8))

    # 3  per-file conversion ----------------------------------------------
    x3, t3, w3 = 0.0, 1.8, W
    rh, rgap = 0.33, 0.04
    h3 = 0.72 + len(FAMILIES) * (rh + rgap) + 0.42
    box(ax, x3, t3, w3, h3, fc=GFC_BG, ec=GFC, lw=0.7)
    title(ax, x3 + pad, t3 + pad, "3", "Per-file conversion: one converter module, one parsing backend")
    txt(ax, x3 + pad, t3 + 0.33,
        "for f in iter_input_files(in_dir, exts, pattern):  # sorted\n"
        "    _convert_file(f, out_dir / (stem(f) + ext))", size=5.6, mono=True)
    arrow(ax, (x2 + w2 / 2, t2 + h2 + 0.02), (x2 + w2 / 2, t3 - 0.02))
    t = t3 + 0.72
    cx = [x3 + 0.18, x3 + 1.75, x3 + 3.75, W - 0.18]
    for (fam, backend, cmds), n in zip(FAMILIES, COUNTS):
        box(ax, x3 + pad, t, w3 - 2 * pad, rh, fc="white", ec=GRID, r=0.04, z=1)
        mid = t + rh / 2
        txt(ax, cx[0], mid, fam, size=5.9, bold=True, va="center")
        txt(ax, cx[1], mid, backend, size=5.5, color=GFC, va="center")
        txt(ax, cx[2], mid, cmds, size=5.2, mono=True, color=INK_2, va="center")
        txt(ax, cx[3], mid, str(n), size=6.0, color=INK_2, ha="right", va="center")
        t += rh + rgap
    txt(ax, x3 + pad, t + 0.04,
        "Converters share helpers, not a data model: _common (file listing, output "
        "directories, logging) and, for the\npopulation-genetics VCF family, "
        "_eigenstrat_common (sample and chromosome maps, polarisation, filters).\n"
        "annotate-tree and orthogroups-to-fasta take named input files rather than "
        "a directory; fasta-gff-to-gbk pairs files by stem.",
        size=5.5, color=INK_2)

    # 4  outputs ---------------------------------------------------------
    t4 = t3 + h3 + 0.25
    x4, w4, h4 = 0.0, 2.45, H - t4
    box(ax, x4, t4, w4, h4)
    title(ax, x4 + pad, t4 + pad, "4", "Outputs")
    txt(ax, x4 + pad, t4 + 0.34,
        "One output per input in --output-dir,\n"
        "named from the input stem\n"
        "(sample.vcf.gz \u2192 sample.geno ...),\n"
        "or the exact --output path.\n"
        "File sets: EIGENSTRAT .geno/.snp/.ind,\n"
        "PLINK .bed/.bim/.fam, FASTA + QUAL.\n"
        "Logs on stderr; errors exit non-zero.",
        size=5.7, color=INK_2)
    arrow(ax, (x4 + w4 / 2, t3 + h3 + 0.02), (x4 + w4 / 2, t4 - 0.02))

    # 5  --threads -----------------------------------------------------------
    x5 = 2.65
    w5, h5 = W - x5, h4
    box(ax, x5, t4, w5, h5)
    title(ax, x5 + pad, t4 + pad, "5", "Parallel batch (--threads N)")
    txt(ax, x5 + pad, t4 + 0.34,
        "The parent lists the input files and runs one\n"
        "worker per file, at most N at a time. Each worker\n"
        "is the same command on one file with stdin closed,\n"
        "so outputs are identical to a sequential run.\n"
        "Pseudohaploid draws are seeded per file:\n"
        "seed XOR crc32(file name).",
        size=5.7, color=INK_2)
    wx, ww = x5 + 2.45, W - x5 - 2.45 - pad
    for k, name in enumerate(["a.vcf.gz", "b.vcf.gz", "c.vcf.gz"]):
        tt = t4 + 0.12 + k * 0.33
        box(ax, wx, tt, ww, 0.29, fc="white", ec=GRID, r=0.04, z=1)
        txt(ax, wx + 0.07, tt + 0.145, f"gfc ... --pattern {name}\n  --threads 1 --force",
            size=5.0, mono=True, va="center")

    ax.set_xlim(-0.02, W + 0.02)
    ax.set_ylim(-0.02, H + 0.02)
    ax.set_aspect("equal")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    a = ap.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(W, H))
    fig.subplots_adjust(0, 0, 1, 1)
    draw(ax)
    for ext in ("pdf", "png"):
        fig.savefig(a.out / f"fig_architecture.{ext}", bbox_inches="tight",
                    dpi=600 if ext == "png" else None)
    plt.close(fig)
    print(f"[fig ] {a.out / 'fig_architecture'}.pdf", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
