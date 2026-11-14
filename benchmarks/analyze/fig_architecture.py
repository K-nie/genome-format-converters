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
# (family, parsing backend, example subcommands, count)
FAMILIES = [
    ("Sequences,\nalignments, trees", "Biopython\nSeqIO / AlignIO / Phylo",
     "fasta-to-fastq\nconvert-alignment\ntree-convert", 11),
    ("Annotation,\nstreaming", "two-pass line\nparser",
     "gff3-to-gtf\ngff3-to-bed12\ngtf-to-gff3", 3),
    ("Annotation,\nfeature model", "BCBio.GFF +\nBiopython SeqFeature",
     "gff3-to-protein\ngenbank-to-gff3\nfasta-gff-to-gbk", 5),
    ("Reads and\nvariant tables", "pysam\n(htslib)",
     "bam-to-bed\nvcf-to-table\nvcf-to-consensus", 5),
    ("Population-\ngenetics VCF", "numpy GT decode;\ncyvcf2 fallback",
     "vcf-to-eigenstrat\nvcf-to-pseudohaploid\nvcf-to-plink", 4),
    ("Tool-output\ntables", "Python standard\nlibrary",
     "blast-to-links\nhmmer-tblout-to-tsv\northogroups-to-fasta", 4),
]
assert sum(f[3] for f in FAMILIES) == 32

# Axes units are inches at print size; y is measured down from the top.
W, H = 7.2, 5.47


def Y(t):
    return H - t


def box(ax, x, t, w, h, fc="white", ec=GRID, lw=0.7, r=0.06, ls="-", z=1):
    ax.add_patch(FancyBboxPatch((x, Y(t) - h), w, h, fc=fc, ec=ec, lw=lw, ls=ls,
                                boxstyle=f"round,pad=0,rounding_size={r}",
                                zorder=z))


def txt(ax, x, t, s, size=6.0, color=INK, mono=False, bold=False,
        ha="center", va="center", z=3, style=None):
    ax.text(x, Y(t), s, ha=ha, va=va, fontsize=size, color=color,
            family=MONO if mono else None, fontweight="bold" if bold else None,
            linespacing=1.25, zorder=z, style=style)


def arrow(ax, p, q, color=INK_2, lw=1.0, ls="-", conn="arc3", label=None,
          lpos=None, lha="left"):
    ax.add_patch(FancyArrowPatch((p[0], Y(p[1])), (q[0], Y(q[1])),
                                 arrowstyle="-|>", mutation_scale=9, lw=lw,
                                 color=color, ls=ls, zorder=2, shrinkA=0,
                                 shrinkB=0, connectionstyle=conn))
    if label:
        lx, lt = lpos
        txt(ax, lx, lt, label, size=5.6, color=INK_2, ha=lha, style="italic")


def stage(ax, x, t, w, h, head, body, accent=False):
    box(ax, x, t, w, h, fc=GFC_BG if accent else "white",
        ec=GFC if accent else GRID, lw=0.8 if accent else 0.7)
    txt(ax, x + w / 2, t + 0.17, head, size=6.8, bold=True)
    txt(ax, x + w / 2, t + 0.19 + (h - 0.2) / 2, body, size=5.7, color=INK_2)


def draw(ax) -> None:
    ax.set_axis_off()
    cx = 3.0                      # centre line of the main flow
    mw = 3.0                      # main-flow box width
    mx = cx - mw / 2

    # input files
    t_in = 0.0
    for k in range(3):            # a small stack of file cards
        box(ax, cx - 0.62 + 0.05 * k, t_in + 0.05 * (2 - k), 1.24, 0.42,
            fc="white", ec=INK_2, lw=0.6, r=0.04, z=1 + k)
    txt(ax, cx + 0.05, t_in + 0.2, "Input file or directory", size=6.6, bold=True, z=6)
    txt(ax, cx + 0.05, t_in + 0.36, "FASTA, GFF3, VCF, BAM, ...", size=5.6,
        color=INK_2, z=6)

    # command line
    t_cli, h_cli = 0.85, 0.62
    arrow(ax, (cx, t_in + 0.53), (cx, t_cli))
    stage(ax, mx, t_cli, mw, h_cli, "gfc <subcommand>",
          "argparse: 32 subcommands\n"
          "--input F --output O  or  --input-dir D --output-dir O")

    # dispatch + file loop
    t_dis, h_dis = 1.82, 0.62
    arrow(ax, (cx, t_cli + h_cli), (cx, t_dis),
          label="--threads 1", lpos=(cx + 0.08, (t_cli + h_cli + t_dis) / 2))
    stage(ax, mx, t_dis, mw, h_dis, "Load one converter, list its files",
          "lazy import: only this converter's libraries load\n"
          "for each input file (sorted): _convert_file(in, out)", accent=True)

    # --threads branch
    tx, tw = 5.05, 2.1
    t_thr, h_thr = t_cli, 1.72
    box(ax, tx, t_thr, tw, h_thr, fc="white", ec=GRID)
    txt(ax, tx + tw / 2, t_thr + 0.17, "--threads N", size=6.8, bold=True)
    txt(ax, tx + tw / 2, t_thr + 0.42,
        "parent lists the files and\nstarts one worker per file,\nN at a time",
        size=5.7, color=INK_2)
    for k in range(3):
        wt = t_thr + 0.78 + k * 0.25
        box(ax, tx + 0.12, wt, tw - 0.24, 0.2, fc=BOX, ec=GRID, r=0.03, z=2)
        txt(ax, tx + tw / 2, wt + 0.1, f"worker {k + 1}: --pattern file{k + 1}",
            size=5.2, mono=True, z=3)
    arrow(ax, (mx + mw, t_cli + h_cli / 2), (tx, t_cli + h_cli / 2),
          label="N > 1", lpos=(mx + mw + 0.12, t_cli + h_cli / 2 - 0.1))
    arrow(ax, (tx, t_thr + 1.3), (mx + mw, t_dis + h_dis / 2),
          conn="arc3,rad=-0.15")
    txt(ax, tx + tw / 2, t_thr + h_thr - 0.1,
        "\u2190 each worker re-enters with --threads 1", size=5.4,
        color=INK_2, style="italic")

    # backends (fan out)
    t_be, h_be = 3.05, 1.25
    n = len(FAMILIES)
    gap = 0.07
    bw = (W - gap * (n - 1)) / n
    bottom_dis = t_dis + h_dis
    for i, (fam, backend, examples, cnt) in enumerate(FAMILIES):
        bx = i * (bw + gap)
        box(ax, bx, t_be, bw, h_be, fc="white", ec=GFC, lw=0.7)
        txt(ax, bx + bw / 2, t_be + 0.2, fam, size=6.0, bold=True)
        txt(ax, bx + bw / 2, t_be + 0.53, backend, size=5.3, color=GFC)
        txt(ax, bx + bw / 2, t_be + 0.92, examples, size=4.7, mono=True, color=INK_2)
        txt(ax, bx + bw - 0.06, t_be + 0.06, str(cnt), size=5.6, color=INK_2,
            ha="right", va="top")
        sx = cx - mw / 2 + mw * (i + 0.5) / n
        arrow(ax, (sx, bottom_dis), (bx + bw / 2, t_be), lw=0.8)

    # outputs (converge)
    t_out, h_out = 4.85, 0.6
    for i in range(n):
        bx = i * (bw + gap)
        ex = cx - mw / 2 + mw * (i + 0.5) / n
        arrow(ax, (bx + bw / 2, t_be + h_be), (ex, t_out), lw=0.8)
    stage(ax, mx, t_out, mw, h_out, "Output files in --output-dir",
          "one per input, named from its stem (or the exact --output);\n"
          "file sets for EIGENSTRAT, PLINK and FASTA + QUAL")
    txt(ax, mx + mw + 0.15, t_out + 0.3,
        "logs on stderr;\nerrors exit non-zero", size=5.6, color=INK_2, ha="left")

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
