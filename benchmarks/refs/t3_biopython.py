#!/usr/bin/env python3
"""Handwritten Biopython reference for T3 (paired FASTA + GFF3 -> GenBank).

The "what a careful bioinformatician would write in an hour" baseline.
The gfc -> ref delta measures the CLI wrapper / argument-validation
overhead, not algorithmic differences.

Pairs each ``foo.fasta`` with ``foo.gff3`` from the input dir; emits one
``foo.gb`` per pair. Files without a paired GFF3 are skipped.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from Bio import SeqIO
from BCBio import GFF


def convert_one(fasta_path: Path, gff_path: Path, out_path: Path) -> None:
    fa_records = SeqIO.to_dict(SeqIO.parse(str(fasta_path), "fasta"))
    with open(gff_path) as gff_fh:
        # base_dict attaches each GFF entry's features to the corresponding
        # FASTA SeqRecord so the output GenBank carries both sequence and
        # annotation in one file.
        gff_records = list(GFF.parse(gff_fh, base_dict=fa_records))
    # Force an alphabet-free molecule_type since Biopython 1.78+ requires
    # one for GenBank output and bcbio-gff doesn't always set it.
    for rec in gff_records:
        rec.annotations.setdefault("molecule_type", "DNA")
    SeqIO.write(gff_records, str(out_path), "genbank")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-dir", required=True)
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()

    in_dir = Path(args.input_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pairs = 0
    for fasta in sorted(in_dir.glob("*.fasta")):
        gff = fasta.with_suffix(".gff3")
        if not gff.exists():
            continue
        out_path = out_dir / f"{fasta.stem}.gb"
        try:
            convert_one(fasta, gff, out_path)
            pairs += 1
        except Exception as exc:
            print(f"[skip-T3-ref] {fasta.name}: {exc}", file=sys.stderr)
    if pairs == 0:
        print("[warn-T3-ref] no fasta+gff3 pairs found", file=sys.stderr)


if __name__ == "__main__":
    main()
