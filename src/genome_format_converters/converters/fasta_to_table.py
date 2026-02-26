#!/usr/bin/env python3
"""
Convert FASTA to two‑column TSV (id, sequence).
"""

import csv
from pathlib import Path
from Bio import SeqIO

def _convert_file(in_file: Path, out_file: Path) -> None:
    """Convert a single FASTA file to TSV."""
    with open(in_file) as in_handle, open(out_file, 'w', newline='') as out_handle:
        writer = csv.writer(out_handle, delimiter='\t')
        writer.writerow(["id", "sequence"])
        for rec in SeqIO.parse(in_handle, "fasta"):
            writer.writerow([rec.id, str(rec.seq)])

def batch_convert(input_dir: str, output_dir: str) -> None:
    """Convert all FASTA files in input_dir to TSV tables in output_dir."""
    in_path = Path(input_dir)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    exts = [".fasta", ".fa", ".fna", ".fas"]
    for ext in exts:
        for fa in in_path.glob(f"*{ext}"):
            out_file = out_path / fa.with_suffix(".tsv").name
            print(f"Converting {fa.name} -> {out_file.name}")
            _convert_file(fa, out_file)