#!/usr/bin/env python3
"""
Convert FASTA to FASTQ with a default quality value.
"""

from pathlib import Path
from Bio import SeqIO

def _convert_file(in_file: Path, out_file: Path, qual_char: str = 'I') -> None:
    """Convert a single FASTA file to FASTQ with constant quality."""
    with open(in_file) as in_handle, open(out_file, 'w') as out_handle:
        for rec in SeqIO.parse(in_handle, "fasta"):
            rec.letter_annotations["phred_quality"] = [ord(qual_char)-33] * len(rec)
            SeqIO.write(rec, out_handle, "fastq")

def batch_convert(input_dir: str, output_dir: str, qual_char: str = 'I') -> None:
    """Convert all FASTA files in input_dir to FASTQ in output_dir."""
    in_path = Path(input_dir)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    exts = [".fasta", ".fa", ".fna", ".fas"]
    for ext in exts:
        for fa in in_path.glob(f"*{ext}"):
            out_file = out_path / fa.with_suffix(".fastq").name
            print(f"Converting {fa.name} -> {out_file.name}")
            _convert_file(fa, out_file, qual_char)