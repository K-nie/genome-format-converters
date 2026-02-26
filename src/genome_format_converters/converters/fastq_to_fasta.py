#!/usr/bin/env python3
"""
Convert FASTQ to FASTA (drop qualities).
"""

from pathlib import Path
from Bio import SeqIO

def batch_convert(input_dir: str, output_dir: str) -> None:
    """Convert all FASTQ files in input_dir to FASTA in output_dir."""
    in_path = Path(input_dir)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    exts = [".fastq", ".fq"]
    for ext in exts:
        for fq in in_path.glob(f"*{ext}"):
            out_file = out_path / fq.with_suffix(".fasta").name
            print(f"Converting {fq.name} -> {out_file.name}")
            SeqIO.convert(fq, "fastq", out_file, "fasta")