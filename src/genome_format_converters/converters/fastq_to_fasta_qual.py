#!/usr/bin/env python3
"""
Split a FASTQ file into FASTA and QUAL files.
"""

from pathlib import Path
from Bio import SeqIO

def _split_file(in_file: Path, out_fasta: Path, out_qual: Path) -> None:
    """Split one FASTQ file into FASTA and QUAL."""
    with open(in_file) as in_handle, open(out_fasta, 'w') as fa_handle, open(out_qual, 'w') as qu_handle:
        for rec in SeqIO.parse(in_handle, "fastq"):
            SeqIO.write(rec, fa_handle, "fasta")
            qual_str = " ".join(str(q) for q in rec.letter_annotations["phred_quality"])
            qu_handle.write(f">{rec.id}\n{qual_str}\n")

def batch_convert(input_dir: str, output_dir: str) -> None:
    """Split all FASTQ files in input_dir into FASTA+QUAL pairs in output_dir."""
    in_path = Path(input_dir)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    exts = [".fastq", ".fq"]
    for ext in exts:
        for fq in in_path.glob(f"*{ext}"):
            stem = fq.stem
            out_fasta = out_path / f"{stem}.fasta"
            out_qual  = out_path / f"{stem}.qual"
            print(f"Splitting {fq.name} -> {out_fasta.name}, {out_qual.name}")
            _split_file(fq, out_fasta, out_qual)