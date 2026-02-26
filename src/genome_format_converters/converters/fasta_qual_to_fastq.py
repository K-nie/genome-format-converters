#!/usr/bin/env python3
"""
Combine FASTA and QUAL files into FASTQ.
"""

from pathlib import Path
from Bio import SeqIO
from Bio.SeqIO.QualityIO import PairedFastaQualIterator

def _combine_pair(fasta_file: Path, qual_file: Path, out_file: Path) -> None:
    """Combine one FASTA and one QUAL file into FASTQ."""
    with open(fasta_file) as f_handle, open(qual_file) as q_handle, open(out_file, 'w') as out_handle:
        records = PairedFastaQualIterator(f_handle, q_handle)
        SeqIO.write(records, out_handle, "fastq")

def batch_convert(input_dir: str, output_dir: str) -> None:
    """For each FASTA/QUAL pair in input_dir, produce a FASTQ file."""
    in_path = Path(input_dir)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    fasta_exts = [".fasta", ".fa", ".fna", ".fas"]
    fastas = {f.stem: f for ext in fasta_exts for f in in_path.glob(f"*{ext}")}
    quals   = {f.stem: f for f in in_path.glob("*.qual")}

    common = set(fastas) & set(quals)
    for stem in common:
        out_file = out_path / f"{stem}.fastq"
        print(f"Combining {fastas[stem].name} + {quals[stem].name} -> {out_file.name}")
        _combine_pair(fastas[stem], quals[stem], out_file)