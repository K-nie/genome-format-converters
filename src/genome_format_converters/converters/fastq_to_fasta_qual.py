#!/usr/bin/env python3
"""Split a FASTQ file into parallel FASTA and QUAL files."""

from pathlib import Path
from typing import Optional

from Bio import SeqIO

from ._common import iter_input_files, log_info, prepare_output_dir

_FASTQ_EXTS = [".fastq", ".fq"]


def _split_file(in_file: Path, out_fasta: Path, out_qual: Path) -> None:
    with open(in_file) as in_handle, open(out_fasta, "w") as fa_handle, \
         open(out_qual, "w") as qu_handle:
        for rec in SeqIO.parse(in_handle, "fastq"):
            SeqIO.write(rec, fa_handle, "fasta")
            qual_str = " ".join(str(q) for q in rec.letter_annotations["phred_quality"])
            qu_handle.write(f">{rec.id}\n{qual_str}\n")


def batch_convert(input_dir: str, output_dir: str,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)
    for fq in iter_input_files(in_path, _FASTQ_EXTS, pattern=pattern):
        stem = fq.stem
        out_fasta = out_path / f"{stem}.fasta"
        out_qual = out_path / f"{stem}.qual"
        log_info(f"Splitting {fq.name} -> {out_fasta.name}, {out_qual.name}")
        _split_file(fq, out_fasta, out_qual)
