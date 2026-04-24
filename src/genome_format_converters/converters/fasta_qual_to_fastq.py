#!/usr/bin/env python3
"""Combine FASTA + QUAL files into FASTQ. Pairs files by stem."""

from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.SeqIO.QualityIO import PairedFastaQualIterator

from ._common import log_info, log_warn, prepare_output_dir

_FASTA_EXTS = [".fasta", ".fa", ".fna", ".fas"]


def _combine_pair(fasta_file: Path, qual_file: Path, out_file: Path) -> None:
    with open(fasta_file) as f_handle, open(qual_file) as q_handle, \
         open(out_file, "w") as out_handle:
        records = PairedFastaQualIterator(f_handle, q_handle)
        SeqIO.write(records, out_handle, "fastq")


def batch_convert(input_dir: str, output_dir: str,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)

    fastas = {f.stem: f for ext in _FASTA_EXTS for f in in_path.glob(f"*{ext}")}
    quals = {f.stem: f for f in in_path.glob("*.qual")}

    common = sorted(set(fastas) & set(quals))
    if not common:
        log_warn("No matching FASTA / QUAL pairs found (matched by stem).")
        return

    for stem in common:
        out_file = out_path / f"{stem}.fastq"
        log_info(f"Combining {fastas[stem].name} + {quals[stem].name} -> {out_file.name}")
        _combine_pair(fastas[stem], quals[stem], out_file)
