#!/usr/bin/env python3
"""Convert FASTA to FASTQ with a constant per-base quality value."""

from pathlib import Path
from typing import Optional

from Bio import SeqIO

from ._common import iter_input_files, log_info, prepare_output_dir

_FASTA_EXTS = [".fasta", ".fa", ".fna", ".fas"]


def _convert_file(in_file: Path, out_file: Path, qual_char: str = "I") -> None:
    with open(in_file) as in_handle, open(out_file, "w") as out_handle:
        q = ord(qual_char) - 33
        for rec in SeqIO.parse(in_handle, "fasta"):
            rec.letter_annotations["phred_quality"] = [q] * len(rec)
            SeqIO.write(rec, out_handle, "fastq")


def convert_one(in_path: Path, out_path: Path, qual_char: str = "I") -> None:
    """Single-file entry point used by the `--input FILE --output FILE` mode."""
    log_info(f"Converting {in_path.name} -> {out_path.name}")
    _convert_file(in_path, out_path, qual_char)


def batch_convert(input_dir: str, output_dir: str,
                  qual_char: str = "I",
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)
    for fa in iter_input_files(in_path, _FASTA_EXTS, pattern=pattern):
        out_file = out_path / (fa.stem + ".fastq")
        log_info(f"Converting {fa.name} -> {out_file.name}")
        _convert_file(fa, out_file, qual_char)
