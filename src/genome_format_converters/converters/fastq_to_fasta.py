#!/usr/bin/env python3
"""Convert FASTQ to FASTA (drop qualities)."""

from pathlib import Path
from typing import Optional

from Bio import SeqIO

from ._common import iter_input_files, log_info, prepare_output_dir

_FASTQ_EXTS = [".fastq", ".fq"]


def _convert_file(in_file: Path, out_file: Path) -> None:
    SeqIO.convert(str(in_file), "fastq", str(out_file), "fasta")


def convert_one(in_path: Path, out_path: Path) -> None:
    log_info(f"Converting {in_path.name} -> {out_path.name}")
    _convert_file(in_path, out_path)


def batch_convert(input_dir: str, output_dir: str,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)
    for fq in iter_input_files(in_path, _FASTQ_EXTS, pattern=pattern):
        out_file = out_path / (fq.stem + ".fasta")
        log_info(f"Converting {fq.name} -> {out_file.name}")
        _convert_file(fq, out_file)
