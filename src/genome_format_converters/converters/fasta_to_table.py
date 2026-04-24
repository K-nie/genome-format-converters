#!/usr/bin/env python3
"""Convert FASTA to two-column TSV (id, sequence)."""

import csv
from pathlib import Path
from typing import Optional

from Bio import SeqIO

from ._common import iter_input_files, log_info, prepare_output_dir

_FASTA_EXTS = [".fasta", ".fa", ".fna", ".fas"]


def _convert_file(in_file: Path, out_file: Path) -> None:
    with open(in_file) as in_handle, open(out_file, "w", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t")
        writer.writerow(["id", "sequence"])
        for rec in SeqIO.parse(in_handle, "fasta"):
            writer.writerow([rec.id, str(rec.seq)])


def convert_one(in_path: Path, out_path: Path) -> None:
    log_info(f"Converting {in_path.name} -> {out_path.name}")
    _convert_file(in_path, out_path)


def batch_convert(input_dir: str, output_dir: str,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)
    for fa in iter_input_files(in_path, _FASTA_EXTS, pattern=pattern):
        out_file = out_path / (fa.stem + ".tsv")
        log_info(f"Converting {fa.name} -> {out_file.name}")
        _convert_file(fa, out_file)
