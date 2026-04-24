#!/usr/bin/env python3
"""GenBank to GFF3 conversion."""

from pathlib import Path
from typing import Optional

from BCBio import GFF
from Bio import SeqIO

from ._common import iter_input_files, log_info, prepare_output_dir

_GBK_EXTS = [".gbk", ".gb", ".gbff"]


def _convert_file(in_file: Path, out_file: Path) -> None:
    with open(in_file) as in_handle, open(out_file, "w") as out_handle:
        records = SeqIO.parse(in_handle, "genbank")
        GFF.write(records, out_handle)


def batch_convert(input_dir: str, output_dir: str,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)
    for gbk in iter_input_files(in_path, _GBK_EXTS, pattern=pattern):
        out_file = out_path / (gbk.stem + ".gff3")
        log_info(f"Converting {gbk.name} -> {out_file.name}")
        _convert_file(gbk, out_file)
