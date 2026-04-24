#!/usr/bin/env python3
"""Convert alignment files between formats (FASTA / PHYLIP / NEXUS / CLUSTAL)."""

from pathlib import Path
from typing import Optional

from Bio import AlignIO

from ._common import log_info, prepare_output_dir

# Biopython format name -> common file extensions for input discovery.
_FORMAT_EXTS = {
    "fasta": [".fasta", ".fa", ".fna", ".fas", ".mfa", ".aln"],
    "phylip": [".phy", ".phylip"],
    "nexus": [".nex", ".nexus"],
    "clustal": [".aln", ".clustal", ".clw"],
}


def _convert_file(in_file: Path, out_file: Path, in_fmt: str, out_fmt: str) -> None:
    with open(in_file) as in_handle, open(out_file, "w") as out_handle:
        aln = AlignIO.read(in_handle, in_fmt)
        AlignIO.write(aln, out_handle, out_fmt)


def batch_convert(input_dir: str, output_dir: str,
                  in_format: str, out_format: str,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)

    exts = {e.lower() for e in _FORMAT_EXTS.get(in_format, [])}
    if pattern:
        candidates = sorted(in_path.glob(pattern))
    else:
        candidates = [f for f in sorted(in_path.iterdir()) if f.is_file()]

    for in_file in candidates:
        # When --pattern is supplied, trust the user's glob; otherwise require
        # a recognised extension for the declared input format so we don't
        # silently try to parse unrelated files.
        if not pattern and exts and in_file.suffix.lower() not in exts:
            continue
        out_file = out_path / f"{in_file.stem}.{out_format}"
        log_info(f"Converting {in_file.name} -> {out_file.name}")
        _convert_file(in_file, out_file, in_format, out_format)
