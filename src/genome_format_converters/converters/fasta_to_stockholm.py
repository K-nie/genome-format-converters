#!/usr/bin/env python3
"""Convert FASTA alignments to Stockholm.

Each input FASTA file is expected to contain an aligned set of sequences
(same length). The output is a single Stockholm block per input file, with
no per-column annotation (no `#=GC SS_cons`). Callers who need secondary-
structure or consensus lines should edit the output, or use AlignIO from
Biopython directly.
"""

from pathlib import Path
from typing import Optional

from Bio import AlignIO

from ._common import iter_input_files, log_info, prepare_output_dir

_FASTA_EXTS = [".fasta", ".fa", ".fna", ".fas", ".aln", ".mfa"]


def _convert_file(in_file: Path, out_file: Path) -> None:
    with open(in_file) as handle:
        aln = AlignIO.read(handle, "fasta")
    with open(out_file, "w") as fh:
        AlignIO.write([aln], fh, "stockholm")


def batch_convert(input_dir: str, output_dir: str,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)
    for fa in iter_input_files(in_path, _FASTA_EXTS, pattern=pattern):
        out_file = out_path / (fa.stem + ".sto")
        log_info(f"Converting {fa.name} -> {out_file.name}")
        _convert_file(fa, out_file)
