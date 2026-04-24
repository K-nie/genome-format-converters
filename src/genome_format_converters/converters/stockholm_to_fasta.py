#!/usr/bin/env python3
"""Convert Stockholm-format alignments (Rfam / Pfam style) to FASTA.

Stockholm can hold multiple alignments per file; this converter preserves
that by writing one FASTA per input alignment when more than one is
present (``<stem>.fasta`` for a single alignment, ``<stem>.<n>.fasta`` for
batched alignments).
"""

from pathlib import Path
from typing import Optional

from Bio import AlignIO

from ._common import iter_input_files, log_info, log_warn, prepare_output_dir

_STOCKHOLM_EXTS = [".sto", ".stk", ".stockholm"]


def _convert_file(in_file: Path, out_dir: Path, stem: str) -> int:
    count = 0
    with open(in_file) as handle:
        alignments = list(AlignIO.parse(handle, "stockholm"))
    if not alignments:
        log_warn(f"{in_file.name}: no Stockholm alignments found.")
        return 0
    for idx, aln in enumerate(alignments, 1):
        if len(alignments) == 1:
            out = out_dir / f"{stem}.fasta"
        else:
            out = out_dir / f"{stem}.{idx}.fasta"
        with open(out, "w") as fh:
            AlignIO.write([aln], fh, "fasta")
        count += 1
    return count


def batch_convert(input_dir: str, output_dir: str,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)
    for sto in iter_input_files(in_path, _STOCKHOLM_EXTS, pattern=pattern):
        log_info(f"Converting {sto.name} -> FASTA")
        _convert_file(sto, out_path, sto.stem)
