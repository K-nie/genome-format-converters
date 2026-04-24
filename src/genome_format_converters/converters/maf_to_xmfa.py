#!/usr/bin/env python3
"""Convert MAF to XMFA (progressiveMauve format).

XMFA preserves alignment columns (gaps are retained as `-`); only the MAF
wrapper lines are dropped. Each alignment block is terminated by `=`.
"""

from pathlib import Path
from typing import Optional

from Bio import AlignIO

from ._common import iter_input_files, log_info, prepare_output_dir


def _convert_file(in_file: Path, out_file: Path) -> None:
    with open(in_file) as in_handle, open(out_file, "w") as out_handle:
        for aln in AlignIO.parse(in_handle, "maf"):
            for seq in aln:
                # XMFA header: `>seqid:start-end strand`. When annotations are
                # missing (synthetic test MAFs) fall back to bare `>seqid`.
                start = getattr(seq, "annotations", {}).get("start")
                size = getattr(seq, "annotations", {}).get("size")
                strand = getattr(seq, "annotations", {}).get("strand")
                if start is not None and size is not None:
                    end = start + size
                    strand_ch = "+" if strand in (1, "+", None) else "-"
                    out_handle.write(f"> {seq.id}:{start + 1}-{end} {strand_ch}\n")
                else:
                    out_handle.write(f"> {seq.id}\n")
                out_handle.write(str(seq.seq) + "\n")
            out_handle.write("=\n")


def batch_convert(input_dir: str, output_dir: str,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)
    for maf in iter_input_files(in_path, [".maf"], pattern=pattern):
        out_file = out_path / (maf.stem + ".xmfa")
        log_info(f"Converting {maf.name} -> {out_file.name}")
        _convert_file(maf, out_file)
