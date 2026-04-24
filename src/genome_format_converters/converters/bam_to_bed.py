#!/usr/bin/env python3
"""Convert BAM / SAM / CRAM to BED6.

Uses pysam. Chooses the open mode from the file extension so SAM input works
without an index, and iterates with `until_eof=True` so unsorted/unindexed
BAMs don't need a `.bai`.
"""

from pathlib import Path
from typing import Optional

try:
    import pysam
except ImportError:
    pysam = None

from ._common import iter_input_files, log_info, prepare_output_dir, require_pysam

_BAM_EXTS = [".bam", ".sam", ".cram"]


def _open_mode(path: Path) -> str:
    ext = path.suffix.lower()
    if ext == ".sam":
        return "r"
    if ext == ".cram":
        return "rc"
    return "rb"


def _convert_file(in_file: Path, out_file: Path) -> None:
    mode = _open_mode(in_file)
    with pysam.AlignmentFile(str(in_file), mode) as bam, open(out_file, "w") as bed:
        # until_eof handles unsorted/unindexed BAMs and SAM streams.
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            chrom = bam.get_reference_name(read.reference_id)
            start = read.reference_start
            end = read.reference_end
            if start is None or end is None:
                continue
            name = read.query_name or "."
            score = read.mapping_quality
            strand = "-" if read.is_reverse else "+"
            bed.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")


def batch_convert(input_dir: str, output_dir: str,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    if pysam is None:
        require_pysam()
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)
    for bam in iter_input_files(in_path, _BAM_EXTS, pattern=pattern):
        out_file = out_path / (bam.stem + ".bed")
        log_info(f"Converting {bam.name} -> {out_file.name}")
        _convert_file(bam, out_file)
