#!/usr/bin/env python3
"""
Convert BAM/SAM to BED6.
"""

from pathlib import Path
import pysam

def _convert_file(in_file: Path, out_file: Path) -> None:
    """Convert a single BAM/SAM file to BED6."""
    # pysam handles both BAM and SAM transparently
    with pysam.AlignmentFile(in_file, "rb") as bam, open(out_file, 'w') as bed:
        for read in bam.fetch():
            if read.is_unmapped:
                continue
            chrom = bam.get_reference_name(read.reference_id)
            start = read.reference_start
            end = read.reference_end
            name = read.query_name
            score = read.mapping_quality
            strand = '+' if not read.is_reverse else '-'
            bed.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")

def batch_convert(input_dir: str, output_dir: str) -> None:
    """Convert all BAM/SAM files in input_dir to BED in output_dir."""
    in_path = Path(input_dir)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    exts = [".bam", ".sam"]
    for ext in exts:
        for bam in in_path.glob(f"*{ext}"):
            out_file = out_path / bam.with_suffix(".bed").name
            print(f"Converting {bam.name} -> {out_file.name}")
            _convert_file(bam, out_file)