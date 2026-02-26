#!/usr/bin/env python3
"""
Convert alignment files between formats (fasta, phylip, nexus, clustal).
Only processes files containing 'aln' in the filename.
"""

from pathlib import Path
from Bio import AlignIO

def _convert_file(in_file: Path, out_file: Path, in_fmt: str, out_fmt: str) -> None:
    """Convert a single alignment file."""
    with open(in_file) as in_handle, open(out_file, 'w') as out_handle:
        aln = AlignIO.read(in_handle, in_fmt)
        AlignIO.write(aln, out_handle, out_fmt)

def batch_convert(input_dir: str, output_dir: str, in_format: str, out_format: str) -> None:
    """Convert alignment files in input_dir (only those with 'aln' in name)."""
    in_path = Path(input_dir)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # Only process files with 'aln' in the name
    files = [f for f in in_path.iterdir() if f.is_file() and 'aln' in f.name.lower()]

    for in_file in files:
        out_file = out_path / f"{in_file.stem}.{out_format}"
        print(f"Converting {in_file.name} -> {out_file.name}")
        _convert_file(in_file, out_file, in_format, out_format)