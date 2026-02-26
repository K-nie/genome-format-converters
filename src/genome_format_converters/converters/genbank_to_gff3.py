#!/usr/bin/env python3
"""
GenBank to GFF3 conversion.
"""

from pathlib import Path
from BCBio import GFF
from Bio import SeqIO

def _convert_file(in_file: Path, out_file: Path) -> None:
    """Convert a single GenBank file to GFF3."""
    with open(in_file) as in_handle, open(out_file, 'w') as out_handle:
        records = SeqIO.parse(in_handle, "genbank")
        GFF.write(records, out_handle)

def batch_convert(input_dir: str, output_dir: str) -> None:
    """Convert all GenBank files in input_dir to GFF3 in output_dir."""
    in_path = Path(input_dir)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    exts = [".gbk", ".gb", ".gbff"]
    for ext in exts:
        for gbk in in_path.glob(f"*{ext}"):
            out_file = out_path / gbk.with_suffix(".gff3").name
            print(f"Converting {gbk.name} -> {out_file.name}")
            _convert_file(gbk, out_file)