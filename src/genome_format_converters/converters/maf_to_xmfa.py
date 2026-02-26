#!/usr/bin/env python3
"""
Convert MAF to XMFA (progressiveMauve format).
"""

from pathlib import Path
from Bio import AlignIO

def _convert_file(in_file: Path, out_file: Path) -> None:
    """Convert a single MAF file to XMFA."""
    with open(in_file) as in_handle, open(out_file, 'w') as out_handle:
        for aln in AlignIO.parse(in_handle, "maf"):
            out_handle.write("#\n")
            for seq in aln:
                out_handle.write(f">{seq.id}\n{str(seq.seq).replace('-', '')}\n")

def batch_convert(input_dir: str, output_dir: str) -> None:
    """Convert all MAF files in input_dir to XMFA in output_dir."""
    in_path = Path(input_dir)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    for maf in in_path.glob("*.maf"):
        out_file = out_path / maf.with_suffix(".xmfa").name
        print(f"Converting {maf.name} -> {out_file.name}")
        _convert_file(maf, out_file)