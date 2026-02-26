#!/usr/bin/env python3
"""
Convert phylogenetic tree formats (newick, nexus, phyloxml).
"""
from pathlib import Path
from Bio import Phylo

def _convert_file(in_file: Path, out_file: Path, in_fmt: str, out_fmt: str) -> None:
    """Convert a single tree file."""
    tree = Phylo.read(in_file, in_fmt)
    Phylo.write(tree, out_file, out_fmt)

def batch_convert(input_dir: str, output_dir: str, in_format: str, out_format: str) -> None:
    """Convert all tree files in input_dir to the target format."""
    in_path = Path(input_dir)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # Define extensions for each input format
    ext_map = {
        'newick': ['.nwk', '.tree', '.tre', '.newick'],
        'nexus': ['.nex', '.nexus'],
        'phyloxml': ['.xml', '.phyloxml']
    }
    exts = ext_map.get(in_format, [])
    files = [f for f in in_path.iterdir() if f.is_file() and f.suffix.lower() in exts]
    for in_file in files:
        out_file = out_path / f"{in_file.stem}.{out_format}"
        print(f"Converting {in_file.name} -> {out_file.name}")
        _convert_file(in_file, out_file, in_format, out_format)