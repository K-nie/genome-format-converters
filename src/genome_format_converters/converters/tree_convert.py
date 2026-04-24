#!/usr/bin/env python3
"""Convert phylogenetic tree formats (Newick, NEXUS, PhyloXML)."""

from pathlib import Path
from typing import Optional

from Bio import Phylo

from ._common import log_info, log_warn, prepare_output_dir

_EXT_MAP = {
    "newick": [".nwk", ".tree", ".tre", ".newick"],
    "nexus": [".nex", ".nexus"],
    "phyloxml": [".xml", ".phyloxml"],
}


def _convert_file(in_file: Path, out_file: Path, in_fmt: str, out_fmt: str) -> None:
    tree = Phylo.read(str(in_file), in_fmt)
    Phylo.write(tree, str(out_file), out_fmt)


def batch_convert(input_dir: str, output_dir: str,
                  in_format: str, out_format: str,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)
    exts = _EXT_MAP.get(in_format, [])
    if not exts:
        log_warn(f"No file-extension list registered for input format {in_format!r}.")
    if pattern:
        candidates = sorted(in_path.glob(pattern))
    else:
        candidates = [f for f in sorted(in_path.iterdir()) if f.is_file()]
    for in_file in candidates:
        if in_file.suffix.lower() not in exts:
            continue
        out_file = out_path / f"{in_file.stem}.{out_format}"
        log_info(f"Converting {in_file.name} -> {out_file.name}")
        _convert_file(in_file, out_file, in_format, out_format)
