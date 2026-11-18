#!/usr/bin/env python3
"""Convert phylogenetic tree formats (Newick, NEXUS, PhyloXML).

Branch lengths and support values are written with up to 12 significant
digits, so values from the input survive unrounded (Biopython's defaults
would cut supports to two decimals). NEXUS output is written here rather
than by Biopython, whose writer leaves taxon names with punctuation or
spaces unquoted: such names (for example Pfam's `NAME/start-end`) are put
in single quotes in both the taxa block and the tree.
"""

import re
from io import StringIO
from pathlib import Path
from typing import List, Optional

from Bio import Phylo

from ._common import log_info, log_warn, prepare_output_dir
from .annotate_tree import _nexus_name

_NUM = {"format_branch_length": "%.12g", "format_confidence": "%.12g"}
_TOKEN = re.compile(r"gfcTAXON(\d+)x")

_EXT_MAP = {
    "newick": [".nwk", ".tree", ".tre", ".newick"],
    "nexus": [".nex", ".nexus"],
    "phyloxml": [".xml", ".phyloxml"],
}


def _write_nexus(trees: List, out_file: Path) -> None:
    taxa: List[str] = []
    lines = []
    for n, tree in enumerate(trees, 1):
        names: List[str] = []
        for clade in tree.find_clades():
            if clade.name:
                names.append(clade.name)
                clade.name = f"gfcTAXON{len(names) - 1}x"
        buf = StringIO()
        Phylo.write(tree, buf, "newick", **_NUM)
        text = _TOKEN.sub(lambda m: _nexus_name(names[int(m.group(1))]), buf.getvalue().strip())
        lines.append(f" Tree tree{n}={text}")
        for t in tree.get_terminals():
            name = names[int(_TOKEN.fullmatch(t.name).group(1))] if t.name else None
            if name and name not in taxa:
                taxa.append(name)
    with open(out_file, "w") as out:
        out.write("#NEXUS\nBegin Taxa;\n")
        out.write(f" Dimensions NTax={len(taxa)};\n")
        out.write(" TaxLabels " + " ".join(_nexus_name(t) for t in taxa) + ";\nEnd;\n")
        out.write("Begin Trees;\n" + "\n".join(lines) + "\nEnd;\n")


def _convert_file(in_file: Path, out_file: Path, in_fmt: str, out_fmt: str) -> None:
    trees = list(Phylo.parse(str(in_file), in_fmt))
    if in_fmt == "nexus":
        # Bio.Phylo keeps NEXUS quotes on names; unwrap them so they are
        # not quoted twice on output.
        for tree in trees:
            for clade in tree.find_clades():
                n = clade.name
                if n and len(n) >= 2 and n[0] == n[-1] == "'":
                    clade.name = n[1:-1].replace("''", "'")
    if out_fmt == "nexus":
        _write_nexus(trees, out_file)
    elif out_fmt == "newick":
        Phylo.write(trees, str(out_file), out_fmt, **_NUM)
    else:
        Phylo.write(trees, str(out_file), out_fmt)


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
