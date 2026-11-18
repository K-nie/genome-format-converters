#!/usr/bin/env python3
"""Convert MAF to XMFA (the Mauve / progressiveMauve alignment format).

Each MAF block becomes one XMFA block, terminated by `=`. Every distinct
MAF source name gets a sequence index (1, 2, ... in order of first
appearance), listed in `#SequenceNFile` header lines. Sequence headers
follow Mauve: `> index:start-end strand name`, with 1-based inclusive
coordinates on the forward strand of the source. MAF gives minus-strand
starts on the reverse complement, so those are converted:
`start = srcSize - mafStart - size + 1`, `end = srcSize - mafStart`.
Alignment columns, including gaps, are kept as they are.
"""

from pathlib import Path
from typing import Dict, Optional

from Bio import AlignIO

from ._common import iter_input_files, log_info, prepare_output_dir

_WRAP = 80


def _source_names(in_file: Path) -> Dict[str, int]:
    """Map each MAF source name to its XMFA index, in file order."""
    index: Dict[str, int] = {}
    with open(in_file) as fh:
        for line in fh:
            if line.startswith("s "):
                name = line.split(None, 2)[1]
                index.setdefault(name, len(index) + 1)
    return index


def _convert_file(in_file: Path, out_file: Path) -> None:
    index = _source_names(in_file)
    with open(in_file) as in_handle, open(out_file, "w") as out_handle:
        out_handle.write("#FormatVersion Mauve1\n")
        for name, i in index.items():
            out_handle.write(f"#Sequence{i}File\t{name}\n")
        for aln in AlignIO.parse(in_handle, "maf"):
            for seq in aln:
                ann = seq.annotations
                start, size = ann["start"], ann["size"]
                if ann.get("strand") in (-1, "-"):
                    lo, hi, strand = ann["srcSize"] - start - size + 1, ann["srcSize"] - start, "-"
                else:
                    lo, hi, strand = start + 1, start + size, "+"
                out_handle.write(f"> {index[seq.id]}:{lo}-{hi} {strand} {seq.id}\n")
                text = str(seq.seq)
                for k in range(0, len(text), _WRAP):
                    out_handle.write(text[k:k + _WRAP] + "\n")
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
