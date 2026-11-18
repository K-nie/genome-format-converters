#!/usr/bin/env python3
"""Convert alignment files between formats (FASTA / PHYLIP / NEXUS / CLUSTAL).

`phylip` is strict PHYLIP (names cut to 10 characters); `phylip-relaxed`
keeps full names, separated from the sequence by whitespace. PHYLIP has no
`.` gap symbol, so `.` (as in Stockholm-derived aligned FASTA) is written
as `-`. When cutting names to 10 characters would make two rows share a
name, strict PHYLIP output is refused with a pointer to `phylip-relaxed`.
"""

import os
from pathlib import Path
from typing import Optional

from Bio import AlignIO
from Bio.Seq import Seq

from ._common import die, log_info, prepare_output_dir

# Biopython format name -> common file extensions for input discovery.
_FORMAT_EXTS = {
    "fasta": [".fasta", ".fa", ".fna", ".fas", ".mfa", ".aln"],
    "phylip": [".phy", ".phylip"],
    "phylip-relaxed": [".phy", ".phylip"],
    "nexus": [".nex", ".nexus"],
    "clustal": [".aln", ".clustal", ".clw"],
}


def _convert_file(in_file: Path, out_file: Path, in_fmt: str, out_fmt: str) -> None:
    with open(in_file) as in_handle:
        aln = AlignIO.read(in_handle, in_fmt)
    if out_fmt.startswith("phylip"):
        for rec in aln:
            if "." in rec.seq:
                rec.seq = Seq(str(rec.seq).replace(".", "-"))
        if out_fmt == "phylip":
            seen = {}
            for rec in aln:
                short = rec.id[:10]
                if short in seen and seen[short] != rec.id:
                    die(f"{in_file.name}: {seen[short]!r} and {rec.id!r} both become "
                        f"{short!r} in strict PHYLIP; use --out-format phylip-relaxed "
                        "to keep full names")
                seen[short] = rec.id
    # Write to a temporary name so a failure never leaves a partial file.
    tmp = out_file.with_name(out_file.name + ".tmp")
    with open(tmp, "w") as out_handle:
        AlignIO.write(aln, out_handle, out_fmt)
    os.replace(tmp, out_file)


def batch_convert(input_dir: str, output_dir: str,
                  in_format: str, out_format: str,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)

    exts = {e.lower() for e in _FORMAT_EXTS.get(in_format, [])}
    if pattern:
        candidates = sorted(in_path.glob(pattern))
    else:
        candidates = [f for f in sorted(in_path.iterdir()) if f.is_file()]

    for in_file in candidates:
        # When --pattern is supplied, trust the user's glob; otherwise require
        # a recognised extension for the declared input format so we don't
        # silently try to parse unrelated files.
        if not pattern and exts and in_file.suffix.lower() not in exts:
            continue
        ext = "phylip" if out_format.startswith("phylip") else out_format
        out_file = out_path / f"{in_file.stem}.{ext}"
        log_info(f"Converting {in_file.name} -> {out_file.name}")
        _convert_file(in_file, out_file, in_format, out_format)
