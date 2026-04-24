#!/usr/bin/env python3
"""GFF3 to 6-column BED.

Recurses into sub-features so nested `gene / mRNA / exon / CDS` hierarchies
are all emitted. Strand is written as `+`, `-`, or `.` (unknown) per the BED
spec — `0` / `None` strands are explicitly mapped to `.`, not silently to
`+`.
"""

from pathlib import Path
from typing import Optional

from BCBio import GFF

from ._common import iter_input_files, log_info, prepare_output_dir


def _strand_char(strand) -> str:
    if strand == 1:
        return "+"
    if strand == -1:
        return "-"
    return "."


def _process_feature(rec_id, feature, out_handle):
    name = feature.qualifiers.get("ID", [feature.type])[0]
    score = feature.qualifiers.get("score", ["0"])[0]
    strand = _strand_char(feature.location.strand)
    start = int(feature.location.start)
    end = int(feature.location.end)
    out_handle.write(f"{rec_id}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")
    for sub in feature.sub_features:
        _process_feature(rec_id, sub, out_handle)


def _convert_file(in_file: Path, out_file: Path) -> None:
    with open(in_file) as in_handle, open(out_file, "w") as out_handle:
        for rec in GFF.parse(in_handle):
            for feature in rec.features:
                _process_feature(rec.id, feature, out_handle)


def batch_convert(input_dir: str, output_dir: str,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)
    for gff in iter_input_files(in_path, [".gff3", ".gff"], pattern=pattern):
        out_file = out_path / (gff.stem + ".bed")
        log_info(f"Converting {gff.name} -> {out_file.name}")
        _convert_file(gff, out_file)
