#!/usr/bin/env python3
"""GFF3 to a tab-separated feature table (includes nested features)."""

import csv
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


def _write_feature(rec_id, feature, writer):
    row = [
        rec_id,
        feature.qualifiers.get("source", [""])[0],
        feature.type,
        int(feature.location.start) + 1,
        int(feature.location.end),
        feature.qualifiers.get("score", ["."])[0],
        _strand_char(feature.location.strand),
        feature.qualifiers.get("phase", ["."])[0],
        feature.qualifiers.get("ID", [""])[0],
        feature.qualifiers.get("Name", [""])[0],
        ";".join(feature.qualifiers.get("Parent", [])),
        feature.qualifiers.get("product", [""])[0],
    ]
    writer.writerow(row)


def _process_features(rec_id, feature, writer):
    _write_feature(rec_id, feature, writer)
    for sub in feature.sub_features:
        _process_features(rec_id, sub, writer)


def _convert_file(in_file: Path, out_file: Path) -> None:
    with open(in_file) as in_handle, open(out_file, "w", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t")
        writer.writerow([
            "seqid", "source", "type", "start", "end", "score",
            "strand", "phase", "ID", "Name", "Parent", "product",
        ])
        for rec in GFF.parse(in_handle):
            for feature in rec.features:
                _process_features(rec.id, feature, writer)


def batch_convert(input_dir: str, output_dir: str,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)
    for gff in iter_input_files(in_path, [".gff3", ".gff"], pattern=pattern):
        out_file = out_path / (gff.stem + ".tsv")
        log_info(f"Converting {gff.name} -> {out_file.name}")
        _convert_file(gff, out_file)
