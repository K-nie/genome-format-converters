#!/usr/bin/env python3
"""
GFF3 to tab‑separated feature table (TSV).
Includes all nested features via recursion.
"""

import csv
from pathlib import Path
from BCBio import GFF

def _write_feature(rec_id, feature, writer):
    """Write a single feature as a TSV row."""
    row = [
        rec_id,
        feature.qualifiers.get('source', [''])[0],
        feature.type,
        feature.location.start + 1,
        feature.location.end,
        feature.qualifiers.get('score', ['.'])[0],
        '+' if feature.location.strand == 1 else '-',
        feature.qualifiers.get('phase', ['.'])[0],
        feature.qualifiers.get('ID', [''])[0],
        feature.qualifiers.get('Name', [''])[0],
        ';'.join(feature.qualifiers.get('Parent', [])),
        feature.qualifiers.get('product', [''])[0],
    ]
    writer.writerow(row)

def _process_features(rec_id, feature, writer):
    """Recursively write a feature and its subfeatures."""
    _write_feature(rec_id, feature, writer)
    for sub in feature.sub_features:
        _process_features(rec_id, sub, writer)

def _convert_file(in_file: Path, out_file: Path) -> None:
    """Convert a single GFF3 file to TSV table."""
    with open(in_file) as in_handle, open(out_file, 'w', newline='') as out_handle:
        writer = csv.writer(out_handle, delimiter='\t')
        writer.writerow(["seqid", "source", "type", "start", "end", "score",
                         "strand", "phase", "ID", "Name", "Parent", "product"])
        for rec in GFF.parse(in_handle):
            for feature in rec.features:
                _process_features(rec.id, feature, writer)

def batch_convert(input_dir: str, output_dir: str) -> None:
    """Convert all GFF3 files in input_dir to TSV tables in output_dir."""
    in_path = Path(input_dir)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    for gff in in_path.glob("*.gff*"):
        out_file = out_path / gff.with_suffix(".tsv").name
        print(f"Converting {gff.name} -> {out_file.name}")
        _convert_file(gff, out_file)