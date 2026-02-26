#!/usr/bin/env python3
"""
GFF3 to 6‑column BED conversion.
"""

from pathlib import Path
from BCBio import GFF

def _process_feature(rec_id, feature, out_handle):
    """Recursively write a feature and its subfeatures in BED format."""
    name = feature.qualifiers.get('ID', [feature.type])[0]
    score = feature.qualifiers.get('score', [0])[0]
    strand = '+' if feature.location.strand == 1 else '-'
    start = feature.location.start
    end = feature.location.end
    out_handle.write(f"{rec_id}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")
    for sub in feature.sub_features:
        _process_feature(rec_id, sub, out_handle)

def _convert_file(in_file: Path, out_file: Path) -> None:
    """Convert a single GFF3 file to BED."""
    with open(in_file) as in_handle, open(out_file, 'w') as out_handle:
        for rec in GFF.parse(in_handle):
            for feature in rec.features:
                _process_feature(rec.id, feature, out_handle)

def batch_convert(input_dir: str, output_dir: str) -> None:
    """Convert all GFF3 files in input_dir to BED in output_dir."""
    in_path = Path(input_dir)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    for gff in in_path.glob("*.gff*"):
        out_file = out_path / gff.with_suffix(".bed").name
        print(f"Converting {gff.name} -> {out_file.name}")
        _convert_file(gff, out_file)