#!/usr/bin/env python3
"""
GFF3 to GTF conversion.
"""

from pathlib import Path
from BCBio import GFF

def _process_feature(rec_id, feature, gene_id, transcript_id, out_handle):
    """
    Recursively write a feature and its subfeatures in GTF format.
    - gene_id: the gene ID inherited from the top-level gene.
    - transcript_id: the transcript ID inherited from the parent mRNA (None for non-transcript features).
    """
    # Determine the transcript_id for this feature
    if feature.type in ["mRNA", "transcript"]:
        # This feature defines a new transcript; its own ID becomes the transcript_id for itself and its children
        current_transcript_id = feature.qualifiers.get('ID', [''])[0]
    else:
        # Use the transcript_id passed from the parent (if any)
        current_transcript_id = transcript_id

    # For non-transcript features without a transcript context, we may still need a transcript_id.
    # If current_transcript_id is None, fall back to the feature's own ID (or leave empty)
    if current_transcript_id is None:
        current_transcript_id = feature.qualifiers.get('ID', [''])[0]

    # Build attributes
    attrs = f'gene_id "{gene_id}"; transcript_id "{current_transcript_id}";'
    strand_char = '-' if feature.location.strand == -1 else '+'

    out_handle.write(
        f"{rec_id}\t{feature.qualifiers.get('source', ['GFF'])[0]}\t"
        f"{feature.type}\t{feature.location.start+1}\t{feature.location.end}\t"
        f"{feature.qualifiers.get('score', ['.'])[0]}\t"
        f"{strand_char}\t"
        f"{feature.qualifiers.get('phase', ['.'])[0]}\t{attrs}\n"
    )

    # Recursively process subfeatures, passing the same gene_id and the current transcript_id
    for sub in feature.sub_features:
        _process_feature(rec_id, sub, gene_id, current_transcript_id, out_handle)

def _convert_file(in_file: Path, out_file: Path) -> None:
    """Convert a single GFF3 file to GTF."""
    with open(in_file) as in_handle, open(out_file, 'w') as out_handle:
        for rec in GFF.parse(in_handle):
            for feature in rec.features:
                # For top-level features, gene_id is derived from the feature's own ID
                gene_id = feature.qualifiers.get('ID', [''])[0].split(':')[0]
                _process_feature(rec.id, feature, gene_id, None, out_handle)

def batch_convert(input_dir: str, output_dir: str) -> None:
    """Convert all GFF3 files in input_dir to GTF in output_dir."""
    in_path = Path(input_dir)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    for gff in in_path.glob("*.gff*"):
        out_file = out_path / gff.with_suffix(".gtf").name
        print(f"Converting {gff.name} -> {out_file.name}")
        _convert_file(gff, out_file)