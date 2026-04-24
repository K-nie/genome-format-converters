#!/usr/bin/env python3
"""GFF3 → GTF (Ensembl-style).

Emits `gene_id` on every row and `transcript_id` on every row that lives
under an `mRNA` / `transcript` parent. Feature IDs are preserved verbatim —
we don't split on `:` because Ensembl-style IDs like `gene:YAL001C` legitimately
contain colons and splitting produces duplicate IDs across genes.
"""

from pathlib import Path
from typing import Optional

from BCBio import GFF

from ._common import iter_input_files, log_info, prepare_output_dir

_TRANSCRIPT_TYPES = {"mRNA", "transcript", "ncRNA", "lnc_RNA", "tRNA", "rRNA", "snoRNA", "snRNA"}


def _strand_char(strand) -> str:
    if strand == 1:
        return "+"
    if strand == -1:
        return "-"
    return "."


def _clean_id(raw: str) -> str:
    # Preserve the ID. GTF has no convention for URL-encoded or colon-containing
    # IDs; `feature.qualifiers["ID"]` already arrives unescaped from bcbio-gff.
    return raw.strip()


def _process_feature(rec_id, feature, gene_id, transcript_id, out_handle):
    feat_id = _clean_id(feature.qualifiers.get("ID", [""])[0])

    if feature.type in _TRANSCRIPT_TYPES:
        current_transcript = feat_id or transcript_id
    else:
        current_transcript = transcript_id

    # GTF attribute string: always emit gene_id; emit transcript_id when we
    # know it. For bare `gene` rows with no transcript context we leave
    # transcript_id off rather than setting it equal to gene_id (which no
    # downstream tool expects).
    attr_parts = [f'gene_id "{gene_id}"']
    if current_transcript:
        attr_parts.append(f'transcript_id "{current_transcript}"')
    attrs = "; ".join(attr_parts) + ";"

    source = feature.qualifiers.get("source", ["GFF"])[0]
    score = feature.qualifiers.get("score", ["."])[0]
    phase = feature.qualifiers.get("phase", ["."])[0]
    strand = _strand_char(feature.location.strand)

    out_handle.write(
        f"{rec_id}\t{source}\t{feature.type}\t"
        f"{int(feature.location.start) + 1}\t{int(feature.location.end)}\t"
        f"{score}\t{strand}\t{phase}\t{attrs}\n"
    )

    for sub in feature.sub_features:
        _process_feature(rec_id, sub, gene_id, current_transcript, out_handle)


def _convert_file(in_file: Path, out_file: Path) -> None:
    with open(in_file) as in_handle, open(out_file, "w") as out_handle:
        for rec in GFF.parse(in_handle):
            for feature in rec.features:
                gene_id = _clean_id(feature.qualifiers.get("ID", [feature.type])[0])
                _process_feature(rec.id, feature, gene_id, None, out_handle)


def batch_convert(input_dir: str, output_dir: str,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)
    for gff in iter_input_files(in_path, [".gff3", ".gff"], pattern=pattern):
        out_file = out_path / (gff.stem + ".gtf")
        log_info(f"Converting {gff.name} -> {out_file.name}")
        _convert_file(gff, out_file)
