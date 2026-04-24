#!/usr/bin/env python3
"""GFF3 → BED12.

Emits one BED12 line per transcript (``mRNA`` or similar), with exon
children collapsed into the ``blockCount`` / ``blockSizes`` /
``blockStarts`` fields. The thick region (CDS) is inferred from CDS
children when present; otherwise ``thickStart == thickEnd == chromStart``
so IGV renders the whole transcript as UTR.

For GFF3 inputs that contain only ``gene`` → ``CDS`` hierarchies (no
explicit ``exon`` features, common in some yeast annotations), CDS spans
are used as exon blocks instead.
"""

from pathlib import Path
from typing import List, Optional, Tuple

from BCBio import GFF

from ._common import iter_input_files, log_info, prepare_output_dir

_TRANSCRIPT_TYPES = {"mRNA", "transcript", "ncRNA", "lnc_RNA", "tRNA", "rRNA",
                     "snoRNA", "snRNA"}


def _strand_char(strand) -> str:
    if strand == 1:
        return "+"
    if strand == -1:
        return "-"
    return "."


def _collect_by_type(feature, target_types):
    out = []
    stack = [feature]
    while stack:
        f = stack.pop()
        if f.type in target_types:
            out.append(f)
        stack.extend(f.sub_features)
    return out


def _emit_bed12(rec_id: str, transcript, out_handle) -> None:
    exons = _collect_by_type(transcript, {"exon"})
    if not exons:
        # Fall back to CDS as exon blocks when no exon child is declared.
        exons = _collect_by_type(transcript, {"CDS"})
    if not exons:
        # Not a splice-bearing feature — emit a single-block BED12 line.
        start = int(transcript.location.start)
        end = int(transcript.location.end)
        name = transcript.qualifiers.get("ID", [transcript.type])[0]
        score = transcript.qualifiers.get("score", ["0"])[0]
        strand = _strand_char(transcript.location.strand)
        out_handle.write(
            f"{rec_id}\t{start}\t{end}\t{name}\t{score}\t{strand}\t"
            f"{start}\t{end}\t0\t1\t{end - start},\t0,\n"
        )
        return

    exons = sorted(exons, key=lambda e: int(e.location.start))
    chrom_start = int(exons[0].location.start)
    chrom_end = int(exons[-1].location.end)

    cds = _collect_by_type(transcript, {"CDS"})
    if cds:
        cds = sorted(cds, key=lambda c: int(c.location.start))
        thick_start = int(cds[0].location.start)
        thick_end = int(cds[-1].location.end)
    else:
        thick_start = chrom_start
        thick_end = chrom_start

    block_sizes = [int(e.location.end) - int(e.location.start) for e in exons]
    block_starts = [int(e.location.start) - chrom_start for e in exons]

    name = transcript.qualifiers.get("ID", [transcript.type])[0]
    score = transcript.qualifiers.get("score", ["0"])[0]
    strand = _strand_char(transcript.location.strand)

    out_handle.write(
        f"{rec_id}\t{chrom_start}\t{chrom_end}\t{name}\t{score}\t{strand}\t"
        f"{thick_start}\t{thick_end}\t0\t{len(exons)}\t"
        f"{','.join(str(s) for s in block_sizes)},\t"
        f"{','.join(str(s) for s in block_starts)},\n"
    )


def _walk(rec_id: str, feature, out_handle) -> None:
    if feature.type in _TRANSCRIPT_TYPES:
        _emit_bed12(rec_id, feature, out_handle)
        return
    # Gene-level features with no transcript children: emit directly.
    if feature.type == "gene" and not any(
        sub.type in _TRANSCRIPT_TYPES for sub in feature.sub_features
    ):
        _emit_bed12(rec_id, feature, out_handle)
        return
    for sub in feature.sub_features:
        _walk(rec_id, sub, out_handle)


def _convert_file(in_file: Path, out_file: Path) -> None:
    with open(in_file) as in_handle, open(out_file, "w") as out_handle:
        for rec in GFF.parse(in_handle):
            for feature in rec.features:
                _walk(rec.id, feature, out_handle)


def batch_convert(input_dir: str, output_dir: str,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)
    for gff in iter_input_files(in_path, [".gff3", ".gff"], pattern=pattern):
        out_file = out_path / (gff.stem + ".bed12")
        log_info(f"Converting {gff.name} -> {out_file.name}")
        _convert_file(gff, out_file)
