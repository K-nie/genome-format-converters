#!/usr/bin/env python3
"""GFF3 -> BED12.

Author: Benjamin Narh-Madey

Emits one BED12 line per transcript (``mRNA`` or similar), with exon
children collapsed into the ``blockCount`` / ``blockSizes`` /
``blockStarts`` fields. The thick region (CDS) is inferred from CDS
children when present; otherwise ``thickStart == thickEnd == chromStart``
so IGV renders the whole transcript as UTR.

For GFF3 inputs that contain only ``gene`` -> ``CDS`` hierarchies (no
explicit ``exon`` features, common in some yeast annotations), CDS spans
are used as exon blocks instead.

Implementation note (perf): we replaced the BCBio.GFF parser with a
two-pass streaming text parser. Pass 1 records every transcript and its
exon/CDS children indexed by Parent ID; pass 2 emits one BED12 row per
transcript directly from the indexed spans. This avoids per-feature
SeqFeature object construction and recursive sub_features walking, which
on real Y1000+ annotations dominate runtime.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple

from ._common import iter_input_files, log_info, prepare_output_dir

_TRANSCRIPT_TYPES = frozenset({
    "mRNA", "transcript", "ncRNA", "lnc_RNA",
    "tRNA", "rRNA", "snoRNA", "snRNA",
})

# A "row" we keep per transcript and per exon/CDS child. Tuples are far
# cheaper than dataclasses or named-tuples in tight inner loops.
# Transcript record: (rec_id, start1, end, strand, name, score, fid)
# Block record: (start1, end)


def _strand_char(s: str) -> str:
    if s == "+" or s == "-":
        return s
    return "."


def _parse_attrs_min(attr_field: str) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """Extract ID, Parent (first), and score (rare on attrs) from column 9."""
    if not attr_field or attr_field == ".":
        return None, None, None
    fid: Optional[str] = None
    parent: Optional[str] = None
    for item in attr_field.split(";"):
        if not item:
            continue
        if item[0] == " ":
            item = item.lstrip()
            if not item:
                continue
        eq = item.find("=")
        if eq < 1:
            continue
        key = item[:eq]
        if key == "ID":
            fid = item[eq + 1:].rstrip()
        elif key == "Parent":
            v = item[eq + 1:].rstrip()
            comma = v.find(",")
            parent = v[:comma] if comma >= 0 else v
    return fid, parent, None


def _emit_bed12(
    rec_id: str, start1: int, end: int, strand: str, name: str, score: str,
    exons: List[Tuple[int, int]], cdses: List[Tuple[int, int]],
    write,
) -> None:
    """Write one BED12 row.

    ``start1`` / ``end`` are GFF3 1-based-inclusive; we convert to 0-based
    half-open BED here. ``exons`` / ``cdses`` likewise.
    """
    if exons:
        # Sort once; ascending by start. ``list.sort`` is in-place and stable.
        exons.sort()
        chrom_start = exons[0][0] - 1
        chrom_end = exons[-1][1]
    else:
        # Single-block: span the whole transcript.
        chrom_start = start1 - 1
        chrom_end = end
        exons = [(start1, end)]

    if cdses:
        cdses.sort()
        thick_start = cdses[0][0] - 1
        thick_end = cdses[-1][1]
    else:
        thick_start = chrom_start
        thick_end = chrom_start

    n = len(exons)
    if n == 1:
        s0, e0 = exons[0]
        block_sizes = f"{e0 - s0 + 1},"
        block_starts = "0,"
    else:
        # Local accumulators - one f-string per block, joined once.
        sizes_parts = []
        starts_parts = []
        for s, e in exons:
            sizes_parts.append(str(e - s + 1))
            starts_parts.append(str(s - 1 - chrom_start))
        block_sizes = ",".join(sizes_parts) + ","
        block_starts = ",".join(starts_parts) + ","

    write(
        f"{rec_id}\t{chrom_start}\t{chrom_end}\t{name}\t{score}\t{strand}\t"
        f"{thick_start}\t{thick_end}\t0\t{n}\t{block_sizes}\t{block_starts}\n"
    )


def _convert_file(in_file: Path, out_file: Path) -> None:
    # Maps:
    #   transcripts[fid] = (rec_id, start1, end, strand, name, score)
    #   exons_by_parent[pid] = list of (start1, end)
    #   cds_by_parent[pid] = list of (start1, end)
    #   gene_singleton[gid] = transcript-shaped tuple, used only when a
    #       gene has no transcript children.
    #   gene_has_transcript[gid] = True if any transcript declares this gene
    #       as parent. We track this so we can decide at emit time whether
    #       to emit the gene as a single-block row.
    transcripts: Dict[str, Tuple[str, int, int, str, str, str]] = {}
    transcript_order: List[str] = []
    exons_by_parent: Dict[str, List[Tuple[int, int]]] = {}
    cds_by_parent: Dict[str, List[Tuple[int, int]]] = {}
    gene_singleton: Dict[str, Tuple[str, int, int, str, str, str]] = {}
    gene_order: List[str] = []
    gene_has_transcript: Dict[str, bool] = {}

    with open(in_file, "r", encoding="utf-8") as in_handle:
        for line in in_handle:
            if not line:
                continue
            first = line[0]
            if first == "#" or first == "\n":
                continue
            if line.endswith("\n"):
                line = line[:-1]
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue
            ftype = cols[2]
            try:
                start1 = int(cols[3])
                end = int(cols[4])
            except ValueError:
                continue
            strand = _strand_char(cols[6])
            fid, parent, _ = _parse_attrs_min(cols[8])

            if ftype == "gene":
                if fid is not None:
                    gene_singleton[fid] = (
                        cols[0], start1, end, strand, fid,
                        cols[5] if cols[5] != "." else "0",
                    )
                    if fid not in gene_has_transcript:
                        gene_has_transcript[fid] = False
                    gene_order.append(fid)
                continue

            if ftype in _TRANSCRIPT_TYPES:
                if fid is not None:
                    transcripts[fid] = (
                        cols[0], start1, end, strand, fid,
                        cols[5] if cols[5] != "." else "0",
                    )
                    transcript_order.append(fid)
                if parent is not None:
                    gene_has_transcript[parent] = True
                continue

            if ftype == "exon" and parent is not None:
                lst = exons_by_parent.get(parent)
                if lst is None:
                    exons_by_parent[parent] = [(start1, end)]
                else:
                    lst.append((start1, end))
                continue

            if ftype == "CDS" and parent is not None:
                lst = cds_by_parent.get(parent)
                if lst is None:
                    cds_by_parent[parent] = [(start1, end)]
                else:
                    lst.append((start1, end))
                continue

    with open(out_file, "w", encoding="utf-8") as out_handle:
        write = out_handle.write
        # 1. Emit one BED12 per transcript, in encounter order.
        for tid in transcript_order:
            rec_id, start1, end, strand, name, score = transcripts[tid]
            exons = exons_by_parent.get(tid, [])
            cdses = cds_by_parent.get(tid, [])
            if not exons and cdses:
                # GFF3 inputs with only gene -> CDS (no exon child) - use
                # CDS as the block list. This matches the bcbio-backed
                # behaviour exactly.
                exons = list(cdses)
            _emit_bed12(rec_id, start1, end, strand, name, score,
                        exons, cdses, write)

        # 2. Emit gene-level singletons for genes with no transcript children
        #    (mirrors the recursive walker that fell through to ``gene``).
        seen_gene = set()
        for gid in gene_order:
            if gid in seen_gene:
                continue
            seen_gene.add(gid)
            if gene_has_transcript.get(gid):
                continue
            rec_id, start1, end, strand, name, score = gene_singleton[gid]
            exons = exons_by_parent.get(gid, [])
            cdses = cds_by_parent.get(gid, [])
            if not exons and cdses:
                exons = list(cdses)
            _emit_bed12(rec_id, start1, end, strand, name, score,
                        exons, cdses, write)


def batch_convert(input_dir: str, output_dir: str,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)
    for gff in iter_input_files(in_path, [".gff3", ".gff"], pattern=pattern):
        out_file = out_path / (gff.stem + ".bed12")
        log_info(f"Converting {gff.name} -> {out_file.name}")
        _convert_file(gff, out_file)
