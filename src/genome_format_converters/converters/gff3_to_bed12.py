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


def _is_gene_type(ftype: str) -> bool:
    return ftype == "gene" or ftype == "pseudogene" or ftype.endswith("_gene")


def _convert_file(in_file: Path, out_file: Path) -> None:
    """One BED12 row per transcript.

    A transcript is a feature of a known transcript type, or any feature
    that is the direct Parent of an exon or CDS, whatever its type: Ensembl
    uses pseudogenic_transcript, miRNA, unconfirmed_transcript,
    V_gene_segment and others, and tRNAscan-SE hangs exons straight on tRNA
    ``pseudogene`` rows (gffread follows the same rule). A gene-like row
    with no transcript child is written as a single row of its own.
    """
    # feats[fid] = (ftype, rec_id, start1, end, strand, score, parent)
    feats: Dict[str, Tuple[str, str, int, int, str, str, Optional[str]]] = {}
    order: List[str] = []
    exons_by_parent: Dict[str, List[Tuple[int, int]]] = {}
    cds_by_parent: Dict[str, List[Tuple[int, int]]] = {}

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

            if ftype == "exon" or ftype == "CDS":
                if parent is not None:
                    d = exons_by_parent if ftype == "exon" else cds_by_parent
                    lst = d.get(parent)
                    if lst is None:
                        d[parent] = [(start1, end)]
                    else:
                        lst.append((start1, end))
                continue
            if fid is not None and fid not in feats:
                feats[fid] = (ftype, cols[0], start1, end, strand,
                              cols[5] if cols[5] != "." else "0", parent)
                order.append(fid)

    def is_transcript(fid: str) -> bool:
        ftype = feats[fid][0]
        return ftype in _TRANSCRIPT_TYPES or fid in exons_by_parent or fid in cds_by_parent

    tx_ids = [f for f in order if is_transcript(f)]
    tx_set = set(tx_ids)
    has_tx_child = {feats[t][6] for t in tx_ids if feats[t][6] is not None}

    with open(out_file, "w", encoding="utf-8") as out_handle:
        write = out_handle.write
        # 1. One BED12 per transcript, in encounter order.
        for tid in tx_ids:
            _ftype, rec_id, start1, end, strand, score, _parent = feats[tid]
            exons = exons_by_parent.get(tid, [])
            cdses = cds_by_parent.get(tid, [])
            if not exons and cdses:
                # gene -> CDS without exon rows: the CDS parts are the blocks.
                exons = list(cdses)
            _emit_bed12(rec_id, start1, end, strand, tid, score,
                        exons, cdses, write)

        # 2. Gene-like rows with no transcript child, as one block each.
        for gid in order:
            ftype, rec_id, start1, end, strand, score, _parent = feats[gid]
            if gid in tx_set or gid in has_tx_child or not _is_gene_type(ftype):
                continue
            _emit_bed12(rec_id, start1, end, strand, gid, score, [], [], write)


def batch_convert(input_dir: str, output_dir: str,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)
    for gff in iter_input_files(in_path, [".gff3", ".gff"], pattern=pattern):
        out_file = out_path / (gff.stem + ".bed12")
        log_info(f"Converting {gff.name} -> {out_file.name}")
        _convert_file(gff, out_file)
