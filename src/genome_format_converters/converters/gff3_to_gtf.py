#!/usr/bin/env python3
"""GFF3 -> GTF (Ensembl-style).

Author: Benjamin Narh-Madey

Emits ``gene_id`` on every row and ``transcript_id`` on every row that lives
under an ``mRNA`` / ``transcript`` parent. Feature IDs are preserved verbatim
- we don't split on ``:`` because Ensembl-style IDs like ``gene:YAL001C``
legitimately contain colons and splitting produces duplicate IDs across genes.

Implementation note (perf): the GFF3 -> GTF conversion is structurally a
column-9 rewrite plus a parent-of-parent walk to find each row's gene_id /
transcript_id. We do this with a streaming text parser instead of going
through ``BCBio.GFF`` / ``gffutils`` because:

    * those libraries build a full SeqFeature graph with per-row dict
      allocation and recursive ``sub_features`` walking, which dominates
      runtime on large eukaryotic GFFs;
    * a streaming pass needs only an ``ID -> Parent`` map and an
      ``ID -> type`` map of ancestors, which fits in memory for any
      genome-scale annotation.

Two passes over the file are still cheaper than one bcbio-gff pass.
"""

from pathlib import Path
from typing import Dict, Optional, Tuple

from ._common import iter_input_files, log_info, prepare_output_dir

_TRANSCRIPT_TYPES = frozenset({
    "mRNA", "transcript", "ncRNA", "lnc_RNA",
    "tRNA", "rRNA", "snoRNA", "snRNA",
})


def _parse_attrs_min(attr_field: str) -> Tuple[Optional[str], Optional[str]]:
    """Extract just (ID, Parent) from a GFF3 column-9 string.

    GFF3 attributes follow ``key=value;key=value`` (RFC 3986 percent-encoded
    on values). For the ID/Parent extraction we only need the literal token,
    so we skip percent-decoding here - GTF emitters don't decode either, and
    round-trip identity is what downstream tools (Ensembl, Cufflinks, IGV)
    expect on these particular keys.
    """
    if not attr_field or attr_field == ".":
        return None, None
    fid: Optional[str] = None
    parent: Optional[str] = None
    # split() is faster than a regex here and the format guarantees ``;`` as
    # separator and ``=`` as key/value delimiter.
    for item in attr_field.split(";"):
        # Skip leading whitespace cheaply without a full strip().
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
            # Multi-parent (``Parent=a,b``) is legal in GFF3 but rare; GTF
            # has no representation for it, so we take the first parent.
            v = item[eq + 1:].rstrip()
            comma = v.find(",")
            parent = v[:comma] if comma >= 0 else v
    return fid, parent


def _resolve_gene_transcript(
    fid: Optional[str],
    ftype: str,
    parent: Optional[str],
    id_to_type: Dict[str, str],
    id_to_parent: Dict[str, Optional[str]],
) -> Tuple[Optional[str], Optional[str]]:
    """Walk the Parent chain to find ``gene_id`` and ``transcript_id``.

    We climb at most a few hops (gene -> mRNA -> exon/CDS), so the linear
    walk is cheap compared to an actual graph build.
    """
    gene_id: Optional[str] = None
    transcript_id: Optional[str] = None

    if ftype == "gene":
        return fid, None
    if ftype in _TRANSCRIPT_TYPES:
        transcript_id = fid
        # gene is the parent if it's a gene; else the row itself acts as gene.
        if parent and id_to_type.get(parent) == "gene":
            gene_id = parent
        else:
            gene_id = parent or fid
        return gene_id, transcript_id

    # Non-transcript, non-gene row (exon, CDS, UTR, ...). Climb the chain.
    cur = parent
    hops = 0
    while cur is not None and hops < 8:
        ctype = id_to_type.get(cur)
        if ctype is None:
            # Unknown ancestor; treat the immediate parent as transcript and
            # stop. This mirrors what gffread does on partial files.
            if transcript_id is None:
                transcript_id = cur
            break
        if ctype in _TRANSCRIPT_TYPES and transcript_id is None:
            transcript_id = cur
        if ctype == "gene":
            gene_id = cur
            break
        cur = id_to_parent.get(cur)
        hops += 1

    if gene_id is None:
        # No gene ancestor found - fall back to the transcript ID, then
        # whatever we have. Matches biopython-based reference behaviour.
        gene_id = transcript_id or parent or fid
    return gene_id, transcript_id


def _convert_file(in_file: Path, out_file: Path) -> None:
    # ---- Pass 1: build ID -> type and ID -> Parent maps. -------------------
    # Both maps are bounded by feature-instance count (a few 100k for a
    # eukaryotic annotation), so a plain dict is fine. We store None for
    # missing parents so a single ``.get`` answers both "exists" and "is root".
    id_to_type: Dict[str, str] = {}
    id_to_parent: Dict[str, Optional[str]] = {}

    with open(in_file, "r", encoding="utf-8") as in_handle:
        for line in in_handle:
            if not line or line[0] == "#":
                continue
            # rstrip just the newline; preserving any trailing tabs/spaces
            # in the attribute column doesn't matter because we only read
            # ID/Parent here.
            if line.endswith("\n"):
                line = line[:-1]
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue
            ftype = cols[2]
            fid, parent = _parse_attrs_min(cols[8])
            if fid is not None:
                id_to_type[fid] = ftype
                id_to_parent[fid] = parent

    # ---- Pass 2: stream-rewrite each row's column 9. -----------------------
    # We build a small attr buffer per row and write into a buffered text
    # stream. ``write`` calls into a buffered file handle batch into a single
    # OS write, so we don't need a manual list-and-join.
    with open(in_file, "r", encoding="utf-8") as in_handle, \
            open(out_file, "w", encoding="utf-8") as out_handle:
        # Local aliases - ~15% inner-loop speedup on CPython 3.11+ because
        # the bytecode skips the LOAD_GLOBAL.
        _write = out_handle.write
        _parse = _parse_attrs_min
        _resolve = _resolve_gene_transcript
        _types = _TRANSCRIPT_TYPES  # noqa: F841 - referenced in _resolve

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
            fid, parent = _parse(cols[8])
            gene_id, transcript_id = _resolve(
                fid, ftype, parent, id_to_type, id_to_parent
            )

            # Ensembl/GTF attribute string. Always emit gene_id; emit
            # transcript_id when known; bare ``gene`` rows omit transcript_id.
            if transcript_id is not None:
                attrs = (
                    f'gene_id "{gene_id}"; '
                    f'transcript_id "{transcript_id}";'
                )
            else:
                attrs = f'gene_id "{gene_id}";'

            # Columns 1-8 are unchanged except: phase for non-CDS rows in
            # GFF3 is often ``.`` already; we leave it. score may also be ``.``.
            _write(
                cols[0]
            )
            _write("\t")
            _write(cols[1])
            _write("\t")
            _write(ftype)
            _write("\t")
            _write(cols[3])
            _write("\t")
            _write(cols[4])
            _write("\t")
            _write(cols[5])
            _write("\t")
            _write(cols[6])
            _write("\t")
            _write(cols[7])
            _write("\t")
            _write(attrs)
            _write("\n")


def batch_convert(input_dir: str, output_dir: str,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)
    for gff in iter_input_files(in_path, [".gff3", ".gff"], pattern=pattern):
        out_file = out_path / (gff.stem + ".gtf")
        log_info(f"Converting {gff.name} -> {out_file.name}")
        _convert_file(gff, out_file)
