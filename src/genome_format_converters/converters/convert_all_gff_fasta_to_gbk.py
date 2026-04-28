#!/usr/bin/env python3
"""Convert paired FASTA + GFF3 files to GenBank format.

Author: Benjamin Narh-Madey

Files are paired by stem, with common yeast-suffix strippers (`.final`,
`_final`, `final`). Strand is mapped to `+1` / `-1` / `0`, matching BioPython's
SeqFeature convention - `0` / `.` strands are preserved rather than silently
forced to `+`.

Implementation note (perf): the GFF -> GenBank pipeline is dominated by
two costs - per-row attribute parsing and per-row SeqFeature allocation.
This module keeps biopython's SeqIO.write for the GenBank emit step
(rewriting the GenBank serialiser would be a much larger change), but
trims the inner loop to the minimum:

    * single-pass line.split("\t", 8) avoids the strip() + split allocation
      pair and the parts[:9] slice;
    * _parse_gff3_attributes uses a single bytes-level find()/lstrip()
      pattern instead of strip() per item;
    * a per-seqid feature-list local hoists the records[seqid].features
      attribute lookup out of the inner loop.
"""

from pathlib import Path
from collections import defaultdict
from typing import List, Optional, Tuple

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

from ._common import log_info, log_warn, prepare_output_dir

_FASTA_EXTS = [".fasta", ".fa", ".fna", ".fas"]
_GFF_EXTS = [".gff3", ".gff"]
_SUFFIX_STRIPS = [".final", "_final", "final"]
# Frozen set membership check is O(1); used in the inner loop for the
# "exclude these keys from the note bucket" filter.
_NOTE_EXCLUDE = frozenset({"ID", "Parent", "Name", "product"})
# Pre-encoded codon_start strings - the only valid GFF3 phases.
_CODON_START = {"0": ["1"], "1": ["2"], "2": ["3"]}


def _strip_suffixes(stem: str, suffixes: List[str]) -> str:
    for suf in suffixes:
        if stem.endswith(suf):
            return stem[: -len(suf)]
    return stem


def _find_pairs(input_dir: Path) -> List[Tuple[str, Path, Path]]:
    fasta_by_stem = defaultdict(list)
    for ext in _FASTA_EXTS:
        for f in input_dir.glob(f"*{ext}"):
            fasta_by_stem[f.stem].append(f)

    gff_by_stem = defaultdict(list)
    for ext in _GFF_EXTS:
        for f in input_dir.glob(f"*{ext}"):
            gff_by_stem[_strip_suffixes(f.stem, _SUFFIX_STRIPS)].append(f)

    pairs = []
    for stem, gff_list in gff_by_stem.items():
        if stem in fasta_by_stem:
            pairs.append((stem, fasta_by_stem[stem][0], gff_list[0]))
    return pairs


def _strand_int(strand: str) -> int:
    # Used in the per-row inner loop. The two equality checks land before
    # the implicit "else 0" path, which is the common case for explicit
    # strand columns.
    if strand == "+":
        return 1
    if strand == "-":
        return -1
    return 0


def _gff3_to_gb_qualifiers(ftype: str, attrs: dict) -> dict:
    quals = {}
    has_id = "ID" in attrs
    has_parent = "Parent" in attrs
    if ftype != "gene" and has_parent:
        quals["locus_tag"] = attrs["Parent"]
    elif has_id:
        quals["locus_tag"] = attrs["ID"]
    has_name = "Name" in attrs
    if has_name:
        if ftype == "gene":
            quals["gene"] = attrs["Name"]
        elif ftype == "mRNA":
            quals["transcript_id"] = attrs["Name"]
        elif ftype == "CDS":
            quals["protein_id"] = attrs["Name"]
            if "product" not in attrs:
                quals["product"] = attrs["Name"]
    if ftype == "CDS" and "product" in attrs:
        quals["product"] = attrs["product"]
    # Note bucket: anything not already promoted to a structured qualifier.
    # We build it lazily - most rows have no extra attributes.
    note_parts = None
    for k, v in attrs.items():
        if k not in _NOTE_EXCLUDE:
            if note_parts is None:
                note_parts = [f"{k}:{','.join(v)}"]
            else:
                note_parts.append(f"{k}:{','.join(v)}")
    if note_parts is not None:
        quals["note"] = note_parts
    return quals


def _create_feature(ftype: str, start: int, end: int, strand: str, phase: str, attrs: dict) -> SeqFeature:
    start0 = start - 1
    if start0 < 0:
        start0 = 0
    qualifiers = _gff3_to_gb_qualifiers(ftype, attrs)
    if ftype == "CDS":
        cs = _CODON_START.get(phase)
        if cs is not None:
            qualifiers["codon_start"] = cs
    location = FeatureLocation(start0, end, strand=_strand_int(strand))
    return SeqFeature(location, type=ftype, qualifiers=qualifiers)


def _parse_gff3_attributes(attr_field: str) -> dict:
    attrs: dict = {}
    if not attr_field or attr_field == ".":
        return attrs
    # Single split pass; lstrip each item only when it actually starts with
    # whitespace (cheap leading-byte check, no allocation in the common case).
    for item in attr_field.split(";"):
        if not item:
            continue
        if item[0] == " " or item[0] == "\t":
            item = item.lstrip()
            if not item:
                continue
        eq = item.find("=")
        if eq < 1:
            continue
        key = item[:eq]
        val = item[eq + 1:]
        # Strip a trailing newline / CR / whitespace cheaply (rare but possible).
        if val and val[-1] <= " ":
            val = val.rstrip()
        # Most values are a single token. Avoid the comprehension allocation
        # for the single-token path.
        if "," in val:
            vals = [v for v in (s.strip() for s in val.split(",")) if v]
            if not vals:
                continue
        else:
            vals = [val] if val else []
            if not vals:
                continue
        existing = attrs.get(key)
        if existing is None:
            attrs[key] = vals
        else:
            existing.extend(vals)
    return attrs


def _convert_pair(fasta_file: Path, gff_file: Path, out_file: Path) -> None:
    records = {}
    # Per-record feature-list aliases. Two attribute lookups per row
    # (records[seqid].features.append) is non-trivial across millions of
    # rows; cache the list reference once per seqid.
    feat_lists: dict = {}
    with open(fasta_file) as f:
        for rec in SeqIO.parse(f, "fasta"):
            rec.id = rec.id.split()[0]
            rec.name = rec.id
            rec.description = rec.id
            rec.annotations["molecule_type"] = "DNA"
            records[rec.id] = rec
            feat_lists[rec.id] = rec.features

    skipped_contigs = set()
    malformed = 0
    parse_attrs = _parse_gff3_attributes
    create_feature = _create_feature

    with open(gff_file) as f:
        for line in f:
            if not line:
                continue
            first = line[0]
            if first == "#" or first == "\n":
                continue
            # Single rstrip + bounded split. The bound (8) means split stops
            # after 9 columns so any embedded tabs in column 9 (rare but
            # legal in some GFF dialects) are preserved verbatim.
            if line.endswith("\n"):
                line = line[:-1]
            if not line:
                continue
            parts = line.split("\t", 8)
            if len(parts) < 9:
                malformed += 1
                continue
            seqid = parts[0]
            target = feat_lists.get(seqid)
            if target is None:
                skipped_contigs.add(seqid)
                continue
            try:
                start_i = int(parts[3])
                end_i = int(parts[4])
            except ValueError:
                malformed += 1
                continue
            attrs = parse_attrs(parts[8])
            target.append(
                create_feature(parts[2], start_i, end_i,
                               parts[6], parts[7], attrs)
            )

    if skipped_contigs:
        log_warn(
            f"{gff_file.name}: GFF references contigs not in FASTA, skipped: "
            f"{', '.join(sorted(skipped_contigs))}"
        )
    if malformed:
        log_warn(f"{gff_file.name}: skipped {malformed} malformed GFF lines.")

    with open(out_file, "w") as out_handle:
        SeqIO.write(list(records.values()), out_handle, "genbank")


def batch_convert(input_dir: str, output_dir: str,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)

    pairs = _find_pairs(in_path)
    if not pairs:
        log_warn("No matching FASTA / GFF pairs found (matched by stem).")
        return

    for stem, fasta, gff in pairs:
        out_file = out_path / f"{stem}.gbk"
        log_info(f"Converting {fasta.name} + {gff.name} -> {out_file.name}")
        _convert_pair(fasta, gff, out_file)
