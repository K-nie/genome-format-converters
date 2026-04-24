#!/usr/bin/env python3
"""Convert paired FASTA + GFF3 files to GenBank format.

Files are paired by stem, with common yeast-suffix strippers (`.final`,
`_final`, `final`). Strand is mapped to `+1` / `-1` / `0`, matching BioPython's
SeqFeature convention — `0` / `.` strands are preserved rather than silently
forced to `+`.
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
    if strand == "+":
        return 1
    if strand == "-":
        return -1
    return 0


def _gff3_to_gb_qualifiers(ftype: str, attrs: dict) -> dict:
    quals = {}
    if ftype != "gene" and "Parent" in attrs:
        quals["locus_tag"] = attrs["Parent"]
    elif "ID" in attrs:
        quals["locus_tag"] = attrs["ID"]
    if ftype == "gene" and "Name" in attrs:
        quals["gene"] = attrs["Name"]
    if ftype == "mRNA" and "Name" in attrs:
        quals["transcript_id"] = attrs["Name"]
    if ftype == "CDS" and "Name" in attrs:
        quals["protein_id"] = attrs["Name"]
        if "product" not in attrs:
            quals["product"] = attrs["Name"]
    if ftype == "CDS" and "product" in attrs:
        quals["product"] = attrs["product"]
    note_parts = []
    for k, v in attrs.items():
        if k not in ("ID", "Parent", "Name", "product"):
            note_parts.append(f"{k}:{','.join(v)}")
    if note_parts:
        quals["note"] = note_parts
    return quals


def _create_feature(ftype: str, start: int, end: int, strand: str, phase: str, attrs: dict) -> SeqFeature:
    start0 = max(0, start - 1)
    end0 = end
    strand_int = _strand_int(strand)
    qualifiers = _gff3_to_gb_qualifiers(ftype, attrs)
    if ftype == "CDS" and phase in ("0", "1", "2"):
        qualifiers["codon_start"] = [str(int(phase) + 1)]
    location = FeatureLocation(start0, end0, strand=strand_int)
    return SeqFeature(location, type=ftype, qualifiers=qualifiers)


def _parse_gff3_attributes(attr_field: str) -> dict:
    attrs = {}
    if not attr_field or attr_field == ".":
        return attrs
    for item in attr_field.split(";"):
        item = item.strip()
        if not item or "=" not in item:
            continue
        key, val = item.split("=", 1)
        vals = [v.strip() for v in val.split(",") if v.strip()]
        attrs.setdefault(key.strip(), []).extend(vals)
    return attrs


def _convert_pair(fasta_file: Path, gff_file: Path, out_file: Path) -> None:
    records = {}
    with open(fasta_file) as f:
        for rec in SeqIO.parse(f, "fasta"):
            rec.id = rec.id.split()[0]
            rec.name = rec.id
            rec.description = rec.id
            rec.annotations["molecule_type"] = "DNA"
            records[rec.id] = rec

    skipped_contigs = set()
    malformed = 0
    with open(gff_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                malformed += 1
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs_raw = parts[:9]
            if seqid not in records:
                skipped_contigs.add(seqid)
                continue
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                malformed += 1
                continue
            attrs = _parse_gff3_attributes(attrs_raw)
            feature = _create_feature(ftype, start_i, end_i, strand, phase, attrs)
            records[seqid].features.append(feature)

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
