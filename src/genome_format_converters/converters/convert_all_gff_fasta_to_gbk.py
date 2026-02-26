#!/usr/bin/env python3
"""
Convert paired FASTA and GFF3 files to GenBank format.
"""

from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

FASTA_EXTS = [".fasta", ".fa", ".fna", ".fas"]
GFF_EXTS = [".gff3", ".gff"]
SUFFIX_STRIPS = [".final", "_final", "final"]

def _strip_suffixes(stem: str, suffixes: list) -> str:
    for suf in suffixes:
        if stem.endswith(suf):
            stem = stem[: -len(suf)]
            break
    return stem

def _find_pairs(input_dir: Path):
    fastas = [f for ext in FASTA_EXTS for f in input_dir.glob(f"*{ext}")]
    gffs   = [f for ext in GFF_EXTS for f in input_dir.glob(f"*{ext}")]

    fasta_by_stem = defaultdict(list)
    for f in fastas:
        fasta_by_stem[f.stem].append(f)

    gff_by_stem = defaultdict(list)
    for g in gffs:
        stripped = _strip_suffixes(g.stem, SUFFIX_STRIPS)
        gff_by_stem[stripped].append(g)

    pairs = []
    for stripped, gff_list in gff_by_stem.items():
        if stripped in fasta_by_stem:
            fa_list = fasta_by_stem[stripped]
            fasta = fa_list[0]   # take first if multiple
            gff   = gff_list[0]
            pairs.append((stripped, fasta, gff))
    return pairs

def _gff3_to_genbank_qualifiers(ftype: str, attrs: dict) -> dict:
    """Convert GFF3 attributes to GenBank qualifiers."""
    quals = {}
    # locus_tag: from Parent if present, else from ID
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
    # Add note for other attributes (simplified)
    note_parts = []
    for k, v in attrs.items():
        if k not in ("ID", "Parent", "Name", "product"):
            note_parts.append(f"{k}:{','.join(v)}")
    if note_parts:
        quals["note"] = note_parts
    return quals

def _create_feature(ftype: str, start: int, end: int, strand: str, phase: str, attrs: dict) -> SeqFeature:
    """Create a Biopython SeqFeature from GFF3 fields."""
    start0 = max(0, start - 1)
    end0 = end
    strand_int = 1 if strand == "+" else (-1 if strand == "-" else 0)
    qualifiers = _gff3_to_genbank_qualifiers(ftype, attrs)
    if ftype == "CDS" and phase in ("0", "1", "2"):
        qualifiers["codon_start"] = [str(int(phase) + 1)]
    location = FeatureLocation(start0, end0, strand=strand_int)
    return SeqFeature(location, type=ftype, qualifiers=qualifiers)

def _parse_gff3_attributes(attr_field: str) -> dict:
    """Parse the 9th column of GFF3 into a dict."""
    attrs = {}
    if not attr_field or attr_field == ".":
        return attrs
    for item in attr_field.split(";"):
        item = item.strip()
        if not item or "=" not in item:
            continue
        key, val = item.split("=", 1)
        key = key.strip()
        vals = [v.strip() for v in val.split(",") if v.strip()]
        attrs.setdefault(key, []).extend(vals)
    return attrs

def _convert_pair(fasta_file: Path, gff_file: Path, out_file: Path) -> None:
    """Convert one FASTA+GFF pair to GenBank."""
    # Load FASTA
    records = {}
    with open(fasta_file) as f:
        for rec in SeqIO.parse(f, "fasta"):
            rec.id = rec.id.split()[0]   # take first token
            rec.name = rec.id
            rec.description = rec.id
            rec.annotations["molecule_type"] = "DNA"
            records[rec.id] = rec

    # Parse GFF and add features
    with open(gff_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs_raw = parts[:9]
            if seqid not in records:
                continue
            try:
                start_i = int(start)
                end_i   = int(end)
            except ValueError:
                continue
            attrs = _parse_gff3_attributes(attrs_raw)
            feature = _create_feature(ftype, start_i, end_i, strand, phase, attrs)
            records[seqid].features.append(feature)

    # Write GenBank
    with open(out_file, 'w') as out_handle:
        SeqIO.write(list(records.values()), out_handle, "genbank")

def batch_convert(input_dir: str, output_dir: str) -> None:
    """For each FASTA/GFF pair in input_dir, produce a GenBank file."""
    in_path = Path(input_dir)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    pairs = _find_pairs(in_path)
    if not pairs:
        print("No matching FASTA/GFF pairs found.")
        return

    for stem, fasta, gff in pairs:
        out_file = out_path / f"{stem}.gbk"
        print(f"Converting {fasta.name} + {gff.name} -> {out_file.name}")
        _convert_pair(fasta, gff, out_file)