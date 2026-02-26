#!/usr/bin/env python3
"""
Extract protein sequences from GFF3 + FASTA.
Files are paired by stem (common prefix) with optional suffix stripping.
"""

import sys
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from BCBio import GFF

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
    fastas = {}
    for ext in FASTA_EXTS:
        for f in input_dir.glob(f"*{ext}"):
            fastas[f.stem] = f
    gffs = {}
    for ext in GFF_EXTS:
        for f in input_dir.glob(f"*{ext}"):
            gffs[f.stem] = f

    stripped_map = {}
    for stem, f in gffs.items():
        stripped = _strip_suffixes(stem, SUFFIX_STRIPS)
        stripped_map[stripped] = f

    pairs = []
    for stem, fasta in fastas.items():
        if stem in stripped_map:
            pairs.append((stem, fasta, stripped_map[stem]))
        elif stem in gffs:
            pairs.append((stem, fasta, gffs[stem]))
    return pairs

def _process_cds(rec, feature, out_handle):
    """Extract and translate a single CDS feature."""
    if feature.type != "CDS":
        return
    # Extract the nucleotide sequence
    seq = feature.extract(rec.seq)
    if not seq:
        # This can happen if coordinates lie outside the available sequence
        sys.stderr.write(
            f"Warning: Empty sequence for CDS {feature.qualifiers.get('ID', ['?'])[0]} "
            f"at {rec.id}:{feature.location.start+1}-{feature.location.end}. "
            f"Skipping.\n"
        )
        return
    # Apply phase (0,1,2) – trim from start
    phase = int(feature.qualifiers.get('phase', [0])[0])
    if phase:
        seq = seq[phase:]
    # Translate to protein (do not stop at internal stops)
    prot = seq.translate(to_stop=False)
    # Write to output
    out_handle.write(f">{feature.qualifiers.get('ID', ['CDS'])[0]}\n")
    out_handle.write(f"{str(prot)}\n")

def _process_features(rec, feature, out_handle):
    """Recursively search for CDS features."""
    _process_cds(rec, feature, out_handle)
    for sub in feature.sub_features:
        _process_features(rec, sub, out_handle)

def _convert_pair(fasta_file: Path, gff_file: Path, out_handle):
    """Extract and translate all CDS from one pair, writing to out_handle."""
    # Load genome
    genome = SeqIO.to_dict(SeqIO.parse(str(fasta_file), "fasta"))
    # Parse GFF with the genome dictionary so that rec.seq is available
    with open(gff_file) as in_handle:
        for rec in GFF.parse(in_handle, base_dict=genome):
            for feature in rec.features:
                _process_features(rec, feature, out_handle)

def batch_convert(input_dir: str, output_dir: str) -> None:
    """For each GFF3/FASTA pair, produce a protein FASTA file."""
    in_path = Path(input_dir)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    pairs = _find_pairs(in_path)
    if not pairs:
        print("No matching GFF3/FASTA pairs found.")
        return

    for stem, fasta, gff in pairs:
        out_file = out_path / f"{stem}.faa"
        print(f"Processing {stem}: {fasta.name} + {gff.name} -> {out_file.name}")
        with open(out_file, 'w') as out_handle:
            _convert_pair(fasta, gff, out_handle)