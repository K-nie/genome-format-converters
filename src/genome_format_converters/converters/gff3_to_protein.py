#!/usr/bin/env python3
"""Extract protein sequences from paired GFF3 + FASTA input.

For each mRNA (or transcript-like feature) we collect its CDS children, sort
them in biological order (ascending by start on the + strand, descending by
end on the - strand), concatenate the nucleotide spans, trim the phase of
the leading CDS, and translate once. Single-exon / yeast-style genes fall
out of this as a degenerate case (one CDS per mRNA).
"""

from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from BCBio import GFF

from ._common import log_info, log_warn, prepare_output_dir

_FASTA_EXTS = [".fasta", ".fa", ".fna", ".fas"]
_GFF_EXTS = [".gff3", ".gff"]
_SUFFIX_STRIPS = [".final", "_final", "final"]
_TRANSCRIPT_TYPES = {"mRNA", "transcript"}


def _strip_suffixes(stem: str, suffixes: List[str]) -> str:
    for suf in suffixes:
        if stem.endswith(suf):
            return stem[: -len(suf)]
    return stem


def _find_pairs(input_dir: Path) -> List[Tuple[str, Path, Path]]:
    fastas: Dict[str, Path] = {}
    for ext in _FASTA_EXTS:
        for f in input_dir.glob(f"*{ext}"):
            fastas[f.stem] = f

    gffs: Dict[str, Path] = {}
    for ext in _GFF_EXTS:
        for f in input_dir.glob(f"*{ext}"):
            gffs[f.stem] = f

    stripped_map: Dict[str, Path] = {}
    for stem, f in gffs.items():
        stripped_map[_strip_suffixes(stem, _SUFFIX_STRIPS)] = f

    pairs: List[Tuple[str, Path, Path]] = []
    for stem, fasta in fastas.items():
        if stem in stripped_map:
            pairs.append((stem, fasta, stripped_map[stem]))
        elif stem in gffs:
            pairs.append((stem, fasta, gffs[stem]))
    return pairs


def _collect_cds(feature) -> List:
    """Flatten a feature's sub-tree and return only CDS children."""
    out = []
    stack = [feature]
    while stack:
        f = stack.pop()
        if f.type == "CDS":
            out.append(f)
        stack.extend(f.sub_features)
    return out


def _translate_transcript(rec, transcript_feature, out_handle) -> None:
    """Concatenate CDS children of one transcript, translate, write FASTA."""
    cds_list = _collect_cds(transcript_feature)
    if not cds_list:
        return

    strand = transcript_feature.location.strand
    # Sort + strand by start ascending, - strand by start descending. A CDS
    # on the '.' strand is unusual for a coding feature; default to + ordering.
    cds_list.sort(key=lambda f: int(f.location.start), reverse=(strand == -1))

    pieces = []
    for cds in cds_list:
        piece = cds.extract(rec.seq)
        if piece:
            pieces.append(str(piece))
    if not pieces:
        log_warn(
            f"Empty CDS concatenation for "
            f"{transcript_feature.qualifiers.get('ID', ['?'])[0]} on {rec.id}; "
            "skipping."
        )
        return

    joined = Seq("".join(pieces))
    # Phase of the first CDS (already in 5'→3' order after the sort above).
    leading_phase_raw = cds_list[0].qualifiers.get("phase", ["0"])[0]
    try:
        leading_phase = int(leading_phase_raw)
    except (TypeError, ValueError):
        leading_phase = 0
    if leading_phase in (1, 2):
        joined = joined[leading_phase:]

    # Trim trailing partial codon so Biopython doesn't warn.
    overflow = len(joined) % 3
    if overflow:
        joined = joined[:-overflow]

    protein = joined.translate(to_stop=False)
    fid = transcript_feature.qualifiers.get("ID", ["CDS"])[0]
    out_handle.write(f">{fid}\n{str(protein)}\n")


def _walk(rec, feature, out_handle) -> None:
    if feature.type in _TRANSCRIPT_TYPES:
        _translate_transcript(rec, feature, out_handle)
        return
    # Some annotations attach CDS directly under `gene`; translate that
    # as if the gene were a single-transcript feature.
    if feature.type == "gene" and not any(
        sub.type in _TRANSCRIPT_TYPES for sub in feature.sub_features
    ):
        _translate_transcript(rec, feature, out_handle)
        return
    for sub in feature.sub_features:
        _walk(rec, sub, out_handle)


def _convert_pair(fasta_file: Path, gff_file: Path, out_handle) -> None:
    genome = SeqIO.to_dict(SeqIO.parse(str(fasta_file), "fasta"))
    with open(gff_file) as in_handle:
        for rec in GFF.parse(in_handle, base_dict=genome):
            for feature in rec.features:
                _walk(rec, feature, out_handle)


def batch_convert(input_dir: str, output_dir: str,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)

    pairs = _find_pairs(in_path)
    if not pairs:
        log_warn("No matching GFF3 / FASTA pairs found (matched by stem).")
        return

    for stem, fasta, gff in pairs:
        out_file = out_path / f"{stem}.faa"
        log_info(f"Processing {stem}: {fasta.name} + {gff.name} -> {out_file.name}")
        with open(out_file, "w") as out_handle:
            _convert_pair(fasta, gff, out_handle)
