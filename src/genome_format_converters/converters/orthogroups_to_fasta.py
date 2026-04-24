#!/usr/bin/env python3
"""OrthoFinder ``Orthogroups.tsv`` + per-species FASTAs → one FASTA per OG.

The input TSV is the one OrthoFinder emits alongside its results: a header
row with species names followed by rows that list, per species, the
comma-separated gene IDs assigned to each orthogroup.

The per-species FASTAs must match the column names (``<species>.fa``,
``<species>.fasta``, etc.). Output FASTA records are labelled
``>species|gene`` so downstream tools can tell who came from where.
"""

from pathlib import Path
from typing import Dict, List, Optional

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ._common import log_info, log_warn, prepare_output_dir

_FASTA_EXTS = [".fasta", ".fa", ".fna", ".fas", ".faa", ".pep"]


def _load_species_fastas(fasta_dir: Path,
                         species_cols: List[str]) -> Dict[str, Dict[str, SeqRecord]]:
    """For each column name, find the FASTA matching that stem and load it."""
    out: Dict[str, Dict[str, SeqRecord]] = {}
    for species in species_cols:
        found = None
        for ext in _FASTA_EXTS:
            candidate = fasta_dir / f"{species}{ext}"
            if candidate.exists():
                found = candidate
                break
        if found is None:
            log_warn(f"No FASTA found for species {species!r} in {fasta_dir}; skipping.")
            out[species] = {}
            continue
        out[species] = SeqIO.to_dict(SeqIO.parse(str(found), "fasta"))
    return out


def convert(orthogroups_tsv: str, fasta_dir: str, output_dir: str,
            force: bool = False) -> None:
    og_path = Path(orthogroups_tsv)
    fa_path = Path(fasta_dir)
    out_path = prepare_output_dir(output_dir, force=force)

    with open(og_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if not header or header[0] != "Orthogroup":
            log_warn(
                f"{og_path.name}: first column is {header[0]!r}, not 'Orthogroup'. "
                "Proceeding, but double-check the input."
            )
        species_cols = header[1:]

        fastas = _load_species_fastas(fa_path, species_cols)

        og_count = 0
        empty = 0
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0]:
                continue
            og_name = parts[0]
            gene_lists = parts[1:]
            records: List[SeqRecord] = []
            for species, blob in zip(species_cols, gene_lists):
                for gene in [g.strip() for g in blob.split(",") if g.strip()]:
                    rec = fastas.get(species, {}).get(gene)
                    if rec is None:
                        continue
                    label = f"{species}|{gene}"
                    records.append(SeqRecord(rec.seq, id=label, description=""))
            if not records:
                empty += 1
                continue
            out_file = out_path / f"{og_name}.fa"
            SeqIO.write(records, str(out_file), "fasta")
            og_count += 1
        log_info(f"{og_path.name}: wrote {og_count} orthogroup FASTA(s) "
                 f"({empty} empty / unresolved groups skipped).")
