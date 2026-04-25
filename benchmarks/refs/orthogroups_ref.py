#!/usr/bin/env python3
"""Reference handwritten implementation of OrthoFinder Orthogroups.tsv → per-OG FASTA.

This is the "what a careful bioinformatician would write in an hour"
baseline used by benchmarks/tasks/T7_orthogroups.sh to compare against
`gfc orthogroups-to-fasta`. The logic is deliberately straightforward so
the comparison measures the overhead of gfc's uniform-CLI layer (argparse
+ dispatch + stderr logging + loud-failure wrapping), not any
fundamental algorithmic difference.

Usage:
    python orthogroups_ref.py --orthogroups Orthogroups.tsv \\
        --fasta-dir ./proteomes --output-dir ./per_og
"""
from __future__ import annotations

import argparse
from pathlib import Path

from Bio import SeqIO

_FASTA_EXTS = [".fasta", ".fa", ".fna", ".fas", ".faa", ".pep"]


def find_fasta(fasta_dir: Path, species: str) -> Path | None:
    for ext in _FASTA_EXTS:
        p = fasta_dir / f"{species}{ext}"
        if p.exists():
            return p
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--orthogroups", required=True)
    ap.add_argument("--fasta-dir", required=True)
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()

    og_path = Path(args.orthogroups)
    fa_dir = Path(args.fasta_dir)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    with open(og_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        species_cols = header[1:]

        # Pre-load every species' FASTA as a dict {gene_id: SeqRecord}.
        fastas: dict[str, dict[str, object]] = {}
        for sp in species_cols:
            p = find_fasta(fa_dir, sp)
            if p is None:
                fastas[sp] = {}
                continue
            fastas[sp] = SeqIO.to_dict(SeqIO.parse(str(p), "fasta"))

        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0]:
                continue
            og_name = parts[0]
            gene_lists = parts[1:]
            records = []
            for species, blob in zip(species_cols, gene_lists):
                for gene in [g.strip() for g in blob.split(",") if g.strip()]:
                    rec = fastas.get(species, {}).get(gene)
                    if rec is None:
                        continue
                    # Rewrite id/description to match gfc's >species|gene format.
                    rec = rec.__class__(rec.seq, id=f"{species}|{gene}", description="")
                    records.append(rec)
            if not records:
                continue
            SeqIO.write(records, str(out / f"{og_name}.fa"), "fasta")


if __name__ == "__main__":
    main()
