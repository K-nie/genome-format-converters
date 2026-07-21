#!/usr/bin/env python3
"""OrthoFinder ``Orthogroups.tsv`` + per-species FASTAs -> one FASTA per OG.

Author: Benjamin Narh-Madey

The input TSV is the one OrthoFinder emits alongside its results: a header
row with species names followed by rows that list, per species, the
comma-separated gene IDs assigned to each orthogroup.

The per-species FASTAs must match the column names (``<species>.fa``,
``<species>.fasta``, etc.). Output FASTA records are labelled
``>species|gene`` so downstream tools can tell who came from where.

Implementation note (perf): for typical OrthoFinder runs we read each
species' FASTA once (read-only) and then emit thousands of small per-OG
FASTAs. Going through SeqIO.to_dict / SeqIO.write builds a SeqRecord tree
per gene and re-runs biopython's per-write encoder on every OG, which
dominates runtime. We instead:

    * parse each species FASTA with a streaming text parser into a
      ``dict[gene_id, str]`` of joined sequence bytes;
    * emit each per-OG FASTA in 60-char-wrapped form using a single
      buffered write per OG.

Output formatting (60 cols, no description, single trailing newline)
matches what ``SeqIO.write(records, ..., "fasta")`` produces, so existing
downstream pipelines are unaffected.
"""

from pathlib import Path
from typing import Dict, Optional

from ._common import log_info, log_warn, prepare_output_dir

_FASTA_EXTS = [".fasta", ".fa", ".fna", ".fas", ".faa", ".pep"]
# Match biopython's default FASTA wrap width.
_FASTA_WRAP = 60


def _parse_fasta_to_seqs(path: Path) -> Dict[str, str]:
    """Stream a FASTA into ``{id: sequence}`` with no SeqRecord wrapping.

    The id is the first whitespace-delimited token of the description line,
    matching what ``SeqIO.to_dict`` keys on by default.
    """
    seqs: Dict[str, str] = {}
    cur_id: Optional[str] = None
    cur_parts: list = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line:
                continue
            if line[0] == ">":
                if cur_id is not None:
                    seqs[cur_id] = "".join(cur_parts)
                # Strip leading '>' and trailing newline; take first token.
                desc = line[1:]
                if desc.endswith("\n"):
                    desc = desc[:-1]
                # Equivalent to ``desc.split(None, 1)[0]`` but avoids the
                # split allocation when there's no whitespace.
                sp = desc.find(" ")
                tab = desc.find("\t")
                if sp < 0:
                    cut = tab
                elif tab < 0:
                    cut = sp
                else:
                    cut = sp if sp < tab else tab
                cur_id = desc[:cut] if cut >= 0 else desc
                cur_parts = []
            else:
                # Strip just the trailing newline; keep any internal
                # whitespace out of the joined sequence.
                if line.endswith("\n"):
                    line = line[:-1]
                if line:
                    cur_parts.append(line)
    if cur_id is not None:
        seqs[cur_id] = "".join(cur_parts)
    return seqs


def _wrap_seq(seq: str, width: int = _FASTA_WRAP) -> str:
    """Wrap a sequence at ``width`` characters with a trailing newline.

    Reuses Python's slice machinery in a tight loop. For genome-scale
    sequences this is materially faster than the biopython per-record
    formatter.
    """
    n = len(seq)
    if n <= width:
        return seq + "\n"
    parts = []
    i = 0
    while i < n:
        parts.append(seq[i:i + width])
        i += width
    parts.append("")  # trailing newline after final block
    return "\n".join(parts)


def _load_species_fastas(fasta_dir: Path,
                         species_cols) -> Dict[str, Dict[str, str]]:
    """For each column name, find the FASTA matching that stem and load it."""
    out: Dict[str, Dict[str, str]] = {}
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
        out[species] = _parse_fasta_to_seqs(found)
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
        n_species = len(species_cols)

        fastas = _load_species_fastas(fa_path, species_cols)

        og_count = 0
        empty = 0
        wrap = _wrap_seq
        for line in fh:
            if not line:
                continue
            if line.endswith("\n"):
                line = line[:-1]
            if not line:
                continue
            parts = line.split("\t")
            if not parts[0]:
                continue
            og_name = parts[0]
            # Build the per-OG content in a single string buffer, then write
            # once. One file open + one write per OG instead of per record.
            chunks = []
            for i in range(n_species):
                if i + 1 >= len(parts):
                    break
                blob = parts[i + 1]
                if not blob:
                    continue
                species = species_cols[i]
                sp_seqs = fastas.get(species)
                if not sp_seqs:
                    continue
                # Inline the gene split: avoid the comprehension allocation
                # when there's only one gene (the common case for 1:1 OGs).
                if "," in blob:
                    genes = [g.strip() for g in blob.split(",") if g.strip()]
                else:
                    g = blob.strip()
                    genes = [g] if g else []
                for gene in genes:
                    seq = sp_seqs.get(gene)
                    if seq is None:
                        continue
                    chunks.append(">")
                    chunks.append(species)
                    chunks.append("|")
                    chunks.append(gene)
                    chunks.append("\n")
                    chunks.append(wrap(seq))
            if not chunks:
                empty += 1
                continue
            out_file = out_path / f"{og_name}.fa"
            with open(out_file, "w", encoding="utf-8") as out_fh:
                out_fh.write("".join(chunks))
            og_count += 1
        log_info(f"{og_path.name}: wrote {og_count} orthogroup FASTA(s) "
                 f"({empty} empty / unresolved groups skipped).")
