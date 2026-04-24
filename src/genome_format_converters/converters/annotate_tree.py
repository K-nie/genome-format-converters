#!/usr/bin/env python3
"""Annotate a phylogenetic tree with alignment sequences.

Produces a NEXUS file with `taxa`, `trees` (with a translate block), and
`characters` blocks. Tree leaves that don't appear in the alignment are
left with their original labels and flagged via a warning; the tree itself
is deep-copied before relabelling so callers that use this as a library
don't see their input mutated.
"""

import copy
from pathlib import Path

from Bio import AlignIO, Phylo

from ._common import log_warn


def annotate_tree(tree_file: str, aln_file: str, out_file: str) -> None:
    tree = Phylo.read(tree_file, "newick")
    aln = AlignIO.read(aln_file, "fasta")
    seq_dict = {rec.id: str(rec.seq) for rec in aln}
    ordered_ids = list(seq_dict.keys())
    id_to_index = {tid: i + 1 for i, tid in enumerate(ordered_ids)}

    tree_copy = copy.deepcopy(tree)
    missing = set()

    def rename_leaves(clade):
        if clade.name:
            idx = id_to_index.get(clade.name)
            if idx is not None:
                clade.name = str(idx)
            else:
                missing.add(clade.name)
        for child in clade.clades:
            rename_leaves(child)

    rename_leaves(tree_copy.root)

    if missing:
        log_warn(
            "Tree leaves not present in alignment (left unrenamed): "
            + ", ".join(sorted(missing))
        )

    newick = tree_copy.format("newick").strip()

    out_path = Path(out_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as out:
        out.write("#NEXUS\nbegin taxa;\n")
        out.write(f"  dimensions ntax={len(seq_dict)};\n")
        out.write("  taxlabels\n")
        for tid in ordered_ids:
            out.write(f"    {tid}\n")
        out.write("  ;\nend;\n\n")

        out.write("begin trees;\n  translate\n")
        last = len(ordered_ids)
        for i, tid in enumerate(ordered_ids, 1):
            sep = "," if i < last else ""
            out.write(f"    {i} {tid}{sep}\n")
        out.write("  ;\n")
        out.write(f"  tree annotated = [&U] {newick}\n")
        out.write("end;\n\n")

        out.write("begin characters;\n")
        out.write(f"  dimensions nchar={aln.get_alignment_length()};\n")
        out.write("  format datatype=dna missing=? gap=-;\n")
        out.write("  matrix\n")
        for tid, seq in seq_dict.items():
            out.write(f"    {tid}\t{seq}\n")
        out.write("  ;\nend;\n")
