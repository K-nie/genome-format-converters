#!/usr/bin/env python3
"""Annotate a phylogenetic tree with alignment sequences.

Produces a NEXUS file with `taxa`, `trees` (with a translate block), and
`characters` blocks. Tree leaves that don't appear in the alignment are
left with their original labels and flagged via a warning; the tree itself
is deep-copied before relabelling so callers that use this as a library
don't see their input mutated.

Taxon names that NEXUS would read as punctuation or whitespace (for
example Pfam's `NAME/start-end`) are written in single quotes, and the
character block's datatype is `dna` or `protein` from the alignment's
residues. NEXUS has no symbol for selenocysteine, so in a protein matrix
`U` is written as `X` (unknown residue), with a warning giving the count.
"""

import copy
import re
from io import StringIO
from pathlib import Path

from Bio import AlignIO, Phylo

from ._common import log_warn


_NEXUS_PLAIN = re.compile(r"^[^\s()\[\]{}/\\,;:=*'\"`+<>-]+$")
_NUCLEOTIDES = set("ACGTUNRYSWKMBDHV-?.")


def _nexus_name(name: str) -> str:
    """Quote a NEXUS token when it holds punctuation or whitespace."""
    if _NEXUS_PLAIN.match(name):
        return name
    return "'" + name.replace("'", "''") + "'"


def annotate_tree(tree_file: str, aln_file: str, out_file: str) -> None:
    tree = Phylo.read(tree_file, "newick")
    aln = AlignIO.read(aln_file, "fasta")
    # NEXUS reads '.' as a match character, not a gap.
    seq_dict = {rec.id: str(rec.seq).replace(".", "-") for rec in aln}
    ordered_ids = list(seq_dict.keys())
    id_to_index = {tid: i + 1 for i, tid in enumerate(ordered_ids)}

    tree_copy = copy.deepcopy(tree)
    missing = set()

    # Iterative walk: deep (caterpillar) trees would overflow recursion.
    for clade in tree_copy.get_terminals():
        if clade.name:
            idx = id_to_index.get(clade.name)
            if idx is not None:
                clade.name = str(idx)
            else:
                missing.add(clade.name)

    if missing:
        log_warn(
            "Tree leaves not present in alignment (left unrenamed): "
            + ", ".join(sorted(missing))
        )

    buf = StringIO()
    Phylo.write(tree_copy, buf, "newick", format_branch_length="%.12g",
                format_confidence="%.12g")
    newick = buf.getvalue().strip()
    residues = set("".join(seq_dict.values()).upper())
    datatype = "dna" if residues <= _NUCLEOTIDES else "protein"
    if datatype == "protein":
        n_sec = sum(seq.upper().count("U") for seq in seq_dict.values())
        if n_sec:
            seq_dict = {k: v.replace("U", "X").replace("u", "x") for k, v in seq_dict.items()}
            log_warn(f"{Path(aln_file).name}: wrote {n_sec} selenocysteine (U) "
                     "residue(s) as X, since NEXUS protein data has no U")

    out_path = Path(out_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as out:
        out.write("#NEXUS\nbegin taxa;\n")
        out.write(f"  dimensions ntax={len(seq_dict)};\n")
        out.write("  taxlabels\n")
        for tid in ordered_ids:
            out.write(f"    {_nexus_name(tid)}\n")
        out.write("  ;\nend;\n\n")

        out.write("begin trees;\n  translate\n")
        last = len(ordered_ids)
        for i, tid in enumerate(ordered_ids, 1):
            sep = "," if i < last else ""
            out.write(f"    {i} {_nexus_name(tid)}{sep}\n")
        out.write("  ;\n")
        out.write(f"  tree annotated = [&U] {newick}\n")
        out.write("end;\n\n")

        out.write("begin characters;\n")
        out.write(f"  dimensions nchar={aln.get_alignment_length()};\n")
        out.write(f"  format datatype={datatype} missing=? gap=-;\n")
        out.write("  matrix\n")
        for tid, seq in seq_dict.items():
            out.write(f"    {_nexus_name(tid)}\t{seq}\n")
        out.write("  ;\nend;\n")
