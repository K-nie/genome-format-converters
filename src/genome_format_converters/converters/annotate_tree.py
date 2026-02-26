#!/usr/bin/env python3
"""
Annotate a phylogenetic tree with sequences from an alignment.
Outputs a NEXUS file with a translate block and the alignment.
"""

from pathlib import Path
from Bio import Phylo, AlignIO, SeqIO

def annotate_tree(tree_file: str, aln_file: str, out_file: str) -> None:
    """Read a Newick tree and a FASTA alignment, produce a NEXUS file."""
    tree = Phylo.read(tree_file, "newick")
    aln = AlignIO.read(aln_file, "fasta")
    seq_dict = {rec.id: str(rec.seq) for rec in aln}

    with open(out_file, 'w') as out:
        out.write("#NEXUS\nbegin taxa;\n")
        out.write(f"  dimensions ntax={len(seq_dict)};\n")
        out.write("  taxlabels\n")
        for tid in seq_dict:
            out.write(f"    {tid}\n")
        out.write("  ;\nend;\n\n")
        out.write("begin trees;\n")
        out.write("  translate\n")
        for i, tid in enumerate(seq_dict, 1):
            out.write(f"    {i} {tid},\n")
        out.write("  ;\n")
        # Rename tree leaves to numbers
        def rename_leaves(clade):
            if clade.name:
                idx = list(seq_dict.keys()).index(clade.name) + 1
                clade.name = str(idx)
            for child in clade.clades:
                rename_leaves(child)
        tree_copy = tree
        rename_leaves(tree_copy.root)
        newick = tree_copy.format("newick")
        out.write(f"  tree annotated = [&U] {newick}\n")
        out.write("end;\n\n")
        out.write("begin characters;\n")
        out.write(f"  dimensions nchar={aln.get_alignment_length()};\n")
        out.write("  format datatype=dna missing=? gap=-;\n")
        out.write("  matrix\n")
        for tid, seq in seq_dict.items():
            out.write(f"    {tid}\t{seq}\n")
        out.write("  ;\nend;\n")