#!/usr/bin/env python3
"""Reference for E03: split a FASTQ into FASTA + QUAL with Biopython's writers.

Usage: fastq_to_fasta_qual.py in.fastq out_prefix
"""
import sys

from Bio import SeqIO

src, prefix = sys.argv[1], sys.argv[2]
with open(prefix + ".fasta", "w") as fa, open(prefix + ".qual", "w") as qu:
    for rec in SeqIO.parse(src, "fastq"):
        SeqIO.write(rec, fa, "fasta")
        SeqIO.write(rec, qu, "qual")
