#!/usr/bin/env python3
"""Reference for E04: join a FASTA and a QUAL file into FASTQ with Biopython.

Usage: fasta_qual_to_fastq.py in.fasta in.qual out.fastq
"""
import sys

from Bio import SeqIO
from Bio.SeqIO.QualityIO import PairedFastaQualIterator

fa, qu, out = sys.argv[1:4]
with open(fa) as f, open(qu) as q:
    SeqIO.write(PairedFastaQualIterator(f, q), out, "fastq")
