#!/usr/bin/env python3
"""
Create consensus FASTA from reference FASTA + VCF.
Pairs FASTA and VCF files by stem.
"""

from pathlib import Path
from collections import defaultdict
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def _find_pairs(input_dir: Path):
    fastas = {f.stem: f for f in input_dir.glob("*.fasta")}
    vcfs   = {f.stem: f for f in input_dir.glob("*.vcf")}
    vcfs.update({f.stem: f for f in input_dir.glob("*.vcf.gz")})
    common = set(fastas) & set(vcfs)
    return [(stem, fastas[stem], vcfs[stem]) for stem in common]

def _create_consensus(fasta_file: Path, vcf_file: Path, out_file: Path) -> None:
    """Generate consensus sequences for all samples in the VCF."""
    # Load reference
    ref_rec = next(SeqIO.parse(fasta_file, "fasta"))
    ref_seq = list(str(ref_rec.seq).upper())
    # Load VCF
    vcf = pysam.VariantFile(vcf_file)
    samples = list(vcf.header.samples)
    sample_seqs = {sample: ref_seq[:] for sample in samples}

    for rec in vcf:
        pos = rec.pos - 1  # 0‑based
        for sample in samples:
            gt = rec.samples[sample].get('GT', (None, None))
            # Simple: if first allele is 1 (ALT), substitute; else keep reference
            if gt[0] == 1:
                alt = rec.alts[0]
                sample_seqs[sample][pos] = alt

    with open(out_file, 'w') as out_handle:
        for sample in samples:
            seq_str = ''.join(sample_seqs[sample])
            rec = SeqRecord(Seq(seq_str), id=sample, description="consensus")
            SeqIO.write(rec, out_handle, "fasta")

def batch_convert(input_dir: str, output_dir: str) -> None:
    """For each FASTA/VCF pair, produce a consensus FASTA file."""
    in_path = Path(input_dir)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    pairs = _find_pairs(in_path)
    if not pairs:
        print("No matching FASTA/VCF pairs found.")
        return

    for stem, fasta, vcf in pairs:
        out_file = out_path / f"{stem}_consensus.fa"
        print(f"Processing {stem}: {fasta.name} + {vcf.name} -> {out_file.name}")
        _create_consensus(fasta, vcf, out_file)