#!/usr/bin/env python3
"""
Convert VCF to BED intervals.
For SNPs: 1‑bp interval (pos-1, pos). For indels: covers the entire variant.
"""

from pathlib import Path
import pysam

def _convert_file(in_file: Path, out_file: Path) -> None:
    """Convert a single VCF/BCF file to BED."""
    with pysam.VariantFile(in_file) as vcf, open(out_file, 'w') as bed:
        for rec in vcf:
            chrom = rec.chrom
            pos = rec.pos
            ref = rec.ref
            for alt in rec.alts:
                name = f"{rec.id or '.'}_{alt}"
                end = pos + max(len(ref), len(alt)) - 1
                bed.write(f"{chrom}\t{pos-1}\t{end}\t{name}\t.\t.\n")

def batch_convert(input_dir: str, output_dir: str) -> None:
    """Convert all VCF/BCF files in input_dir to BED in output_dir."""
    in_path = Path(input_dir)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    exts = [".vcf", ".vcf.gz", ".bcf"]
    for ext in exts:
        for vcf in in_path.glob(f"*{ext}"):
            out_file = out_path / vcf.with_suffix(".bed").name
            print(f"Converting {vcf.name} -> {out_file.name}")
            _convert_file(vcf, out_file)