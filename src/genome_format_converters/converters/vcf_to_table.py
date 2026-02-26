#!/usr/bin/env python3
"""
Convert VCF to tab‑separated table (TSV) with all INFO fields as columns.
"""

import csv
from pathlib import Path
import pysam

def _convert_file(in_file: Path, out_file: Path) -> None:
    """Convert a single VCF/BCF file to TSV."""
    with pysam.VariantFile(in_file) as vcf, open(out_file, 'w', newline='') as tsv:
        # Get all INFO keys from the first record
        first_rec = next(vcf)
        info_keys = list(first_rec.info.keys())
        vcf.reset()
        writer = csv.writer(tsv, delimiter='\t')
        writer.writerow(["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"] + info_keys)

        for rec in vcf:
            alt_str = ','.join(rec.alts) if rec.alts else '.'
            row = [rec.chrom, rec.pos, rec.id or '.', rec.ref, alt_str, rec.qual or '.', ','.join(rec.filter.keys())]
            for key in info_keys:
                val = rec.info.get(key, '.')
                if isinstance(val, tuple):
                    val = ','.join(str(v) for v in val)
                row.append(val)
            writer.writerow(row)

def batch_convert(input_dir: str, output_dir: str) -> None:
    """Convert all VCF/BCF files in input_dir to TSV in output_dir."""
    in_path = Path(input_dir)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    exts = [".vcf", ".vcf.gz", ".bcf"]
    for ext in exts:
        for vcf in in_path.glob(f"*{ext}"):
            out_file = out_path / vcf.with_suffix(".tsv").name
            print(f"Converting {vcf.name} -> {out_file.name}")
            _convert_file(vcf, out_file)