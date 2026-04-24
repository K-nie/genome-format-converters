#!/usr/bin/env python3
"""Convert VCF / BCF to a tab-separated table.

INFO columns come from the VCF header (not just the first record), so
fields that appear only in later records don't get silently dropped.
"""

import csv
from pathlib import Path
from typing import Optional

try:
    import pysam
except ImportError:
    pysam = None

from ._common import (
    iter_input_files,
    log_info,
    prepare_output_dir,
    require_pysam,
    strip_compound_suffix,
)

_VCF_EXTS = [".vcf", ".vcf.gz", ".bcf"]


def _format_info_value(val) -> str:
    if val is None:
        return "."
    if isinstance(val, (tuple, list)):
        return ",".join(str(v) for v in val)
    if isinstance(val, bool):
        return "1" if val else "0"
    return str(val)


def _convert_file(in_file: Path, out_file: Path) -> None:
    with pysam.VariantFile(str(in_file)) as vcf, open(out_file, "w", newline="") as tsv:
        info_keys = list(vcf.header.info.keys())
        writer = csv.writer(tsv, delimiter="\t")
        writer.writerow(
            ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"] + info_keys
        )
        for rec in vcf:
            alt_str = ",".join(rec.alts) if rec.alts else "."
            qual = "." if rec.qual is None else rec.qual
            filt = ",".join(rec.filter.keys()) if rec.filter.keys() else "."
            row = [rec.chrom, rec.pos, rec.id or ".", rec.ref, alt_str, qual, filt]
            for key in info_keys:
                row.append(_format_info_value(rec.info.get(key, None)))
            writer.writerow(row)


def batch_convert(input_dir: str, output_dir: str,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    if pysam is None:
        require_pysam()
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)
    for vcf in iter_input_files(in_path, _VCF_EXTS, pattern=pattern):
        stem = strip_compound_suffix(vcf, _VCF_EXTS)
        out_file = out_path / (stem + ".tsv")
        log_info(f"Converting {vcf.name} -> {out_file.name}")
        _convert_file(vcf, out_file)
