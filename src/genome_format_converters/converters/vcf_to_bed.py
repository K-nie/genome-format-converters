#!/usr/bin/env python3
"""Convert VCF / BCF to BED intervals.

SNPs emit a 1-bp interval (`pos-1`, `pos`). Indels use the span of the
longest allele. Symbolic alts (`<DEL>`, `<INS>`, `<DUP>`, …) and sites with
no ALT (gVCF non-variant blocks) are skipped with a warning rather than
crashing.
"""

from pathlib import Path
from typing import Optional

try:
    import pysam
except ImportError:
    pysam = None  # Windows path; require_pysam() in batch_convert emits the help

from ._common import (
    iter_input_files,
    log_info,
    log_warn,
    prepare_output_dir,
    require_pysam,
    strip_compound_suffix,
)

_VCF_EXTS = [".vcf", ".vcf.gz", ".bcf"]


def _is_symbolic(alt: str) -> bool:
    return alt.startswith("<") or alt.startswith(".") or "[" in alt or "]" in alt


def _convert_file(in_file: Path, out_file: Path) -> None:
    symbolic = 0
    with pysam.VariantFile(str(in_file)) as vcf, open(out_file, "w") as bed:
        for rec in vcf:
            if not rec.alts:
                # gVCF non-variant block — nothing to emit.
                continue
            ref = rec.ref
            pos = rec.pos
            for alt in rec.alts:
                if _is_symbolic(alt):
                    symbolic += 1
                    continue
                span = max(len(ref), len(alt))
                start = pos - 1
                end = start + span
                name = f"{rec.id or '.'}_{alt}"
                bed.write(f"{rec.chrom}\t{start}\t{end}\t{name}\t.\t.\n")
    if symbolic:
        log_warn(f"{in_file.name}: skipped {symbolic} symbolic / breakend records.")


def batch_convert(input_dir: str, output_dir: str,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    if pysam is None:
        require_pysam()
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)
    for vcf in iter_input_files(in_path, _VCF_EXTS, pattern=pattern):
        stem = strip_compound_suffix(vcf, _VCF_EXTS)
        out_file = out_path / (stem + ".bed")
        log_info(f"Converting {vcf.name} -> {out_file.name}")
        _convert_file(vcf, out_file)
