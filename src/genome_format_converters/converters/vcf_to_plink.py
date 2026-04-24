#!/usr/bin/env python3
"""Convert VCF / BCF to PLINK binary (``.bed`` + ``.bim`` + ``.fam``).

Shares the filter / polarisation pipeline with the EIGENSTRAT converters.
The binary ``.bed`` format encoding (per the PLINK 1.9 spec):

    magic bytes (3):  0x6C 0x1B 0x01        (little-endian, SNP-major)
    then, per SNP:    ceil(N/4) bytes of packed 2-bit genotypes
                      where the 2-bit codes are:
                        00 = hom A1 (minor-allele homozygote)
                        01 = missing
                        10 = heterozygous
                        11 = hom A2 (major-allele homozygote)
"""

from pathlib import Path
from typing import Optional, Tuple

try:
    import pysam
except ImportError:
    pysam = None

from ._common import (
    iter_input_files,
    log_info,
    log_warn,
    prepare_output_dir,
    require_pysam,
    strip_compound_suffix,
    load_two_col_map,
)
from ._eigenstrat_common import (
    BASES,
    AncestralProvider,
    is_transition,
    load_chrom_map,
    load_genetic_map,
    morgans_for,
    passes_maf,
    passes_missing,
    polarise,
    resolve_chrom,
    validate_map_coverage,
)

_VCF_EXTS = [".vcf", ".vcf.gz", ".bcf"]
_PLINK_MAGIC = bytes([0x6C, 0x1B, 0x01])


def _eig_code(alleles, flipped: bool) -> str:
    """Same logic as vcf_to_eigenstrat._genotype_code: counts copies of the
    major / reference allele, with flipping when polarisation swapped
    ref/alt."""
    if alleles is None or any(a is None for a in alleles):
        return "9"
    major_idx = 1 if flipped else 0
    major_count = sum(1 for a in alleles if a == major_idx)
    if major_count == 2:
        return "2"
    if major_count == 1:
        return "1"
    return "0"


_EIG_TO_PLINK_BITS = {"0": 0b00, "9": 0b01, "1": 0b10, "2": 0b11}


def _pack_snp_row(row: str) -> bytes:
    """Pack an EIGENSTRAT-encoded genotype row into PLINK .bed bytes.

    Sample 0 lands in the least-significant 2 bits of byte 0, sample 1 in
    bits 2-3, and so on. Trailing samples in the last byte are padded with
    `01` (missing).
    """
    n = len(row)
    n_bytes = (n + 3) // 4
    buf = bytearray(n_bytes)
    for i, code_char in enumerate(row):
        bit_pair = _EIG_TO_PLINK_BITS[code_char]
        byte_idx = i // 4
        bit_offset = (i % 4) * 2
        buf[byte_idx] |= bit_pair << bit_offset
    # Pad remaining sample slots in the last byte with 01 (missing).
    remainder = n % 4
    if remainder:
        last = n_bytes - 1
        for slot in range(remainder, 4):
            buf[last] |= 0b01 << (slot * 2)
    return bytes(buf)


def _write_fam(fam_path: Path, samples, pop_map, sex_map) -> None:
    sex_code = {"M": "1", "F": "2"}
    with open(fam_path, "w") as fh:
        for s in samples:
            fam = pop_map.get(s, s)
            sex = sex_code.get(sex_map.get(s, "U"), "0")
            fh.write(f"{fam} {s} 0 0 {sex} -9\n")


def _convert_file(in_file: Path, out_prefix: Path,
                  pop_map, sex_map, chrom_map,
                  default_chrom: Optional[int],
                  gmap,
                  ancestral: AncestralProvider,
                  info_aa: bool,
                  transversions_only: bool,
                  min_maf: float,
                  max_missing: float,
                  strict_maps: bool) -> None:
    bed_path = out_prefix.with_suffix(".bed")
    bim_path = out_prefix.with_suffix(".bim")
    fam_path = out_prefix.with_suffix(".fam")

    kept = 0
    skipped_indel = 0
    skipped_multiallelic = 0
    skipped_non_acgt = 0
    skipped_transition = 0
    skipped_maf = 0
    skipped_missing = 0

    with pysam.VariantFile(str(in_file)) as vcf:
        samples = list(vcf.header.samples)
        contigs = set(vcf.header.contigs.keys())

        validate_map_coverage("pop-map", pop_map, samples, strict_maps)
        validate_map_coverage("sex-map", sex_map, samples, strict_maps)
        validate_map_coverage("chrom-map", chrom_map, contigs, strict_maps)

        _write_fam(fam_path, samples, pop_map, sex_map)

        with open(bed_path, "wb") as bed_fh, open(bim_path, "w") as bim_fh:
            bed_fh.write(_PLINK_MAGIC)
            for rec in vcf:
                if not rec.alts or len(rec.alts) != 1:
                    skipped_multiallelic += 1
                    continue
                ref_a = rec.ref.upper()
                alt_a = rec.alts[0].upper()
                if len(ref_a) != 1 or len(alt_a) != 1:
                    skipped_indel += 1
                    continue
                if ref_a not in BASES or alt_a not in BASES:
                    skipped_non_acgt += 1
                    continue

                info_aa_val = None
                if info_aa:
                    aa = rec.info.get("AA")
                    if isinstance(aa, (list, tuple)):
                        aa = aa[0] if aa else None
                    info_aa_val = str(aa) if aa else None

                major, minor, flipped = polarise(
                    ref_a, alt_a, rec.chrom, rec.pos, ancestral, info_aa_val
                )

                if transversions_only and is_transition(major, minor):
                    skipped_transition += 1
                    continue

                row = "".join(
                    _eig_code(rec.samples[s].allele_indices, flipped)
                    for s in samples
                )

                if len(row) != len(samples):
                    log_warn(
                        f"{in_file.name}: row-length mismatch at {rec.chrom}:{rec.pos} "
                        f"({len(row)} vs {len(samples)} samples). Skipping."
                    )
                    continue

                if not passes_missing(row, max_missing):
                    skipped_missing += 1
                    continue
                if not passes_maf(row, min_maf):
                    skipped_maf += 1
                    continue

                snp_id = rec.id if rec.id not in (None, ".") else f"{rec.chrom}_{rec.pos}"
                emitted_chrom = resolve_chrom(rec.chrom, chrom_map, default_chrom)
                morgans = morgans_for(rec.chrom, rec.pos, gmap)

                # .bim format: chrom snp_id cM bp allele1(minor) allele2(major)
                bim_fh.write(
                    f"{emitted_chrom}\t{snp_id}\t{morgans * 100.0:.6f}\t"
                    f"{rec.pos}\t{minor}\t{major}\n"
                )
                bed_fh.write(_pack_snp_row(row))
                kept += 1

    log_info(
        f"{in_file.name}: {kept} biallelic SNPs written (PLINK binary); "
        f"skipped indels={skipped_indel}, multiallelic={skipped_multiallelic}, "
        f"non-ACGT={skipped_non_acgt}, transitions={skipped_transition}, "
        f"maf={skipped_maf}, missing={skipped_missing}."
    )


def batch_convert(input_dir: str, output_dir: str,
                  pop_map: Optional[str] = None,
                  sex_map: Optional[str] = None,
                  chrom_map: Optional[str] = None,
                  default_chrom: Optional[int] = None,
                  genetic_map: Optional[str] = None,
                  ancestral_fasta: Optional[str] = None,
                  info_aa: bool = False,
                  transversions_only: bool = False,
                  min_maf: float = 0.0,
                  max_missing: float = 1.0,
                  strict_maps: bool = False,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    if pysam is None:
        require_pysam()
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)

    pop = load_two_col_map(pop_map)
    sex = load_two_col_map(sex_map)
    cmap = load_chrom_map(chrom_map)
    gmap = load_genetic_map(genetic_map)

    with AncestralProvider(ancestral_fasta) as ancestral:
        for vcf_file in iter_input_files(in_path, _VCF_EXTS, pattern=pattern):
            stem = strip_compound_suffix(vcf_file, _VCF_EXTS)
            out_prefix = out_path / stem
            log_info(f"Converting {vcf_file.name} -> {stem}.bed / .bim / .fam")
            _convert_file(
                vcf_file, out_prefix,
                pop, sex, cmap, default_chrom, gmap, ancestral, info_aa,
                transversions_only, min_maf, max_missing, strict_maps,
            )
