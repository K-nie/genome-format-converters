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

import numpy as np

try:
    import cyvcf2
except ImportError:
    cyvcf2 = None

try:
    import pysam  # only needed when --ancestral-fasta is supplied
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
    ref/alt. Kept for reference / external callers; the chr22 hot path goes
    through the bulk numpy encoder below."""
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

# cyvcf2.Variant.gt_types -> EIGENSTRAT byte ('0' / '1' / '2' / '9'), used
# only for the existing string-based passes_missing / passes_maf filters.
# Same lookups as in vcf_to_eigenstrat.py — should be consolidated in
# _eigenstrat_common.py in a follow-up.
_EIG_LOOKUP_UNFLIPPED = np.array(
    [ord("2"), ord("1"), ord("9"), ord("0")], dtype=np.uint8
)
_EIG_LOOKUP_FLIPPED = np.array(
    [ord("0"), ord("1"), ord("9"), ord("2")], dtype=np.uint8
)

# cyvcf2.Variant.gt_types -> PLINK 2-bit code, used for the .bed payload.
# PLINK 1.9 encoding (per the docstring at the top of this module):
#   00 = hom A1 (minor),  01 = missing,  10 = het,  11 = hom A2 (major).
# cyvcf2 gt_types: 0 = HOM_REF, 1 = HET, 2 = MISSING, 3 = HOM_ALT.
# Unflipped (REF is major):
#   HOM_REF -> hom-major -> 11 (3); HET -> 10 (2); MISSING -> 01 (1); HOM_ALT -> hom-minor -> 00 (0).
# Flipped (ALT is major):
#   HOM_REF -> hom-minor -> 00 (0); HET -> 10 (2); MISSING -> 01 (1); HOM_ALT -> hom-major -> 11 (3).
_PLINK_BITS_UNFLIPPED = np.array([3, 2, 1, 0], dtype=np.uint8)
_PLINK_BITS_FLIPPED = np.array([0, 2, 1, 3], dtype=np.uint8)

# Per-byte bit shifts: sample 0 lands in bits 0-1, sample 1 in bits 2-3,
# sample 2 in bits 4-5, sample 3 in bits 6-7. Hoisted to module scope so
# we don't re-allocate it per record on the chr22 hot path.
_BYTE_SHIFTS = np.array([0, 2, 4, 6], dtype=np.uint8)


def _pack_snp_row(row: str) -> bytes:
    """Pack an EIGENSTRAT-encoded genotype row into PLINK .bed bytes.

    Reference implementation kept for tests and external callers. The chr22
    hot path uses `_pack_bits_bulk` below, which packs straight from a
    numpy uint8 array of 2-bit codes without a per-sample Python loop.

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


def _pack_bits_bulk(bits: np.ndarray) -> bytes:
    """Pack n_samples 2-bit PLINK codes into ceil(n/4) bytes via numpy.

    Produces output byte-equivalent to `_pack_snp_row` when fed the same
    sample sequence: sample i lands in byte (i // 4), bit-offset 2 * (i %
    4); trailing slots in the final byte are padded with 01 (missing).
    """
    n = bits.shape[0]
    n_bytes = (n + 3) // 4
    padded = np.full(n_bytes * 4, 0b01, dtype=np.uint8)
    padded[:n] = bits
    packed = (padded.reshape(n_bytes, 4) << _BYTE_SHIFTS).sum(axis=1).astype(np.uint8)
    return packed.tobytes()


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

    with cyvcf2.VCF(str(in_file)) as vcf, \
            open(bed_path, "wb") as bed_fh, \
            open(bim_path, "w") as bim_fh:
        samples = list(vcf.samples)
        contigs = set(vcf.seqnames)

        validate_map_coverage("pop-map", pop_map, samples, strict_maps)
        validate_map_coverage("sex-map", sex_map, samples, strict_maps)
        validate_map_coverage("chrom-map", chrom_map, contigs, strict_maps)

        _write_fam(fam_path, samples, pop_map, sex_map)

        bed_fh.write(_PLINK_MAGIC)
        n_samples = len(samples)

        for var in vcf:
            if not var.ALT or len(var.ALT) != 1:
                skipped_multiallelic += 1
                continue
            ref_a = var.REF.upper()
            alt_a = var.ALT[0].upper()
            if len(ref_a) != 1 or len(alt_a) != 1:
                skipped_indel += 1
                continue
            if ref_a not in BASES or alt_a not in BASES:
                skipped_non_acgt += 1
                continue

            info_aa_val = None
            if info_aa:
                aa = var.INFO.get("AA")
                if isinstance(aa, (list, tuple)):
                    aa = aa[0] if aa else None
                info_aa_val = str(aa) if aa else None

            major, minor, flipped = polarise(
                ref_a, alt_a, var.CHROM, var.POS, ancestral, info_aa_val
            )

            if transversions_only and is_transition(major, minor):
                skipped_transition += 1
                continue

            # Bulk encode: gt_types is the per-record numpy int32 array of
            # {0,1,2,3}. One numpy gather builds the EIGENSTRAT row (used
            # only for the existing string-based filters); a second gather
            # builds the PLINK 2-bit codes that get packed for the .bed.
            gt = var.gt_types
            eig_codes = (_EIG_LOOKUP_FLIPPED if flipped else _EIG_LOOKUP_UNFLIPPED)[gt]
            row = eig_codes.tobytes().decode("ascii")

            if len(row) != n_samples:
                log_warn(
                    f"{in_file.name}: row-length mismatch at {var.CHROM}:{var.POS} "
                    f"({len(row)} vs {n_samples} samples). Skipping."
                )
                continue

            if not passes_missing(row, max_missing):
                skipped_missing += 1
                continue
            if not passes_maf(row, min_maf):
                skipped_maf += 1
                continue

            snp_id = var.ID if var.ID not in (None, ".") else f"{var.CHROM}_{var.POS}"
            emitted_chrom = resolve_chrom(var.CHROM, chrom_map, default_chrom)
            morgans = morgans_for(var.CHROM, var.POS, gmap)

            # .bim format: chrom snp_id cM bp allele1(minor) allele2(major)
            bim_fh.write(
                f"{emitted_chrom}\t{snp_id}\t{morgans * 100.0:.6f}\t"
                f"{var.POS}\t{minor}\t{major}\n"
            )
            bits = (_PLINK_BITS_FLIPPED if flipped else _PLINK_BITS_UNFLIPPED)[gt]
            bed_fh.write(_pack_bits_bulk(bits))
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
    if cyvcf2 is None:
        raise ImportError(
            "cyvcf2 is required for vcf_to_plink. "
            "Install via `pip install cyvcf2` or `mamba install -c bioconda cyvcf2`."
        )
    if ancestral_fasta is not None and pysam is None:
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
