#!/usr/bin/env python3
"""Convert VCF / BCF to EIGENSTRAT (`.geno` / `.snp` / `.ind`).

Biallelic SNPs only. Genotype encoding (count of major/reference alleles):

    2 = homozygous for the major allele
    1 = heterozygous
    0 = homozygous for the minor allele
    9 = missing

The `.snp` file's chromosome column can be re-mapped to integers via
`--chrom-map` / `--default-chrom` for legacy `smartpca`. Genetic positions
are 0.0 unless a PLINK-style `.map` is supplied via `--genetic-map`. Allele
polarisation can be driven from an outgroup `--ancestral-fasta` or the
`--info-aa` field of the VCF.
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
    write_ind,
    write_plink_sidecar,
)

_VCF_EXTS = [".vcf", ".vcf.gz", ".bcf"]


def _genotype_code(alleles, flipped: bool) -> str:
    """Count copies of the major allele; returns '9' if any allele missing.

    When `flipped=True` the original ref/alt were swapped during polarisation,
    so we count copies of the VCF ALT instead (major in post-polarisation
    terms).
    """
    if alleles is None or any(a is None for a in alleles):
        return "9"
    # a == 0 -> VCF REF; a == 1 -> VCF ALT.
    # Major allele in post-polarisation terms: flipped=False -> REF, else ALT.
    major_idx = 1 if flipped else 0
    major_count = sum(1 for a in alleles if a == major_idx)
    if major_count == 2:
        return "2"
    if major_count == 1:
        return "1"
    return "0"


def _convert_file(in_file: Path, out_prefix: Path,
                  pop_map, sex_map, chrom_map,
                  default_chrom: Optional[int],
                  gmap,
                  ancestral: AncestralProvider,
                  info_aa: bool,
                  transversions_only: bool,
                  min_maf: float,
                  max_missing: float,
                  strict_maps: bool,
                  also_plink: bool) -> None:
    geno_path = out_prefix.with_suffix(".geno")
    snp_path = out_prefix.with_suffix(".snp")
    ind_path = out_prefix.with_suffix(".ind")

    kept = 0
    skipped_indel = 0
    skipped_multiallelic = 0
    skipped_non_acgt = 0
    skipped_transition = 0
    skipped_maf = 0
    skipped_missing = 0

    # Stream `.geno` and `.snp` writes inside the read loop so chr22-scale
    # inputs don't accumulate in RAM. Only the `also_plink` path needs the
    # rows kept in memory (write_plink_sidecar consumes them after the loop).
    snp_rows: list = [] if also_plink else None
    geno_rows: list = [] if also_plink else None

    with pysam.VariantFile(str(in_file)) as vcf, \
            open(geno_path, "w") as geno_fh, \
            open(snp_path, "w") as snp_fh:
        samples = list(vcf.header.samples)
        contigs = set(vcf.header.contigs.keys())

        validate_map_coverage("pop-map", pop_map, samples, strict_maps)
        validate_map_coverage("sex-map", sex_map, samples, strict_maps)
        validate_map_coverage("chrom-map", chrom_map, contigs, strict_maps)

        write_ind(ind_path, samples, pop_map, sex_map)

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
                _genotype_code(rec.samples[s].allele_indices, flipped)
                for s in samples
            )

            if len(row) != len(samples):
                log_warn(
                    f"{in_file.name}: row-length mismatch at {rec.chrom}:{rec.pos} "
                    f"({len(row)} vs {len(samples)} samples). Skipping site."
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

            snp_fh.write(
                f"{snp_id}\t{emitted_chrom}\t{morgans}\t{rec.pos}\t{major}\t{minor}\n"
            )
            geno_fh.write(row + "\n")

            if also_plink:
                snp_rows.append(
                    (snp_id, emitted_chrom, morgans, rec.pos, major, minor)
                )
                geno_rows.append(row)
            kept += 1

    log_info(
        f"{in_file.name}: {kept} biallelic SNPs written; skipped "
        f"indels={skipped_indel}, multiallelic={skipped_multiallelic}, "
        f"non-ACGT={skipped_non_acgt}, transitions={skipped_transition}, "
        f"maf={skipped_maf}, missing={skipped_missing}."
    )

    if also_plink:
        write_plink_sidecar(out_prefix, samples, pop_map, sex_map, snp_rows, geno_rows)


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
                  also_plink: bool = False,
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
            log_info(f"Converting {vcf_file.name} -> {stem}.geno / .snp / .ind")
            _convert_file(
                vcf_file, out_prefix,
                pop, sex, cmap, default_chrom, gmap, ancestral, info_aa,
                transversions_only, min_maf, max_missing, strict_maps,
                also_plink,
            )
