#!/usr/bin/env python3
"""Convert VCF / BCF to pseudohaploid EIGENSTRAT.

For each biallelic SNP and each sample, one of the two genotype alleles is
drawn uniformly at random and emitted as a homozygous call. The resulting
`.geno` file contains only `0`, `2`, and `9` (missing) — never `1` — which
is the encoding expected by `admixtools` / `smartpca` when ancient or
low-coverage samples are represented as pseudohaploid pulls (cf.
`pileupCaller`'s random-call mode).

`--seed` controls reproducibility. The RNG is re-seeded per input file
(`seed XOR stable_hash(filename)`) so that the per-file output does not
depend on the order in which files are processed.

All polarisation / filter / sidecar flags from `vcf-to-eigenstrat` are
supported.
"""

import random
import zlib
from pathlib import Path
from typing import Optional

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


def _pseudohaploid_code(alleles, flipped: bool, rng: random.Random) -> str:
    """Randomly pick one of the genotype alleles; emit '2' for major, '0'
    for minor, '9' if any allele is missing."""
    if alleles is None or any(a is None for a in alleles):
        return "9"
    picked = rng.choice(list(alleles))
    major_idx = 1 if flipped else 0
    return "2" if picked == major_idx else "0"


def _file_seed(base_seed: Optional[int], file_name: str) -> Optional[int]:
    if base_seed is None:
        return None
    # zlib.crc32 gives a deterministic 32-bit hash that doesn't vary by
    # PYTHONHASHSEED (unlike the built-in `hash()`).
    file_token = zlib.crc32(file_name.encode("utf-8"))
    return base_seed ^ file_token


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
                  also_plink: bool,
                  rng: random.Random) -> None:
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

    with pysam.VariantFile(str(in_file)) as vcf:
        samples = list(vcf.header.samples)
        contigs = set(vcf.header.contigs.keys())

        validate_map_coverage("pop-map", pop_map, samples, strict_maps)
        validate_map_coverage("sex-map", sex_map, samples, strict_maps)
        validate_map_coverage("chrom-map", chrom_map, contigs, strict_maps)

        write_ind(ind_path, samples, pop_map, sex_map)

        snp_rows = []
        geno_rows = []

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
                _pseudohaploid_code(rec.samples[s].allele_indices, flipped, rng)
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

            snp_rows.append(
                (snp_id, emitted_chrom, morgans, rec.pos, major, minor)
            )
            geno_rows.append(row)
            kept += 1

        with open(geno_path, "w") as geno_fh, open(snp_path, "w") as snp_fh:
            for (snp_id, chrom, morgans, bp, major, minor), row in zip(snp_rows, geno_rows):
                snp_fh.write(f"{snp_id}\t{chrom}\t{morgans}\t{bp}\t{major}\t{minor}\n")
                geno_fh.write(row + "\n")

    log_info(
        f"{in_file.name}: {kept} biallelic SNPs written (pseudohaploid); "
        f"skipped indels={skipped_indel}, multiallelic={skipped_multiallelic}, "
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
                  seed: Optional[int] = None,
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
            log_info(f"Converting {vcf_file.name} -> {stem}.geno / .snp / .ind (pseudohaploid)")
            rng = random.Random(_file_seed(seed, vcf_file.name))
            _convert_file(
                vcf_file, out_prefix,
                pop, sex, cmap, default_chrom, gmap, ancestral, info_aa,
                transversions_only, min_maf, max_missing, strict_maps,
                also_plink, rng,
            )
