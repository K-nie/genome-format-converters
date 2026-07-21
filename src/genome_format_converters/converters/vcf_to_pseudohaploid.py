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

import zlib
from pathlib import Path
from typing import Optional

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
    write_ind,
    write_plink_sidecar,
)

_VCF_EXTS = [".vcf", ".vcf.gz", ".bcf"]

# Reference per-sample encoder kept for tests / external callers. Produces
# byte-identical output to the bulk numpy path below when fed the same
# allele tuples and a Python `random.Random` seeded identically — but the
# bulk path uses numpy's RNG, so cross-implementation byte-equivalence is
# NOT preserved (the reproducibility contract is "same seed -> same output
# in this implementation", not across RNG backends; the test suite only
# checks the within-implementation form).
def _pseudohaploid_code(alleles, flipped: bool, rng) -> str:
    """Randomly pick one of the genotype alleles; emit '2' for major, '0'
    for minor, '9' if any allele is missing.

    `rng` may be either a `random.Random` or a `numpy.random.Generator`.
    """
    if alleles is None or any(a is None for a in alleles):
        return "9"
    if hasattr(rng, "choice") and not hasattr(rng, "integers"):
        # stdlib random.Random.choice
        picked = rng.choice(list(alleles))
    else:
        # numpy.random.Generator
        picked = alleles[int(rng.integers(0, len(alleles)))]
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
                  rng: np.random.Generator) -> None:
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

    # Streaming write inside the read loop so chr22-scale inputs don't sit
    # in RAM. The PLINK sidecar (also_plink=True) still needs the rows
    # materialised, so accumulate only when that flag is set.
    snp_rows: list = [] if also_plink else None
    geno_rows: list = [] if also_plink else None

    _BYTE_0 = ord("0")
    _BYTE_2 = ord("2")
    _BYTE_9 = ord("9")

    with cyvcf2.VCF(str(in_file)) as vcf, \
            open(geno_path, "w") as geno_fh, \
            open(snp_path, "w") as snp_fh:
        samples = list(vcf.samples)
        contigs = set(vcf.seqnames)

        validate_map_coverage("pop-map", pop_map, samples, strict_maps)
        validate_map_coverage("sex-map", sex_map, samples, strict_maps)
        validate_map_coverage("chrom-map", chrom_map, contigs, strict_maps)

        write_ind(ind_path, samples, pop_map, sex_map)

        n_samples = len(samples)
        flipped_int = 1  # placeholder, set per record below

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

            # Bulk pseudohaploid encoding via numpy. cyvcf2's gt_types is the
            # int32 array of {0=HOM_REF, 1=HET, 2=MISSING, 3=HOM_ALT}. We
            # draw one random {0,1} per sample (only the HET draws actually
            # affect output, but generating in bulk is cheap and avoids
            # branching). Then "picked" is the allele index per sample
            # (0=REF, 1=ALT). Output is '2' when picked==major, '0' otherwise,
            # overridden to '9' for missing.
            gt = var.gt_types
            random_picks = rng.integers(0, 2, size=n_samples, dtype=np.int32)
            picked = np.where(gt == 1, random_picks,
                              np.where(gt == 3, 1, 0)).astype(np.int32)
            flipped_int = 1 if flipped else 0
            is_major = picked == flipped_int
            codes = np.where(is_major, _BYTE_2, _BYTE_0).astype(np.uint8)
            codes = np.where(gt == 2, _BYTE_9, codes).astype(np.uint8)
            row = codes.tobytes().decode("ascii")

            if len(row) != n_samples:
                log_warn(
                    f"{in_file.name}: row-length mismatch at {var.CHROM}:{var.POS} "
                    f"({len(row)} vs {n_samples} samples). Skipping site."
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

            snp_fh.write(
                f"{snp_id}\t{emitted_chrom}\t{morgans}\t{var.POS}\t{major}\t{minor}\n"
            )
            geno_fh.write(row + "\n")

            if also_plink:
                snp_rows.append(
                    (snp_id, emitted_chrom, morgans, var.POS, major, minor)
                )
                geno_rows.append(row)
            kept += 1

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
    if cyvcf2 is None:
        raise ImportError(
            "cyvcf2 is required for vcf_to_pseudohaploid. "
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
            log_info(f"Converting {vcf_file.name} -> {stem}.geno / .snp / .ind (pseudohaploid)")
            # numpy.random.Generator with the same per-file seed produces the
            # same integer sequence on every run (the test suite verifies
            # `seed -> deterministic output` directly).
            rng = np.random.default_rng(_file_seed(seed, vcf_file.name))
            _convert_file(
                vcf_file, out_prefix,
                pop, sex, cmap, default_chrom, gmap, ancestral, info_aa,
                transversions_only, min_maf, max_missing, strict_maps,
                also_plink, rng,
            )
