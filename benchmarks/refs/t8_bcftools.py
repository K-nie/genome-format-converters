#!/usr/bin/env python3
"""Handwritten bcftools-query + Python reference for T8 (pseudohaploid).

The conventional ad-hoc pipeline: pipe ``bcftools view -m2 -M2 -v snps |
bcftools query`` output into Python, randomly pick one of the two GT
alleles per sample, encode '2' / '0' / '9'. The gfc -> ref delta
measures the streaming-cyvcf2-to-numpy implementation versus the more
typical "shell glue" approach that papers like to compare against.

This reference does NOT polarise via INFO/AA or an ancestral FASTA --
it treats REF as the major allele unconditionally (matches the default
of gfc's vcf-to-pseudohaploid when no ancestral source is given).

Reproducibility: same ``--seed`` -> same output, file by file.
"""
from __future__ import annotations

import argparse
import random
import subprocess
import sys
import zlib
from pathlib import Path


_BIALLELIC_BASES = {"A", "C", "G", "T"}


def _file_seed(base_seed: int, file_name: str) -> int:
    return base_seed ^ zlib.crc32(file_name.encode("utf-8"))


def _samples_from_header(vcf_path: Path) -> list[str]:
    proc = subprocess.run(
        ["bcftools", "view", "-h", str(vcf_path)],
        capture_output=True, check=True, text=True,
    )
    for line in proc.stdout.splitlines():
        if line.startswith("#CHROM"):
            return line.split("\t")[9:]
    return []


def convert_one(vcf_path: Path, out_prefix: Path, seed: int) -> None:
    rng = random.Random(_file_seed(seed, vcf_path.name))
    samples = _samples_from_header(vcf_path)

    with open(f"{out_prefix}.ind", "w") as ind_fh:
        for s in samples:
            ind_fh.write(f"{s}\tU\t{s}\n")

    # Two-stage pipe: filter biallelic SNPs, then extract REF/ALT/GT cols.
    view = subprocess.Popen(
        ["bcftools", "view", "-m2", "-M2", "-v", "snps", "-Ou", str(vcf_path)],
        stdout=subprocess.PIPE,
    )
    query = subprocess.Popen(
        ["bcftools", "query", "-f", "%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n"],
        stdin=view.stdout, stdout=subprocess.PIPE, text=True,
    )
    if view.stdout is not None:
        view.stdout.close()

    geno_fh = open(f"{out_prefix}.geno", "w")
    snp_fh = open(f"{out_prefix}.snp", "w")
    try:
        for line in query.stdout:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            chrom, pos, snp_id, ref, alt = cols[0], cols[1], cols[2], cols[3], cols[4]
            gts = cols[5:]
            if len(ref) != 1 or len(alt) != 1:
                continue
            if ref not in _BIALLELIC_BASES or alt not in _BIALLELIC_BASES:
                continue
            row_chars: list[str] = []
            for gt in gts:
                # Normalise phasing separators; missing -> '9'.
                normed = gt.replace("|", "/")
                if normed in (".", "./.", "."):
                    row_chars.append("9")
                    continue
                alleles = normed.split("/")
                if len(alleles) != 2 or "." in alleles:
                    row_chars.append("9")
                    continue
                picked = rng.choice([alleles[0], alleles[1]])
                # REF is index "0"; treat REF as major (no polarisation).
                row_chars.append("2" if picked == "0" else "0")
            geno_fh.write("".join(row_chars) + "\n")
            sid = snp_id if snp_id and snp_id != "." else f"{chrom}_{pos}"
            snp_fh.write(f"{sid}\t{chrom}\t0.0\t{pos}\t{ref}\t{alt}\n")
    finally:
        geno_fh.close()
        snp_fh.close()
    query.wait()
    view.wait()
    if query.returncode not in (0, None):
        raise RuntimeError(f"bcftools query exited {query.returncode}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-dir", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--pattern", default="*.vcf*")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    in_dir = Path(args.input_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    found = 0
    for vcf in sorted(in_dir.glob(args.pattern)):
        # Strip compound suffixes the same way gfc does so output files line up.
        stem = vcf.name
        for ext in (".vcf.gz", ".vcf", ".bcf"):
            if stem.endswith(ext):
                stem = stem[: -len(ext)]
                break
        try:
            convert_one(vcf, out_dir / stem, args.seed)
            found += 1
        except Exception as exc:
            print(f"[skip-T8-ref] {vcf.name}: {exc}", file=sys.stderr)
    if found == 0:
        print("[warn-T8-ref] no VCF inputs matched the pattern", file=sys.stderr)


if __name__ == "__main__":
    main()
