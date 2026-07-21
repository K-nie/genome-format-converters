#!/usr/bin/env bash
# Download 1000 Genomes Phase 3 chr22, biallelic SNPs only.
# Target: benchmarks/data/1kg/chr22.biallelic.snps.vcf.gz (+ .tbi)
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
dest="$here/1kg"
mkdir -p "$dest"

vcf="$dest/chr22.biallelic.snps.vcf.gz"
if [[ -f "$vcf" ]]; then
    echo "[skip] $vcf already present" >&2
    exit 0
fi

# URL target. EBI's 1000G FTP path 404s as of 2026; NCBI mirror still serves
# the identical phase3 chr22 VCF. Update if the upstream moves again.
src="https://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

echo "[get ] $src" >&2
curl -fSL --retry 3 -o "$vcf" "$src"
curl -fSL --retry 3 -o "${vcf}.tbi" "${src}.tbi" || {
    echo "[warn] index download failed; regenerating"
    command -v tabix >/dev/null && tabix -p vcf "$vcf"
}

# Subset to biallelic SNPs for benchmark consistency (T1/T2/T8 assume biallelic).
if command -v bcftools >/dev/null 2>&1; then
    bcftools view -m2 -M2 -v snps "$vcf" -Oz -o "$dest/chr22.biallelic.snps.filtered.vcf.gz"
    command -v tabix >/dev/null && tabix -p vcf "$dest/chr22.biallelic.snps.filtered.vcf.gz"
fi
