#!/usr/bin/env bash
# Real gene annotations for T4 (GFF3 -> GTF) and T5 (GFF3 -> BED12).
#
# The Stage 1 Y1000+ slice is tRNA-scan output: single-exon features with
# no CDS or UTR structure, which does not exercise transcript handling.
# These two Ensembl annotations do:
#
#   sc_r64/    S. cerevisiae R64-1-1, whole genome (~7k genes, few introns)
#   hs_chr22/  GRCh38 chromosome 22 (~1.4k protein-coding genes, multi-exon,
#              CDS + UTRs, many isoforms); same chromosome as T1/T2/T8.
#
# Pinned to Ensembl release 112 so the paper's numbers stay reproducible.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
rel="${GFC_BENCH_ENSEMBL_RELEASE:-112}"
base="https://ftp.ensembl.org/pub/release-${rel}/gff3"

fetch() {  # fetch <url> <dest.gff3>
    local url="$1" dest="$2"
    if [[ -s "$dest" ]]; then
        echo "[skip] $dest already present" >&2
        return
    fi
    mkdir -p "$(dirname "$dest")"
    echo "[get ] $url" >&2
    curl -fSL --retry 3 "$url" | gunzip -c > "$dest.tmp"
    mv "$dest.tmp" "$dest"
}

fetch "$base/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.${rel}.gff3.gz" \
      "$here/annotations/sc_r64/Saccharomyces_cerevisiae.R64-1-1.${rel}.gff3"
fetch "$base/homo_sapiens/Homo_sapiens.GRCh38.${rel}.chromosome.22.gff3.gz" \
      "$here/annotations/hs_chr22/Homo_sapiens.GRCh38.${rel}.chromosome.22.gff3"

for d in "$here"/annotations/*/; do
    n_mrna=$(awk -F'\t' '$3=="mRNA"' "$d"/*.gff3 | wc -l | tr -d ' ')
    n_exon=$(awk -F'\t' '$3=="exon"' "$d"/*.gff3 | wc -l | tr -d ' ')
    echo "[ok  ] $(basename "$d"): $n_mrna mRNA, $n_exon exon rows" >&2
done
