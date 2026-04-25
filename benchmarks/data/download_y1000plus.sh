#!/usr/bin/env bash
# Download a 20-species Y1000+ slice (Shen 2018 / Opulente 2024).
# For the 1154-species full run, set GFC_BENCH_Y1000_FULL=1.
#
# Target: benchmarks/data/y1000plus/<species>.{fasta,final.gff3}
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
dest="$here/y1000plus"
mkdir -p "$dest"

# TODO: pin the exact figshare / NCBI URLs from the Shen 2018 paper and
# the Opulente 2024 update. The 20-species slice should be representative
# across major Saccharomycotina clades (Saccharomyces, Candida, Yarrowia,
# Lipomyces, Ascoidea, Dipodascaceae, etc.) to avoid clade bias.
echo "[TODO] populate benchmarks/data/y1000plus/ with 20 species" >&2
echo "       see benchmarks/data/README.md for dataset provenance"   >&2

# Placeholder: so downstream scripts don't crash if nothing is there yet.
[[ -d "$dest" ]] || mkdir -p "$dest"
