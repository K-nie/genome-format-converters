#!/usr/bin/env bash
# Y1000+ data provisioning for the gfc benchmark suite.
#
# Stage 1 (current): a 20-species tRNA-scan FASTA + GFF slice ships
# committed in the repo at benchmarks/data/y1000plus/. No download
# required for T3 / T4 / T5 input.
#
# Stage 2 (planned): pull a full-genome 20-species slice from the
# Shen 2018 / Opulente 2024 figshare deposits for paper-grade T3 numbers.
# Set GFC_BENCH_Y1000_FULL=1 to opt in once the URLs are pinned. Tracked
# in project_gfc_bench_figures.md → "Stage 2".
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
dest="$here/y1000plus"
mkdir -p "$dest"

slice_count=$(ls "$dest"/*.fasta 2>/dev/null | wc -l | tr -d ' ')
if [[ "$slice_count" -ge 20 ]]; then
    echo "[skip] Y1000+ 20-species slice already present ($slice_count fasta files)"
    exit 0
fi

cat >&2 <<'EOM'
[warn] Y1000+ slice missing from benchmarks/data/y1000plus/
       The slice is committed in the repo and should be present after a
       fresh clone or rsync of the full repository. If this directory is
       empty, the slice was either pruned or the .gitignore whitelist for
       data/y1000plus/*.fasta + *.gff was reverted.

       To restore from a clean checkout:
         git checkout -- benchmarks/data/y1000plus/

       To pull the Stage 2 full-genome slice (not yet wired):
         GFC_BENCH_Y1000_FULL=1 bash benchmarks/data/download_y1000plus.sh
EOM
exit 0
