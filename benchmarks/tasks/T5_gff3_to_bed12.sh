#!/usr/bin/env bash
# T5 — GFF3 → BED12.  Competitors: UCSC gtfToGenePred + genePredToBed, AGAT.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
task="T5_gff3_to_bed12"
out="$repo/benchmarks/results/raw/${task}.tsv"
mkdir -p "$(dirname "$out")"

# Stage 1: prefer the committed Y1000+ tRNA-scan slice when present.
y1000_dir="$repo/benchmarks/data/y1000plus"
if compgen -G "$y1000_dir/*.gff3" >/dev/null 2>&1; then
    GFC_BENCH_INPUT_DIR="$y1000_dir"
else
    GFC_BENCH_INPUT_DIR="${GFC_BENCH_INPUT_DIR:-$repo/tests/test_data}"
fi
: "${GFC_BENCH_REPLICATES:=5}"

bench_dir="/tmp/gfc_bench_${task}"
rm -rf "$bench_dir" && mkdir -p "$bench_dir"

python "$repo/benchmarks/bench_one.py" --header > "$out"

gfc_version="$(gfc --version | cut -d' ' -f2)"
for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
    out_dir="$bench_dir/gfc_rep${rep}"
    python "$repo/benchmarks/bench_one.py" \
        --task "T5" --tool "gfc" --version "$gfc_version" --replicate "$rep" \
        --cmd "gfc gff3-to-bed12 --input-dir '$GFC_BENCH_INPUT_DIR' \
               --output-dir '$out_dir' --pattern '*.gff3' --force" \
        >> "$out"
done
# TODO: UCSC chain (gff3 → gtf → genepred → bed12).
echo "[done] T5 rows written to $out" >&2
