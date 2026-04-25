#!/usr/bin/env bash
# T3 — Paired FASTA + GFF3 → GenBank.
# Competitors: EMBOSS seqret, handwritten Biopython, gff3toembl.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
task="T3_fasta_gff_to_gbk"
out="$repo/benchmarks/results/raw/${task}.tsv"
mkdir -p "$(dirname "$out")"

input_dir="$repo/tests/test_data"
: "${GFC_BENCH_INPUT_DIR:=$input_dir}"
: "${GFC_BENCH_REPLICATES:=5}"

bench_dir="/tmp/gfc_bench_${task}"
rm -rf "$bench_dir" && mkdir -p "$bench_dir"

python "$repo/benchmarks/bench_one.py" --header > "$out"

gfc_version="$(gfc --version | cut -d' ' -f2)"
for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
    out_dir="$bench_dir/gfc_rep${rep}"
    python "$repo/benchmarks/bench_one.py" \
        --task "T3" --tool "gfc" --version "$gfc_version" --replicate "$rep" \
        --cmd "gfc fasta-gff-to-gbk --input-dir '$GFC_BENCH_INPUT_DIR' \
               --output-dir '$out_dir' --force" \
        --notes "20-species Y1000+ slice when available" \
        >> "$out"
done

# TODO: EMBOSS seqret -sequence <fa> -feature -fformat gff -osformat genbank
# TODO: handwritten biopython reference (scripts/ref/t3_biopython.py)
echo "[done] T3 rows written to $out" >&2
