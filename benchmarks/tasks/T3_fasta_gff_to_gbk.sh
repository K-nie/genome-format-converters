#!/usr/bin/env bash
# T3 — Paired FASTA + GFF3 → GenBank.
# Competitors: EMBOSS seqret, handwritten Biopython, gff3toembl.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
task="T3_fasta_gff_to_gbk"
out="$repo/benchmarks/results/raw/${task}.tsv"
mkdir -p "$(dirname "$out")"

# Stage 1: prefer the committed Y1000+ tRNA-scan slice over the tiny
# fixture when present, regardless of what run_on_condor.sh exported into
# GFC_BENCH_INPUT_DIR (which it sets to the 1kg VCF dir for T1/T2/T8).
y1000_dir="$repo/benchmarks/data/y1000plus"
if compgen -G "$y1000_dir/*.fasta" >/dev/null 2>&1; then
    GFC_BENCH_INPUT_DIR="$y1000_dir"
    notes_default="20-species Y1000+ tRNA-scan slice"
else
    GFC_BENCH_INPUT_DIR="${GFC_BENCH_INPUT_DIR:-$repo/tests/test_data}"
    notes_default="tiny fixture (no Y1000+ slice present)"
fi
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
        --notes "$notes_default" \
        >> "$out"
done

# TODO: EMBOSS seqret -sequence <fa> -feature -fformat gff -osformat genbank
# TODO: handwritten biopython reference (scripts/ref/t3_biopython.py)
echo "[done] T3 rows written to $out" >&2
