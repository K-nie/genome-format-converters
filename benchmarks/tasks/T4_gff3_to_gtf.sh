#!/usr/bin/env bash
# T4 — GFF3 → GTF.  Competitors: AGAT agat_sp_gff2gtf.pl, gffread.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
task="T4_gff3_to_gtf"
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
        --task "T4" --tool "gfc" --version "$gfc_version" --replicate "$rep" \
        --cmd "gfc gff3-to-gtf --input-dir '$GFC_BENCH_INPUT_DIR' \
               --output-dir '$out_dir' --pattern '*.gff3' --force" \
        >> "$out"
done
# ---------- gffread ---------------------------------------------------------
# gffread -T emits GTF on stdout from a GFF3 input; iterate over files in the
# input dir so the benchmark still reflects batch behaviour.
if command -v gffread >/dev/null 2>&1; then
    gffread_version="$(gffread --version 2>&1 | head -1)"
    for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
        out_dir="$bench_dir/gffread_rep${rep}"
        mkdir -p "$out_dir"
        python "$repo/benchmarks/bench_one.py" \
            --task "T4" --tool "gffread" --version "$gffread_version" \
            --replicate "$rep" \
            --cmd "for f in '$GFC_BENCH_INPUT_DIR'/*.gff3; do \
                      stem=\$(basename \"\$f\" .gff3); \
                      gffread -T \"\$f\" -o \"$out_dir/\$stem.gtf\"; \
                   done" \
            >> "$out"
    done
else
    echo "[skip] T4 gffread (gffread not on PATH)" >&2
fi

# ---------- AGAT ------------------------------------------------------------
# AGAT ships a dedicated agat_sp_gff2gtf.pl. It is a Perl tool that does
# deep validation/conversion of the GFF hierarchy, so it is substantially
# slower than gffread or gfc — this is by design and worth showing.
if command -v agat_sp_gff2gtf.pl >/dev/null 2>&1; then
    agat_version="$(agat_sp_gff2gtf.pl --help 2>&1 | grep -i 'Version' | head -1 | tr -s ' ' | cut -d':' -f2- | tr -d ' ' || echo unknown)"
    for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
        out_dir="$bench_dir/agat_rep${rep}"
        mkdir -p "$out_dir"
        python "$repo/benchmarks/bench_one.py" \
            --task "T4" --tool "AGAT" --version "$agat_version" \
            --replicate "$rep" \
            --cmd "for f in '$GFC_BENCH_INPUT_DIR'/*.gff3; do \
                      stem=\$(basename \"\$f\" .gff3); \
                      agat_sp_gff2gtf.pl --gff \"\$f\" -o \"$out_dir/\$stem.gtf\"; \
                   done" \
            >> "$out"
    done
else
    echo "[skip] T4 AGAT (agat_sp_gff2gtf.pl not on PATH)" >&2
fi

echo "[done] T4 rows written to $out" >&2
