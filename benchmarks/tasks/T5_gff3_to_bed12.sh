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

# ---------- UCSC gff3ToGenePred + genePredToBed chain -----------------------
# Two-step pipeline that's the canonical reference for genePred-aware tooling.
# Same input dir as gfc; per-file output is `.gp` (intermediate) then `.bed`.
# UCSC tools don't expose a clean `--version`, so the version string is
# pulled from the active conda env's package metadata.
if command -v gff3ToGenePred >/dev/null 2>&1 && command -v genePredToBed >/dev/null 2>&1; then
    ucsc_version="$(conda list ucsc-gff3togenepred 2>/dev/null \
        | awk '$1=="ucsc-gff3togenepred" {print "v"$2}' | head -1)"
    : "${ucsc_version:=unknown}"
    for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
        out_dir="$bench_dir/ucsc_rep${rep}"
        mkdir -p "$out_dir"
        # Per-file `&&` chains so genePredToBed only runs when its input
        # exists; trailing `true` makes the whole loop exit 0 even when
        # one file's GFF3 is rejected by gff3ToGenePred (UCSC's parser is
        # stricter than the GFF3 spec). Per-file failures are captured by
        # bench_one.py's per-run stderr log.
        python "$repo/benchmarks/bench_one.py" \
            --task "T5" --tool "ucsc-chain" --version "$ucsc_version" \
            --replicate "$rep" \
            --cmd "for f in '$GFC_BENCH_INPUT_DIR'/*.gff3; do \
                      stem=\$(basename \"\$f\" .gff3); \
                      gff3ToGenePred \"\$f\" '$out_dir/'\$stem.gp \
                        && genePredToBed '$out_dir/'\$stem.gp '$out_dir/'\$stem.bed; \
                   done; true" \
            --notes "gff3ToGenePred + genePredToBed chain" \
            >> "$out"
    done
else
    echo "[skip] T5 UCSC chain (gff3ToGenePred / genePredToBed not on PATH)" >&2
fi

echo "[done] T5 rows written to $out" >&2
