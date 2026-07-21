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

# ---------- AGAT agat_convert_sp_gff2bed.pl ---------------------------------
# AGAT's bed converter is the "deep validation" reference for T5. Same
# trailing-pipe version-string fixup as T4. Output goes to .bed (not BED12)
# but the timing comparison is what matters here.
if command -v agat_convert_sp_gff2bed.pl >/dev/null 2>&1; then
    agat_version="$(agat_convert_sp_gff2bed.pl --help 2>&1 | grep -i 'Version' | head -1 | tr -s ' ' | cut -d':' -f2- | tr -d ' |' || echo unknown)"
    : "${agat_version:=unknown}"
    for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
        out_dir="$bench_dir/agat_rep${rep}"
        mkdir -p "$out_dir"
        # AGAT 1.4 unconditionally writes <input>.agat.log next to its
        # CWD; the cd confines that clutter to the per-rep output dir.
        # Same fix as T4. Trailing `; true` keeps the loop's exit clean
        # even when one file fails — bench_one.py captures stderr.
        python "$repo/benchmarks/bench_one.py" \
            --task "T5" --tool "AGAT" --version "$agat_version" \
            --replicate "$rep" \
            --cmd "cd '$out_dir' && \
                   for f in '$GFC_BENCH_INPUT_DIR'/*.gff3; do \
                      stem=\$(basename \"\$f\" .gff3); \
                      agat_convert_sp_gff2bed.pl --gff \"\$f\" -o \"\$stem.bed\"; \
                   done; true" \
            >> "$out"
    done
else
    echo "[skip] T5 AGAT (agat_convert_sp_gff2bed.pl not on PATH)" >&2
fi

echo "[done] T5 rows written to $out" >&2

# Sanity-check the per-tool rep1 dirs BEFORE invoking the correctness
# comparator so the bench log explains a blank `correct` cell at the
# earliest possible point. The Apr 28 scarcity-16 run produced blank
# correctness for both gfc and AGAT because gff3ToGenePred (UCSC's
# stricter-than-spec GFF3 parser) silently emitted zero .bed files on
# the Y1000+ tRNA-scan slice — the `; true` in the UCSC loop above
# swallowed the per-file failures and check_correctness.py had nothing
# to compare against.
for tool_subdir in gfc_rep1 ucsc_rep1 agat_rep1; do
    d="$bench_dir/$tool_subdir"
    if [[ -d "$d" ]]; then
        n_bed=$(find "$d" -maxdepth 1 \( -name "*.bed" -o -name "*.bed12" \) | wc -l | tr -d ' ')
        echo "[corr-precheck] $tool_subdir: $n_bed .bed/.bed12 file(s)" >&2
        if [[ "$n_bed" -eq 0 ]]; then
            echo "[corr-precheck] WARNING: $tool_subdir wrote zero .bed/.bed12 files — correctness will be blank for this tool" >&2
        fi
    fi
done

python "$repo/benchmarks/check_correctness.py" --task T5 --bench-dir "$bench_dir" --tsv "$out" 2>&1 | head -20 || echo "[warn] T5 correctness check failed (non-fatal)" >&2
