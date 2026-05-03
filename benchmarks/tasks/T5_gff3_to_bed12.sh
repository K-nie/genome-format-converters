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

# ---------- bedops convert2bed (secondary) ---------------------------------
# bedops convert2bed (Neph 2012, Bioinformatics 28:1919) is the canonical
# BED-conversion utility from the BED-native toolkit. `convert2bed
# --input=gff` emits BED6 (not BED12 - bedops doesn't reconstruct
# transcript blockCounts/blockSizes), so it is a fairer comparison
# against AGAT's BED6 output than against the UCSC chain's BED12.
if command -v convert2bed >/dev/null 2>&1; then
    bedops_version="$(convert2bed --version 2>&1 | head -1 | tr -s ' ' | awk '{print $NF}' | tr -d '|')"
    : "${bedops_version:=unknown}"
    for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
        out_dir="$bench_dir/bedops_rep${rep}"
        mkdir -p "$out_dir"
        # convert2bed reads stdin and writes stdout. Per-file loop with
        # `< $f > $out_dir/$stem.bed` mirrors the UCSC chain.
        # `--input=gff` accepts gff/gff3 alike. Trailing `; true` keeps
        # one bad input from voiding the whole rep.
        python "$repo/benchmarks/bench_one.py" \
            --task "T5" --tool "bedops" --version "$bedops_version" \
            --replicate "$rep" \
            --cmd "for f in '$GFC_BENCH_INPUT_DIR'/*.gff3; do \
                      stem=\$(basename \"\$f\" .gff3); \
                      convert2bed --input=gff < \"\$f\" > '$out_dir/'\$stem.bed; \
                   done; true" \
            --notes "bedops convert2bed --input=gff (Neph 2012); BED6 output" \
            >> "$out"
    done
else
    echo "[skip] T5 bedops convert2bed (convert2bed not on PATH; bioconda: bedops)" >&2
fi

# ---------- BioConvert (intentionally NOT benchmarked on T5) --------------
# BioConvert 1.2.0 has no GFF3-to-BED converter. Verified against the
# `bioconvert --help` subcommand list on 2026-05-02; available BED-targeted
# paths are vcf2bed, bigbed2bed, and wig2bed. There is no gff32bed,
# gff3:bed, or any equivalent converter. A prior version of this script
# invoked `bioconvert gff3:bed`, which exited non-zero (unknown converter)
# and produced no output - the failure was masked by the trailing `; true`
# in the per-file loop. UCSC kent (gff3ToGenePred + genePredToBed), AGAT
# (agat_convert_sp_gff2bed.pl), and bedops (convert2bed --input=gff)
# remain the comparators for T5.
#
# If a future BioConvert release adds a GFF3-to-BED converter, restore
# the block from git history (commit prior to the Phase 4 BioConvert
# audit) and re-add the ("T5", "bioconvert") entry to _TOOL_TO_DIR in
# benchmarks/check_correctness.py.

echo "[done] T5 rows written to $out" >&2
python "$repo/benchmarks/check_correctness.py" --task T5 --bench-dir "$bench_dir" --tsv "$out" 2>&1 | head -20 || echo "[warn] T5 correctness check failed (non-fatal)" >&2
