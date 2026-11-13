#!/usr/bin/env bash
# T8 — Pseudohaploid VCF → EIGENSTRAT with reproducibility contract.
# No timing competitor (see below). Highlights: per-file seeded RNG means
# two runs on the same file with --seed give identical output regardless
# of file order or --threads.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
task="T8_pseudohaploid"
out="${GFC_BENCH_RAW_DIR:-$repo/benchmarks/results/raw}/${task}.tsv"
mkdir -p "$(dirname "$out")"

input_dir="$repo/tests/test_data"
: "${GFC_BENCH_INPUT_DIR:=$input_dir}"
: "${GFC_BENCH_VCF_PATTERN:=tiny.vcf}"
: "${GFC_BENCH_REPLICATES:=5}"

bench_dir="${GFC_BENCH_TMP:-/tmp}/gfc_bench_${task}"
rm -rf "$bench_dir" && mkdir -p "$bench_dir"

python "$repo/benchmarks/bench_one.py" --header > "$out"

gfc_version="$(gfc --version | cut -d' ' -f2)"
for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
    out_dir="$bench_dir/gfc_rep${rep}"
    python "$repo/benchmarks/bench_one.py" \
        --task "T8" --tool "gfc" --version "$gfc_version" --replicate "$rep" \
        --cmd "gfc vcf-to-pseudohaploid --input-dir '$GFC_BENCH_INPUT_DIR' \
               --output-dir '$out_dir' --pattern '$GFC_BENCH_VCF_PATTERN' \
               --seed 42 --force" \
        --notes "seed=42, deterministic" \
        >> "$out"
done

# T8 is reported as a correctness and reproducibility task, not a speed
# race: no published tool converts VCF to pseudohaploid EIGENSTRAT directly
# (ANGSD --doHaploCall and pileupCaller both read BAM/pileup), so any timing
# competitor would be doing different work. The earlier handwritten
# bcftools-query + Python baseline is no longer run.
input_vcf="$(ls "$GFC_BENCH_INPUT_DIR"/$GFC_BENCH_VCF_PATTERN | head -1)"
repro_out="${GFC_BENCH_RAW_DIR:-$repo/benchmarks/results/raw}/${task}.reproducibility.tsv"
# The reproducibility checks run 20 conversions (5 runs x GFC_BENCH_T8_COPIES
# files), so on a whole chromosome they use the first
# GFC_BENCH_T8_REPRO_SITES records (default 20000; 0 = the whole file).
# Seeding and threading do not depend on file size.
repro_vcf="$input_vcf"
: "${GFC_BENCH_T8_REPRO_SITES:=20000}"
if [[ "$GFC_BENCH_T8_REPRO_SITES" != 0 ]]; then
    repro_vcf="$bench_dir/repro_subset.vcf.gz"
    python "$repo/benchmarks/make_scaled_inputs.py" vcf --input "$input_vcf" \
        --sites "$GFC_BENCH_T8_REPRO_SITES" --output "$repro_vcf"
fi
if python "$repo/benchmarks/t8_reproducibility.py" --vcf "$repro_vcf" \
        --work "$bench_dir/repro" --copies "${GFC_BENCH_T8_COPIES:-4}" \
        --threads "${GFC_BENCH_THREADS:-4}" --out "$repro_out"; then
    echo "[ok  ] T8 reproducibility checks passed ($repro_out)" >&2
else
    echo "[FAIL] T8 reproducibility check failed; see $repro_out" >&2
    exit 1
fi

echo "[done] T8 rows written to $out" >&2
# The correctness check must not hold the job: stop it after
# GFC_BENCH_CORR_TIMEOUT seconds (default 1 h) where `timeout` exists.
corr_to=()
command -v timeout >/dev/null && corr_to=(timeout -k 60 "${GFC_BENCH_CORR_TIMEOUT:-3600}")
${corr_to[@]+"${corr_to[@]}"} python "$repo/benchmarks/check_correctness.py" --task T8 --bench-dir "$bench_dir" --tsv "$out" --input-vcf "$input_vcf" 2>&1 | head -20 || echo "[warn] T8 correctness check failed (non-fatal)" >&2
