#!/usr/bin/env bash
# T8 — Pseudohaploid VCF → EIGENSTRAT with reproducibility contract.
# Competitor: pileupCaller (Haskell). Highlights: per-file seeded RNG in
# gfc means two runs on the same file with --seed give identical output
# regardless of file order.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
task="T8_pseudohaploid"
out="$repo/benchmarks/results/raw/${task}.tsv"
mkdir -p "$(dirname "$out")"

input_dir="$repo/tests/test_data"
: "${GFC_BENCH_INPUT_DIR:=$input_dir}"
: "${GFC_BENCH_VCF_PATTERN:=tiny.vcf}"
: "${GFC_BENCH_REPLICATES:=5}"

bench_dir="/tmp/gfc_bench_${task}"
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

# Reproducibility assertion: same seed → identical .geno across replicates.
# This is the honest selling point for T8 — not speed.
if [[ -f "$bench_dir/gfc_rep1/tiny.geno" && -f "$bench_dir/gfc_rep2/tiny.geno" ]]; then
    if diff -q "$bench_dir/gfc_rep1/tiny.geno" "$bench_dir/gfc_rep2/tiny.geno" >/dev/null; then
        echo "[ok  ] T8 reproducibility: rep1 .geno == rep2 .geno with --seed 42" >&2
    else
        echo "[FAIL] T8 reproducibility broken — rep1 != rep2 with same seed" >&2
        exit 1
    fi
fi

# ---------- handwritten bcftools query + python reference ------------------
# The conventional "shell glue" approach to pseudohaploid: stream
# bcftools query | python and randomly pick one of the two GT alleles.
# benchmarks/refs/t8_bcftools.py is that script.
ref_script="$repo/benchmarks/refs/t8_bcftools.py"
if [[ -f "$ref_script" ]] && command -v bcftools >/dev/null 2>&1; then
    py_version="$(python --version 2>&1 | cut -d' ' -f2)"
    for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
        out_dir="$bench_dir/bcftoolsref_rep${rep}"
        mkdir -p "$out_dir"
        python "$repo/benchmarks/bench_one.py" \
            --task "T8" --tool "bcftools-pyref" --version "$py_version" \
            --replicate "$rep" \
            --cmd "python '$ref_script' --input-dir '$GFC_BENCH_INPUT_DIR' \
                   --output-dir '$out_dir' --pattern '$GFC_BENCH_VCF_PATTERN' \
                   --seed 42" \
            --notes "handwritten bcftools query + python; seed=42" \
            >> "$out"
    done
else
    echo "[skip] T8 bcftools-pyref (script or bcftools missing)" >&2
fi

# Skipped competitor: pileupCaller. Operates on samtools mpileup, not VCF —
# a fair head-to-head needs a VCF -> pileup feeder pipeline (samtools
# mpileup against the original BAMs), which the bench harness doesn't
# stage. Documented for paper discussion; not implemented in Stage 1.

echo "[done] T8 rows written to $out" >&2
