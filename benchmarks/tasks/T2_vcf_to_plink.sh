#!/usr/bin/env bash
# T2 — VCF → PLINK binary (.bed/.bim/.fam).
# Competitors: plink2 --vcf, plink1.9 --vcf --make-bed.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
task="T2_vcf_to_plink"
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
        --task "T2" --tool "gfc" --version "$gfc_version" --replicate "$rep" \
        --cmd "gfc vcf-to-plink --input-dir '$GFC_BENCH_INPUT_DIR' \
               --output-dir '$out_dir' --pattern '$GFC_BENCH_VCF_PATTERN' --force" \
        >> "$out"
done

# TODO: plink2 --vcf <in> --make-bed --out <stem>
# TODO: plink1.9 equivalent
# TODO: correctness check — byte-diff gfc .bed vs plink2 .bed (both should be
# SNP-major binary format; differences will come from allele ordering
# conventions; document in the paper).
echo "[done] T2 rows written to $out" >&2
