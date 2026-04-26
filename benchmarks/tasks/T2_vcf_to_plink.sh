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

# ---------- plink2 ----------------------------------------------------------
if command -v plink2 >/dev/null 2>&1; then
    plink2_version="$(plink2 --version 2>&1 | head -1 | tr -s ' ' | cut -d' ' -f2)"
    vcf_input="$(ls "$GFC_BENCH_INPUT_DIR"/$GFC_BENCH_VCF_PATTERN 2>/dev/null | head -1)"
    # Strip compound suffix for the --out stem; plink will write stem.bed/bim/fam.
    stem="$(basename "$vcf_input")"
    stem="${stem%.vcf.gz}"; stem="${stem%.vcf}"; stem="${stem%.bcf}"
    for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
        out_dir="$bench_dir/plink2_rep${rep}"
        mkdir -p "$out_dir"
        python "$repo/benchmarks/bench_one.py" \
            --task "T2" --tool "plink2" --version "$plink2_version" \
            --replicate "$rep" \
            --cmd "plink2 --vcf '$vcf_input' --threads 1 --allow-extra-chr --make-bed \
                   --out '$out_dir/$stem' --silent" \
            >> "$out"
    done
else
    echo "[skip] T2 plink2 (plink2 not on PATH)" >&2
fi

# ---------- plink 1.9 -------------------------------------------------------
if command -v plink >/dev/null 2>&1 && plink --version 2>&1 | grep -q 'PLINK v1'; then
    plink_version="$(plink --version 2>&1 | head -1 | tr -s ' ' | cut -d' ' -f2)"
    vcf_input="$(ls "$GFC_BENCH_INPUT_DIR"/$GFC_BENCH_VCF_PATTERN 2>/dev/null | head -1)"
    stem="$(basename "$vcf_input")"
    stem="${stem%.vcf.gz}"; stem="${stem%.vcf}"; stem="${stem%.bcf}"
    for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
        out_dir="$bench_dir/plink1_rep${rep}"
        mkdir -p "$out_dir"
        python "$repo/benchmarks/bench_one.py" \
            --task "T2" --tool "plink1.9" --version "$plink_version" \
            --replicate "$rep" \
            --cmd "plink --vcf '$vcf_input' --threads 1 --allow-extra-chr --make-bed \
                   --out '$out_dir/$stem' --silent" \
            >> "$out"
    done
else
    echo "[skip] T2 plink1.9 (plink v1 not on PATH)" >&2
fi

# TODO (correctness check, manual): byte-diff between gfc .bed and plink2 .bed.
# Both are SNP-major 2-bit packed; diffs typically come from allele-ordering
# conventions (a1/a2 swap) — plink2 sometimes emits major as a1 and minor as
# a2, which flips every 00/11 pair. Document the exact conversion rule in the
# paper rather than asserting byte-identity here.
echo "[done] T2 rows written to $out" >&2
