#!/usr/bin/env bash
# T1 — VCF → EIGENSTRAT triplet (.geno/.snp/.ind).
#
# Competitors referenced in docs/BENCHMARK_PLAN.md:
#   - EIGENSOFT convertf (install from source; see data/README.md)
#   - handwritten bcftools + awk + python (included inline, TODO)
#   - pileupCaller (pseudohaploid variant) — covered by T8
#
# This script runs gfc end-to-end on the bundled tiny.vcf fixture so the
# harness is immediately executable. For paper-grade numbers, point it at
# the 1000G chr22 or Peter-2018 yeast-VCF datasets downloaded by
# data/download_1kg_chr22.sh.

set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"

task="T1_vcf_to_eigenstrat"
out="$repo/benchmarks/results/raw/${task}.tsv"
mkdir -p "$(dirname "$out")"

# Smoke dataset (committed fixture). Swap for a downloaded real dataset
# when generating paper numbers.
input_dir="$repo/tests/test_data"
: "${GFC_BENCH_INPUT_DIR:=$input_dir}"
: "${GFC_BENCH_VCF_PATTERN:=tiny.vcf}"
: "${GFC_BENCH_REPLICATES:=5}"

bench_dir="/tmp/gfc_bench_${task}"
rm -rf "$bench_dir"
mkdir -p "$bench_dir"

# Header row
python "$repo/benchmarks/bench_one.py" --header > "$out"

# ---------- tool 1: gfc -----------------------------------------------------
gfc_version="$(gfc --version | cut -d' ' -f2)"
for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
    out_dir="$bench_dir/gfc_rep${rep}"
    python "$repo/benchmarks/bench_one.py" \
        --task "T1" --tool "gfc" --version "$gfc_version" --replicate "$rep" \
        --cmd "gfc vcf-to-eigenstrat --input-dir '$GFC_BENCH_INPUT_DIR' \
               --output-dir '$out_dir' --pattern '$GFC_BENCH_VCF_PATTERN' --force" \
        --notes "fixture=${GFC_BENCH_VCF_PATTERN}" \
        >> "$out"
done

# ---------- tool 2: EIGENSOFT convertf (Stage 3 reference) ------------------
# convertf is the canonical reference for VCF → EIGENSTRAT. It lives in
# DReichLab/EIG, is not on bioconda, and ships via the Dockerfile builder
# stage (see benchmarks/Dockerfile). When it is on PATH, run it through a
# par-file generated from the input VCF; plink2 stages the VCF as PED first.
if command -v convertf >/dev/null 2>&1; then
    # convertf has no `--version`; called bare it prints `fatalx:\nparameter
    # p compulsory` and exits non-zero. Older code parsed `head -1 | cut`
    # which captured `fatalx:` and leaked a newline into the TSV (run
    # 136834: row split across two physical lines). Pull the version from
    # the active conda env's eigensoft package metadata instead — single
    # line, no shell errors, robust to convertf's lack of a version flag.
    cf_version="$(conda list eigensoft 2>/dev/null \
        | awk '$1=="eigensoft" {print "v"$2; exit}')"
    : "${cf_version:=unknown}"
    vcf_input="$(ls "$GFC_BENCH_INPUT_DIR"/$GFC_BENCH_VCF_PATTERN 2>/dev/null | head -1)"
    stem="$(basename "$vcf_input")"
    stem="${stem%.vcf.gz}"; stem="${stem%.vcf}"; stem="${stem%.bcf}"
    for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
        out_dir="$bench_dir/convertf_rep${rep}"
        mkdir -p "$out_dir"
        # Single shell command so bench_one.py times the canonical pipeline,
        # not just convertf in isolation.
        python "$repo/benchmarks/bench_one.py" \
            --task "T1" --tool "convertf" --version "$cf_version" \
            --replicate "$rep" \
            --cmd "plink2 --vcf '$vcf_input' --threads 1 --allow-extra-chr \
                          --recode --out '$out_dir/$stem' --silent && \
                   printf 'genotypename:    %s.ped\nsnpname:         %s.map\nindivname:       %s.ped\noutputformat:    EIGENSTRAT\ngenotypeoutname: %s.geno\nsnpoutname:      %s.snp\nindivoutname:    %s.ind\n' \
                          '$out_dir/$stem' '$out_dir/$stem' '$out_dir/$stem' \
                          '$out_dir/$stem' '$out_dir/$stem' '$out_dir/$stem' \
                          > '$out_dir/par.PED.EIGENSTRAT' && \
                   convertf -p '$out_dir/par.PED.EIGENSTRAT'" \
            --notes "via plink2 PED → convertf EIGENSTRAT" \
            >> "$out"
    done
else
    echo "[skip] T1 convertf (EIGENSOFT convertf not on PATH; build via Dockerfile)" >&2
fi

# Removed in Stage 1 (bench/stage1-fairness): the bcftools+awk hand-script
# baseline previously here. Its awk encoded ref-allele dosage instead of
# variant-allele dosage, so the .geno output was inverted relative to the
# EIGENSTRAT spec — the comparison was therefore meaningless. convertf is
# the correct reference.

echo "[done] T1 rows written to $out" >&2
