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
# DReichLab/EIG, ships via bioconda (`eigensoft`), and reads PACKEDPED
# (plink1.9 binary .bed/.bim/.fam) much more reliably than text PED at
# 1KG scale.
#
# Run 136834 history: the previous form used `plink2 --recode` (text
# PED) + `inputformat: PED`. On chr22 biallelic SNPs (~1.1 M sites x
# 2504 samples) that produces a ~5 GB .ped file, and convertf bailed
# at 0.13 s wall (exit 8) — diagnostically that is convertf's
# `fatalx` for malformed/oversized PED, almost certainly the
# combination of (a) text-PED size and (b) plink2's 0/M/F sex column
# convention not matching what convertf's PED reader expects (plink2
# emits "0" for unknown sex; convertf prefers "U" or numeric 9, and
# its strict PED parser sometimes treats "0" as a parse failure).
#
# The fix: write plink1.9 binary PACKEDPED via `plink2 --make-bed`,
# then point convertf at the .bed/.bim/.fam triplet with
# `inputformat: PACKEDPED`. This is the standard convertf input for
# large datasets, sidesteps the text-PED size and parser fragility,
# and is the path EIGENSOFT itself recommends in its README.
#
# Defensive choices retained:
#   - `--silent` keeps plink2's per-step banners off stderr.
#   - convertf has no `--version`; pull it from the active conda env.
#   - On exit 8 / convertf failure, bench_one.py's per-run stderr log
#     captures the convertf message; this lets us iterate without
#     re-running the whole 4-min plink2 stage.
if command -v convertf >/dev/null 2>&1; then
    cf_version="$(conda list eigensoft 2>/dev/null \
        | awk '$1=="eigensoft" {print "v"$2; exit}')"
    : "${cf_version:=unknown}"
    vcf_input="$(ls "$GFC_BENCH_INPUT_DIR"/$GFC_BENCH_VCF_PATTERN 2>/dev/null | head -1)"
    stem="$(basename "$vcf_input")"
    stem="${stem%.vcf.gz}"; stem="${stem%.vcf}"; stem="${stem%.bcf}"
    for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
        out_dir="$bench_dir/convertf_rep${rep}"
        mkdir -p "$out_dir"
        # Single shell command so bench_one.py times the canonical
        # pipeline (VCF -> PACKEDPED -> EIGENSTRAT), not just convertf
        # in isolation. We pre-create the par file so convertf gets a
        # complete spec on stdin-equivalent. `numchrom: 90` is the
        # EIGENSOFT convention for "accept any chromosome code"; chr22
        # is well within range, and Y1000+ scaffolds (if anyone re-uses
        # this for non-human) are accommodated up to that ceiling.
        python "$repo/benchmarks/bench_one.py" \
            --task "T1" --tool "convertf" --version "$cf_version" \
            --replicate "$rep" \
            --cmd "plink2 --vcf '$vcf_input' --threads 1 --allow-extra-chr \
                          --make-bed --out '$out_dir/$stem' --silent && \
                   printf 'genotypename:    %s.bed\nsnpname:         %s.bim\nindivname:       %s.fam\ninputformat:     PACKEDPED\noutputformat:    EIGENSTRAT\ngenotypeoutname: %s.geno\nsnpoutname:      %s.snp\nindivoutname:    %s.ind\nfamilynames:     NO\nnumchrom:        90\n' \
                          '$out_dir/$stem' '$out_dir/$stem' '$out_dir/$stem' \
                          '$out_dir/$stem' '$out_dir/$stem' '$out_dir/$stem' \
                          > '$out_dir/par.PACKEDPED.EIGENSTRAT' && \
                   convertf -p '$out_dir/par.PACKEDPED.EIGENSTRAT'" \
            --notes "via plink2 --make-bed (PACKEDPED) → convertf EIGENSTRAT" \
            >> "$out"
    done
else
    echo "[skip] T1 convertf (EIGENSOFT convertf not on PATH; install eigensoft from bioconda)" >&2
fi

# Removed in Stage 1 (bench/stage1-fairness): the bcftools+awk hand-script
# baseline previously here. Its awk encoded ref-allele dosage instead of
# variant-allele dosage, so the .geno output was inverted relative to the
# EIGENSTRAT spec — the comparison was therefore meaningless. convertf is
# the correct reference.

echo "[done] T1 rows written to $out" >&2
python "$repo/benchmarks/check_correctness.py" --task T1 --bench-dir "$bench_dir" --tsv "$out" 2>&1 | head -20 || echo "[warn] T1 correctness check failed (non-fatal)" >&2
