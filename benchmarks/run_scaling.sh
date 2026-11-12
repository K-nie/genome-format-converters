#!/usr/bin/env bash
# Size-scaling series: separates start-up cost from per-record throughput.
#
#   VCF tasks (T1, T2):   first N sites of chr22, N in GFC_SCALE_SITES
#   Annotation (T4, T5):  a GFF3 set replicated k times, k in GFC_SCALE_FACTORS
#
# Each (task, size) run calls the normal task script with GFC_BENCH_INPUT_DIR
# pointed at the scaled input, then copies the task TSV (and its
# <task>.verdicts.tsv sidecar) to results/raw/scaling/<task>__<size>.tsv. Merge and plot
# with benchmarks/analyze/ as usual; the size is in the file name.
#
# Inputs:
#   GFC_SCALE_VCF        chr22 VCF (default: data/1kg/chr22.biallelic.snps.filtered.vcf.gz)
#   GFC_SCALE_GFF3_DIR   GFF3 directory to replicate (default: data/annotations/hs_chr22)
#   GFC_SCALE_SITES      default "10000 100000 1000000"
#   GFC_SCALE_FACTORS    default "1 10 100"
#   GFC_BENCH_REPLICATES default 5 (scaling rows need fewer than the main panel)
#
# Usage:  bash benchmarks/run_scaling.sh            # all
#         bash benchmarks/run_scaling.sh vcf        # VCF tasks only
#         bash benchmarks/run_scaling.sh gff3       # annotation tasks only
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(dirname "$here")"
cd "$repo"

: "${GFC_SCALE_VCF:=$here/data/1kg/chr22.biallelic.snps.filtered.vcf.gz}"
: "${GFC_SCALE_GFF3_DIR:=$here/data/annotations/hs_chr22}"
: "${GFC_SCALE_SITES:=10000 100000 1000000}"
: "${GFC_SCALE_FACTORS:=1 10 100}"
: "${GFC_BENCH_REPLICATES:=5}"
export GFC_BENCH_REPLICATES

scratch="${GFC_SCALE_SCRATCH:-${TMPDIR:-/tmp}/gfc_scaling_inputs}"
dest="$here/results/raw/scaling"
mkdir -p "$scratch" "$dest"
what="${1:-all}"

save() {  # save <task_script_name> <size_label>
    local task="$1" size="$2" raw="$here/results/raw/$1"
    cp "$raw.tsv" "$dest/${task}__${size}.tsv"
    if [[ -f "$raw.verdicts.tsv" ]]; then
        cp "$raw.verdicts.tsv" "$dest/${task}__${size}.verdicts.tsv"
    fi
    echo "[ok  ] $task @ $size -> $dest/${task}__${size}.tsv" >&2
}

if [[ "$what" == "all" || "$what" == "vcf" ]]; then
    [[ -f "$GFC_SCALE_VCF" ]] || { echo "[err ] missing $GFC_SCALE_VCF (run data/download_1kg_chr22.sh)" >&2; exit 1; }
    for n in $GFC_SCALE_SITES; do
        d="$scratch/vcf_${n}"
        if [[ ! -f "$d/chr22_${n}.vcf.gz" ]]; then
            python "$here/make_scaled_inputs.py" vcf --input "$GFC_SCALE_VCF" \
                --sites "$n" --output "$d/chr22_${n}.vcf.gz"
        fi
        for task in T1_vcf_to_eigenstrat T2_vcf_to_plink; do
            GFC_BENCH_INPUT_DIR="$d" GFC_BENCH_VCF_PATTERN="chr22_${n}.vcf.gz" \
                bash "$here/tasks/${task}.sh"
            save "$task" "sites${n}"
        done
    done
fi

if [[ "$what" == "all" || "$what" == "gff3" ]]; then
    [[ -d "$GFC_SCALE_GFF3_DIR" ]] || { echo "[err ] missing $GFC_SCALE_GFF3_DIR (run data/download_gene_annotations.sh)" >&2; exit 1; }
    label="$(basename "$GFC_SCALE_GFF3_DIR")"
    for k in $GFC_SCALE_FACTORS; do
        d="$scratch/gff3_${label}_x${k}"
        if [[ ! -d "$d" ]]; then
            python "$here/make_scaled_inputs.py" gff3 --input-dir "$GFC_SCALE_GFF3_DIR" \
                --factor "$k" --output-dir "$d"
        fi
        for task in T4_gff3_to_gtf T5_gff3_to_bed12; do
            GFC_BENCH_INPUT_DIR="$d" bash "$here/tasks/${task}.sh"
            save "$task" "${label}_x${k}"
        done
    done
fi

echo "[done] scaling TSVs in $dest" >&2
