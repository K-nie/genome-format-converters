#!/usr/bin/env bash
# Top-level driver: run every benchmark task and collect results.
#
# Usage:
#   bash benchmarks/run_bench.sh              # run everything
#   bash benchmarks/run_bench.sh T1 T3        # run only listed tasks
#
# Each task writes its own TSV into benchmarks/results/raw/T*.tsv; the
# merged file is written by `benchmarks/analyze/collect_results.py`.

set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$(dirname "$here")"

tasks=(
    T1_vcf_to_eigenstrat
    T2_vcf_to_plink
    T3_fasta_gff_to_gbk
    T4_gff3_to_gtf
    T5_gff3_to_bed12
    T6_hmmer_tblout
    T7_orthogroups
    T8_pseudohaploid
)

if [[ $# -gt 0 ]]; then
    # User passed a selection like `T1 T3`; filter to those.
    selected=()
    for pick in "$@"; do
        for t in "${tasks[@]}"; do
            if [[ "$t" == "$pick"* ]]; then
                selected+=("$t")
            fi
        done
    done
    tasks=("${selected[@]}")
fi

mkdir -p benchmarks/results/raw benchmarks/results/figures

for task in "${tasks[@]}"; do
    script="benchmarks/tasks/${task}.sh"
    if [[ ! -x "$script" ]]; then
        echo "[skip] $task (no executable script at $script)" >&2
        continue
    fi
    echo "[run ] $task" >&2
    bash "$script" || echo "[warn] $task exited non-zero" >&2
done

python benchmarks/analyze/collect_results.py
echo "[done] results merged at benchmarks/results/raw/all.tsv" >&2
