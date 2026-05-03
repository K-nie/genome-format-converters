#!/usr/bin/env bash
# Author: Benjamin Narh-Madey
# HTCondor execute-host driver for the Phase 4 BioConvert + T4 smoke.
#
# Mirrors run_on_condor.sh's conda activation + gfc reinstall steps
# but runs ONLY benchmarks/tasks/T4_gff3_to_gtf.sh and pins
# GFC_BENCH_REPLICATES=2 (smoke tier). Used to validate the Phase 4
# BioConvert wrapper on the cluster before launching the full
# multi-comparator production run.

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"

# ---- conda activation (same robust path as run_on_condor.sh) ----------
conda_base="$(command -v conda >/dev/null && conda info --base 2>/dev/null || true)"
if [[ -z "$conda_base" && -d "$HOME/miniconda3" ]]; then
    conda_base="$HOME/miniconda3"
fi
if [[ -z "$conda_base" && -d "$HOME/miniforge3" ]]; then
    conda_base="$HOME/miniforge3"
fi
# shellcheck disable=SC1091
source "$conda_base/etc/profile.d/conda.sh"

: "${LD_LIBRARY_PATH:=}"
: "${PYTHONPATH:=}"
export LD_LIBRARY_PATH PYTHONPATH

set +u
conda activate gfc-bench
set -u

# Reinstall gfc from this checkout so the smoke runs against the
# committed source, not the conda env's pinned 0.1.5.
pip install -e . --no-deps --quiet --force-reinstall 2>&1 | tail -3 || \
    echo "[warn] pip install -e . failed; running env's pinned gfc" >&2

echo "[info] PATH=$PATH"
echo "[info] python=$(command -v python)"
echo "[info] gfc=$(command -v gfc)"
echo "[info] bioconvert=$(command -v bioconvert || echo NOT_INSTALLED)"
gfc --version

# Smoke = n=2 replicates per tool.
export GFC_BENCH_REPLICATES=2

mkdir -p benchmarks/results/raw
echo "[info] running ONLY T4 with n=2 replicates" >&2
bash benchmarks/tasks/T4_gff3_to_gtf.sh

echo "[done] smoke finished. T4 TSV:"
column -t -s $'\t' benchmarks/results/raw/T4_gff3_to_gtf.tsv | tail -25
