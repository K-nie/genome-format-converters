#!/usr/bin/env bash
# One-time setup on the HTCondor submit host (scarcity-ap-1.glbrc.org).
# Call this from the repo root, e.g.:
#
#   cd ~/gfc-bench  (or wherever rsync landed the repo)
#   bash benchmarks/condor/bootstrap.sh
#
# Effects:
#   1. Creates / updates the `gfc-bench` conda env from environment.yml.
#   2. Installs gfc itself in editable mode into that env.
#   3. Pre-creates benchmarks/condor/logs/ so HTCondor has somewhere to
#      write stdout/stderr/log files.
#   4. Prints the condor_submit one-liner for the next step.

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"

echo "[info] repo_root=$repo_root"

# ---- conda env ---------------------------------------------------------
# Use mamba if available (faster solver); fall back to conda.
solver="mamba"
if ! command -v mamba >/dev/null 2>&1; then
    solver="conda"
fi

if conda env list 2>/dev/null | awk '{print $1}' | grep -qx gfc-bench; then
    echo "[info] gfc-bench env exists; updating from environment.yml"
    "$solver" env update -n gfc-bench -f benchmarks/environment.yml
else
    echo "[info] creating gfc-bench env from environment.yml (this can take 5-15 minutes)"
    "$solver" env create -n gfc-bench -f benchmarks/environment.yml
fi

# ---- extra envs for the extended phase (E01-E24) ------------------------
# Separate envs so the pinned gfc-bench solve stays as it is. A failure
# here only affects the extended phase, so warn and carry on.
for spec in gfc-ext:benchmarks/extended/environment_extended.yml \
            gfc-bioconvert:benchmarks/extended/environment_bioconvert.yml; do
    env_name="${spec%%:*}"; env_file="${spec#*:}"
    if conda env list 2>/dev/null | awk '{print $1}' | grep -qx "$env_name"; then
        "$solver" env update -n "$env_name" -f "$env_file" || \
            echo "[warn] could not update $env_name" >&2
    else
        "$solver" env create -n "$env_name" -f "$env_file" || \
            echo "[warn] could not create $env_name (extended phase will skip its tools)" >&2
    fi
done

# ---- install gfc into the env -----------------------------------------
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gfc-bench

pip install -e "$repo_root"

echo "[info] gfc version in the env:"
gfc --version

# ---- HTCondor log dir -------------------------------------------------
mkdir -p "$repo_root/benchmarks/condor/logs"

conda deactivate

cat <<EOM

============================================================
Bootstrap complete.

Submit the benchmark job with:

    cd $repo_root/benchmarks/condor && condor_submit bench.sub

Watch its status with:

    condor_q -nobatch

When it finishes, the result tarball appears at:

    $repo_root/benchmarks/results/results.YYYYMMDD-HHMMSS.tar.gz

and the raw TSVs + figures are also in benchmarks/results/.
============================================================
EOM
