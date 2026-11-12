#!/usr/bin/env bash
# Install footprint: package count and disk size of clean environments.
#
# Environments (one fresh prefix each, files copied rather than hard-linked
# so `du` counts them):
#   gfc          python + gfc from this checkout (pip), its declared deps only
#   stack        the domain tools a cross-domain pipeline otherwise installs:
#                plink2, bcftools, EIGENSOFT (convertf), AGAT, gffread, EMBOSS
#   bioconvert   BioConvert from bioconda
#   biopython.convert  python + biopython.convert (pip)
#
# Requires micromamba (or set MAMBA=conda/mamba). Writes
# results/raw/install_footprint.tsv with columns:
#   env, n_packages, disk_mb, python, spec
#
# Usage: bash benchmarks/install_footprint.sh [env ...]
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(dirname "$here")"
MAMBA="${MAMBA:-micromamba}"
PY="${GFC_FOOTPRINT_PYTHON:-3.11}"
root="${GFC_FOOTPRINT_ROOT:-${TMPDIR:-/tmp}/gfc_footprint}"
out="${GFC_BENCH_RAW_DIR:-$here/results/raw}/install_footprint.tsv"
mkdir -p "$root" "$(dirname "$out")"

declare -A SPEC=(
    [gfc]="python=$PY pip"
    [stack]="python=$PY plink2 bcftools eigensoft agat gffread emboss"
    [bioconvert]="python=$PY bioconvert"
    [biopython.convert]="python=$PY pip"
)
declare -A PIP=(
    [gfc]="$repo"
    [biopython.convert]="biopython.convert"
)
envs=("$@")
[[ ${#envs[@]} -eq 0 ]] && envs=(gfc stack bioconvert biopython.convert)

[[ -f "$out" ]] || printf "env\tn_packages\tdisk_mb\tpython\tspec\n" > "$out"
for e in "${envs[@]}"; do
    prefix="$root/$e"
    rm -rf "$prefix"
    echo "[env ] $e: ${SPEC[$e]}" >&2
    copy_flag="--always-copy"
    [[ "$(basename "$MAMBA")" == "conda" ]] && copy_flag="--copy"
    "$MAMBA" create -y -q "$copy_flag" -p "$prefix" -c conda-forge -c bioconda \
        ${SPEC[$e]} >/dev/null
    if [[ -n "${PIP[$e]:-}" ]]; then
        "$prefix/bin/python" -m pip install -q --no-cache-dir "${PIP[$e]}"
    fi
    # Conda packages plus pip-only packages (pip list includes conda's
    # python packages too, so count pip entries not installed by conda).
    n_conda=$("$MAMBA" list -p "$prefix" --json | python3 -c "import json,sys; print(len(json.load(sys.stdin)))")
    n_pip=$("$prefix/bin/python" -m pip list --format=json 2>/dev/null | python3 -c "
import json, sys, pathlib
pkgs = json.load(sys.stdin)
meta = pathlib.Path('$prefix/conda-meta')
conda = {p.name.rsplit('-', 2)[0].lower().replace('_', '-') for p in meta.glob('*.json')}
print(sum(1 for p in pkgs if p['name'].lower().replace('_', '-') not in conda))")
    disk=$(du -sm "$prefix" | cut -f1)
    pyver=$("$prefix/bin/python" --version 2>&1 | cut -d' ' -f2)
    printf "%s\t%s\t%s\t%s\t%s\n" "$e" "$((n_conda + n_pip))" "$disk" "$pyver" \
        "${SPEC[$e]}${PIP[$e]:+ + pip ${PIP[$e]##*/}}" >> "$out"
    echo "[ok  ] $e: $((n_conda + n_pip)) packages, ${disk} MB" >&2
    [[ "${GFC_FOOTPRINT_KEEP:-0}" == "1" ]] || rm -rf "$prefix"
done
echo "[done] $out" >&2
