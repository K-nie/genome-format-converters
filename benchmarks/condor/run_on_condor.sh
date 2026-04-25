#!/usr/bin/env bash
# HTCondor execute-host driver. Run by condor_submit via bench.sub.
#
# Contract:
#   - Lives under $repo_root/benchmarks/condor/ alongside bench.sub.
#   - Assumes the gfc-bench conda env is reachable on the execute host
#     (GLBRC's scarcity/execute nodes share the submit host's home FS,
#      so the env built by bootstrap.sh is visible here).
#   - Writes per-task TSVs under $repo_root/benchmarks/results/raw/ and
#     figures under $repo_root/benchmarks/results/figures/.
#   - Tars the whole results/ dir into a timestamped archive at the end.

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"

# ---- conda activation (robust across login / non-login shells) --------
conda_base="$(command -v conda >/dev/null && conda info --base 2>/dev/null || true)"
if [[ -z "$conda_base" && -d "$HOME/miniconda3" ]]; then
    conda_base="$HOME/miniconda3"
fi
if [[ -z "$conda_base" && -d "$HOME/miniforge3" ]]; then
    conda_base="$HOME/miniforge3"
fi
# shellcheck disable=SC1091
source "$conda_base/etc/profile.d/conda.sh"
conda activate gfc-bench

echo "[info] PATH=$PATH"
echo "[info] python=$(command -v python)"
echo "[info] gfc=$(command -v gfc)"
gfc --version

# ---- hardware fingerprint for the paper Methods section ---------------
mkdir -p benchmarks/results
{
    echo "=== run start ==="
    date -u +"%Y-%m-%dT%H:%M:%SZ"
    echo
    echo "=== host ==="
    hostname
    uname -a
    echo
    echo "=== cpu ==="
    if command -v lscpu >/dev/null; then
        lscpu | head -25
    else
        sysctl -a 2>/dev/null | grep -E "machdep.cpu|hw.ncpu|hw.memsize" | head -20 || true
    fi
    echo
    echo "=== memory ==="
    if command -v free >/dev/null; then
        free -h
    fi
    echo
    echo "=== disk (benchmarks/) ==="
    df -h benchmarks/ 2>/dev/null || true
    echo
    echo "=== tool versions ==="
    gfc --version
    python --version
    bcftools --version 2>/dev/null | head -1 || echo "bcftools not installed"
    samtools --version 2>/dev/null | head -1 || echo "samtools not installed"
    plink2 --version 2>/dev/null | head -1 || echo "plink2 not installed"
    plink --version 2>/dev/null | head -1 || echo "plink 1.9 not installed"
    gffread --version 2>/dev/null || echo "gffread not installed"
    agat_sp_gff2gtf.pl --help 2>&1 | head -1 || echo "AGAT not installed"
    seqret -help 2>&1 | head -3 || echo "EMBOSS seqret not installed"
    hmmsearch -h 2>&1 | head -1 || echo "HMMER not installed"
    echo
    echo "=== conda package versions ==="
    conda list --export
} > benchmarks/results/hardware.log 2>&1

# ---- download real datasets (idempotent; skip if already present) -----
bash benchmarks/data/download_1kg_chr22.sh     || echo "[warn] 1kg download failed"
bash benchmarks/data/download_y1000plus.sh    || echo "[warn] y1000+ download failed / URLs are TODO"
bash benchmarks/data/download_pfam_scan.sh    || echo "[warn] pfam scan download failed / URLs are TODO"

# ---- point tasks at real data where available -------------------------
# Default: use the bundled tiny fixtures (harness smoke-test).
# Override: point VCF tasks at 1000G chr22 when the download succeeded.
if [[ -f benchmarks/data/1kg/chr22.biallelic.snps.filtered.vcf.gz ]]; then
    export GFC_BENCH_INPUT_DIR="$repo_root/benchmarks/data/1kg"
    export GFC_BENCH_VCF_PATTERN="chr22.biallelic.snps.filtered.vcf.gz"
    echo "[info] VCF tasks will run on 1000G chr22 biallelic SNPs"
else
    echo "[info] VCF tasks will run on bundled tiny.vcf fixture"
fi

# 5 replicates by default; set GFC_BENCH_REPLICATES in env to override.
export GFC_BENCH_REPLICATES="${GFC_BENCH_REPLICATES:-5}"
echo "[info] replicates per tool per task: $GFC_BENCH_REPLICATES"

# ---- run all tasks ----------------------------------------------------
bash benchmarks/run_bench.sh

# ---- render figures ---------------------------------------------------
python benchmarks/analyze/plot_wall_time.py  || echo "[warn] wall-time figure failed"
python benchmarks/analyze/plot_memory.py     || echo "[warn] memory figure failed"
python benchmarks/analyze/feature_matrix.py  || echo "[warn] feature matrix failed"

# ---- package results for transfer back to submit host -----------------
timestamp="$(date +%Y%m%d-%H%M%S)"
archive="benchmarks/results/results.$timestamp.tar.gz"
tar czf "$archive" -C benchmarks results/
echo "[done] benchmark run packaged at $archive"
