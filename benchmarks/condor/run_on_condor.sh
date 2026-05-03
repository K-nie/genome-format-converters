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

# Initialise vars that some conda activate.d hooks (notably aster, mkl,
# certain bioconda packages) extend without first checking they are set.
# Without this, `set -u` blows up at activation with messages like
# `aster_activate.sh: line 1: LD_LIBRARY_PATH: unbound variable` and the
# whole job exits in 26 seconds before any task can run (run 136836).
: "${LD_LIBRARY_PATH:=}"
: "${PYTHONPATH:=}"
export LD_LIBRARY_PATH PYTHONPATH

# Bracket `conda activate` with `set +u` because conda's hook ecosystem
# is allowed to be loose about defined-ness even after the init above —
# safer to drop nounset for the duration of activation than to chase
# every package's hook.
set +u
conda activate gfc-bench
set -u

# Install gfc from THIS checkout in editable mode so the cluster runs the
# converters as currently committed on bench/stage1-fairness. Without this
# the env's pip-installed gfc 0.1.5 (frozen at env-creation time) wins on
# PATH and every src/ optimisation since 0.1.5 is silently bypassed —
# exactly what made T3/T4/T5/T7 wall numbers stuck across runs 136838 and
# 136846 even after commits 7c0c71d, 1d1fd14, ef14752, 8494f14 landed.
# `--no-deps` keeps us from re-resolving the conda env's already-installed
# pysam / cyvcf2 / biopython etc.
echo "[info] reinstalling gfc from local checkout (editable, no deps)" >&2
pip install -e . --no-deps --quiet --force-reinstall 2>&1 | tail -3 || \
    echo "[warn] pip install -e . failed; cluster will run conda's pinned gfc" >&2

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
    agat_convert_sp_gff2gtf.pl --help 2>&1 | head -1 || echo "AGAT not installed"
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

# 10 replicates by default (Stage 1 audit: n=3 was too few for SD-based
# claims given the ~1% CV on T1/T2). Set GFC_BENCH_REPLICATES in env to
# override. Tasks > 30s use n=10; short tasks (T3/T4/T5) absorb startup
# jitter at higher n but n=10 is the floor.
export GFC_BENCH_REPLICATES="${GFC_BENCH_REPLICATES:-10}"
echo "[info] replicates per tool per task: $GFC_BENCH_REPLICATES"

# ---- run all tasks ----------------------------------------------------
bash benchmarks/run_bench.sh

# ---- render figures ---------------------------------------------------
python benchmarks/analyze/plot_wall_time.py  || echo "[warn] wall-time figure failed"
python benchmarks/analyze/plot_peak_memory.py || echo "[warn] memory figure failed"
python benchmarks/analyze/plot_feature_matrix.py || echo "[warn] feature matrix failed"

# ---- package results for transfer back to submit host -----------------
timestamp="$(date +%Y%m%d-%H%M%S)"
archive="benchmarks/results/results.$timestamp.tar.gz"
tar czf "$archive" -C benchmarks results/
echo "[done] benchmark run packaged at $archive"
