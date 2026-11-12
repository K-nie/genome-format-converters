#!/usr/bin/env bash
# Execute-host driver for the paper rerun. One HTCondor job per phase
# (see paper_rerun.sub); the phase is $1:
#
#   main         full 8-task panel, n=10, four-state verdicts, T8 reproducibility
#   annotations  T4/T5 on real Ensembl gene annotations (S. cerevisiae R64,
#                GRCh38 chr22), n=10
#   scaling      T1/T2 at 10k/100k/1M sites, T4/T5 on GRCh38 chr22 at x1/x5/x25
#   cold         the main panel again with input files evicted from the page
#                cache before every replicate (Supplementary Table S2), n=5
#   import       CLI start-up time and RSS per subcommand (bench_import.py)
#   footprint    install footprint of gfc, the multi-tool stack, BioConvert
#                and biopython.convert in fresh conda environments
#
# Results land in benchmarks/results/<phase>/ and are tarred to
# benchmarks/results/<phase>.<timestamp>.tar.gz.
set -euo pipefail

phase="${1:?usage: run_paper_rerun.sh <main|annotations|scaling|cold|import|footprint>}"
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"
bench="$repo_root/benchmarks"
dest="$bench/results/$phase"
# Per-phase raw TSVs, scratch and command logs, so phases can run as
# parallel HTCondor jobs without overwriting each other.
export GFC_BENCH_RAW_DIR="${GFC_BENCH_RAW_DIR:-$bench/results/raw_$phase}"
export GFC_BENCH_TMP="${GFC_BENCH_TMP:-${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/gfc_$phase}"
export GFC_BENCH_LOG_DIR="${GFC_BENCH_LOG_DIR:-$GFC_BENCH_TMP/cmd_logs}"
raw="$GFC_BENCH_RAW_DIR"

# ---- conda activation (same as run_on_condor.sh) ----------------------
conda_base="$(command -v conda >/dev/null && conda info --base 2>/dev/null || true)"
[[ -z "$conda_base" && -d "$HOME/miniconda3" ]] && conda_base="$HOME/miniconda3"
[[ -z "$conda_base" && -d "$HOME/miniforge3" ]] && conda_base="$HOME/miniforge3"
# shellcheck disable=SC1091
source "$conda_base/etc/profile.d/conda.sh"
: "${LD_LIBRARY_PATH:=}"
: "${PYTHONPATH:=}"
export LD_LIBRARY_PATH PYTHONPATH
set +u
conda activate gfc-bench
set -u
# bootstrap.sh installs gfc from this checkout in editable mode. Reinstall
# only if the env's gfc points somewhere else (parallel jobs would race on
# a shared reinstall).
gfc_src="$(python -c 'import genome_format_converters as g; print(g.__file__)' 2>/dev/null || true)"
if [[ "$gfc_src" != "$repo_root"/* ]]; then
    echo "[info] gfc in env is not this checkout ($gfc_src); reinstalling" >&2
    pip install -e . --no-deps --quiet --force-reinstall 2>&1 | tail -3
fi
echo "[info] phase=$phase gfc=$(gfc --version) python=$(command -v python)" >&2

mkdir -p "$raw" "$dest" "$GFC_BENCH_TMP"
{
    date -u +"%Y-%m-%dT%H:%M:%SZ"; hostname; uname -a
    command -v lscpu >/dev/null && lscpu | head -25
    command -v free >/dev/null && free -h
    git -C "$repo_root" rev-parse HEAD 2>/dev/null || true
    gfc --version; python --version
    for t in bcftools plink2 plink gffread; do "$t" --version 2>&1 | head -1 || true; done
    conda list --export
} > "$dest/hardware.log" 2>&1

# ---- data (idempotent) -----------------------------------------------
bash "$bench/data/download_1kg_chr22.sh" || echo "[warn] 1kg download failed" >&2
bash "$bench/data/download_gene_annotations.sh" || echo "[warn] annotation download failed" >&2
vcf_dir="$bench/data/1kg"
vcf_name="chr22.biallelic.snps.filtered.vcf.gz"

use_chr22() {
    if [[ -f "$vcf_dir/$vcf_name" ]]; then
        export GFC_BENCH_INPUT_DIR="$vcf_dir" GFC_BENCH_VCF_PATTERN="$vcf_name"
    else
        echo "[err ] $vcf_dir/$vcf_name missing" >&2; exit 1
    fi
}

save_raw() {  # copy every per-task TSV (+ sidecars) and all.tsv into $dest
    cp "$raw"/T*.tsv "$dest"/ 2>/dev/null || true
    [[ -f "$raw/all.tsv" ]] && cp "$raw/all.tsv" "$dest/"
}

case "$phase" in
    main)
        use_chr22
        export GFC_BENCH_REPLICATES="${GFC_BENCH_REPLICATES:-10}"
        # T4/T5 read GFC_BENCH_INPUT_DIR too, so run the VCF tasks with chr22
        # and the annotation tasks on their default Y1000+ slice.
        bash "$bench/run_bench.sh" T1 T2 T8
        (unset GFC_BENCH_INPUT_DIR GFC_BENCH_VCF_PATTERN; bash "$bench/run_bench.sh" T3 T4 T5 T6 T7)
        python "$bench/analyze/collect_results.py"
        save_raw
        ;;
    annotations)
        export GFC_BENCH_REPLICATES="${GFC_BENCH_REPLICATES:-10}"
        for d in "$bench"/data/annotations/*/; do
            label="$(basename "$d")"
            for task in T4_gff3_to_gtf T5_gff3_to_bed12; do
                GFC_BENCH_INPUT_DIR="$d" bash "$bench/tasks/${task}.sh"
                cp "$raw/${task}.tsv" "$dest/${task}__${label}.tsv"
                [[ -f "$raw/${task}.verdicts.tsv" ]] && \
                    cp "$raw/${task}.verdicts.tsv" "$dest/${task}__${label}.verdicts.tsv"
            done
        done
        ;;
    scaling)
        export GFC_BENCH_REPLICATES="${GFC_BENCH_REPLICATES:-3}"
        export GFC_SCALE_VCF="$vcf_dir/$vcf_name"
        export GFC_SCALE_GFF3_DIR="$bench/data/annotations/hs_chr22"
        export GFC_SCALE_FACTORS="${GFC_SCALE_FACTORS:-1 5 25}"
        bash "$bench/run_scaling.sh"
        cp -r "$raw/scaling/." "$dest/"
        ;;
    cold)
        use_chr22
        export GFC_BENCH_REPLICATES="${GFC_BENCH_REPLICATES:-5}"
        export GFC_BENCH_EVICT_DIRS="$vcf_dir:$bench/data/y1000plus:$bench/data/hmmer:$bench/data/orthofinder"
        bash "$bench/run_bench.sh" T1 T2 T8
        (unset GFC_BENCH_INPUT_DIR GFC_BENCH_VCF_PATTERN; bash "$bench/run_bench.sh" T3 T4 T5 T6 T7)
        python "$bench/analyze/collect_results.py"
        save_raw
        ;;
    import)
        python "$bench/bench_import.py" --trials "${GFC_IMPORT_TRIALS:-10}" \
            --out "$dest/import_times.tsv"
        ;;
    footprint)
        export MAMBA="${MAMBA:-$(command -v mamba >/dev/null && echo mamba || echo conda)}"
        GFC_FOOTPRINT_ROOT="$GFC_BENCH_TMP/footprint" bash "$bench/install_footprint.sh"
        cp "$raw/install_footprint.tsv" "$dest/"
        ;;
    *)
        echo "[err ] unknown phase $phase" >&2; exit 2 ;;
esac

# Render a preview of the paper figures and tables from every phase finished
# so far (the final set is rebuilt from all tarballs with the same script).
python "$bench/analyze/paper_figures.py" --results "$bench/results" \
    --out "$dest/paper_preview" || echo "[warn] figure rendering failed" >&2

stamp="$(date +%Y%m%d-%H%M%S)"
tar czf "$bench/results/${phase}.${stamp}.tar.gz" -C "$bench/results" "$phase"
echo "[done] $phase -> benchmarks/results/${phase}.${stamp}.tar.gz" >&2
