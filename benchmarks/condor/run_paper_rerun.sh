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
#   mutreal      mutation test on the real inputs (mutation_test_real.py)
#   footprint    install footprint of gfc, the multi-tool stack, BioConvert
#                and biopython.convert in fresh conda environments
#   extended     the 24 subcommands outside T1-T8 (E01-E24): builds the
#                inputs (benchmarks/extended/prepare_inputs.sh), then wall
#                time, peak RSS and verdicts against other tools, n=10
#
# Results land in benchmarks/results/<phase>/ and are tarred to
# benchmarks/results/<phase>.<timestamp>.tar.gz.
set -euo pipefail

phase="${1:?usage: run_paper_rerun.sh <main|annotations|scaling|cold|import|footprint|t1|verdicts|mutreal|extended>}"
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
# Run gfc from THIS checkout without touching the shared conda env: the
# env's own install may point at another checkout, and reinstalling it from
# parallel jobs fails (run 150657: OSError on .../envs/gfc-bench/bin/gfc).
# PYTHONPATH makes `import genome_format_converters` resolve to src/, and a
# per-job shim puts a matching `gfc` first on PATH.
export PYTHONPATH="$repo_root/src${PYTHONPATH:+:$PYTHONPATH}"
shim_dir="${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/gfc_shim_${phase}_$$"
mkdir -p "$shim_dir"
cat > "$shim_dir/gfc" <<SHIM
#!/usr/bin/env bash
exec "$(command -v python)" -m genome_format_converters.cli "\$@"
SHIM
chmod +x "$shim_dir/gfc"
export PATH="$shim_dir:$PATH"
gfc_src="$(python -c 'import genome_format_converters as g; print(g.__file__)')"
want="$(cd "$repo_root/src" && pwd -P)"
if [[ "$(cd "$(dirname "$gfc_src")/.." && pwd -P)" != "$want" ]]; then
    echo "[err ] gfc resolves to $gfc_src, not $want" >&2
    exit 1
fi
echo "[info] gfc from $gfc_src via $(command -v gfc)" >&2
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
    t1)
        # T1 alone (gfc, convertf, convertf-only, plink2) on one node, so the
        # plink2 --export eig competitor is timed beside the other two.
        use_chr22
        export GFC_BENCH_REPLICATES="${GFC_BENCH_REPLICATES:-10}"
        bash "$bench/run_bench.sh" T1
        save_raw
        ;;
    verdicts)
        # Correctness only (one replicate; timings not used): T3 on R64 and
        # T4/T5 on every annotation set.
        export GFC_BENCH_REPLICATES=1
        (unset GFC_BENCH_INPUT_DIR GFC_BENCH_VCF_PATTERN; bash "$bench/run_bench.sh" T3)
        save_raw
        for d in "$bench"/data/annotations/*/; do
            label="$(basename "$d")"
            for task in T4_gff3_to_gtf T5_gff3_to_bed12; do
                GFC_BENCH_INPUT_DIR="$d" bash "$bench/tasks/${task}.sh" || \
                    echo "[warn] $task on $label exited non-zero" >&2
                [[ -f "$raw/${task}.verdicts.tsv" ]] && \
                    cp "$raw/${task}.verdicts.tsv" "$dest/${task}__${label}.verdicts.tsv"
            done
        done
        ;;
    mutreal)
        # Mutation test on the real benchmark inputs: defects planted at
        # random sites, samples and features of gfc's real outputs.
        python "$bench/mutation_test_real.py" \
            --vcf "$vcf_dir/$vcf_name" --sites "${GFC_MUTREAL_SITES:-20000}" \
            --annotation-dir "$bench/data/annotations/hs_chr22" \
            --annotation-dir "$bench/data/annotations/sc_r64" \
            --annotation-dir "$bench/data/y1000plus" \
            --genome-dir "$bench/data/annotations/sc_r64" \
            --reps "${GFC_MUTREAL_REPS:-5}" --seed 1 \
            --out "$dest/mutation_real.tsv" 2> "$dest/mutation_real.log" || \
            echo "[warn] mutation_test_real reported unexpected cases (see mutation_real.log)" >&2
        tail -n 20 "$dest/mutation_real.log" >&2
        ;;
    annotations)
        export GFC_BENCH_REPLICATES="${GFC_BENCH_REPLICATES:-10}"
        for d in "$bench"/data/annotations/*/; do
            label="$(basename "$d")"
            for task in T4_gff3_to_gtf T5_gff3_to_bed12; do
                GFC_BENCH_INPUT_DIR="$d" bash "$bench/tasks/${task}.sh" || \
                    echo "[warn] $task on $label exited non-zero" >&2
                [[ -f "$raw/${task}.tsv" ]] || continue
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
    extended)
        # Extra tools live in the gfc-ext env and BioConvert in
        # gfc-bioconvert (bootstrap.sh builds both). gfc-ext goes after
        # gfc-bench on PATH so shared tools keep the gfc-bench versions.
        envs="$(conda info --base)/envs"
        [[ -d "$envs/gfc-ext/bin" ]] && export PATH="$PATH:$envs/gfc-ext/bin" || \
            echo "[warn] gfc-ext env missing: run bootstrap.sh" >&2
        if [[ -x "$envs/gfc-bioconvert/bin/bioconvert" ]]; then
            export GFC_BENCH_BIOCONVERT="$envs/gfc-bioconvert/bin/bioconvert"
        else
            echo "[warn] gfc-bioconvert env missing: BioConvert rows will be skipped" >&2
        fi
        for t in art_illumina minimap2 nucmer show-coords lastal FastTree esl-reformat \
                 diamond plink2 bcftools seqkit seqret gffread agat_convert_sp_gxf2gxf.pl bedtools; do
            command -v "$t" >/dev/null && echo "[info] $t: $(command -v "$t")" >&2 || \
                echo "[warn] $t not on PATH" >&2
        done
        GFC_BENCH_PREP_THREADS="${GFC_BENCH_PREP_THREADS:-4}" \
            bash "$bench/extended/prepare_inputs.sh"
        # GFC_EXT_TASKS="E07 E16" reruns only those tasks; give it its own
        # GFC_EXT_OUT so the full run's summaries are not overwritten.
        # shellcheck disable=SC2086
        python "$bench/extended/run_extended.py" --reps "${GFC_BENCH_REPLICATES:-10}" \
            ${GFC_EXT_TASKS:+--tasks $GFC_EXT_TASKS} \
            --out "${GFC_EXT_OUT:-$dest}" --tmp "$GFC_BENCH_TMP/extended"
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
