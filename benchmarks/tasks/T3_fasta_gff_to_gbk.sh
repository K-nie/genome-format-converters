#!/usr/bin/env bash
# T3 — Paired FASTA + GFF3 → GenBank.
# Competitors: EMBOSS seqret, handwritten Biopython, gff3toembl.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
task="T3_fasta_gff_to_gbk"
out="${GFC_BENCH_RAW_DIR:-$repo/benchmarks/results/raw}/${task}.tsv"
mkdir -p "$(dirname "$out")"

# T3 needs a genome FASTA whose sequence IDs are the GFF3 seqids. The
# Y1000+ slice does not qualify: its FASTA holds the tRNA sequences while its
# GFF3 places them on contigs that are not in the FASTA, so every tool was
# asked to annotate sequences it did not have (13 Nov rerun: py-ref wrote 358
# GenBank records with no sequence). Use S. cerevisiae R64 (Ensembl genome
# FASTA + GFF3, from data/download_gene_annotations.sh) and fall back to the
# tiny fixture. GFC_BENCH_T3_DIR overrides; GFC_BENCH_INPUT_DIR is ignored
# because run_paper_rerun.sh points it at the chr22 VCF.
sc_dir="$repo/benchmarks/data/annotations/sc_r64"
if [[ -n "${GFC_BENCH_T3_DIR:-}" ]]; then
    GFC_BENCH_INPUT_DIR="$GFC_BENCH_T3_DIR"
    notes_default="$(basename "$GFC_BENCH_T3_DIR")"
elif compgen -G "$sc_dir/*.fasta" >/dev/null 2>&1; then
    GFC_BENCH_INPUT_DIR="$sc_dir"
    notes_default="S. cerevisiae R64-1-1 genome + Ensembl GFF3"
else
    GFC_BENCH_INPUT_DIR="$repo/tests/test_data"
    notes_default="tiny fixture (no sc_r64 genome present)"
fi
: "${GFC_BENCH_REPLICATES:=5}"

bench_dir="${GFC_BENCH_TMP:-/tmp}/gfc_bench_${task}"
rm -rf "$bench_dir" && mkdir -p "$bench_dir"

python "$repo/benchmarks/bench_one.py" --header > "$out"

gfc_version="$(gfc --version | cut -d' ' -f2)"
for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
    out_dir="$bench_dir/gfc_rep${rep}"
    python "$repo/benchmarks/bench_one.py" \
        --task "T3" --tool "gfc" --version "$gfc_version" --replicate "$rep" \
        --cmd "gfc fasta-gff-to-gbk --input-dir '$GFC_BENCH_INPUT_DIR' \
               --output-dir '$out_dir' --force" \
        --notes "$notes_default" \
        >> "$out"
done

# ---------- EMBOSS seqret ---------------------------------------------------
# seqret reads paired FASTA + feature table and emits GenBank. Per-file loop
# so that one malformed input doesn't abort the whole task; trailing `; true`
# guarantees a clean exit code so bench_one.py records a row regardless.
if command -v seqret >/dev/null 2>&1; then
    seqret_version="$(seqret -version 2>&1 | head -1 | tr -s ' ' | cut -d' ' -f2 | tr -d '|')"
    : "${seqret_version:=unknown}"
    for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
        out_dir="$bench_dir/seqret_rep${rep}"
        mkdir -p "$out_dir"
        python "$repo/benchmarks/bench_one.py" \
            --task "T3" --tool "EMBOSS-seqret" --version "$seqret_version" \
            --replicate "$rep" \
            --cmd "for fa in '$GFC_BENCH_INPUT_DIR'/*.fasta; do \
                      stem=\$(basename \"\$fa\" .fasta); \
                      gff='$GFC_BENCH_INPUT_DIR'/\$stem.gff3; \
                      [[ -f \"\$gff\" ]] && seqret -sequence \"\$fa\" -feature -fformat gff -ufo \"\$gff\" -osformat genbank -outseq '$out_dir/'\$stem.gb -auto; \
                   done; true" \
            --notes "$notes_default" \
            >> "$out"
    done
else
    echo "[skip] T3 EMBOSS seqret (seqret not on PATH)" >&2
fi

# ---------- handwritten biopython reference --------------------------------
# benchmarks/refs/t3_biopython.py is the "what a careful bioinformatician
# would write in an hour" baseline using SeqIO + bcbio-gff.
ref_script="$repo/benchmarks/refs/t3_biopython.py"
if [[ -f "$ref_script" ]]; then
    py_version="$(python --version 2>&1 | cut -d' ' -f2)"
    for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
        out_dir="$bench_dir/pyref_rep${rep}"
        mkdir -p "$out_dir"
        python "$repo/benchmarks/bench_one.py" \
            --task "T3" --tool "py-ref" --version "$py_version" \
            --replicate "$rep" \
            --cmd "python '$ref_script' --input-dir '$GFC_BENCH_INPUT_DIR' --output-dir '$out_dir'" \
            --notes "handwritten Biopython baseline; $notes_default" \
            >> "$out"
    done
fi

# Skipped competitor: gff3toembl. It's a Python-2-era tool not on bioconda
# and pip-install requires legacy deps. Out of scope for Stage 1.
echo "[done] T3 rows written to $out" >&2
# The correctness check must not hold the job: stop it after
# GFC_BENCH_CORR_TIMEOUT seconds (default 1 h) where `timeout` exists.
corr_to=()
command -v timeout >/dev/null && corr_to=(timeout -k 60 "${GFC_BENCH_CORR_TIMEOUT:-3600}")
${corr_to[@]+"${corr_to[@]}"} python "$repo/benchmarks/check_correctness.py" --task T3 --bench-dir "$bench_dir" --tsv "$out" 2>&1 | head -20 || echo "[warn] T3 correctness check failed (non-fatal)" >&2
