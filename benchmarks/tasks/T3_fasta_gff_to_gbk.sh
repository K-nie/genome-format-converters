#!/usr/bin/env bash
# T3 — Paired FASTA + GFF3 → GenBank.
# Competitors: EMBOSS seqret, handwritten Biopython, gff3toembl.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
task="T3_fasta_gff_to_gbk"
out="$repo/benchmarks/results/raw/${task}.tsv"
mkdir -p "$(dirname "$out")"

# Stage 1: prefer the committed Y1000+ tRNA-scan slice over the tiny
# fixture when present, regardless of what run_on_condor.sh exported into
# GFC_BENCH_INPUT_DIR (which it sets to the 1kg VCF dir for T1/T2/T8).
y1000_dir="$repo/benchmarks/data/y1000plus"
if compgen -G "$y1000_dir/*.fasta" >/dev/null 2>&1; then
    GFC_BENCH_INPUT_DIR="$y1000_dir"
    notes_default="20-species Y1000+ tRNA-scan slice"
else
    GFC_BENCH_INPUT_DIR="${GFC_BENCH_INPUT_DIR:-$repo/tests/test_data}"
    notes_default="tiny fixture (no Y1000+ slice present)"
fi
: "${GFC_BENCH_REPLICATES:=5}"

bench_dir="/tmp/gfc_bench_${task}"
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
python "$repo/benchmarks/check_correctness.py" --task T3 --bench-dir "$bench_dir" --tsv "$out" 2>&1 | head -20 || echo "[warn] T3 correctness check failed (non-fatal)" >&2
