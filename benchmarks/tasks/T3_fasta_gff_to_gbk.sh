#!/usr/bin/env bash
# T3 — Paired FASTA + GFF3 → GenBank.
# Phase 4 primary competitors: NCBI table2asn (Sayers 2023), EMBLmyGFF3
# (Norling 2018), EMBOSS seqret (Rice 2000). The handwritten Biopython
# baseline (t3_biopython.py) is retained as a labelled naive baseline
# but no longer carries primary-comparator weight — see
# docs/LITERATURE_COMPARATORS_2026-05-03.md.
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

# ---------- table2asn (NCBI; PRIMARY) --------------------------------------
# table2asn is the official NCBI GenBank-submission tool, mandatory for
# every GenBank deposit since 2024-06-01 (replacement for tbl2asn). The
# closest thing to a primary citation is Sayers 2023 (NAR D141) and the
# NCBI tool homepage. table2asn reads paired FASTA + GFF3 and emits
# .sqn (binary ASN.1) and .gbf (GenBank flat-file). For our comparison
# the .gbf is the apples-to-apples output. Flags: `-M n` selects the
# normal-genome workflow; `-J` enables FASTA validation (harmless).
if command -v table2asn >/dev/null 2>&1; then
    t2a_version="$(table2asn -version 2>&1 | head -1 | tr -s ' ' | awk '{print $NF}' | tr -d '|')"
    : "${t2a_version:=unknown}"
    for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
        out_dir="$bench_dir/table2asn_rep${rep}"
        mkdir -p "$out_dir"
        # table2asn reads -i <fasta> and looks for a paired .gff next
        # to it. Stage symlinks per-file in $out_dir so its outputs
        # land beside its inputs. Trailing `; true` lets one bad
        # input not void the whole rep.
        python "$repo/benchmarks/bench_one.py" \
            --task "T3" --tool "table2asn" --version "$t2a_version" \
            --replicate "$rep" \
            --cmd "for fa in '$GFC_BENCH_INPUT_DIR'/*.fasta; do \
                      stem=\$(basename \"\$fa\" .fasta); \
                      gff='$GFC_BENCH_INPUT_DIR'/\$stem.gff3; \
                      [[ -f \"\$gff\" ]] || continue; \
                      ln -sf \"\$fa\" '$out_dir/'\$stem.fasta; \
                      ln -sf \"\$gff\" '$out_dir/'\$stem.gff; \
                      table2asn -i '$out_dir/'\$stem.fasta -f '$out_dir/'\$stem.gff -M n -J -outdir '$out_dir' 2>/dev/null || true; \
                   done; true" \
            --notes "NCBI table2asn -M n -J (Sayers 2023); $notes_default" \
            >> "$out"
    done
else
    echo "[skip] T3 table2asn (download from https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/)" >&2
fi

# ---------- EMBLmyGFF3 (NBIS; PRIMARY) -------------------------------------
# EMBLmyGFF3 (Norling 2018, BMC Research Notes 11:584) is a peer-reviewed
# converter validated by ENA. Outputs EMBL (not GenBank), but the
# format-conversion problem is identical and the same code path
# generates the feature tables that GenBank uses. Citable replacement
# for the demoted handwritten py-ref.
if command -v EMBLmyGFF3 >/dev/null 2>&1; then
    embl_version="$(EMBLmyGFF3 --version 2>&1 | head -1 | awk '{print $NF}' | tr -d '|')"
    : "${embl_version:=unknown}"
    for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
        out_dir="$bench_dir/emblmygff3_rep${rep}"
        mkdir -p "$out_dir"
        # EMBLmyGFF3 wants gff + fasta as separate args plus required
        # ENA project metadata flags. Placeholder values for the
        # metadata since we are timing the conversion, not staging a
        # real submission.
        python "$repo/benchmarks/bench_one.py" \
            --task "T3" --tool "EMBLmyGFF3" --version "$embl_version" \
            --replicate "$rep" \
            --cmd "for fa in '$GFC_BENCH_INPUT_DIR'/*.fasta; do \
                      stem=\$(basename \"\$fa\" .fasta); \
                      gff='$GFC_BENCH_INPUT_DIR'/\$stem.gff3; \
                      [[ -f \"\$gff\" ]] || continue; \
                      EMBLmyGFF3 \"\$gff\" \"\$fa\" \
                          --topology linear --molecule_type 'genomic DNA' \
                          --transl_table 1 --species 'Y1000+ slice' \
                          --locus_tag GFCBENCH --project_id PRJ00000 \
                          --output '$out_dir/'\$stem.embl 2>/dev/null || true; \
                   done; true" \
            --notes "EMBLmyGFF3 (Norling 2018); EMBL output, structurally equivalent to GenBank; $notes_default" \
            >> "$out"
    done
else
    echo "[skip] T3 EMBLmyGFF3 (pip install EMBLmyGFF3)" >&2
fi

# ---------- handwritten biopython baseline (DEMOTED) -----------------------
# Per docs/LITERATURE_COMPARATORS_2026-05-03.md: t3_biopython.py is
# uncitable (BCBio.GFF has no peer-reviewed publication, only the
# chapmanb/bcbb GitHub repo). Phase 4 reframes it as a labelled
# "naive in-script baseline" while table2asn becomes the primary
# citable comparator above. The block is retained so reviewers can
# see how a one-hour Biopython script compares to gfc and to the
# published tools — that is itself instructive — but it no longer
# carries primary-comparator weight.
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
            --notes "naive Biopython baseline (unpublished; pedagogical); $notes_default" \
            >> "$out"
    done
fi

# Skipped competitor: gff3toembl. It's a Python-2-era tool not on bioconda
# and pip-install requires legacy deps. Out of scope for Stage 1.
echo "[done] T3 rows written to $out" >&2
python "$repo/benchmarks/check_correctness.py" --task T3 --bench-dir "$bench_dir" --tsv "$out" 2>&1 | head -20 || echo "[warn] T3 correctness check failed (non-fatal)" >&2
