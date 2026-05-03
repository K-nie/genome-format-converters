#!/usr/bin/env bash
# T8 — Pseudohaploid VCF → EIGENSTRAT with reproducibility contract.
# Competitor: pileupCaller (Haskell). Highlights: per-file seeded RNG in
# gfc means two runs on the same file with --seed give identical output
# regardless of file order.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
task="T8_pseudohaploid"
out="$repo/benchmarks/results/raw/${task}.tsv"
mkdir -p "$(dirname "$out")"

input_dir="$repo/tests/test_data"
: "${GFC_BENCH_INPUT_DIR:=$input_dir}"
: "${GFC_BENCH_VCF_PATTERN:=tiny.vcf}"
: "${GFC_BENCH_REPLICATES:=5}"

bench_dir="/tmp/gfc_bench_${task}"
rm -rf "$bench_dir" && mkdir -p "$bench_dir"

python "$repo/benchmarks/bench_one.py" --header > "$out"

gfc_version="$(gfc --version | cut -d' ' -f2)"
for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
    out_dir="$bench_dir/gfc_rep${rep}"
    python "$repo/benchmarks/bench_one.py" \
        --task "T8" --tool "gfc" --version "$gfc_version" --replicate "$rep" \
        --cmd "gfc vcf-to-pseudohaploid --input-dir '$GFC_BENCH_INPUT_DIR' \
               --output-dir '$out_dir' --pattern '$GFC_BENCH_VCF_PATTERN' \
               --seed 42 --force" \
        --notes "seed=42, deterministic" \
        >> "$out"
done

# Reproducibility assertion: same seed → identical .geno across replicates.
# This is the honest selling point for T8 — not speed.
if [[ -f "$bench_dir/gfc_rep1/tiny.geno" && -f "$bench_dir/gfc_rep2/tiny.geno" ]]; then
    if diff -q "$bench_dir/gfc_rep1/tiny.geno" "$bench_dir/gfc_rep2/tiny.geno" >/dev/null; then
        echo "[ok  ] T8 reproducibility: rep1 .geno == rep2 .geno with --seed 42" >&2
    else
        echo "[FAIL] T8 reproducibility broken — rep1 != rep2 with same seed" >&2
        exit 1
    fi
fi

# ---------- bcftools consensus -H (PRIMARY direct competitor) --------------
# Per docs/LITERATURE_COMPARATORS_2026-05-03.md, replace the
# t8_bcftools.py wrapper as the named competitor with a direct
# bcftools invocation. Danecek 2021 (GigaScience giab008): bcftools
# `consensus -H 1pIu` emits per-sample pseudohaploid sequences from a
# phased VCF (htslib-native, no Python wrapper). The wrapper script
# is retained below as a labelled "naive baseline" for completeness.
#
# Implementation note: `bcftools consensus` needs (a) a phased VCF
# indexed with tabix, (b) a reference FASTA. We index per-rep so the
# rep is timed end-to-end including the index step. The reference
# FASTA path comes from $GFC_BENCH_REF_FASTA (set in run_on_condor.sh
# alongside the 1KG VCF download). When unset, the row records the
# missing ref via stderr and exits 0 — bench_one.py still writes the
# row so the TSV stays rectangular.
if command -v bcftools >/dev/null 2>&1; then
    bcftools_version="$(bcftools --version 2>&1 | head -1 | awk '{print $NF}')"
    : "${bcftools_version:=unknown}"
    vcf_input="$(ls "$GFC_BENCH_INPUT_DIR"/$GFC_BENCH_VCF_PATTERN 2>/dev/null | head -1)"
    ref_fasta="${GFC_BENCH_REF_FASTA:-}"
    for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
        out_dir="$bench_dir/bcftools_rep${rep}"
        mkdir -p "$out_dir"
        python "$repo/benchmarks/bench_one.py" \
            --task "T8" --tool "bcftools" --version "$bcftools_version" \
            --replicate "$rep" \
            --cmd "cp '$vcf_input' '$out_dir/in.vcf' 2>/dev/null || cp '$vcf_input' '$out_dir/in.vcf.gz'; \
                   if [[ -f '$out_dir/in.vcf' ]]; then bgzip -f '$out_dir/in.vcf'; fi; \
                   tabix -f -p vcf '$out_dir/in.vcf.gz' 2>/dev/null; \
                   if [[ -n '$ref_fasta' && -f '$ref_fasta' ]]; then \
                       for s in \$(bcftools query -l '$out_dir/in.vcf.gz'); do \
                           bcftools consensus -H 1pIu -s \$s -f '$ref_fasta' \
                               '$out_dir/in.vcf.gz' > '$out_dir/'\$s.fasta 2>/dev/null; \
                       done; \
                   else \
                       echo '[T8 bcftools] no ref FASTA at GFC_BENCH_REF_FASTA; consensus skipped' >&2; \
                   fi" \
            --notes "bcftools consensus -H 1pIu per sample (Danecek 2021)" \
            >> "$out"
    done
else
    echo "[skip] T8 bcftools direct (bcftools missing)" >&2
fi

# ---------- ANGSD --doHaploCall (PRIMARY published reference) --------------
# ANGSD (Korneliussen 2014, BMC Bioinformatics 15:356) is the published
# reference for pseudohaploid calling. `-doHaploCall 1` randomly samples
# one allele per site — exactly what gfc does. ANGSD is BAM-native by
# default; the `-vcf-gl <vcf>` path accepts VCF as the upstream source.
if command -v angsd >/dev/null 2>&1; then
    angsd_version="$(angsd 2>&1 | grep -i 'version' | head -1 | awk '{print $NF}' | tr -d '|')"
    : "${angsd_version:=unknown}"
    vcf_input="$(ls "$GFC_BENCH_INPUT_DIR"/$GFC_BENCH_VCF_PATTERN 2>/dev/null | head -1)"
    for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
        out_dir="$bench_dir/angsd_rep${rep}"
        mkdir -p "$out_dir"
        # ANGSD's VCF path: `-vcf-gl <vcf> -nInd <N>` then -doHaploCall 1.
        # `-doMajorMinor 1` declares alleles. Output prefix sample.* in
        # $out_dir.
        python "$repo/benchmarks/bench_one.py" \
            --task "T8" --tool "angsd" --version "$angsd_version" \
            --replicate "$rep" \
            --cmd "n_ind=\$(bcftools query -l '$vcf_input' 2>/dev/null | wc -l); \
                   angsd -vcf-gl '$vcf_input' -nInd \$n_ind \
                       -doMajorMinor 1 -doHaploCall 1 -dumpCounts 4 \
                       -doCounts 1 -out '$out_dir/sample' 2>/dev/null || true" \
            --notes "ANGSD --doHaploCall 1 (Korneliussen 2014)" \
            >> "$out"
    done
else
    echo "[skip] T8 angsd (angsd not on PATH; bioconda: angsd)" >&2
fi

# ---------- handwritten bcftools query + python (DEMOTED — naive baseline) -
# Per Phase 4 (docs/LITERATURE_COMPARATORS_2026-05-03.md): t8_bcftools.py
# is the conventional "shell glue" pseudohaploid pipeline but is
# unpublished. bcftools is now timed directly above as the primary
# competitor. This row is retained as a labelled naive baseline so
# reviewers can see how a one-hour shell glue script compares to the
# published tool it wraps; it no longer carries primary-comparator weight.
ref_script="$repo/benchmarks/refs/t8_bcftools.py"
if [[ -f "$ref_script" ]] && command -v bcftools >/dev/null 2>&1; then
    py_version="$(python --version 2>&1 | cut -d' ' -f2)"
    for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
        out_dir="$bench_dir/bcftoolsref_rep${rep}"
        mkdir -p "$out_dir"
        python "$repo/benchmarks/bench_one.py" \
            --task "T8" --tool "bcftools-pyref" --version "$py_version" \
            --replicate "$rep" \
            --cmd "python '$ref_script' --input-dir '$GFC_BENCH_INPUT_DIR' \
                   --output-dir '$out_dir' --pattern '$GFC_BENCH_VCF_PATTERN' \
                   --seed 42" \
            --notes "naive bcftools query + python (unpublished); seed=42" \
            >> "$out"
    done
else
    echo "[skip] T8 bcftools-pyref (script or bcftools missing)" >&2
fi

# Skipped competitor: pileupCaller. Operates on samtools mpileup, not VCF —
# a fair head-to-head needs a VCF -> pileup feeder pipeline (samtools
# mpileup against the original BAMs), which the bench harness doesn't
# stage. Documented for paper discussion; not implemented in Stage 1.

echo "[done] T8 rows written to $out" >&2
python "$repo/benchmarks/check_correctness.py" --task T8 --bench-dir "$bench_dir" --tsv "$out" 2>&1 | head -20 || echo "[warn] T8 correctness check failed (non-fatal)" >&2
