#!/usr/bin/env bash
# T7 — OrthoFinder Orthogroups.tsv → per-OG FASTA.  Competitor: handwritten python.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
task="T7_orthogroups"
out="$repo/benchmarks/results/raw/${task}.tsv"
mkdir -p "$(dirname "$out")"

# Synthesise a 2-species, 2-OG fixture so the harness is immediately runnable.
fixture_dir="${GFC_BENCH_INPUT_DIR:-$repo/benchmarks/data/orthofinder}"
mkdir -p "$fixture_dir"
if [[ ! -f "$fixture_dir/Orthogroups.tsv" ]]; then
    cat > "$fixture_dir/spA.fasta" <<'EOF'
>a1
MKTALVSLL
>a2
MRKSTTVAA
EOF
    cat > "$fixture_dir/spB.fasta" <<'EOF'
>b1
MTKQVLLAA
>b2
MQQSTALVV
EOF
    cat > "$fixture_dir/Orthogroups.tsv" <<'EOF'
Orthogroup	spA	spB
OG0000000	a1, a2	b1
OG0000001		b2
EOF
fi

: "${GFC_BENCH_REPLICATES:=5}"
bench_dir="/tmp/gfc_bench_${task}"
rm -rf "$bench_dir" && mkdir -p "$bench_dir"

python "$repo/benchmarks/bench_one.py" --header > "$out"

gfc_version="$(gfc --version | cut -d' ' -f2)"
for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
    out_dir="$bench_dir/gfc_rep${rep}"
    python "$repo/benchmarks/bench_one.py" \
        --task "T7" --tool "gfc" --version "$gfc_version" --replicate "$rep" \
        --cmd "gfc orthogroups-to-fasta --orthogroups '$fixture_dir/Orthogroups.tsv' \
               --fasta-dir '$fixture_dir' --output-dir '$out_dir' --force" \
        >> "$out"
done

# ---------- handwritten Python reference ------------------------------------
# benchmarks/refs/orthogroups_ref.py is a minimal biopython-based
# reimplementation. It's the "what a careful bioinformatician would write
# in an hour" baseline — the gfc→ref delta measures the CLI wrapper
# overhead, not algorithmic differences.
ref_script="$repo/benchmarks/refs/orthogroups_ref.py"
if [[ -x "$ref_script" || -f "$ref_script" ]]; then
    py_version="$(python --version 2>&1 | cut -d' ' -f2)"
    for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
        out_dir="$bench_dir/ref_rep${rep}"
        mkdir -p "$out_dir"
        python "$repo/benchmarks/bench_one.py" \
            --task "T7" --tool "py-ref" --version "$py_version" \
            --replicate "$rep" \
            --cmd "python '$ref_script' --orthogroups '$fixture_dir/Orthogroups.tsv' \
                   --fasta-dir '$fixture_dir' --output-dir '$out_dir'" \
            --notes "handwritten biopython baseline" \
            >> "$out"
    done
fi

# ---------- OrthoFinder cold start (PRIMARY citable comparator) ------------
# OrthoFinder (Emms & Kelly 2019, Genome Biology 20:238) is the
# canonical orthology-inference tool. Its standard run emits per-OG
# FASTAs under `Orthogroup_Sequences/` — exactly the artefact gfc
# orthogroups-to-fasta produces from a saved Orthogroups.tsv.
#
# Framing per docs/LITERATURE_COMPARATORS_2026-05-03.md:
#   gfc reads a *saved* Orthogroups.tsv (cheap; ms-scale).
#   OrthoFinder runs the full BLAST/DIAMOND + MCL pipeline cold (slow;
#   minutes to hours). The comparison is honest because it answers the
#   actual use case: "post-hoc per-OG FASTA reconstruction without
#   re-running orthology inference." Reviewers should see both axes:
#   gfc is faster *because* it skips the inference step that
#   OrthoFinder unavoidably re-runs.
#
# `-og` flag: OrthoFinder shorthand for "stop after orthogroups, don't
# run gene-tree inference." This matches the artefact gfc produces
# and trims OrthoFinder's wall-time to the minimum-fair value.
# `-S diamond` keeps the alignment step on its modern default.
# `-t 1 -a 1` pin single-thread for fairness with gfc; remove for the
# scaling figure.
if command -v orthofinder >/dev/null 2>&1; then
    of_version="$(orthofinder -h 2>&1 | grep -i 'OrthoFinder version' | head -1 | awk '{print $3}' | tr -d '|')"
    : "${of_version:=unknown}"
    for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
        out_dir="$bench_dir/orthofinder_rep${rep}"
        # OrthoFinder writes to a sibling dir of -f; stage a per-rep
        # input dir holding the proteome FASTAs (spA.fasta, spB.fasta)
        # and let OrthoFinder write its results tree inside -o.
        # Per-rep so reps don't trip on each other's cached
        # BLAST/DIAMOND DB files.
        in_dir="$bench_dir/orthofinder_rep${rep}_in"
        mkdir -p "$in_dir"
        for fa in "$fixture_dir"/*.fasta; do cp -f "$fa" "$in_dir/"; done
        python "$repo/benchmarks/bench_one.py" \
            --task "T7" --tool "orthofinder" --version "$of_version" \
            --replicate "$rep" \
            --cmd "orthofinder -f '$in_dir' -og -S diamond -t 1 -a 1 \
                   -o '$out_dir' 2>/dev/null" \
            --notes "OrthoFinder cold start -og -S diamond (Emms 2019); per-OG FASTAs in Results_*/Orthogroup_Sequences/" \
            >> "$out"
    done
else
    echo "[skip] T7 orthofinder (orthofinder not on PATH; bioconda: orthofinder)" >&2
fi

# ---------- Proteinortho (secondary) ---------------------------------------
# Proteinortho (Lechner 2011, BMC Bioinformatics 12:124; Proteinortho6
# in Lechner 2023) uses a different orthology-graph algorithm and ships
# its own per-OG FASTA splitter via `proteinortho_grab_proteins.pl`.
# Skipped automatically when not on PATH.
if command -v proteinortho >/dev/null 2>&1 || command -v proteinortho6 >/dev/null 2>&1; then
    po_bin="$(command -v proteinortho6 || command -v proteinortho)"
    po_version="$($po_bin --version 2>&1 | head -1 | awk '{print $NF}' | tr -d '|')"
    : "${po_version:=unknown}"
    for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
        out_dir="$bench_dir/proteinortho_rep${rep}"
        in_dir="$bench_dir/proteinortho_rep${rep}_in"
        mkdir -p "$in_dir" "$out_dir"
        for fa in "$fixture_dir"/*.fasta; do cp -f "$fa" "$in_dir/"; done
        python "$repo/benchmarks/bench_one.py" \
            --task "T7" --tool "proteinortho" --version "$po_version" \
            --replicate "$rep" \
            --cmd "(cd '$out_dir' && \
                   $po_bin -project=gfcbench -cpus=1 -singles \
                       '$in_dir'/*.fasta 2>/dev/null && \
                   proteinortho_grab_proteins.pl -tofiles=. \
                       gfcbench.proteinortho.tsv '$in_dir'/*.fasta 2>/dev/null) || true" \
            --notes "Proteinortho + proteinortho_grab_proteins.pl (Lechner 2011/2023)" \
            >> "$out"
    done
else
    echo "[skip] T7 proteinortho (not on PATH; bioconda: proteinortho)" >&2
fi

echo "[done] T7 rows written to $out" >&2
python "$repo/benchmarks/check_correctness.py" --task T7 --bench-dir "$bench_dir" --tsv "$out" 2>&1 | head -20 || echo "[warn] T7 correctness check failed (non-fatal)" >&2
