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

echo "[done] T7 rows written to $out" >&2
