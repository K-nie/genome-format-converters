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

# TODO: wire a handwritten Python reference — the point is to show that gfc
# matches "what the user would write themselves" while staying under one CLI.
echo "[done] T7 rows written to $out" >&2
