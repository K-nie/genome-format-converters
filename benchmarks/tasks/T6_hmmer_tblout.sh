#!/usr/bin/env bash
# T6 — HMMER --tblout parsing.  Competitor: awk one-liner.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
task="T6_hmmer_tblout"
out="$repo/benchmarks/results/raw/${task}.tsv"
mkdir -p "$(dirname "$out")"

# Synthesise a small tblout fixture inline if one isn't provided.
input_dir="${GFC_BENCH_INPUT_DIR:-$repo/benchmarks/data/hmmer}"
mkdir -p "$input_dir"
if [[ ! -f "$input_dir/sample.tblout" ]]; then
    cat > "$input_dir/sample.tblout" <<'EOF'
#                                                                 --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name        accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#------------------- ---------- -------------------- ---------- --------- ------ -----   --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
YAL001C               -          PF12831.10           -          1.2e-45  154.2   0.3   2.1e-45  153.4   0.3   1.1   1   0   0   1   1   1   1 Uncharacterised protein
YAL002W               -          PF00022.25           -          9.8e-30  101.1   0.1   1.2e-29  100.8   0.1   1.0   1   0   0   1   1   1   1 Actin binding protein
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
        --task "T6" --tool "gfc" --version "$gfc_version" --replicate "$rep" \
        --cmd "gfc hmmer-tblout-to-tsv --input-dir '$input_dir' \
               --output-dir '$out_dir' --force" \
        >> "$out"
done

# ---------- pyhmmer in-process baseline -------------------------------------
# pyhmmer parses tblout via HMMER's own C library bindings — apples-to-apples
# competitor for column-aware tblout parsing, available on bioconda.
if python -c "import pyhmmer" 2>/dev/null; then
    pyhmmer_version="$(python -c 'import pyhmmer; print(pyhmmer.__version__)' 2>/dev/null || echo unknown)"
    for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
        out_dir="$bench_dir/pyhmmer_rep${rep}"
        mkdir -p "$out_dir"
        python "$repo/benchmarks/bench_one.py" \
            --task "T6" --tool "pyhmmer" --version "$pyhmmer_version" \
            --replicate "$rep" \
            --cmd "python '$repo/benchmarks/refs/t6_pyhmmer.py' \
                   --tblout '$input_dir/sample.tblout' \
                   --output '$out_dir/sample.tsv'" \
            >> "$out"
    done
else
    echo "[skip] T6 pyhmmer (pyhmmer not importable)" >&2
fi

# Removed in Stage 1 (bench/stage1-fairness): the awk one-liner baseline
# previously here. Its own --notes column already declared the output
# "NOT equivalent" because it drops the description field, so the speed
# comparison wasn't apples-to-apples. pyhmmer above is the correct
# in-process competitor; the reference parser lives in
# benchmarks/refs/t6_pyhmmer.py.

echo "[done] T6 rows written to $out" >&2
