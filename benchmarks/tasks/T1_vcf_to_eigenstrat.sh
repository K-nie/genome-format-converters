#!/usr/bin/env bash
# T1 — VCF → EIGENSTRAT triplet (.geno/.snp/.ind).
#
# Competitors referenced in docs/BENCHMARK_PLAN.md:
#   - EIGENSOFT convertf (install from source; see data/README.md)
#   - handwritten bcftools + awk + python (included inline, TODO)
#   - pileupCaller (pseudohaploid variant) — covered by T8
#
# This script runs gfc end-to-end on the bundled tiny.vcf fixture so the
# harness is immediately executable. For paper-grade numbers, point it at
# the 1000G chr22 or Peter-2018 yeast-VCF datasets downloaded by
# data/download_1kg_chr22.sh.

set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"

task="T1_vcf_to_eigenstrat"
out="$repo/benchmarks/results/raw/${task}.tsv"
mkdir -p "$(dirname "$out")"

# Smoke dataset (committed fixture). Swap for a downloaded real dataset
# when generating paper numbers.
input_dir="$repo/tests/test_data"
: "${GFC_BENCH_INPUT_DIR:=$input_dir}"
: "${GFC_BENCH_VCF_PATTERN:=tiny.vcf}"
: "${GFC_BENCH_REPLICATES:=5}"

bench_dir="/tmp/gfc_bench_${task}"
rm -rf "$bench_dir"
mkdir -p "$bench_dir"

# Header row
python "$repo/benchmarks/bench_one.py" --header > "$out"

# ---------- tool 1: gfc -----------------------------------------------------
gfc_version="$(gfc --version | cut -d' ' -f2)"
for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
    out_dir="$bench_dir/gfc_rep${rep}"
    python "$repo/benchmarks/bench_one.py" \
        --task "T1" --tool "gfc" --version "$gfc_version" --replicate "$rep" \
        --cmd "gfc vcf-to-eigenstrat --input-dir '$GFC_BENCH_INPUT_DIR' \
               --output-dir '$out_dir' --pattern '$GFC_BENCH_VCF_PATTERN' --force" \
        --notes "fixture=${GFC_BENCH_VCF_PATTERN}" \
        >> "$out"
done

# ---------- tool 2: handwritten bcftools + awk + python ---------------------
# Hand-script that reproduces the EIGENSTRAT .geno encoding from a VCF by
# piping bcftools through awk. This is the "what a careful bioinformatician
# would write in an afternoon" comparison.
#
# Enabled when `bcftools` is on PATH. Skipped otherwise — paper numbers
# should come from a run that has it.
if command -v bcftools >/dev/null 2>&1; then
    bcftools_version="$(bcftools --version | head -1 | cut -d' ' -f2)"
    handscript="$repo/benchmarks/tasks/_T1_handscript.sh"
    cat > "$handscript" <<'EOF'
#!/usr/bin/env bash
# Hand-script: bcftools → awk → .geno / .snp / .ind. Not a real EIGENSTRAT
# converter — just a proof-of-speed baseline.
set -euo pipefail
vcf="$1"; outdir="$2"
mkdir -p "$outdir"
stem="$(basename "$vcf" .vcf)"
stem="${stem%.vcf.gz}"
bcftools view -m2 -M2 -v snps "$vcf" \
  | bcftools norm -m+ \
  | awk -F'\t' '
      BEGIN { OFS="\t" }
      /^#CHROM/ {
          for (i = 10; i <= NF; i++) samples[i] = $i
          print "" > "/dev/null"
          next
      }
      /^##/ { next }
      {
          printf "%s_%s\t%s\t0.0\t%s\t%s\t%s\n", $1, $2, $1, $2, $4, $5 > snp
          line = ""
          for (i = 10; i <= NF; i++) {
              split($i, gt, /[\/|]/)
              a = gt[1]; b = gt[2]
              if (a == "." || b == ".") { line = line "9"; continue }
              n = (a == "0") + (b == "0")
              line = line n
          }
          print line > geno
      }
  ' snp="$outdir/$stem.snp" geno="$outdir/$stem.geno"
bcftools query -l "$vcf" | awk '{ print $1"\tU\t"$1 }' > "$outdir/$stem.ind"
EOF
    chmod +x "$handscript"

    for rep in $(seq 1 "$GFC_BENCH_REPLICATES"); do
        out_dir="$bench_dir/hand_rep${rep}"
        vcf_input="$(ls "$GFC_BENCH_INPUT_DIR"/$GFC_BENCH_VCF_PATTERN 2>/dev/null | head -1)"
        python "$repo/benchmarks/bench_one.py" \
            --task "T1" --tool "bcftools+awk" --version "$bcftools_version" \
            --replicate "$rep" \
            --cmd "bash '$handscript' '$vcf_input' '$out_dir'" \
            --notes "hand-script baseline" \
            >> "$out"
    done
else
    echo "[skip] T1 bcftools+awk baseline (bcftools not on PATH)" >&2
fi

# ---------- tool 3: EIGENSOFT convertf --------------------------------------
# TODO: wire up when eigensoft/convertf is available. convertf takes a
# parameter file (par) pointing at .ped/.map or .pedind; construct those
# from the VCF via plink2 --vcf --recode, then invoke convertf.
#
# if command -v convertf >/dev/null 2>&1; then
#     ... populate par file ...
#     python bench_one.py --task T1 --tool convertf --version "$cf_ver" ...
# fi

echo "[done] T1 rows written to $out" >&2
