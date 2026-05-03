# Benchmark harness

Reviewer-replayable benchmark suite for `gfc`, per the plan in
`docs/BENCHMARK_PLAN.md`. Every number that appears in the paper comes out
of this directory.

## Quick start

```bash
# 1. Create the benchmark environment (conda recommended; competitors are
#    easier to install that way than via pip).
mamba env create -f benchmarks/environment.yml
conda activate gfc-bench

# 2. Download datasets. Each script skips if its target already exists.
bash benchmarks/data/download_y1000plus.sh
bash benchmarks/data/download_1kg_chr22.sh
bash benchmarks/data/download_pfam_scan.sh

# 3. Run the full suite.
bash benchmarks/run_bench.sh

# 4. Render figures and tables.
python benchmarks/analyze/collect_results.py
python benchmarks/analyze/plot_wall_time.py
python benchmarks/analyze/plot_peak_memory.py
python benchmarks/analyze/plot_feature_matrix.py
```

All tasks are opt-out, not opt-in: comment out lines in `run_bench.sh` to
skip a task. Each task emits one TSV into `results/raw/T*.tsv`; the
collector merges them into `results/raw/all.tsv`.

## Layout

```
benchmarks/
├── README.md              # this file
├── environment.yml        # conda env (gfc + competitors)
├── Dockerfile             # reproducible runner (optional)
├── bench_one.py           # wall-time + peak-RSS wrapper around any command
├── run_bench.sh           # top-level driver — invokes every tasks/T*.sh
├── data/                  # dataset download scripts (no datasets committed)
├── tasks/                 # one *.sh per benchmark task (T1 .. T7 + T8)
├── analyze/               # result collection + figure generation
└── results/
    ├── raw/               # per-task TSVs + merged all.tsv
    └── figures/           # rendered PNG / PDF / SVG
```

## What each TSV row records

All task scripts emit rows matching this schema:

```
task    tool    version    replicate    wall_s    peak_rss_mb    exit_code    correct    notes
```

- `task` — `T1` .. `T8`.
- `tool` — e.g. `gfc`, `bcftools`, `plink2`, `convertf`, `AGAT`, `seqkit`.
- `version` — captured at run time (e.g. `gfc --version` → `0.1.5`).
- `replicate` — 1 .. N (defaults to 5).
- `wall_s` — wall-clock seconds, measured with `time.perf_counter()` in
  `bench_one.py` — identical accounting across OSes.
- `peak_rss_mb` — peak resident memory of the child process tree (via
  `psutil`). Monotonic, polled at 50 ms intervals.
- `exit_code` — the child's exit code; non-zero marks the replicate as
  failed.
- `correct` — `1` if the output matches the golden-hash fixture for that
  task; `0` if it does not; empty if correctness isn't checked in this
  task.
- `notes` — free-text column for remarks like "cold cache", "unindexed
  input", etc.

## Honest caveats (also in docs/BENCHMARK_PLAN.md)

1. `gfc` uses `pysam` which wraps `htslib`, so correctness vs `bcftools`
   (also htslib) isn't fully independent.
2. Python startup dominates wall time for tasks with <1000 records; the
   dataset sizes in `data/` are chosen so that startup is <5 % of total.
3. macOS runs drop OS page cache imperfectly; prefer Linux for numbers
   that go into the paper.
4. Windows runs the pure-Python subset (see
   `.github/workflows/ci.yml`); VCF / BAM benchmarks require Linux or
   macOS with `pysam` installed.

## Reproducibility

A pinned Docker image is provided (`Dockerfile`) so reviewers can
reproduce the numbers without having to wrestle AGAT or EIGENSOFT onto
their local machine.

```bash
docker build -t gfc-bench -f benchmarks/Dockerfile .
docker run --rm -v $(pwd)/benchmarks/results:/results gfc-bench
```
