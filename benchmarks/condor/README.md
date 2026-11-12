# HTCondor runbook for gfc benchmarks

Target host: **`scarcity-ap-1.glbrc.org`** (GLBRC HTCondor submit host).

Three files drive the run:

| File | Where it runs | What it does |
|---|---|---|
| `bootstrap.sh` | submit host (once) | Builds the `gfc-bench` conda env, `pip install -e`s gfc into it, pre-creates `logs/` |
| `bench.sub` | submit host (`condor_submit`) | Queues one HTCondor job with 4 CPUs / 16 GB RAM / 50 GB disk |
| `run_on_condor.sh` | execute host | Activates the env, downloads datasets, runs all 8 tasks × 5 replicates, renders figures, tars results |

## End-to-end runbook

### 0. Create the remote working directory

Already done by you — assume the repo will live at `~/gfc-bench/` on `scarcity-ap-1`. Substitute your actual path throughout.

### 1. Push the repo to the cluster

From **your Mac**, in this repo's root:

```
rsync -avz --delete --exclude='.git/' --exclude='benchmarks/data/' --exclude='benchmarks/results/' --exclude='dist/' --exclude='build/' --exclude='.pytest_cache/' --exclude='__pycache__' ./ narhmadey@scarcity-ap-1.glbrc.org:~/gfc-bench/
```

One line, no backslashes (GLBRC paste convention). The exclude list skips the git history, datasets that will be redownloaded cluster-side, and Python build artefacts.

### 2. SSH in and bootstrap the env

```
ssh narhmadey@scarcity-ap-1.glbrc.org
```

Then on the cluster:

```
cd ~/gfc-bench && bash benchmarks/condor/bootstrap.sh
```

First run takes 5-15 minutes — `mamba` resolves every bioconda dependency (bcftools / samtools / plink2 / AGAT / EMBOSS / HMMER / UCSC tools) in one solve. Subsequent runs just update the env.

### 3. Submit the benchmark job

```
cd ~/gfc-bench/benchmarks/condor && condor_submit bench.sub
```

HTCondor prints something like:

```
Submitting job(s).
1 job(s) submitted to cluster 123456.
```

Note the **cluster ID** — needed for monitoring.

### 4. Watch it run

Either of these work:

```
condor_q -nobatch
```

```
tail -f ~/gfc-bench/benchmarks/condor/logs/bench.123456.out
```

Expected runtime on the 1000G chr22 dataset: **30-60 minutes** for all 8 tasks × 5 replicates with the competitors present. Longer if the Y1000+ slice is added, or if AGAT is present (AGAT validates the entire GFF hierarchy — it's slow by design).

### 5. Pull results back to your Mac

When the job finishes (`condor_q` shows no rows for your cluster ID):

```
rsync -avz narhmadey@scarcity-ap-1.glbrc.org:~/gfc-bench/benchmarks/results/ ./benchmarks/results/
```

The tarball `results.YYYYMMDD-HHMMSS.tar.gz` also comes back inside that.

### 6. What you get

```
benchmarks/results/
├── hardware.log                        # CPU / mem / OS / software versions — paper Methods §4
├── raw/
│   ├── T1_vcf_to_eigenstrat.tsv        # 5 replicates × each tool
│   ├── T2_vcf_to_plink.tsv
│   ├── T3_fasta_gff_to_gbk.tsv
│   ├── T4_gff3_to_gtf.tsv
│   ├── T5_gff3_to_bed12.tsv
│   ├── T6_hmmer_tblout.tsv
│   ├── T7_orthogroups.tsv
│   ├── T8_pseudohaploid.tsv
│   └── all.tsv                         # merged
├── figures/
│   ├── wall_time.{png,pdf}             # paper Fig 3
│   ├── peak_memory.{png,pdf}           # paper Fig 4
│   ├── feature_matrix.{png,pdf,md}     # paper Fig 2
└── results.YYYYMMDD-HHMMSS.tar.gz      # bundle of the above
```

## Tuning knobs

Set these in the submit host's shell before `condor_submit bench.sub` if you want to override the defaults. `getenv = True` in `bench.sub` propagates them.

| Env var | Default | Effect |
|---|---|---|
| `GFC_BENCH_REPLICATES` | 5 | N replicates per tool per task. Set to 3 for a faster first pass. |
| `GFC_BENCH_INPUT_DIR` | `$repo/tests/test_data` or `$repo/benchmarks/data/1kg` if downloaded | Input directory for VCF / FASTA tasks. |
| `GFC_BENCH_VCF_PATTERN` | `tiny.vcf` or `chr22.biallelic.snps.filtered.vcf.gz` | Glob pattern for VCF input files. |

## Troubleshooting

**`condor_q` shows your job in `H` (held) state:**
Usually a resource-limit mismatch. `condor_q -hold` gives the reason. Most common: `request_memory` too low — bump `bench.sub`'s `request_memory` to `32 GB` and resubmit.

**Job hits the wall-clock ceiling:**
Reduce `GFC_BENCH_REPLICATES` or disable the slowest task (AGAT is usually the culprit). Edit `benchmarks/run_bench.sh` to comment out specific tasks, or pass a subset on the command line:

```
bash benchmarks/run_bench.sh T1 T2 T8
```

Then resubmit by editing `run_on_condor.sh` to pass the subset.

**Conda env broken on execute host:**
Most likely cause is the execute host not sharing home filesystem with submit host. Confirm with `condor_status -long <execute-host>` or ask GLBRC ops. If the filesystems aren't shared, the env needs to be `conda-pack`ed and transferred — out of scope for this first iteration.

**Need to rerun from scratch:**
```
rm -rf ~/gfc-bench/benchmarks/results
cd ~/gfc-bench/benchmarks/condor && condor_submit bench.sub
```

Datasets in `~/gfc-bench/benchmarks/data/` are cached — the download scripts skip if the target already exists, so re-runs are fast.

## What's still TODO before this is paper-grade

- `benchmarks/data/download_y1000plus.sh` has placeholder URLs — pin real figshare / NCBI paths before the Y1000+ tasks produce numbers.
- `benchmarks/data/download_pfam_scan.sh` has placeholder URL — either host the pre-computed tblout on Zenodo or regenerate on the cluster via `hmmsearch`.
- EIGENSOFT `convertf` and `pileupCaller` aren't in bioconda; they need source builds if you want T1 / T8 with proper references (skipped gracefully today).
- Currently one job, one host, sequential. If a single-host run takes too long on the full dataset, split into 8 per-task jobs and use `condor_dagman` for the aggregation step. Out of scope for first pass.

## Paper rerun (`paper_rerun.sub`)

The paper's numbers come from one rerun split into six HTCondor jobs, one per
phase. Each job writes to its own `benchmarks/results/raw_<phase>/` and
`benchmarks/results/<phase>/` and tars the latter to
`benchmarks/results/<phase>.<timestamp>.tar.gz`, so the jobs can run at the
same time.

| Phase | What it measures | Paper use |
|---|---|---|
| `main` | 8 tasks on chr22 / Y1000+ slice, n=10, four-state verdicts, T8 reproducibility | Figs 2 and 4, Table 1 |
| `annotations` | T4/T5 on Ensembl R64 and GRCh38 chr22 gene annotations, n=10 | Results 3.2 |
| `scaling` | T1/T2 at 10k/100k/1M sites, T4/T5 on chr22 at x1/x5/x25, n=3 | Fig 3 |
| `cold` | main panel with input files evicted from the page cache before each replicate, n=5 | Supp. Table S2 |
| `import` | start-up wall time and RSS per subcommand | Results 3.3 |
| `footprint` | package count and disk size of fresh envs: gfc, the multi-tool stack, BioConvert, biopython.convert | Table 2 |

Steps on the submit host, after pulling this branch:

```
cd ~/gfc-bench && bash benchmarks/condor/bootstrap.sh
cd benchmarks/condor && condor_submit paper_rerun.sub
```

To run a single phase, replace the last line of `paper_rerun.sub` with, for
example, `queue phase in (main)`. To keep all timing phases on the same
hardware, uncomment the `requirements` line. Then copy
`benchmarks/results/*.tar.gz` back.

Knobs (set in the environment before `condor_submit`, which passes
`getenv = True`): `GFC_BENCH_REPLICATES`, `GFC_BENCH_PLINK_MEM_MB` (default
512), `GFC_BENCH_PLINK_DEFAULT_MEM=1` (also record uncapped plink rows),
`GFC_SCALE_SITES`, `GFC_SCALE_FACTORS`.

AGAT on GRCh38 chr22 x25 is the slowest cell. If the scaling job runs past the
cluster's wall-time limit, lower `GFC_SCALE_FACTORS` to `1 5 10`.
