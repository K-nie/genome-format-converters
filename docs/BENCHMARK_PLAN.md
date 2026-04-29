# Benchmarking plan — `gfc`

This document defines the benchmark tasks, competitor tools, datasets, metrics, and reproducibility harness that will go into the methods paper. Everything here should land as code under `benchmarks/` so reviewers can replay the numbers.

---

## 1. Benchmark tasks

Seven tasks, chosen to span the tool's feature coverage without being a sales pitch. At least two must include a case where gfc is **slower** than a specialist tool — the honest comparison is what sells a uniform-CLI paper.

| # | Task | `gfc` subcommand | Reference / competitor tools |
|---|---|---|---|
| T1 | VCF → EIGENSTRAT triplet | `vcf-to-eigenstrat` | convertf (EIGENSOFT), a handwritten `bcftools + awk` script, `pileupCaller` (pseudohap variant). |
| T2 | VCF → PLINK binary (`.bed`/`.bim`/`.fam`) | `vcf-to-plink` | `plink2 --vcf`, `plink1.9 --vcf --make-bed`. |
| T3 | FASTA + GFF → GenBank | `fasta-gff-to-gbk` | `EMBOSS seqret`, a handwritten Biopython script, `gff3toembl`. |
| T4 | GFF3 → GTF | `gff3-to-gtf` | `AGAT agat_convert_sp_gff2gtf.pl`, `gffread`. |
| T5 | GFF3 → BED12 | `gff3-to-bed12` | `UCSC gtfToGenePred` + `genePredToBed` chain, `AGAT agat_convert_sp_gff2bed.pl`. |
| T6 | HMMER `--tblout` → TSV | `hmmer-tblout-to-tsv` | `awk`/`cut` one-liner, `ESL-reformat` (not quite equivalent). |
| T7 | OrthoFinder `Orthogroups.tsv` → per-OG FASTA | `orthogroups-to-fasta` | Handwritten Python (typical lab script). |

Optional stretch goal:
| T8 | Pseudohaploid VCF → EIGENSTRAT with `--seed` reproducibility | `vcf-to-pseudohaploid` | `pileupCaller --randomHaploid`, a handwritten `bcftools query + python` script. Demonstrates reproducibility contract that isn't in the alternatives. |

---

## 2. Datasets

All datasets are either open-access or can be regenerated from published inputs. Each gets a line in Table 1 of the paper.

| Dataset | Size | Source | Used for |
|---|---|---|---|
| **Y1000+ slice — 20 representative Saccharomycotina genomes** | ~0.4 GB total (FASTA + GFF3) | Shen 2018 / Opulente 2024, Figshare / NCBI | T3 (fasta-gff-to-gbk), T4 (gff3-to-gtf), T5 (gff3-to-bed12). |
| **Y1000+ full — 1154 genomes** | ~20 GB | Same source | T3 at scale (the 20-genome case is the speed-comparison, the 1154-genome case is the stress-test). |
| **1000 Genomes Phase 3, chr22, 2504 samples, biallelic SNPs** | ~1.5 GB bgzip'd | 1000 Genomes FTP | T1, T2, T8. |
| **Saccharomyces cerevisiae × 100 isolates VCF (Peter 2018)** | ~0.5 GB | Peter et al. Nature 2018 | T1, T2 — yeast-native VCF comparison at population scale. |
| **HMMER scan of Pfam-A against Y1000+ proteomes** | ~2 GB tblout | Run in-house (already in Hittinger-lab pipeline) | T6. |
| **OrthoFinder run on 20 Saccharomycotina species** | ~50 MB TSV + ~0.5 GB FASTAs | Run in-house | T7. |
| **Simulated small fixtures** | <1 MB | `tests/test_data/` | Correctness pinning (not speed). |

---

## 3. Metrics

Three axes per task, three replicates each (increase to 5 if variance is high):

### 3.1 Wall time
- Measured with `/usr/bin/time -v` on Linux, `gtime -v` on macOS (coreutils), or a Python wrapper using `time.perf_counter()` for portability.
- Report mean ± stdev across 5 replicates.
- Cold-cache: drop OS page cache between replicates (`echo 3 > /proc/sys/vm/drop_caches` on Linux, reboot-level on macOS; document if we can't clear cache).

### 3.2 Peak resident memory
- `/usr/bin/time -v` → "Maximum resident set size (kbytes)".
- For Python-internal comparison we can additionally use `tracemalloc` / `psutil`; for comparing against C tools, RSS from `time -v` is the fair metric.

### 3.3 Output correctness
- **Byte-identical** where possible (e.g. `gfc vcf-to-plink` vs `plink2 --vcf` should produce the same `.bed` modulo `.bim` column conventions).
- **Row-count + spot-check** where byte-identical isn't achievable (different tools sort / name SNPs differently; document conversion rules before comparison).
- **Golden files** checked into `tests/test_data/golden/` with a diff tolerance specified per task.

---

## 4. Hardware and environment

Document in the paper's Methods section exactly once:

- **CPU**: <fill in at benchmark time — e.g. AMD EPYC 7452, 32 cores, 2.35 GHz>
- **RAM**: <fill in — e.g. 256 GB DDR4>
- **Disk**: <fill in — NVMe SSD, 2 TB>
- **OS**: Linux (Ubuntu 22.04 LTS, kernel 5.15), macOS 14 Sonoma (Apple Silicon M2 for secondary run).
- **Python**: 3.11.x (pinned — freeze version in `benchmarks/environment.yml`).
- **Package versions**: `gfc 0.1.5`, pysam 0.22.1, biopython 1.83, bcbio-gff 0.7.x, pandas 2.x.
- **Competitor versions**: bcftools 1.19, samtools 1.19, plink 2.00a5, AGAT 1.4, EMBOSS 6.6, pileupCaller 1.5, OrthoFinder 2.5.

---

## 5. Harness

```
benchmarks/
├── README.md                  # one-page runbook
├── environment.yml            # conda env lock
├── Dockerfile                 # reproducible runner
├── data/                      # download scripts (do not commit raw data)
│   ├── download_y1000plus.sh
│   ├── download_1kg_chr22.sh
│   └── download_pfam.sh
├── run_bench.sh               # top-level driver
├── tasks/
│   ├── T1_vcf_to_eigenstrat.sh
│   ├── T2_vcf_to_plink.sh
│   ├── T3_fasta_gff_to_gbk.sh
│   ├── T4_gff3_to_gtf.sh
│   ├── T5_gff3_to_bed12.sh
│   ├── T6_hmmer_tblout.sh
│   ├── T7_orthogroups.sh
│   └── T8_pseudohaploid.sh
├── analyze/
│   ├── collect_results.py     # parse /usr/bin/time output into one TSV
│   ├── plot_wall_time.py      # Fig 3
│   ├── plot_memory.py         # Fig 4
│   └── feature_matrix.py      # Fig 2 (checkmark grid)
└── results/
    ├── raw/                   # one TSV per task, per replicate
    └── figures/               # rendered PNG / SVG / PDF
```

Each task script:
1. Validates inputs are downloaded.
2. Runs each competitor + gfc 5 times.
3. Captures wall time + peak RSS + exit code.
4. Emits one TSV row per replicate: `task, tool, version, replicate, wall_s, rss_kb, exit_code, correct (bool), notes`.

`run_bench.sh` concatenates every task's output into `results/raw/all.tsv`, then hands it to `analyze/`.

---

## 6. Expected outputs (paper-ready)

- **Figure 2** (feature matrix): one PNG / one PDF at 300 dpi. Checkmark grid, gfc column flagged with a subtle highlight. Generated by `analyze/feature_matrix.py`.
- **Figure 3** (wall time): bar chart with grouped bars per task, 7-8 tasks on the x-axis, ~4 tools per task. Log scale on y-axis (tasks span seconds to hours).
- **Figure 4** (memory): same layout as Figure 3, peak RSS in GB.
- **Table 1** (datasets): handwritten once; regenerate row counts from the download scripts.
- **Table 2** (correctness): byte-diff results, with "✓" or "≈ (reason)" per cell.
- **Supplementary Table S1** (exact commands): the command line for every tool × task combination, with flags and version strings.

---

## 7. Definition of "win"

We're **not** claiming gfc is the fastest tool at any one task — it's written in Python, competing with C-based specialists, so that claim is rarely true and would rightly get pushed back.

We **are** claiming:

1. **Coverage**: gfc handles all 30 tasks under one CLI; no other single tool does this.
2. **Parity on correctness**: gfc output matches the reference tools byte-for-byte (or in documented, justified ways differs).
3. **Within-order-of-magnitude speed**: gfc is within 3× of the specialist tool's wall time on every task. If it's not, fix the implementation before publishing.
4. **Lower memory in several cases**: where the specialist tool loads entire tables into memory (e.g. AGAT loading a whole GFF), gfc's streaming approach should win on memory.
5. **Superior reproducibility**: per-file RNG seeding for pseudohaploid calls (T8) is a contract no competitor currently offers.

---

## 8. Threats to the comparison (and mitigations)

| Threat | Mitigation |
|---|---|
| Different tools interpret the same input slightly differently (e.g. coordinate conventions). | Use small hand-curated fixtures for correctness pinning; document every conversion rule. |
| Disk I/O variance on shared HPC. | Run on a quiescent machine; 5 replicates; report median not just mean. |
| Python startup cost inflates wall time on short inputs. | Use datasets sized so Python startup is <5% of total. |
| Tool versions drift during review. | Pin versions in `environment.yml`; build Docker image; archive on Zenodo alongside raw results. |
| Reviewer asks "did you tune each tool?" | For each competitor, use its recommended flags from the upstream docs; commit the exact invocations in `tasks/*.sh`. |

---

## 9. Milestones

| Milestone | Target | Blocker |
|---|---|---|
| Harness skeleton committed | +1 week | — |
| Datasets downloaded + checksums committed | +2 weeks | network / storage |
| T1-T3 run end-to-end | +3 weeks | Y1000+ slice selection |
| T4-T7 run end-to-end | +4 weeks | AGAT install on macOS (known flaky) |
| Figures 2-4 generated | +5 weeks | — |
| Figures committed to `paper/figures/` | +5 weeks | — |
| First draft of paper text citing these numbers | +6 weeks | — |

---

## 10. Honest caveats to disclose in the paper

1. gfc uses pysam internally, which itself uses htslib. We're not an independent implementation of VCF/BAM parsing — correctness vs bcftools / samtools is therefore not fully independent.
2. Python's global interpreter lock means per-file parallelism (`--threads`) is the realistic ceiling; we're not going to out-parallelise `bcftools --threads`.
3. Windows support is subset-only (pure-Python 20/30 subcommands). Flag this explicitly; don't pretend otherwise.
4. The `--also-plink` sidecar on EIGENSTRAT is text `.ped` / `.map`; for binary use the standalone `vcf-to-plink`. Reviewers will ask which is "primary"; document that the binary format is canonical and `.ped` is convenience.

---

## Status (2026-04-29)

Stage 1 has shipped. Production benchmark run completed on HTCondor
cluster job 136866 (10 replicates × 8 tasks × 14 competitors, ~10 hour
wallclock, sd < 2% on every task). Results, figures, and Methods write-up:

- Manuscript: [`docs/manuscript/gfc_application_note_2026-04-29.md`](manuscript/gfc_application_note_2026-04-29.md)
- Per-task TSVs: `benchmarks/results/raw/T*.tsv` (n=10 replicates each)
- Headline tables: `benchmarks/results/figures/summary_table.md`,
  `summary_per_pair.md`
- Plots (PNG + PDF, 300 dpi): `benchmarks/results/figures/`
- Per-run notes: `docs/bench-runs/`

Stage 2 / paper-grade extensions (full Y1000+ scale on T3–T5; real
~10 K-row Pfam-A scan on T6; biallelic-SNP-restricted plink2 baseline
on T2) are tracked under the manuscript's revision plan rather than
this planning document.
