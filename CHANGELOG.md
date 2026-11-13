# Changelog

All notable changes to this project are documented here. The format is
based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] — paper-plan engineering pass

### Changed
- **Lazy converter imports.** The CLI imports a converter module only when
  its subcommand runs, and `converters/__init__.py` no longer imports every
  submodule. `gfc --help` and the annotation/HMMER commands no longer load
  Biopython, BCBio.GFF, pandas, numpy, pysam or cyvcf2 (CLI start-up peak
  RSS 55 MB -> 19 MB locally). `benchmarks/bench_import.py` measures it.

### Benchmarks
- The T1/T2 comparators pair records by chrom, position and allele set,
  so split multi-allelic sites (two biallelic records at one position, as in
  1000 Genomes phase 3) are no longer reported as disagreements.
- The T1/T2 comparators stream both outputs site by site and buffer only
  records that are out of step, so memory stays flat on whole chromosomes
  (300k sites x 2,504 samples: 3.1 GB -> 57 MB for T1). Every task's
  correctness check runs under `timeout` (GFC_BENCH_CORR_TIMEOUT, default
  1 h) and logs each tool it checks.
- The T8 reference-free check streams the VCF with numpy per site instead of
  holding every genotype in memory (about 20 GB on chr22). The T8
  reproducibility runs use the first 20,000 sites, log each step, and time
  out after an hour per call. `--threads` workers get stdin closed.
- T8 is now a reproducibility and correctness task
  (`benchmarks/t8_reproducibility.py`); the handwritten bcftools baseline is
  gone. plink runs with a capped `--memory`; T1 adds convertf-only rows.
- Size-scaling series (`run_scaling.sh`, `make_scaled_inputs.py`), Ensembl
  gene annotations for T4/T5, cold-cache mode (`GFC_BENCH_EVICT_DIRS`),
  install-footprint measurement (`install_footprint.sh`), and a per-phase
  HTCondor rerun (`condor/paper_rerun.sub`).

### Added
- **`eigenstrat-to-vcf`**: EIGENSTRAT triplet back to VCF, diploid or
  `--haploid` (pseudohaploid), with optional `--ref-fasta` re-orientation.
  Rejects packed ANCESTRYMAP/TGENO input with a clear message.
- **`vcf-to-snp-alignment`**: VCF to relaxed PHYLIP / FASTA / NEXUS SNP
  alignment, IUPAC (default), random or missing coding of hets,
  `--min-samples`, `--drop-invariant`.
- **Four-state correctness comparators** (`benchmarks/four_state.py`):
  byte / structural / limitation / disagreement, with the tolerated
  convention recorded per verdict. `check_correctness.py` writes them to
  `<task>.verdicts.tsv`; T8 is now checked against its input VCF.
- **Mutation test** (`benchmarks/mutation_test.py`): 29 seeded defects
  across T1-T8 are all reported as `disagreement`; convention rewrites
  stay `structural`.

### Changed
- **VCF converters** (`vcf-to-eigenstrat`, `vcf-to-plink`,
  `vcf-to-pseudohaploid`) read GT-only VCF text with a numpy decoder and
  fall back to cyvcf2 for anything irregular; MAF/missingness filters and
  the PLINK `.ped` sidecar are vectorised. Output is byte-identical to
  before; 2-4x faster on a 20k-site x 2,504-sample test file.
- **`--threads N`** now runs batch conversions in parallel worker
  processes (previously it only printed a warning). Output is identical
  to a sequential run.

## [Unreleased] — Stage 1 benchmark + manuscript

### Added
- **Stage 1 benchmark suite** (`benchmarks/`) covering 8 head-to-head
  conversion tasks (T1 VCF→EIGENSTRAT through T8 pseudohaploid) against
  EIGENSOFT convertf, PLINK 2, PLINK 1.9, AGAT, gffread, UCSC kent
  tools, EMBOSS seqret, pyhmmer, and handwritten Biopython baselines.
  HTCondor runner + per-task TSVs + 9 publication-grade plots + 3
  summary tables. Production run at n=10 reps per (task, tool).
- **`benchmarks/check_correctness.py`** — per-task format-aware
  comparators (byte-identical for T1/T6/T7, structural for T3/T4/T5,
  format-only for T8). Format-divergence cases between gfc and the
  canonical reference are documented as `_REASON_T*` constants paste-
  ready into a paper Methods section.
- **`benchmarks/bench_one.py`** — single-rep timer using `os.wait4` +
  `getrusage(ru_maxrss)` for per-child wall + peak RSS, with explicit
  `flush=True` on row writes to prevent mid-shutdown buffer-loss under
  cluster SIGTERM.
- **`benchmarks/analyze/`** — collect_results, summary_table (per-task +
  per-pair), correctness-aware plot_wall_time, plot_memory,
  plot_correctness_matrix, plot_pareto, plot_speedup, plot_variance,
  plot_per_task_panels, plot_compute_cost, hardware_table.
- **VCF converters rewritten on cyvcf2 + numpy bulk encoders**
  (`vcf_to_eigenstrat`, `vcf_to_plink`, `vcf_to_pseudohaploid`):
  ~12.3× speedup vs the prior pysam per-sample dict-lookup pattern;
  peak RSS drops from 3.3 GB to 60 MB on chr22 (1KG biallelic SNPs).
- **GFF / FASTA converters tightened** (`gff3_to_gtf`, `gff3_to_bed12`,
  `convert_all_gff_fasta_to_gbk`, `orthogroups_to_fasta`): two-pass
  streaming text parser replacing bcbio-gff `GFF.parse` on T4/T5; raw
  FASTA reader replacing biopython on T7.
- **Application Note manuscript draft** at
  `docs/manuscript/gfc_application_note_2026-04-29.md` (Q1/Q2 target),
  with a plain-text mirror `gfc_application_note_2026-04-29.txt` for
  journal submission portals that don't accept Markdown.

### Changed
- `pyproject.toml` adds `cyvcf2>=0.30` and `numpy>=1.20` as runtime
  deps (Linux + macOS).

## [0.1.5] — 2026-04-24

### Added
- **Windows support** (pure-Python subset). `pysam` is now a
  platform-conditional dependency (`platform_system != 'Windows'`); the
  20/30 subcommands that don't touch VCF/BAM/BCF work natively on Windows
  via `pip install`. VCF/BAM subcommands emit a readable pointer to WSL or
  conda-forge instead of a cryptic `ImportError`. CI matrix now covers
  `windows-latest` × Python 3.10-3.12 for the pure-Python subset.
- **ASCII banner** on bare `gfc` invocation (TTY-guarded so it doesn't
  corrupt piped output).
- **`gfc guide`** — categorised cheat sheet of every subcommand with one-
  line summaries.
- **`gfc examples <subcommand>`** — worked invocations per subcommand.
- **Per-subcommand `--help` epilog** — every subparser now shows 2-3 real
  example invocations at the bottom of its `--help` output.
- **`docs/PAPER_PLAN.md`** — journal shortlist (JOSS first, then BMC
  Bioinformatics / Bioinformatics Advances / Bioinformatics App Note),
  paper outline, selling points, timeline, authorship template, risk
  register.
- **`docs/BENCHMARK_PLAN.md`** — 7 benchmark tasks with competitor tools,
  datasets (Y1000+, 1000G chr22, Peter 2018, Pfam scan), metrics, harness
  layout, honest caveats.

### Fixed
- Bare `gfc` (no subcommand) now prints the banner + a pointer to
  `gfc guide` / `gfc --help` instead of only the argparse usage line.

## [0.1.4] — 2026-04-24

### Added
- **`gff3-to-bed12`** — emit BED12 with `blockCount` / `blockSizes` /
  `blockStarts` built from each transcript's exon children. The existing
  `gff3-to-bed` (BED6) is unchanged. IGV/UCSC-ready.
- **`gtf-to-gff3`** — inverse of `gff3-to-gtf`. Parses Ensembl-style GTF
  attribute strings and synthesises explicit `gene` rows from `gene_id`
  groupings when the GTF omits them, so downstream GFF3 consumers get a
  valid `ID` / `Parent` hierarchy.
- **`stockholm-to-fasta`** and **`fasta-to-stockholm`** — Rfam / Pfam
  alignment round-trip via Biopython AlignIO.
- **`hmmer-tblout-to-tsv`** — parse HMMER `--tblout` (default) or
  `--domtblout` (`--hmmer-format domtblout`) into clean TSV with a fixed
  column schema and a correctly-preserved free-text `description_of_target`.
- **`orthogroups-to-fasta`** — OrthoFinder `Orthogroups.tsv` + a directory
  of per-species protein FASTAs → one FASTA per orthogroup with
  `>species|gene` headers. Takes `--orthogroups`, `--fasta-dir`,
  `--output-dir` instead of the standard `--input-dir` shape.
- **`vcf-to-plink`** — companion to `vcf-to-eigenstrat`, emitting the
  PLINK binary triplet (`.bed` / `.bim` / `.fam`). Shares the polarisation
  / filter / chrom-map / genetic-map plumbing with the EIGENSTRAT
  converters so the CLI surface stays uniform.

## [0.1.3] — 2026-04-24

### Added
- `vcf-to-eigenstrat` — convert VCF/BCF to EIGENSTRAT `.geno` / `.snp` /
  `.ind` triplets for smartpca, ADMIXTOOLS2, qpAdm, and friends. Only
  biallelic SNPs are emitted; indels, multi-allelic sites, and non-ACGT
  alleles are skipped with per-file skip counts on stderr.
- `vcf-to-pseudohaploid` — as above but one allele per genotype is drawn
  uniformly at random, so `.geno` contains only 0, 2, and 9. RNG is seeded
  per file (`seed XOR stable_hash(filename)`) so batch output is
  independent of file order.
- `--pop-map`, `--sex-map`, `--chrom-map`, `--default-chrom`,
  `--genetic-map`, `--ancestral-fasta`, `--info-aa`, `--transversions-only`,
  `--min-maf`, `--max-missing`, `--strict-maps`, and `--also-plink`
  (PLINK `.ped`/`.map` sidecar) on both EIGENSTRAT converters.
- Uniform `--input`/`-i` and `--output`/`-o` single-file mode alongside the
  existing `--input-dir`/`--output-dir` batch mode (1:1 converters only).
- Global `--quiet` / `--verbose` flags and per-subcommand `--pattern`,
  `--force`, `--dry-run`, and `--threads` flags.
- `--het-policy` (`iupac`/`ref`/`alt`) and `--mask-missing` on
  `vcf-to-consensus`.
- pytest suite (`tests/test_converters.py`) with 23 tests covering every
  converter, a golden-file test for EIGENSTRAT, and a reproducibility test
  for pseudohaploid with `--seed`.
- GitHub Actions CI matrix (Ubuntu + macOS × Python 3.9/3.10/3.11/3.12).
- `Dockerfile` (multi-stage, `python:3.11-slim`).
- conda-forge recipe skeleton at `recipes/genome-format-converters/meta.yaml`.
- `CITATION.cff`, `CONTRIBUTING.md`, `CHANGELOG.md`.

### Fixed
- **`vcf-to-consensus`** — now handles multi-contig references (previously
  only read the first FASTA record), indels (insertions/deletions now
  change output length instead of silently corrupting positions),
  heterozygous SNPs (IUPAC ambiguity code by default), and missing
  genotypes (`--mask-missing` for explicit `N`). Also accepts `.fa`,
  `.fna`, `.fas` extensions in addition to `.fasta`.
- **`gff3-to-gtf`** — feature IDs are no longer silently truncated on `:`
  (Ensembl-style `gene:YAL001C` stays intact). `gene` rows emit only
  `gene_id`; transcript / exon / CDS rows emit `gene_id` + `transcript_id`
  with the correct parent transcript.
- **`gff3-to-protein`** — spliced CDS fragments are now concatenated in
  phase order before translation. Previously each CDS fragment was
  translated independently, producing one partial protein per CDS instead
  of one correct protein per transcript.
- **`gff3-to-bed`** and **`gff3-to-table`** — strand is now `+`, `-`, or
  `.`; `None`/`0` strands are no longer silently mapped to `+`.
- **`vcf-to-bed`** — guards against `rec.alts is None` (gVCF non-variant
  blocks) and symbolic alts (`<DEL>`, `<INS>`, breakends). Symbolic count
  is reported on stderr.
- **`vcf-to-table`** — INFO columns are now pulled from the VCF header so
  fields that only appear in later records aren't silently dropped.
- **`bam-to-bed`** — detects SAM vs BAM vs CRAM from the file extension
  (was hard-coded to `"rb"`) and iterates with `until_eof=True` so
  unsorted/unindexed BAMs no longer require a `.bai`.
- **`convert-alignment`** — no longer requires the substring `aln` in the
  filename.
- **`annotate-tree`** — tree is deep-copied before relabelling (library
  callers no longer see their input mutated); tree leaves missing from the
  alignment now emit a warning instead of crashing on `.index()`.
- **`.vcf.gz`** output filename collisions — `sample.vcf.gz` and
  `sample.vcf` no longer write to the same output file.
- All converters route progress / warning messages to **stderr** instead
  of stdout so piped output stays clean.
- `requirements.txt` — rewritten as proper pip-requirement lines (was a
  Python list literal that `pip install -r` couldn't parse); `bcbio-gff`
  added to `pyproject.toml` dependencies (every GFF converter imports it).
- Author email: `narhmadey@wisc.edu` (was `narhmadey@swisc.edu`).

### Changed
- Package `__version__` is now sourced from installed metadata via
  `importlib.metadata`, not hard-coded in three places.
- Minimum Python version bumped from 3.6 (EOL) to 3.9.
- Biopython / pandas / pysam dependency ranges now have upper bounds.
- README rewritten with real markdown tables, fenced code blocks, and a
  complete Testing section. Scripts Overview section actually exists now
  (the TOC previously linked to a non-existent anchor).

### Removed
- `reorganise_for_package.sh` — one-shot script whose purpose has been
  served since the package layout is stable.

## [0.1.2] — 2026-02-27

Initial public release with 21 converters (GFF3, GTF, BED, GenBank, FASTA,
FASTQ, QUAL, alignment formats, BAM, VCF, MAF, MUMmer delta, phylogenetic
trees).
