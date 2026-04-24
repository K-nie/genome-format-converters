# Methods-paper plan — *genome-format-converters* (`gfc`)

Working title: **"gfc: a uniform command-line toolkit for common bioinformatics file-format conversions"**

Target audience: practising bioinformaticians running routine format conversions on mid-to-large genomic datasets (e.g. Y1000+-scale yeast comparative genomics, population-scale VCFs, orthogroup-driven phylogenomics).

Author: **Benjamin Narh-Madey** (lead), Hittinger Lab, University of Wisconsin-Madison.

---

## 1. Ranked journal shortlist

Four journals are realistic; others are either scope-mismatched or too high-bar for a utility tool. **Recommendation: publish in JOSS first, then submit an expanded Applications Note to *Bioinformatics* once benchmarking is mature.**

### Tier A — submit first

| Journal | Why it fits | Page / word limit | Fee (2026) | Time-to-publish | Acceptance odds |
|---|---|---|---|---|---|
| **Journal of Open Source Software (JOSS)** | Built for research software. Review happens on GitHub; reviewers install and run the tool. Uniform CLI + tests + docs + license satisfy almost every JOSS checklist item already. CrossRef DOI, Scopus indexed. | 250-1000 words paper + references | Free | 2-3 months | **Very high (~90 %)**. Main risk: reviewer requests more examples or API-stability clarifications. |
| **BMC Bioinformatics — Software Articles** | Long-tail home for bioinformatics CLIs that don't pitch a new algorithm. Open access, wide readership. | 4-8k words, unlimited figures | USD 2,490 | 4-6 months | **High**. Reviewers more forgiving on novelty than *Bioinformatics* reviewers. |

### Tier B — aim after benchmarking is strong

| Journal | Why it fits | Page / word limit | Fee | Time | Odds |
|---|---|---|---|---|---|
| **Bioinformatics Advances (Oxford)** | Sister journal of *Bioinformatics*, open-access, slightly lower bar. Accepts tool papers with clear benchmarks. | ~3,000 words + 2 figures | USD 2,480 | 4-6 months | Medium-high. |
| **Bioinformatics — Applications Note** | Prestige track (IF ≈ 6.9). 2-page format, very concise. Needs either a genuinely novel subcommand (e.g. the EIGENSTRAT pseudohaploid path) or a convincing benchmark showing speed / memory / correctness wins vs established tools. | ~1,800 words + 1 fig + 1 table | USD 2,570 (OA optional) | 4-6 months + possibly revisions | Medium. Reviewers will ask "what does gfc do that seqkit + bcftools + samtools + AGAT don't?" — the answer must be benchmarked. |

### Explicitly not a fit

- **PLoS Computational Biology, Nature Methods, Genome Biology** — utility converters don't clear the novelty bar. Don't submit here.
- **Molecular Biology and Evolution (MBE)** — scope is evolutionary methods, not conversion utilities. Only viable if the paper is *about* an evolutionary analysis (e.g. a Y1000+ ADMIXTOOLS pipeline using gfc as its spine), not about gfc itself.
- **GigaScience** — possible but they now prefer large datasets + tools bundled together. If the paper doubles as a Y1000+ EIGENSTRAT release, re-consider.

---

## 2. What to say where — audience tailoring

| Venue | Lead with… | Minimise… |
|---|---|---|
| JOSS | Statement of need, ecosystem context (bcftools, seqkit, AGAT), summary of the 30 subcommands, link to docs + tests + CI. | Benchmarks (not required). |
| BMC Bioinformatics | Uniform-CLI thesis + benchmarks on real Y1000+ data + feature matrix vs competitors. | Algorithmic novelty (this isn't a new algorithm). |
| Bioinformatics Advances | Pseudohaploid EIGENSTRAT for yeast + orthogroup-to-FASTA + HMMER tblout-to-TSV as the "glue" that wasn't in any single prior tool. Benchmarks. | Broad scope pitches ("does everything"). |
| Bioinformatics (App Note) | Strongest single benchmark (probably VCF → EIGENSTRAT / PLINK at Y1000+ scale) + one cleanly-illustrated use case (a yeast admixture pipeline end-to-end). | Feature enumeration. |

---

## 3. Paper outline (Applications-Note form, ~1800 words)

### Abstract (200-240 words)
- **Motivation.** Modern bioinformatics pipelines chain together file-format conversions across FASTA/FASTQ, GFF3/GTF, VCF, BAM, EIGENSTRAT, PLINK, Stockholm, HMMER output, OrthoFinder orthogroups, and phylogenetic trees. Each chain is glued with ad-hoc scripts because no single tool covers all of these, and conventions differ across tools (bcftools vs samtools vs seqkit vs AGAT vs PLINK).
- **Results.** `gfc` is a uniform command-line toolkit with 30 subcommands, sharing a consistent `--input-dir`/`--output-dir`/`--pattern`/`--force`/`--threads` interface. It includes first-class support for EIGENSTRAT and pseudohaploid output (including per-file reproducible RNG seeding), PLINK binary `.bed`/`.bim`/`.fam` output, and OrthoFinder orthogroup expansion — which are otherwise written as one-off scripts. Benchmarks on Y1000+ (1154 Saccharomycotina genomes) show gfc matches or beats the speed of existing single-purpose tools while covering broader format coverage under one interface.
- **Availability and implementation.** `gfc` is implemented in Python (≥3.9), tested on Linux, macOS, and a pure-Python subset on Windows. Source under MIT at https://github.com/K-nie/genome-format-converters; install via `pip install genome-format-converters` or `conda install -c bioconda genome-format-converters`.
- **Contact.** narhmadey@wisc.edu

### 1. Introduction (~300 words)
- File-format interconversion is one of the highest-friction parts of bioinformatics pipelines.
- Existing tools are scattered and have different CLI conventions — cite seqkit (Shen et al. 2016), bcftools/samtools (Danecek et al. 2021), PLINK 2.0 (Chang et al. 2015), AGAT (Dainat et al. 2020), EMBOSS (Rice et al. 2000), ADMIXTOOLS/convertf (Patterson et al. 2012).
- Three concrete pain points we target: (1) pseudohaploid VCF → EIGENSTRAT with reproducible RNG, commonly re-implemented from `pileupCaller` (Schiffels et al.); (2) OrthoFinder orthogroup expansion to per-OG FASTA, scripted in every lab; (3) HMMER tblout → TSV, typically an awk one-liner that loses the description-of-target column.
- Our contribution: one uniform CLI across 30 converters, with loud-failure semantics, stderr-only logging, composable flags, and reproducible RNG.

### 2. Implementation (~400 words)
- Pure-Python with biopython, pysam, bcbio-gff, pandas; no external C binaries beyond htslib (via pysam).
- Architecture: one module per converter under `converters/`, uniform `batch_convert(input_dir, output_dir, …)` entry point, shared helpers (`_common.py`, `_eigenstrat_common.py`) for logging, compound-suffix handling, map loading, polarisation, site filtering, PLINK packing.
- Feature matrix (Figure 1 / Table 1): checkbox grid of gfc vs seqkit, samtools/bcftools, AGAT, convertf, plink, EMBOSS, custom-script.
- Windows: pysam isn't on PyPI for Windows; gfc's pure-Python subset (~20/30 subcommands) still works, VCF/BAM paths emit a helpful error pointing at WSL or conda-forge.

### 3. Results / Benchmarks (~500 words)
- Tasks benchmarked (see `docs/BENCHMARK_PLAN.md`).
- Figure 2: wall-time per task, scaled against input size.
- Figure 3: peak memory.
- Table 2: correctness vs reference tools (byte-identical; row-counted; spot-checked).
- Case study: one-screen yeast ADMIXTOOLS pipeline on a Y1000+ slice — `bcftools merge` + `gfc vcf-to-eigenstrat --transversions-only --min-maf 0.05 --chrom-map chroms.tsv --pop-map pops.tsv` + `smartpca` + `qpAdm`. Highlight the single-tool cleanliness vs the equivalent multi-tool sh script.

### 4. Discussion (~250 words)
- **Limitations.** HAL ↔ MAF needs halTools (C binary) — out of scope. BED12 → GFF3 inverse not yet implemented. Tree support-value conversion (IQ-TREE ↔ RAxML ↔ FigTree) deferred.
- **Future work.** Per-subcommand ProcessPoolExecutor parallelism (currently serial). Streaming VCF conversion for datasets that don't fit in memory. Native Windows support when pysam lands.
- **Design philosophy.** Uniform CLI, loud-failure, reproducibility-by-default — these matter for real pipelines more than any single feature.

### 5. Acknowledgements + Funding + Availability + Ethics
- Hittinger Lab colleagues (specific reviewers/testers — fill in).
- Funding placeholder — confirm NIH grant numbers with Chris Hittinger.
- Data availability: test fixtures under `tests/test_data/`. Benchmark datasets (Y1000+ slice) deposited on Zenodo with DOI, linked in paper.
- Code availability: GitHub + PyPI + bioconda (pending). CITATION.cff provides Zenodo DOI anchor.

---

## 4. Figures and tables (planned)

| # | Type | Content |
|---|---|---|
| **Fig 1** | Architecture schematic | gfc sitting between raw genomic files (FASTA/VCF/GFF/BAM) and downstream tools (smartpca / ADMIXTOOLS / IGV / IQ-TREE / iqtree / OrthoFinder / PLINK). Shows the 30 subcommands grouped into 5 categories. |
| **Fig 2** | Feature matrix | Rows: tools (gfc, seqkit, bcftools, samtools, AGAT, convertf, plink, EMBOSS, pysam-based custom). Cols: 20-25 representative tasks. Checkmarks with colour shading for "first-class" vs "possible with scripting". |
| **Fig 3** | Benchmark bars | Wall-time across 6-8 tasks (see benchmark plan). Grouped bars: gfc vs each competitor. Error bars for 5 replicates. |
| **Fig 4** | Benchmark bars | Peak RSS memory, same layout as Fig 3. |
| **Fig 5** | Case study | Yeast admixture pipeline end-to-end, one-screen command sequence + resulting PCA/f4 figure. |
| **Table 1** | Dataset descriptions | Y1000+, 1000G chr22 slice, Pfam-A scan output, Ensembl Fungi GTFs — sizes, record counts, source. |
| **Table 2** | Correctness | Per-task: byte-identical vs ref (✓/✗), row counts, any tolerances, chosen reference tool. |

---

## 5. Selling points (what makes this publishable)

1. **Uniform CLI over 30 format conversions** — nobody else has this breadth under a single interface with consistent flags.
2. **Pseudohaploid EIGENSTRAT with reproducible seeding** — per-file RNG seeding (XOR CRC32) is a real engineering improvement over naive global-seed implementations.
3. **OrthoFinder orthogroup expansion + HMMER tblout parsing** — these two alone will resonate with every comparative-genomics group; they're currently written as disposable awk / Python every time.
4. **Ancestral-allele polarisation via `--ancestral-fasta` / `--info-aa`** — right conceptual primitive for population-genetics work that convertf doesn't ship.
5. **Loud-failure semantics + stderr-only logging** — not sexy but is the correct behaviour for pipelines.

---

## 6. Timeline

| Milestone | Target date | Owner |
|---|---|---|
| JOSS paper draft (`paper.md` + `paper.bib`) in repo | **+3 weeks** from today | BNM |
| Benchmark harness + figures (see benchmark plan) | **+6 weeks** | BNM |
| Internal Hittinger-lab review | **+7 weeks** | Chris + one postdoc |
| JOSS submission | **+8 weeks** | BNM |
| JOSS review → published | JOSS timeline | — |
| App-Note draft building on JOSS paper | **+3 months** | BNM |
| Applications Note submission (Bioinformatics or Bioinformatics Advances) | **+4 months** | BNM |

---

## 7. Authorship / acknowledgements template

- **Benjamin Narh-Madey** — design, implementation, benchmarks, writing.
- **Chris Todd Hittinger** (PI) — supervision, writing contributions.
- Acknowledge: Hittinger-lab members who test / review the tool before submission; list them by name on the acknowledgement line.
- Funding line: confirm active grants with Chris (NIH R01, NSF DEB, etc.).
- Competing interests: none declared.

---

## 8. Risks / likely reviewer pushback

| Risk | Mitigation |
|---|---|
| "This is just a wrapper around biopython/pysam." | Emphasise the uniform CLI contract, loud-failure semantics, reproducible RNG, EIGENSTRAT/PLINK first-class support. Show pseudohaploid, orthogroup, HMMER — these are real engineering contributions. |
| "Why not add support to AGAT / seqkit?" | Scope argument: gfc's thesis is *uniform interface across formats*, not best-in-class for any single format. AGAT is GFF-only; seqkit is FASTA/FASTQ-only. |
| "Your benchmarks cherry-pick tasks that favour gfc." | Publish the benchmark harness + Docker image; include at least one task where gfc is slower but correct. |
| "No Windows support for VCF / BAM." | Document clearly; most bioinformatics users are Linux/macOS (or WSL); pysam limitation is upstream. |
| "Does it scale to population-size VCFs?" | Document memory / time scaling curves. Add a note on the two-pass in-memory approach for EIGENSTRAT. |

---

## 9. Pre-submission checklist (before JOSS)

- [ ] `paper.md` (JOSS format) + `paper.bib` committed to `paper/` directory.
- [ ] `CITATION.cff` pointing at Zenodo DOI.
- [ ] All 30 subcommands documented in README with one example each.
- [ ] CI green on all three OS matrix jobs (Linux/macOS full + Windows subset).
- [ ] pytest suite ≥ 30 tests, coverage report.
- [ ] Bioconda recipe lodged (or at least submitted to bioconda-recipes).
- [ ] Zenodo archived release.
- [ ] License: MIT (already present).
- [ ] Contribution guide (already present).
- [ ] 3 external users who have installed and run gfc, for JOSS reviewer pool diversity.
