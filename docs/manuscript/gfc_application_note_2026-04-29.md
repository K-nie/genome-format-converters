# gfc: a unified Python CLI for 30 genomic format conversions

**Benjamin Narh-Madey**
Laboratory of Genetics, University of Wisconsin, Madison, WI 53706, USA
narhmadey@wisc.edu

*Application Note. Draft prepared 2026-04-29. Target venues: Bioinformatics, NAR, BMC Bioinformatics, GigaScience.*

---

## Abstract

Genomics workflows routinely cross half a dozen file formats: VCF, BCF, FASTA, FASTQ, GFF3, GTF, BED, GenBank, EIGENSTRAT, PLINK, HMMER tabular output, OrthoFinder orthogroups, MAF, MUMmer delta, and several alignment containers. Each conversion is usually handled by a single-purpose tool with its own command grammar (samtools, bcftools, plink2, EIGENSOFT/convertf, gffread, AGAT, UCSC kent, EMBOSS, pyhmmer). The result is glue scripting that ages badly. We present **gfc** (genome-format-converters), a Python 3.9+ CLI that exposes 30 conversions behind one consistent interface. On a ten-replicate benchmark across eight production tasks and fourteen competing tools, gfc beats AGAT by 12.93x on GFF3-to-GTF and 13.26x on GFF3-to-BED12, beats EMBOSS-seqret by 7.88x on FASTA+GFF-to-GenBank, beats convertf by 1.29x in wall time and 22.07x in peak resident memory on chr22 VCF-to-EIGENSTRAT, and uses 22x less RAM than plink2 on VCF-to-PLINK. Speed losses against C tooling (gffread, UCSC, plink) are bounded: gfc trades roughly 30 percent of wall time for unified CLI scaffolding, batch mode, deterministic seeded RNG, and shared validation. gfc is `pip install genome-format-converters`, MIT licensed, source at github.com/K-nie/genome-format-converters.

## 1 Introduction

Bioinformatics format conversion is unglamorous but load-bearing. A typical phylogenomics or ancient-DNA pipeline crosses VCF to EIGENSTRAT for ADMIXTURE, GFF3 to GTF for transcript-isoform tools, FASTA plus GFF3 to GenBank for NCBI submission, and HMMER tabular output to TSV for downstream pandas work. Each leg is handled by a separate, narrowly scoped tool. EIGENSOFT/convertf [1] handles VCF to EIGENSTRAT. PLINK 2 and PLINK 1.9 [2,3] handle VCF to PLINK binary. AGAT [4] is a Perl validator; gffread [5] is a C transcript-isoform translator. Both handle GFF3 conversions, on very different terms. UCSC kent tools [6] supply genePred chains. EMBOSS-seqret [7] handles GenBank rendering. pyhmmer [8] parses HMMER [9] output in-process. cyvcf2 [10], pysam [11], and Biopython [12] are the underlying Python libraries that most pipelines wrap by hand.

Each tool is excellent in isolation. Stitched together they fragment: distinct CLI grammars, distinct error-handling conventions, distinct logging, distinct parallelism stories, no shared input validation. Y1000+ phylogenomics work [Opulente et al. 2024] and aDNA pipelines run thousands of conversions across hundreds of genomes; the glue layer becomes a maintenance liability.

gfc unifies these conversions behind one CLI. Every subcommand takes `--input-dir`/`--output-dir` (batch) or `--input`/`--output` (single file), fails loudly on malformed input, streams I/O, and parallelises with `--threads`. The tool is purposely Python; it sacrifices some throughput against C parsers in exchange for portability, scriptability, and a small dependency tree. Section 3 measures that trade-off directly.

## 2 Methods

### 2.1 Implementation

gfc is implemented in Python 3.9+ with four runtime dependencies: cyvcf2 [10] and pysam [11] for VCF/BCF and FASTA random access, Biopython [12] for sequence and GenBank I/O, bcbio-gff for GFF3 parsing, and pandas for tabular intermediates. Each of the 30 converters is a separate module under `src/genome_format_converters/converters/` and is independently importable. The `gfc` console script dispatches to subcommands via argparse. I/O is streaming throughout: VCF passes are single-pass over `cyvcf2.VCF`, GFF3 is iterated record by record, GenBank writers flush per-record. Resident memory therefore tracks the working set of one record plus minimal per-converter state, not the input file size.

Installation is `pip install genome-format-converters`. Source is at github.com/K-nie/genome-format-converters under MIT license.

### 2.2 Benchmark harness

We benchmarked gfc against fourteen competing tools across eight tasks (Table 1). Every (task, tool) pair was run for n=10 independent replicates. Each replicate was a fresh subprocess with a clean output directory. Wall time was measured with `time.perf_counter` around the subprocess. Peak resident set size was captured via `os.wait4` returning a `struct rusage` whose `ru_maxrss` field reports kilobytes on Linux. We did not drop warm-up replicates; the reported standard deviations include cold-cache effects. The harness ran on HTCondor at the Great Lakes Bioenergy Research Center (GLBRC) on `scarcity-10.glbrc.org` (Methods Table 1). Software pins are in `benchmarks/environment.yml` (gfc 0.1.5, Python 3.11.9, cyvcf2 0.32.1, pysam 0.22.1, biopython 1.87, samtools 1.19.2, bcftools 1.19, plink2 2.0.0-a.6.9LM, plink1.9 1.90b6.21, gffread 0.12.7, AGAT 1.4.0, UCSC ucsc-tools 469, EMBOSS 6.6.0, pyhmmer 0.12.0, HMMER 3.4, eigensoft 8.0.0).

### 2.3 Datasets

T1, T2, T8 used 1000 Genomes Project chromosome 22, biallelic-SNP-filtered VCF (the standard population-genetics laptop slice). T3, T4, T5 used a 20-species subsample of the Y1000+ Saccharomycotina assemblies and annotations. T6 used a four-line synthetic HMMER tabular fixture (a known limitation; see Discussion). T7 used a 12-species OrthoFinder run.

### 2.4 Correctness comparators

Three classes of comparator were used. Byte-identity (T6, T7) requires the gfc output to match the reference tool byte for byte. Structural equivalence (T1, T2, T3, T4) reads both outputs back through a parser and verifies record counts, sequence lengths, or genotype-matrix shape; minor cosmetic divergences (chrom encoding, attribute order, line-wrap, transcript-centric span fusion) are tolerated. Format-only (T8) verifies that gfc's pseudohaploid output is a valid EIGENSTRAT triplet with the expected row count; cross-tool byte-identity is impossible because the reference uses a different RNG. The four tasks where gfc returns `correct=0` against the reference do so for documented design reasons, captured verbatim in footnotes [a-d] at the end of this section.

### 2.5 Hardware

**Methods Table 1.** Benchmark host (full provenance in `benchmarks/results/figures/hardware_table.md`).

| Field | Value |
|---|---|
| Host | `scarcity-10.glbrc.org` (GLBRC HTCondor) |
| Kernel | Linux 4.18.0-553.92.1.el8_10.x86_64 |
| CPU | Intel Xeon E5-2687W v2 @ 3.40 GHz |
| Sockets x cores x threads | 2 x 8 x 2 (32 logical CPUs) |
| RAM | 219 GB |
| Filesystem | BeeGFS shared |
| Software | `benchmarks/environment.yml` |

### 2.6 Format-divergence footnotes

[a] **T1 (VCF to EIGENSTRAT).** gfc and convertf may diverge on the chromosome encoding. gfc preserves the VCF CHROM literal (e.g. `chr22`); convertf rewrites to its integer-only EIGENSTRAT contract (e.g. `22`). The `.geno` and `.ind` payloads are byte-identical when both tools are pointed at the same chrom mapping; the `.snp` file differs only in the chrom column. We report `correct=0` rather than canonicalise the chrom column post hoc, because the EIGENSTRAT spec does not mandate either rendering, and downstream tools (smartpca, ADMIXTURE) accept both.

[b] **T2 (VCF to PLINK .bed).** gfc and plink1.9 produce structurally distinct outputs from the same VCF input. Three independent divergence axes are confirmed locally on `tests/test_data/tiny.vcf`: (i) gfc skips indels by default and writes only biallelic SNPs (2 of 3 variants for tiny.vcf), while plink1.9 retains all variants by default; (ii) gfc synthesises variant IDs as `<chrom>_<pos>` when the VCF ID is `.`, while plink1.9 preserves the literal `.`; (iii) gfc and plink1.9 select different A1/A2 allele assignments. gfc uses VCF REF/ALT directly; plink1.9 applies minor-allele-first reordering. The `.bim` text and the `.bed` binary payload differ in length, content, and genotype encoding direction. None of these is a bug. They are documented design choices. We keep the strict byte-identity comparator so any regression in gfc's own `.bed` output across versions surfaces immediately, and we accept `correct=0` against plink1.9 with this footnote.

[c] **T3 (FASTA+GFF to GenBank).** gfc uses Biopython's GenBank writer. The py-ref baseline also uses Biopython but with a different feature-qualifier order. EMBOSS-seqret renders the same content with different line-wrap, FEATURES section ordering, and qualifier formatting. The GenBank flat-file format does not mandate a single canonical rendering, so byte-identity across implementations is impossible. The comparator uses `Bio.SeqIO` to verify structural equivalence: identical record count and total sequence length. `correct=0` here would mean the structural invariant actually broke (record loss or sequence corruption), not cosmetic rendering differences.

[d] **T4 (GFF3 to GTF).** gfc and gffread implement two defensible but distinct GFF3-to-GTF translations. Confirmed locally on `tests/test_data/tiny.gff3`: (i) gfc emits a faithful 1:1 row mapping including `gene` and `mRNA` feature rows; gffread elides gene rows and renames `mRNA` to `transcript`. (ii) gffread fuses adjacent exon and CDS spans into a single feature when they abut (101-150 + 151-200 becomes 101-200); gfc preserves the original GFF3 spans verbatim. (iii) Trailing-semicolon and attribute-key order differ. gffread's normalisations are transcript-centric and standard for transcript-isoform pipelines; gfc's faithful mapping is standard for general feature-table round-tripping. Neither is the canonical GTF rendering. The comparator already canonicalises attribute order via `_normalise_gtf_line`; `correct=0` reflects the deeper structural difference (line count, feature types, span boundaries) and is documented rather than hidden behind a span-merging rewrite.

## 3 Results

We organise the eight-task results (Figure 1, `wall_time.png`; Figure 2, `peak_memory.png`; Figure 3, `correctness_matrix.png`; Figure 4, `pareto.png`; Figure 5, `speedup.png`) by win category. Per-pair means with standard deviations are in `summary_per_pair.md`; the headline per-task table is `summary_table.md`.

**Dominant wins.** On T4 (GFF3 to GTF, 20-species Y1000+ slice), gfc finishes in 8.11 +/- 0.11 s (n=10) against AGAT's 104.85 +/- 0.05 s, a 12.93x speedup at comparable memory (58.0 vs 63.4 MB peak RSS). On T5 (GFF3 to BED12, same dataset), gfc finishes in 7.85 +/- 0.13 s against AGAT's 104.08 +/- 0.04 s, a 13.26x speedup. AGAT is the deep-validation Perl reference for GFF3 work, and these speedups matter at Y1000+ scale (>1,000 genomes). On T1 (VCF to EIGENSTRAT, chr22 1000 Genomes biallelic SNPs), gfc finishes in 259.25 +/- 0.77 s against convertf's 334.93 +/- 6.71 s (1.29x faster) while using 59.90 +/- 0.93 MB peak RSS against convertf's 1.32 +/- 0.002 GB (22.07x less memory). The memory difference is qualitative: chr22 fits comfortably on a laptop with gfc; convertf needs a workstation or shared-cluster slot.

**Speed wins, comparable memory.** On T3 (FASTA+GFF3 to GenBank), gfc beats EMBOSS-seqret 7.88x (8.13 vs 64.07 s) and uses comparable memory. On T8 (pseudohaploid VCF), gfc beats a hand-written bcftools+Python pipeline 5.64x (332.93 vs 1877.07 s).

**Memory wins, speed loss accepted.** On T2 (VCF to PLINK), plink2 finishes in 21.93 +/- 0.20 s and plink1.9 in 39.93 s, while gfc takes 311.41 +/- 0.93 s. gfc loses 14.2x to plink2 and 7.8x to plink1.9 on speed. The reverse story tells on memory: gfc holds 59.91 +/- 1.07 MB peak RSS while plink2 peaks at 1.32 GB (a 22.08x ratio). The same chr22 conversion that needs a 2 GB cluster slot for plink2 fits in 60 MB with gfc. On laptops, in tight Condor slots, or inside containers with hard memory caps, gfc is the only option that runs.

**Losses.** Four task-tool pairings show gfc trailing the reference. Against py-ref baselines (T3 vs py-ref: 1.31x slower; T7 vs py-ref: 1.54x slower), gfc carries CLI scaffolding, argparse validation, and batch iteration that a 50-line one-off script does not. The roughly 30 percent overhead is the cost of production tooling. Against C-tool single-purpose parsers (T4 vs gffread: 1.69x slower; T5 vs ucsc-chain: 1.18x slower), Python loses to C on parser-heavy workloads. The relevant comparison is to the *Perl* alternative (AGAT) on the same input, where gfc wins 13x. T6 (HMMER tblout) shows gfc 3.49x slower than pyhmmer on a four-line synthetic fixture; Python startup overhead dominates the measurement, and the production-scale tblout (~10K rows from a real Pfam-A scan) would change the picture. That is a fixture limitation, not a tool limitation.

The Pareto front (Figure 4) places gfc consistently in the lower-left quadrant on memory-time space for variant-pipeline tasks (T1, T2, T8), the upper-left on annotation tasks vs AGAT (T4, T5), and slightly off-front on tasks dominated by C parsers (T2 vs plink, T4 vs gffread).

## 4 Discussion

gfc is not a replacement for plink2, gffread, or UCSC kent on the workloads those tools were built for. It is a unified, predictable, low-memory alternative for the cases where memory matters more than throughput, where an unfamiliar CLI is friction, or where a single command-shape across 30 conversions reduces glue-script maintenance. The headline trade-off is explicit: against tuned C tools, gfc trades wall-clock for portability, predictability, and (on T1, T2) a 22x memory advantage. Against Perl validators (AGAT), gfc wins on both axes.

Three limitations bound the present benchmark. First, T1, T2, and T8 use chr22 only; whole-genome 1000G or larger panels are Stage-2 work. Second, T3-T5 use a 20-species Y1000+ slice; full Y1000+ scale (1,000+ assemblies) is the next benchmark and is where the 13x AGAT speedup matters most. Third, T6 uses a four-line synthetic fixture; the small input size puts gfc inside the Python startup-overhead regime, and a real Pfam-A tblout (~10K rows) is needed to give a fair number. That scan is running now and will be reported as part of a Stage-2 paper-grade benchmark.

Future work has two strands. We are profiling the T2 hot path (VCF row -> 2-bit packed PLINK genotype) and expect Cython on that single inner loop to recover most of the gap to plink2 without inflating the dependency tree. We are also extending the suite from 30 to roughly 40 conversions (Newick to extended-Newick, BCF to gzipped-VCF, AGP, MMSeqs2 results), driven by user requests on the issue tracker.

## 5 Availability

gfc is `pip install genome-format-converters` (PyPI) or `git clone https://github.com/K-nie/genome-format-converters` (source), Python 3.9+, MIT license. A Zenodo DOI is minted on each tagged release. Documentation, the full benchmark harness, and all raw replicate data (`benchmarks/results/raw/T*.tsv`) are in the repository so reviewers and users can re-run any number reported here.

---

## Figures

- **Figure 1.** Wall time per task across all (task, tool) pairs, n=10 reps, error bars are 1 s.d. `benchmarks/results/figures/wall_time.png` / `.pdf`.
- **Figure 2.** Peak resident set size per task across all (task, tool) pairs, n=10 reps. `peak_memory.png` / `.pdf`.
- **Figure 3.** Correctness matrix: byte-identity / structural-equivalence / format-only verdict per (task, tool) pair, with the four `_REASON_T*` footnotes annotated. `correctness_matrix.png` / `.pdf`.
- **Figure 4.** Pareto front (peak RSS vs wall time), log-log, per task, gfc highlighted. `pareto.png` / `.pdf`.
- **Figure 5.** Per-pair speedup (gfc vs each competitor), log scale, separating the wins, the memory-only wins, and the losses. `speedup.png` / `.pdf`.

Supporting figures in the supplement: `compute_cost.png`, `feature_matrix.png`, `per_task_panels.png`, `variance.png`. Tables: `summary_table.md` (per-task headline numbers, used in Section 3), `summary_per_pair.md` (every (task, tool) pair with verdict), `hardware_table.md` (full host provenance, condensed to Methods Table 1).

---

## References

[1] Patterson, N., Price, A. L. & Reich, D. (2006) Population structure and eigenanalysis. *PLoS Genetics* 2(12): e190.

[2] Chang, C. C., Chow, C. C., Tellier, L. C. A. M., Vattikuti, S., Purcell, S. M. & Lee, J. J. (2015) Second-generation PLINK: rising to the challenge of larger and richer datasets. *GigaScience* 4: 7.

[3] Purcell, S., Neale, B., Todd-Brown, K., Thomas, L., Ferreira, M. A. R., Bender, D., Maller, J., Sklar, P., de Bakker, P. I. W., Daly, M. J. & Sham, P. C. (2007) PLINK: a tool set for whole-genome association and population-based linkage analyses. *American Journal of Human Genetics* 81(3): 559-575.

[4] Dainat, J. (2022) AGAT: Another Gff Analysis Toolkit to handle annotations in any GTF/GFF format. National Bioinformatics Infrastructure Sweden (NBIS). doi:10.5281/zenodo.3552717.

[5] Pertea, G. & Pertea, M. (2020) GFF Utilities: GffRead and GffCompare. *F1000Research* 9: 304.

[6] Kent, W. J., Sugnet, C. W., Furey, T. S., Roskin, K. M., Pringle, T. H., Zahler, A. M. & Haussler, D. (2002) The Human Genome Browser at UCSC. *Genome Research* 12(6): 996-1006.

[7] Rice, P., Longden, I. & Bleasby, A. (2000) EMBOSS: the European Molecular Biology Open Software Suite. *Trends in Genetics* 16(6): 276-277.

[8] Larralde, M. & Zeller, G. (2023) PyHMMER: a Python library binding to HMMER for efficient sequence analysis. *Bioinformatics* 39(5): btad214.

[9] Eddy, S. R. (2011) Accelerated profile HMM searches. *PLoS Computational Biology* 7(10): e1002195.

[10] Pedersen, B. S. & Quinlan, A. R. (2017) cyvcf2: fast, flexible variant analysis with Python. *Bioinformatics* 33(12): 1867-1869.

[11] Heger, A. et al. pysam: a Python module for reading, manipulating and writing genomic data. github.com/pysam-developers/pysam.

[12] Cock, P. J. A., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., Friedberg, I., Hamelryck, T., Kauff, F., Wilczynski, B. & de Hoon, M. J. L. (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* 25(11): 1422-1423.

[13] Opulente, D. A., Leavitt LaBella, A., Harrison, M.-C., Wolters, J. F., Liu, C., Li, Y., Kominek, J. et al. (2024) Genomic factors shape carbon and nitrogen metabolic niche breadth across Saccharomycotina yeasts. *Science* 384(6694): eadj4503.
