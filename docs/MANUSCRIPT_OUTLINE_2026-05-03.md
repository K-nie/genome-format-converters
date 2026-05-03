Author: Benjamin Narh-Madey
Date: 2026-05-03
Target venue: Bioinformatics (Oxford Academic), Original Paper

# Manuscript outline: gfc-bench

This skeleton is the section-by-section build plan for the regular Original Paper submission. Word budgets follow the Bioinformatics author guidance for full-length papers (~5000 words main text, abstract <=250). Section bullets list what goes in each block; they are not draft prose. Every quantitative claim cites its anchor file or TSV row.

## Working title (three options to pick from)

1. **gfc: a single-process Python suite for cross-domain genome-file conversion with bounded memory and per-conversion correctness reporting**
2. **Cross-domain genome-file conversion in one Python process: gfc, with constant memory and a four-state correctness matrix**
3. **gfc: dependency-light, memory-bounded format conversion across sequence, alignment, variant, and annotation domains**

Title (1) is the safer default. It signals the architectural choice (single-process Python, no shell-outs), the structural property (bounded memory), and the methodological contribution (per-conversion correctness reporting) without leaning on a format-count claim that biopython.convert beats us on (FRAMING_DECISIONS Decision 1).

## Sections and word budgets

### Abstract: 250 words (hard cap)

Five-paragraph structure. Drafted separately in `ABSTRACT_DRAFT_2026-05-03.md`. Anchor numbers: 22.09x plink2 RSS [21.70-22.46]; CV 2.04% across 80 runs; install footprint 163 MB / 32 deps vs BioConvert 1186 MB / 201 deps; speed losses 14.2x to plink2, 3.5x to pyhmmer, 1.7x to gffread.

### Introduction: 600 words

Paragraph 1 (~120 words). Format conversion as the connective tissue of bioinformatics pipelines: every multi-tool workflow crosses at least one format boundary. Cite Mölder et al. 2021 (Snakemake, F1000Research 10:33) and Di Tommaso et al. 2017 (Nextflow, Nat Biotech 35:316) for the workflow context.

Paragraph 2 (~150 words). Two existing tool families. Domain-deep utilities (samtools, bcftools, AGAT, gffread, convertf) are fast and correct within one format family but require users to chain four to six binaries to span a typical workflow. Wrapper frameworks (BioConvert per Caro et al. 2023, NAR Genomics & Bioinformatics 5:lqad074; EMBOSS-seqret per Rice et al. 2000, Trends Genet 16:276; biopython.convert per brinkmanlab GitHub) aggregate external binaries through a unified interface but inherit the dependency closure of the binaries they wrap. State both honestly; do not invent a "third architecture."

Paragraph 3 (~150 words). The operational gap. Embedding a converter into a containerised pipeline pays for the wrapper's dependency closure even when only one conversion is invoked. BioConvert's conda recipe pulls 201 packages and 1186 MB; biopython.convert ships 32 packages and 118 MB but covers only sequence and structural formats (CLAIM_VERIFICATION). Memory profiles across these tools vary by orders of magnitude on the same input, and no tool reports per-conversion correctness against a published reference in a way that distinguishes byte divergence from structural equivalence (STATS_AUDIT issue 7).

Paragraph 4 (~130 words). What we do. gfc is a single Python package that performs ~30 conversions natively across sequence, alignment, variant, and annotation domains. It calls no external binaries at runtime. We benchmark gfc against citable per-task references on eight representative tasks (T1-T8) and against BioConvert as a cross-task wrapper baseline. We report wall time, peak resident set size, install footprint, and a four-state correctness matrix that distinguishes byte equivalence, structural equivalence, comparator-limitation, and disagreement. We disclose the per-task speed losses without apology and frame the contribution around constant memory across heterogeneous inputs.

Paragraph 5 (~50 words). Roadmap of the paper.

### Methods: 1500 words

#### Architecture (~250 words)

- Single-package Python design. Native dispatch table mapping (input format, output format) to a handler function. No `subprocess.run(["samtools", ...])` calls anywhere in the runtime path (verifiable in source).
- Handler implementations build on pysam (htslib bindings), Biopython (sequence I/O), and pure-Python parsers for EIGENSTRAT, PLINK, HMMER tblout, and orthogroup formats.
- Common CLI surface: `gfc <in> <out> --from FMT --to FMT`, with stdin/stdout streaming, magic-byte gzip detection, and uniform exit codes across all 30 conversions.
- Per-output provenance header: gfc version, input SHA-256, conversion direction, command-line arguments, UTC timestamp.
- Cite Cock et al. 2009 (Biopython, Bioinformatics 25:1422) and the pysam GitHub release for the underlying libraries.

#### Benchmark design (~400 words)

- Eight tasks chosen to span the four bioinformatics domains and to admit at least one citable reference per task.
  - T1: VCF -> EIGENSTRAT triplet, reference EIGENSOFT convertf (Patterson et al. 2006, PLoS Genet 2:e190).
  - T2: VCF -> PLINK binary, references PLINK 2.0 and PLINK 1.9 (Chang et al. 2015, GigaScience 4:7).
  - T3: paired FASTA + GFF3 -> GenBank, reference EMBOSS seqret (Rice et al. 2000) and (Phase 4 expansion) table2asn (Sayers et al. 2023, NAR 51:D141) plus EMBLmyGFF3 (Norling et al. 2018, BMC Res Notes 11:584).
  - T4: GFF3 -> GTF, references gffread (Pertea & Pertea 2020, F1000Res 9:304) and AGAT (Dainat 2022 Zenodo doi:10.5281/zenodo.3552717).
  - T5: GFF3 -> BED12, references UCSC `gff3ToGenePred + genePredToBed` (Kent 2002, Genome Res 12:996) and AGAT.
  - T6: HMMER tblout -> TSV, reference pyhmmer (Larralde & Zeller 2023, Bioinformatics 39:btad214).
  - T7: OrthoFinder Orthogroups.tsv -> per-OG FASTA, reference OrthoFinder (Emms & Kelly 2019, Genome Biol 20:238).
  - T8: VCF -> pseudohaploid EIGENSTRAT, reference bcftools consensus -H (Danecek et al. 2021, GigaScience 10:giab008) and (Phase 4) ANGSD --doHaploCall (Korneliussen et al. 2014, BMC Bioinformatics 15:356).
- Datasets per LITERATURE_COMPARATORS table: 1000G phase 3 chr22 for T1/T2/T8, GENCODE v45 chr22 for T4/T5, RefSeq E. coli K-12 for T3, Pfam-A.hmm vs UniRef50 slice for T6, Quest for Orthologs reference proteomes for T7. Cite Frankish et al. 2023 (NAR 51:D942) for GENCODE; Wagner et al. 2022 (Cell Genomics 2:100128) for GIAB HG002 truth set.
- Hardware and harness: single Intel Xeon E5-2687W v2 Condor slot on GLBRC scarcity-10. Wall measured by `time.perf_counter`; peak RSS by `os.wait4` `ru_maxrss` with platform normalisation. Cite Weber et al. 2019 (Genome Biol 20:125) for the benchmark methodology guideline.

#### Statistical analysis (~250 words)

- n=10 replicates per (task, tool), sequential on the same Condor slot. Page cache warm after replicate 1; this caveat appears in the cold-cache footnote.
- Pairwise comparisons by paired Wilcoxon signed-rank test on log-transformed wall and peak RSS, paired by replicate index.
- Multiple-comparison correction: Benjamini-Hochberg at q=0.05 across the 24 (task, competitor) pairs in the primary panel.
- Effect sizes: median ratio of competitor to gfc with 95% bootstrap CI (10 000 resamples per pair) plus Cliff's delta.
- Coefficient of variation reported per (task, tool); n column reported per cell to flag any cell with n<10.
- Pre-registration acknowledgement: protocol finalised at git commit `<hash>`; numbers reported are from the single n=10 production run on `<date>`.
- Implementation: SciPy 1.13 (`scipy.stats.wilcoxon`, `scipy.stats.bootstrap`, `scipy.stats.false_discovery_control(method='bh')`). Code at `benchmarks/analyze/stats_pairs.py`.

#### Four-state correctness matrix (~600 words, drafted separately as `METHODS_4STATE_CORRECTNESS_2026-05-03.md`)

- Frames the matrix as a generalisable benchmarking norm for converter tools. Sections covered in the dedicated draft. This subsection of the manuscript Methods inlines the dedicated draft verbatim or by reference.

### Results: 1200 words

#### Cross-domain coverage (~200 words)

- Figure 1 (capability matrix): gfc spans sequence, alignment, variant, and annotation domains under a single CLI. biopython.convert covers more raw formats but only in sequence and structural domains. BioConvert covers more formats and more domains via external binaries. EMBOSS seqret covers more sequence formats but no modern alignment or variant family.
- Frame as cross-domain coverage in a single Python process, not as raw count.

#### Wall-time and memory (~400 words)

- Figure 2 (wall-time bars per task): show all 24 (task, tool) cells. Annotate the 5 clean wins (T1 vs convertf 1.29x [1.27-1.31]; T3 vs EMBOSS-seqret 7.95x [7.71-8.07]; T4 vs AGAT 12.96x [12.80-13.03]; T5 vs AGAT 13.29x [13.13-13.44]; T8 vs bcftools-pyref 5.58x [5.37-5.82]). Annotate the speed losses without apology (T2 vs plink2 14.2x slower; T6 vs pyhmmer 3.5x; T4 vs gffread 1.7x).
- Figure 3 (peak RSS bars per task, log y): gfc holds 58-61 MB across all 8 tasks. Reference C tools that mmap the input scale to 1322 MB. Annotate the constant-RSS band visibly.
- Headline RAM number: gfc 22.09x lower RSS than plink2 on T2 [95% CI 21.70-22.46]. Anchor: T2_vcf_to_plink.tsv per-rep RSS confirmed in STATS_AUDIT issue 5.
- T1 RAM: gfc 22.0x lower than convertf [21.77-22.48]. T2 plink1.9 RAM: gfc 1.18x lower [1.14-1.20]. T4 vs AGAT RAM: 1.09x lower [1.08-1.11]. T5 vs AGAT RAM: 1.11x lower [1.09-1.11].
- 10 of 24 paired Wilcoxon clean wins at q=0.05 (BH-FDR corrected). State the count plainly; do not round up.

#### RSS stability across tasks (~200 words)

- Figure 4 (Pareto wall vs RSS): gfc clusters in a horizontal band at 58-61 MB across all tasks. Convertf and plink2 sit in the high-RSS quadrant. Annotate the Pareto frontier explicitly.
- The RSS-stability claim: gfc mean 58.78 MB across 80 runs (8 tasks x 10 reps), CV 2.04%, range 56.88-61.31 MB. Anchor: brief, R3 numbers (locked); FRAMING_DECISIONS line 50.
- Frame as a structural property (Python streaming implementation does not mmap the input), not a per-task accident.

#### Install footprint (~150 words)

- Table 2 (R2 footprint): pip install gfc -> 163 MB / 32 deps; conda install bioconvert -> 1186 MB / 201 deps; pip install biopython.convert -> 118 MB / 32 deps.
- biopython.convert is 27% lighter than gfc with the same dependency count but covers only sequence and structural formats.
- BioConvert is 7.3x heavier on disk and 6.3x more dependencies, covers more formats but via external binaries.
- The wrapper-tax measurement at runtime: BioConvert T4 at 139s vs gffread T4 at 3.2s (44x slower for byte-identical output). Anchor: brief, R2 wrapper-tax number.

#### Per-task correctness (~250 words)

- Figure 5 (four-state correctness matrix): byte-equal (green), structurally-equivalent (blue, with footnote), comparator-limitation (grey), disagreement (red, none observed).
- Per-task verdicts: T1 structural (chrom encoding); T2 structural (indel handling, A1/A2 ordering); T3 structural (GenBank rendering); T4 structural (mRNA-to-transcript span fusion); T5 grey pending comparator fix (BED12 vs BED6 column-count mismatch); T6 byte-equal; T7 byte-equal; T8 format-only equal (cross-RNG byte equivalence intentionally out of scope).
- Aggregate: 7 of 8 at structural-equivalence-or-better.
- Anchor: STATS_AUDIT issue 7; AUDIT lines 117-129 for visualisation fix; `benchmarks/check_correctness.py:_REASON_T1` through `_REASON_T4` for per-task definitions.

### Discussion: 700 words

Paragraph 1 (~150 words). What gfc does. Restate the headline numbers (22x RAM advantage on T2; constant 58-61 MB across 8 tasks; 7 of 8 at structural-equivalence-or-better; install footprint 163 MB single command). Frame as a cross-domain coverage claim, not a "third architecture" claim.

Paragraph 2 (~150 words). Who gfc serves. Users who value cross-domain coverage and dependency-light installs over per-task speed. Users embedding conversion into containerised pipelines where worst-case task drives resource allocation. Users for whom system-binary version drift is a known pain point. State this user band explicitly per FRAMING_DECISIONS Decision 2.

Paragraph 3 (~150 words). The four-state correctness matrix as a generalisable norm. Lift the dedicated draft's argument: structural divergence is not failure; binary correctness reporting mislabels convention divergence as defect. Cite Weber et al. 2019 for the benchmarking-guidelines precedent. Frame as a contribution worth adopting across converter benchmarks beyond gfc.

Paragraph 4 (~150 words). Limitations. gfc is 14.2x slower than plink2 on T2, 3.5x slower than pyhmmer on T6, 1.7x slower than gffread on T4. The slower tasks involve C tools with specialised in-memory representations (plink2's 2-bit packing, pyhmmer's HMMER C library binding, gffread's optimised AST). For users running one of those conversions repeatedly at scale, the dedicated tool remains the right choice. T6 and T7 currently use Pfam-A vs UniRef50 and Quest for Orthologs respectively at the dataset scale specified in LITERATURE_COMPARATORS; smaller fixtures from earlier benchmarks have been archived as supplementary. Cold-cache replicate variance is reported in Supplementary Table S2.

Paragraph 5 (~100 words). Future directions. Snakemake `benchmark:` wrapper for direct integration (Mölder et al. 2021). Additional format families (BAM index variants; tabix-indexed VCF). Cross-platform parity test on macOS (currently Linux-only). bio.tools registry entry before submission per CLAIM_VERIFICATION risk #6.

### Acknowledgements (~50 words)

GLBRC computing infrastructure. Hittinger Lab. Y1000+ project for the eukaryotic dataset slice (Shen et al. 2018, Cell 175:1533).

### References: TBD

Estimated 35-45 references. Required citations: Caro et al. 2023 (BioConvert); Patterson et al. 2006 (EIGENSOFT); Chang et al. 2015 (PLINK 2.0); Rice et al. 2000 (EMBOSS); Norling et al. 2018 (EMBLmyGFF3); Sayers et al. 2023 (table2asn); Pertea & Pertea 2020 (gffread); Dainat 2022 (AGAT); Kent 2002 (UCSC); Larralde & Zeller 2023 (pyhmmer); Emms & Kelly 2019 (OrthoFinder); Danecek et al. 2021 (bcftools); Korneliussen et al. 2014 (ANGSD); Cock et al. 2009 (Biopython); Weber et al. 2019 (benchmarking guidelines); Mölder et al. 2021 (Snakemake); Frankish et al. 2023 (GENCODE); Wagner et al. 2022 (GIAB); Shen et al. 2018 (Y1000+).

## Figure plan

| Figure | Content | Source script (branch `phase2-plot-overhaul`) |
|---|---|---|
| Fig 1 | Capability matrix: 24 conversions x 11 tools, 3-state colour, with bedtools / vcftools / pyfaidx / biopython.convert columns added | `analyze/feature_matrix.py` |
| Fig 2 | Wall-time per task, log y, mean +/- sd error bars, asterisks for structural-equivalence cases, divergence notes from `_REASON_T*` (not auto-generated "X disagreed with Y") | `analyze/plot_wall_time.py` |
| Fig 3 | Peak RSS per task, log y, gfc band annotated as "constant 58-61 MB across T1-T8" | `analyze/plot_memory.py` |
| Fig 4 | Pareto wall vs RSS, gfc cluster spread with leader lines, Pareto frontier line drawn explicitly | `analyze/plot_pareto.py` |
| Fig 5 | Four-state correctness matrix (byte-equal / structurally-equivalent / comparator-limit / disagreement), categorical labels not percentages, footnote text from `_REASON_T*` | `analyze/plot_correctness_matrix.py` (re-encoded per AUDIT lines 117-129) |
| Fig 6 | Speedup ratio bars per (task, competitor), log y, dashed reference at y=1, bootstrap 95% CI error bars, correctness annotation per bar | `analyze/plot_speedup.py` |

All figures use one shared Okabe-Ito tool palette and one Tableau-10 task palette via `analyze/_style.py`. All legends move to `bbox_to_anchor=(1.02, 1.0)` outside the plot area. AUDIT Phase 2 plot-infra fixes (1-9) are pre-conditions.

## Table plan

| Table | Content |
|---|---|
| Table 1 | Feature matrix summary: gfc cross-domain coverage vs comparators on (sequence, alignment, variant, annotation) breadth |
| Table 2 | Install footprint: gfc 163 MB / 32 deps; BioConvert 1186 MB / 201 deps; biopython.convert 118 MB / 32 deps. Anchor: R2 evidence numbers in brief |
| Table 3 | Per-task RSS stability: gfc per-task mean RSS, sd, CV, range; aggregate 58.78 MB mean, CV 2.04%, range 56.88-61.31 MB across 80 runs. Anchor: R3 evidence in brief |
| Table 4 | Per-pair Wilcoxon: 24 (task, competitor) pairs, paired W statistic, p-value, BH-adjusted q-value, median ratio with 95% bootstrap CI, Cliff's delta. 10 clean wins at q=0.05 |

## Supplementary plan

- Supplementary Figure S1: replicate variance per (task, tool) as CV bars (per AUDIT recommendation, replaces collapsed boxplot)
- Supplementary Figure S2: per-task small multiples, 2x4 grid, fixed tool->colour mapping across panels
- Supplementary Figure S3: cumulative compute cost, stacked bars, with caption explicitly framing "user pays for N tools" argument
- Supplementary Table S1: CV per (task, tool) for wall and peak RSS
- Supplementary Table S2: cold-cache supplementary run, n=5 per (task, tool), with `drop_caches` between replicates; warm-cache vs cold-cache comparison
- Supplementary Table S3: full Wilcoxon output with W, p-raw, p-BH, ratio CI, Cliff's delta for all 24 pairs
- Supplementary Table S4: BioConvert wrapper-tax per task (T4 BioConvert 139s vs gffread 3.2s for byte-identical output)
- Supplementary Methods: cold-cache protocol; provenance-header schema; pre-registration commit hash; deviations from BENCHMARK_PLAN.md

## Open structural questions

1. Where do the four Phase 4 tasks (T3 table2asn / EMBLmyGFF3, T7 OrthoFinder / Proteinortho / SonicParanoid, T8 ANGSD / pileupCaller, BioConvert as cross-task baseline) land in the figure budget? Current plan uses primary numbers from the n=10 production run; Phase 4 additions go into Figures 2, 3, 6 as additional bars or into supplementary tables.
2. Does the Pareto figure (Fig 4) replace one of the other six, or is it supplementary? Current plan: keep Fig 4 as primary; the constant-RSS band visualisation is the load-bearing image for the "RSS-stable" claim.
3. T6/T7 real-scale data: brief implies the Pfam-A vs UniRef50 and Quest for Orthologs runs are landed, but AUDIT line 587-590 flags both as gated on data downloads. Confirm before drafting Results.
