Author: Benjamin Narh-Madey
Date: 2026-05-03
Target venue: Bioinformatics (Oxford Academic), Original Paper

# Abstract draft: gfc-bench

250 words exactly (at the cap).

---

**Motivation.** Genome-file conversion spans sequence, alignment, variant, and annotation domains in nearly every multi-tool pipeline. Existing tools fall into two families. Domain-deep utilities (samtools, bcftools, AGAT, gffread) are fast within one format family but require chaining four to six binaries to cover a typical workflow. Wrapper frameworks (BioConvert, EMBOSS-seqret, biopython.convert) aggregate external binaries through a unified interface and inherit those binaries' install footprint and version-pinning surface.

**Results.** We present gfc, a single Python package performing ~30 conversions natively across all four domains with no external binary calls. Benchmarking against citable per-task references on eight tasks (10 paired replicates, paired Wilcoxon with BH-FDR at q=0.05) yields four contributions. First, peak memory is bounded at 58.78 MB (CV 2.04%, range 56.88-61.31 MB) across 80 runs over eight heterogeneous tasks. Second, gfc uses 22.09x less RAM than plink2 [95% CI 21.70-22.46] on VCF-to-PLINK. Third, gfc avoids wrapper-tax: the install closure is 7.3x lighter and 6.3x fewer dependencies than BioConvert, and 44x faster than BioConvert on T4 GFF3-to-GTF for byte-identical output. Fourth, a four-state correctness matrix (byte-equivalent, structurally equivalent, comparator-limitation, disagreement) flags structural divergences rather than mislabelling them as failures. gfc is 14.2x slower than plink2 on T2, 3.5x slower than pyhmmer on T6, and 1.7x slower than gffread on T4. We position gfc for users who value cross-domain coverage and dependency-light installs over per-task speed.

**Availability and implementation.** gfc is MIT-licensed at https://github.com/K-nie/genome-format-converters, installable via `pip install gfc`. Benchmark code and raw measurements are in the same repository.

**Contact.** narhmadey@wisc.edu

---

## Section word counts

- Motivation: 64 words
- Results: 163 words
- Availability and implementation: 22 words
- Contact: 1 word
- **Total: 250 words** (at cap)

## Anchor numbers used

- Stats framework: n=10 paired replicates, paired Wilcoxon, BH-FDR q=0.05.
- R2 install footprint: 7.3x lighter, 6.3x fewer deps vs BioConvert.
- T4 wrapper-tax at runtime: 44x faster than BioConvert for byte-identical output.
- R3 RSS stability: 58.78 MB mean, CV 2.04%, range 56.88-61.31 MB across 80 runs.
- Headline RAM ratio: 22.09x less RAM than plink2 [95% CI 21.70-22.46].
- Speed losses: 14.2x slower than plink2 on T2, 3.5x slower than pyhmmer on T6, 1.7x slower than gffread on T4.

## Anchor numbers deferred to body

- T1 wall 1.29x [1.27-1.31] vs convertf — Results figure
- T3 wall 7.95x [7.71-8.07] vs EMBOSS-seqret — Results figure
- T4 wall 12.96x [12.80-13.03] and T5 wall 13.29x [13.13-13.44] vs AGAT — Results figure
- T8 wall 5.58x [5.37-5.82] vs bcftools-pyref — Results figure (refresh after Phase 4 ANGSD comparator)
- T1 RAM 22.0x [21.77-22.48] vs convertf, T2 vs plink1.9 1.18x, T4 vs AGAT 1.09x, T5 vs AGAT 1.11x — Table 4
- 10 of 24 paired Wilcoxon clean wins at q=0.05 — stated in Results, not abstract

## Open ordering question

The four contributions appear in order: wrapper-tax, memory bounding, 22x plink2, four-state matrix. FRAMING_DECISIONS Decision 1 says lead with "memory efficiency + reproducibility." The current order puts the load-bearing architectural-differentiator first (wrapper-tax), then the structural memory property, then the headline RAM number, then the methodological contribution. Both orderings are defensible under Decision 1. To match Decision 1 strictly, swap First and Second.
