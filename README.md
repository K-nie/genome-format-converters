# Genome Format Converters

**Author:** Benjamin Narh-Madey
**Affiliation:** Hittinger Lab, Laboratory of Genetics, University of Wisconsin-Madison
**Contact:** narhmadey@wisc.edu

![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)
![Version 0.1.5](https://img.shields.io/badge/version-0.1.5-green.svg)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

A uniform command-line toolkit for converting common bioinformatics file
formats. Each converter follows the same interface: point it at an input
directory (or single file), point it at an output directory (or single file
/ stdout), and it writes the conversion.

## Table of Contents

- [Features](#features)
- [Benchmark](#benchmark)
- [Installation](#installation)
- [Usage](#usage)
- [Command Reference](#command-reference)
  - [Annotation formats](#annotation-formats)
  - [Sequence formats](#sequence-formats)
  - [Alignment / mapping results](#alignment--mapping-results)
  - [Variant formats (VCF / EIGENSTRAT / PLINK)](#variant-formats-vcf--eigenstrat--plink)
  - [Phylogenetic tree formats](#phylogenetic-tree-formats)
  - [Orthology / homology](#orthology--homology)
  - [HMMER output](#hmmer-output)
- [Scripts Overview](#scripts-overview)
- [Testing](#testing)
- [License](#license)
- [Contributing](#contributing)

## Features

- **Uniform interface** — every subcommand accepts `--input-dir` / `--output-dir` (batch mode) or `--input` / `--output` (single-file mode, with `-` for stdin/stdout where it makes sense).
- **Batch processing** — convert every file of a recognised type in a directory in one call. `--pattern GLOB` restricts the set; `--threads N` converts files in parallel worker processes, with output identical to a sequential run.
- **Fails loudly** — converters raise on malformed input and exit non-zero so pipelines can detect failure. Progress and warnings go to stderr; data goes to stdout or the output path.
- **Lightweight** — pure-Python, depends on Biopython, pysam, cyvcf2, bcbio-gff, numpy, and pandas.
- **Checked correctness** — the benchmark scores every output against a reference tool with executable comparators (byte-identical, structurally equivalent under a named convention, not comparable, or disagreement), and a mutation test shows they catch seeded defects.

## Benchmark

`benchmarks/` holds a head-to-head harness covering 8 tasks against
established single-purpose tools: EIGENSOFT convertf, PLINK 1.9/2,
bcftools, AGAT, gffread, UCSC utilities, bedops, EMBOSS seqret, pyhmmer and
BioConvert. Inputs are 1000 Genomes chr22 biallelic SNPs for the VCF tasks
and a Y1000+ slice for the annotation tasks.

| Task | Conversion |
|---|---|
| T1 | VCF → EIGENSTRAT |
| T2 | VCF → PLINK binary |
| T3 | FASTA + GFF3 → GenBank |
| T4 | GFF3 → GTF |
| T5 | GFF3 → BED12 |
| T6 | HMMER tblout → TSV |
| T7 | OrthoFinder orthogroups → per-OG FASTA |
| T8 | VCF → pseudohaploid EIGENSTRAT |

**Correctness.** `benchmarks/four_state.py` assigns each (task, tool)
output one of four states against the task's reference:

- `byte`: byte-identical after line-ending normalisation.
- `structural`: same content once a named convention difference is allowed for. The verdict records which convention, e.g. `chrom-label`, `allele-order` or `exon-fusion`.
- `limitation`: the outputs cannot be compared on this axis, e.g. a BED6 reference against BED12 output.
- `disagreement`: a content difference that no documented convention explains.

T8 is stochastic, so it is checked against its input VCF instead. Every
call must be consistent with the source genotype, and het draws must be
balanced. `python benchmarks/mutation_test.py` seeds 29 defects across
the 8 tasks and confirms each one is reported as `disagreement`, while
convention-only rewrites stay `structural`.

Updated timings with these checks are being re-measured. Reproduce with
`bash benchmarks/run_bench.sh` after `pip install -e .`.

## Installation

### From PyPI

```bash
pip install genome-format-converters
```

The command it installs is `gfc`. (The unrelated PyPI package named `gfc`
is not this tool.) PyPI currently has an older release (0.1.2); install
from source for everything described here.

### From source

```bash
git clone https://github.com/K-nie/genome-format-converters.git
cd genome-format-converters
pip install -e .
```

For the test suite add the `test` extra:

```bash
pip install -e ".[test]"
```

**External dependencies.** None at runtime: no converter calls an external binary. VCF handling uses `cyvcf2` and `numpy`, BAM and indexed FASTA use `pysam`, all installed as Python packages. On Windows, `pysam`/`cyvcf2` are skipped and the VCF/BAM subcommands print a pointer to WSL or conda-forge. BLAST tabular conversion only needs BLAST+ separately if you also want to *generate* the `-outfmt 6` input — the converter itself is self-contained.

## Usage

After installation the `gfc` command becomes available. The general form is:

```bash
gfc <subcommand> --input-dir INPUT_DIR --output-dir OUTPUT_DIR [options]
# or, single-file:
gfc <subcommand> --input FILE --output FILE [options]
```

The output directory is created if it doesn't exist. Each subcommand processes all files with recognised extensions in the input directory. Use `--pattern '*.vcf.gz'` to restrict further, `--threads N` to parallelise, and `--force` to overwrite existing outputs.

### Getting help

```bash
gfc --help                  # list every subcommand
gfc <subcommand> --help     # detailed options for one subcommand
```

## Command Reference

### Annotation formats

| Subcommand          | Description                                      | Example |
| ------------------- | ------------------------------------------------ | ------- |
| `gff3-to-gtf`       | GFF3 → GTF                                       | `gfc gff3-to-gtf --input-dir ./gff --output-dir ./gtf` |
| `gtf-to-gff3`       | GTF → GFF3 (synthesises gene rows)               | `gfc gtf-to-gff3 --input-dir ./gtf --output-dir ./gff3` |
| `gff3-to-bed`       | GFF3 → 6-column BED                              | `gfc gff3-to-bed --input-dir ./gff --output-dir ./bed` |
| `gff3-to-bed12`     | GFF3 → BED12 (exon blocks per transcript)        | `gfc gff3-to-bed12 --input-dir ./gff --output-dir ./bed12` |
| `genbank-to-gff3`   | GenBank → GFF3                                   | `gfc genbank-to-gff3 --input-dir ./gbk --output-dir ./gff3` |
| `gff3-to-table`     | GFF3 → tab-separated feature table               | `gfc gff3-to-table --input-dir ./gff --output-dir ./tsv` |
| `gff3-to-protein`   | Extract proteins from GFF3 + FASTA               | `gfc gff3-to-protein --input-dir ./data --output-dir ./proteins` |
| `fasta-gff-to-gbk`  | Paired FASTA + GFF3 → GenBank                    | `gfc fasta-gff-to-gbk --input-dir ./data --output-dir ./gbk` |

### Sequence formats

| Subcommand             | Description                                            | Example |
| ---------------------- | ------------------------------------------------------ | ------- |
| `fasta-to-fastq`       | FASTA → FASTQ with a default quality character         | `gfc fasta-to-fastq --input-dir ./fa --output-dir ./fq` |
| `fastq-to-fasta`       | FASTQ → FASTA (drop qualities)                         | `gfc fastq-to-fasta --input-dir ./fq --output-dir ./fa` |
| `fasta-qual-to-fastq`  | Combine FASTA + QUAL into FASTQ                        | `gfc fasta-qual-to-fastq --input-dir ./data --output-dir ./fq` |
| `fastq-to-fasta-qual`  | Split FASTQ into FASTA + QUAL                          | `gfc fastq-to-fasta-qual --input-dir ./fq --output-dir ./split` |
| `fasta-to-table`       | FASTA → two-column TSV (id, sequence)                  | `gfc fasta-to-table --input-dir ./fa --output-dir ./tsv` |
| `convert-alignment`    | Convert between FASTA / PHYLIP / NEXUS / CLUSTAL       | `gfc convert-alignment --input-dir ./aln --output-dir ./phy --in-format fasta --out-format phylip` |
| `stockholm-to-fasta`   | Stockholm (Rfam / Pfam) → FASTA alignment              | `gfc stockholm-to-fasta --input-dir ./sto --output-dir ./fa` |
| `fasta-to-stockholm`   | FASTA alignment → Stockholm                            | `gfc fasta-to-stockholm --input-dir ./fa --output-dir ./sto` |

### Alignment / mapping results

| Subcommand        | Description                                         | Example |
| ----------------- | --------------------------------------------------- | ------- |
| `bam-to-bed`      | BAM / SAM → BED6                                    | `gfc bam-to-bed --input-dir ./bam --output-dir ./bed` |
| `blast-to-links`  | BLAST tabular (outfmt 6) → simplified link TSV      | `gfc blast-to-links --input-dir ./blast --output-dir ./links --min-length 100 --min-identity 30` |
| `delta-to-tab`    | MUMmer `.delta` → tabular alignment coordinates     | `gfc delta-to-tab --input-dir ./delta --output-dir ./tsv` |
| `maf-to-xmfa`     | MAF → XMFA (progressiveMauve format)                | `gfc maf-to-xmfa --input-dir ./maf --output-dir ./xmfa` |

### Variant formats (VCF / EIGENSTRAT / PLINK)

| Subcommand             | Description                                                                    | Example |
| ---------------------- | ------------------------------------------------------------------------------ | ------- |
| `vcf-to-bed`           | VCF → BED intervals                                                            | `gfc vcf-to-bed --input-dir ./vcf --output-dir ./bed` |
| `vcf-to-table`         | VCF → tab-separated table                                                      | `gfc vcf-to-table --input-dir ./vcf --output-dir ./tsv` |
| `vcf-to-consensus`     | VCF + reference FASTA → per-sample consensus FASTA                             | `gfc vcf-to-consensus --input-dir ./data --output-dir ./consensus` |
| `vcf-to-eigenstrat`    | VCF → EIGENSTRAT `.geno` / `.snp` / `.ind` (biallelic SNPs only)               | `gfc vcf-to-eigenstrat --input-dir ./vcf --output-dir ./eig --pop-map pops.tsv` |
| `vcf-to-pseudohaploid` | VCF → pseudohaploid EIGENSTRAT (one random allele per heterozygous site)       | `gfc vcf-to-pseudohaploid --input-dir ./vcf --output-dir ./eig --seed 42` |
| `vcf-to-plink`         | VCF → PLINK binary (`.bed` / `.bim` / `.fam`); shares filters with EIGENSTRAT  | `gfc vcf-to-plink --input-dir ./vcf --output-dir ./plink --min-maf 0.05` |
| `eigenstrat-to-vcf`    | EIGENSTRAT triplet → VCF (diploid, or `--haploid` for pseudohaploid data)      | `gfc eigenstrat-to-vcf --input-dir ./eig --output-dir ./vcf --ref-fasta genome.fa` |
| `vcf-to-snp-alignment` | VCF → SNP alignment for phylogenetics (relaxed PHYLIP, FASTA or NEXUS)         | `gfc vcf-to-snp-alignment --input-dir ./vcf --output-dir ./aln --min-samples 10` |

**EIGENSTRAT note.** Outputs a triplet per input VCF (`<stem>.geno`, `<stem>.snp`, `<stem>.ind`). Only biallelic SNPs are emitted; indels, multi-allelic sites, and non-ACGT alleles are skipped per the `convertf` / `smartpca` convention. Genotype encoding is **2** = hom-ref, **1** = het, **0** = hom-alt, **9** = missing. The pseudohaploid variant draws one allele uniformly at random per site per sample, so its `.geno` file contains only 0, 2, and 9.

Useful flags on both EIGENSTRAT converters:

| Flag                  | Effect |
| --------------------- | ------ |
| `--pop-map FILE`      | Two-column TSV (sample → population) to populate column 3 of `.ind`. |
| `--sex-map FILE`      | Two-column TSV (sample → sex `M`/`F`/`U`) for column 2 of `.ind`. |
| `--chrom-map FILE`    | Two-column TSV (contig name → integer) — legacy smartpca requires integer chromosomes. |
| `--default-chrom N`   | Integer fallback for any contig not in the chrom map. |
| `--genetic-map FILE`  | PLINK-style `.map` (`chrom snp_id cM bp`) to populate the `.snp` genetic-position column. |
| `--ancestral-fasta F` | Polarize against an ancestral reference FASTA (per-site major/minor swap). |
| `--info-aa`           | Polarize using the `AA` INFO field from the VCF. |
| `--transversions-only`| Drop transitions (A↔G, C↔T) — useful against aDNA deamination. |
| `--min-maf F`         | Minimum minor-allele frequency. |
| `--max-missing F`     | Maximum fraction of missing genotypes per site. |
| `--strict-maps`       | Error instead of warn when pop/sex/chrom maps mention samples/contigs not in the VCF. |
| `--also-plink`        | Also emit PLINK-compatible `.ped` / `.map` alongside the EIGENSTRAT triplet. |
| `--seed N`            | (`vcf-to-pseudohaploid` only) reproducible random allele selection. |

**`eigenstrat-to-vcf`.** Writes the `.snp` first allele as REF and the second as ALT (2 → `0/0`, 1 → `0/1`, 0 → `1/1`, 9 → `./.`). EIGENSTRAT alleles follow the producing tool's polarisation, so pass `--ref-fasta` (indexed FASTA) to re-orient sites to the reference genome. Packed ANCESTRYMAP / TGENO `.geno` files are rejected with an explanation.

**`vcf-to-snp-alignment`.** One sequence per sample, one column per biallelic SNP. Heterozygous calls become IUPAC codes by default (`--het random` picks one allele, seeded per file with `--seed`; `--het missing` writes `N`). `--format phylip|fasta|nexus`, `--min-samples N` drops poorly called sites, and `--drop-invariant` removes columns where all called samples share one base (useful before IQ-TREE `+ASC`).

### Phylogenetic tree formats

| Subcommand         | Description                                              | Example |
| ------------------ | -------------------------------------------------------- | ------- |
| `tree-convert`     | Convert tree formats (Newick, NEXUS, PhyloXML)           | `gfc tree-convert --input-dir ./trees --output-dir ./converted --in-format newick --out-format nexus` |
| `annotate-tree`    | Annotate a tree with alignment sequences (NEXUS output)  | `gfc annotate-tree --tree tree.nwk --aln aln.fasta --output annotated.nex` |

### Orthology / homology

| Subcommand              | Description                                                                | Example |
| ----------------------- | -------------------------------------------------------------------------- | ------- |
| `orthogroups-to-fasta`  | OrthoFinder `Orthogroups.tsv` + per-species FASTAs → one FASTA per orthogroup, with `>species\|gene` headers. | `gfc orthogroups-to-fasta --orthogroups Orthogroups.tsv --fasta-dir ./proteomes --output-dir ./ogs` |

### HMMER output

| Subcommand              | Description                                                                | Example |
| ----------------------- | -------------------------------------------------------------------------- | ------- |
| `hmmer-tblout-to-tsv`   | HMMER `--tblout` (default) or `--domtblout` (`--hmmer-format domtblout`) → structured TSV. Description-of-target column preserved. | `gfc hmmer-tblout-to-tsv --input-dir ./hmmer --output-dir ./tsv` |

## Scripts Overview

The underlying Python modules live in `src/genome_format_converters/converters/`. Each module can be imported independently, but the `gfc` CLI is the recommended entry point.

| #  | Module                              | Description                                                        | Input extensions                 | Output extension             |
| -- | ----------------------------------- | ------------------------------------------------------------------ | -------------------------------- | ---------------------------- |
| 1  | `gff3_to_gtf.py`                    | GFF3 → GTF                                                          | `.gff3`, `.gff`                  | `.gtf`                       |
| 2  | `gff3_to_bed.py`                    | GFF3 → 6-column BED                                                 | `.gff3`, `.gff`                  | `.bed`                       |
| 3  | `genbank_to_gff3.py`                | GenBank → GFF3                                                      | `.gbk`, `.gb`, `.gbff`           | `.gff3`                      |
| 4  | `gff3_to_table.py`                  | GFF3 → tab-separated feature table                                  | `.gff3`, `.gff`                  | `.tsv`                       |
| 5  | `gff3_to_protein.py`                | GFF3 + FASTA → protein FASTA (phase-aware CDS concatenation)        | `.gff3`/`.gff` + `.fasta`/`.fa`  | `.faa`                       |
| 6  | `fasta_to_fastq.py`                 | FASTA → FASTQ with default quality                                  | `.fasta`, `.fa`, `.fna`, `.fas`  | `.fastq`                     |
| 7  | `fastq_to_fasta.py`                 | FASTQ → FASTA                                                       | `.fastq`, `.fq`                  | `.fasta`                     |
| 8  | `fasta_qual_to_fastq.py`            | Combine FASTA + QUAL → FASTQ                                        | `.fasta`/`.fa` + `.qual`         | `.fastq`                     |
| 9  | `fastq_to_fasta_qual.py`            | Split FASTQ → FASTA + QUAL                                          | `.fastq`, `.fq`                  | `.fasta`, `.qual`            |
| 10 | `convert_alignment.py`              | FASTA / PHYLIP / NEXUS / CLUSTAL interconversion                    | any alignment                    | user-specified               |
| 11 | `fasta_to_table.py`                 | FASTA → two-column TSV                                              | `.fasta`, `.fa`, `.fna`, `.fas`  | `.tsv`                       |
| 12 | `bam_to_bed.py`                     | BAM / SAM → BED6                                                    | `.bam`, `.sam`                   | `.bed`                       |
| 13 | `blast_tab_to_links.py`             | BLAST tabular (outfmt 6) → link TSV                                 | `.tab`, `.tsv`                   | `.links.tsv`                 |
| 14 | `delta_to_tab.py`                   | MUMmer `.delta` → tabular                                           | `.delta`                         | `.tsv`                       |
| 15 | `maf_to_xmfa.py`                    | MAF → XMFA                                                          | `.maf`                           | `.xmfa`                      |
| 16 | `vcf_to_bed.py`                     | VCF / BCF → BED                                                     | `.vcf`, `.vcf.gz`, `.bcf`        | `.bed`                       |
| 17 | `vcf_to_table.py`                   | VCF / BCF → tab-separated table                                     | `.vcf`, `.vcf.gz`, `.bcf`        | `.tsv`                       |
| 18 | `vcf_to_consensus.py`               | VCF + reference FASTA → per-sample consensus FASTA (multi-contig, indel-aware, IUPAC for hets) | `.vcf`/`.vcf.gz` + `.fasta`/`.fa`/`.fna`/`.fas` | `.fa` |
| 19 | `vcf_to_eigenstrat.py`              | VCF / BCF → EIGENSTRAT triplet                                      | `.vcf`, `.vcf.gz`, `.bcf`        | `.geno` + `.snp` + `.ind`    |
| 20 | `vcf_to_pseudohaploid.py`           | VCF / BCF → pseudohaploid EIGENSTRAT                                | `.vcf`, `.vcf.gz`, `.bcf`        | `.geno` + `.snp` + `.ind`    |
| 21 | `tree_convert.py`                   | Newick ↔ NEXUS ↔ PhyloXML                                           | `.nwk`, `.nex`, `.xml`           | user-specified               |
| 22 | `annotate_tree.py`                  | Tree + alignment → annotated NEXUS                                  | `.nwk` + `.fasta` (aligned)      | `.nex`                       |
| 23 | `convert_all_gff_fasta_to_gbk.py`   | FASTA + GFF3 → GenBank                                              | `.fasta`/`.fa` + `.gff3`/`.gff`  | `.gbk`                       |
| 24 | `gff3_to_bed12.py`                  | GFF3 → BED12 (exon blocks per transcript)                           | `.gff3`, `.gff`                  | `.bed12`                     |
| 25 | `gtf_to_gff3.py`                    | GTF (Ensembl-style) → GFF3 with synthesised gene rows               | `.gtf`, `.gff2`                  | `.gff3`                      |
| 26 | `stockholm_to_fasta.py`             | Stockholm (Rfam/Pfam) alignments → FASTA                            | `.sto`, `.stk`, `.stockholm`     | `.fasta`                     |
| 27 | `fasta_to_stockholm.py`             | FASTA alignment → Stockholm                                         | `.fasta`, `.fa`, `.aln`          | `.sto`                       |
| 28 | `hmmer_tblout_to_tsv.py`            | HMMER `--tblout` / `--domtblout` → TSV (column schema + description) | `.tblout`, `.domtblout`, `.hmmer` | `.tsv`                      |
| 29 | `orthogroups_to_fasta.py`           | OrthoFinder `Orthogroups.tsv` + per-species FASTAs → per-OG FASTAs  | `.tsv` + FASTA dir               | one `.fa` per orthogroup     |
| 30 | `vcf_to_plink.py`                   | VCF / BCF → PLINK binary (`.bed` + `.bim` + `.fam`)                 | `.vcf`, `.vcf.gz`, `.bcf`        | `.bed` + `.bim` + `.fam`     |
| 31 | `eigenstrat_to_vcf.py`              | EIGENSTRAT triplet → VCF                                            | `.geno` + `.snp` + `.ind`        | `.vcf`                       |
| 32 | `vcf_to_snp_alignment.py`           | VCF / BCF → SNP alignment                                           | `.vcf`, `.vcf.gz`, `.bcf`        | `.phy` / `.fasta` / `.nex`   |

## Testing

Each converter has a fixture under `tests/test_data/` and tests under `tests/`. The suite also runs the benchmark's correctness mutation test. To run:

```bash
pip install -e ".[test]"
pytest -v
```

Smoke-testing a single subcommand against the bundled fixtures:

```bash
mkdir -p /tmp/gfc_out
gfc vcf-to-eigenstrat --input-dir tests/test_data --output-dir /tmp/gfc_out
ls /tmp/gfc_out  # tiny.geno, tiny.snp, tiny.ind
```

## License

MIT — see [`LICENSE`](LICENSE).

## Contributing

Contributions are welcome. For new converters or bug fixes, please open an issue or a pull request against [K-nie/genome-format-converters](https://github.com/K-nie/genome-format-converters). See [`CONTRIBUTING.md`](CONTRIBUTING.md) for style and workflow notes.
