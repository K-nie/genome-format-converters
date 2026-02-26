Here's the updated README.md with the new command reference section added after "Usage". I've also updated the Table of Contents to include the new section.

markdown
# Genome Format Converters

![Python 3.6+](https://img.shields.io/badge/python-3.6+-blue.svg)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

A collection of Python scripts for converting common bioinformatics file formats.
Each script follows a simple, uniform interface: you point it to an input directory, and it writes converted files to an output directory.

## Table of Contents
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Command Reference](#command-reference)
  - [Annotation Format Conversions](#annotation-format-conversions)
  - [Sequence Format Conversions](#sequence-format-conversions)
  - [Alignment / Mapping Results](#alignment--mapping-results)
  - [Variant Formats (VCF)](#variant-formats-vcf)
  - [Phylogenetic Tree Formats](#phylogenetic-tree-formats)
- [Scripts Overview](#scripts-overview)
- [License](#license)
- [Contributing](#contributing)

## Features
- **Uniform interface**: all scripts accept `--input-dir` and `--output-dir` arguments.
- **Batch processing**: convert all files of a given type in a directory at once.
- **Lightweight**: only requires a few well‑maintained Python libraries.
- **Well tested**: each script has been tested on small example datasets.

## Installation

Clone the repository:
<<<<<<< HEAD
```bash
=======

bash
>>>>>>> d358144a2a02b65cb03d7e47c3e0728031948e20
git clone https://github.com/K-nie/genome-format-converters.git
cd genome-format-converters
Install the required dependencies:

##bash
pip install -r requirements.txt
## Note: For scripts that work with BAM or VCF files, you also need pysam (included in requirements.txt).
## For BLAST tabular conversion, you need BLAST+ installed separately (optional – only if you generate the input files).

Usage
All scripts are used in the same way:

##bash
gfc <subcommand> --input-dir INPUT_DIR --output-dir OUTPUT_DIR [options]
## The input directory should contain the files you want to convert. The output directory will be created if it doesn’t exist.

## Each subcommand processes all files with recognised extensions in the input directory.

## Run gfc --help to see all available subcommands, or gfc <subcommand> --help for detailed options.

## Command Reference

## Annotation Format Conversions
## Subcommand	Description	Example
gff3-to-gtf	Convert GFF3 to GTF	gfc gff3-to-gtf --input-dir ./gff_files --output-dir ./gtf_output
gff3-to-bed	Convert GFF3 to 6‑column BED	gfc gff3-to-bed --input-dir ./gff_files --output-dir ./bed_output
genbank-to-gff3	Convert GenBank to GFF3	gfc genbank-to-gff3 --input-dir ./gbk_files --output-dir ./gff3_output
gff3-to-table	Convert GFF3 to tab‑separated feature table	gfc gff3-to-table --input-dir ./gff_files --output-dir ./table_output
gff3-to-protein	Extract protein sequences from GFF3 + FASTA	gfc gff3-to-protein --input-dir ./data --output-dir ./proteins
fasta-gff-to-gbk	Convert paired FASTA and GFF3 files to GenBank	gfc fasta-gff-to-gbk --input-dir ./data --output-dir ./gbk_output

## Sequence Format Conversions
## Subcommand	Description	Example
fasta-to-fastq	FASTA → FASTQ with default quality (I)	gfc fasta-to-fastq --input-dir ./fasta --output-dir ./fastq
fastq-to-fasta	FASTQ → FASTA (drop qualities)	gfc fastq-to-fasta --input-dir ./fastq --output-dir ./fasta
fasta-qual-to-fastq	Combine FASTA + QUAL into FASTQ	gfc fasta-qual-to-fastq --input-dir ./data --output-dir ./fastq
fastq-to-fasta-qual	Split FASTQ into FASTA and QUAL	gfc fastq-to-fasta-qual --input-dir ./fastq --output-dir ./split
fasta-to-table	FASTA → two‑column TSV (id, sequence)	gfc fasta-to-table --input-dir ./fasta --output-dir ./tables
convert-alignment	Convert alignment formats (fasta, phylip, nexus, clustal)	gfc convert-alignment --input-dir ./aln --output-dir ./phylip --in-format fasta --out-format phylip

## Alignment / Mapping Results
## Subcommand	Description	Example
bam-to-bed	Convert BAM/SAM to BED6	gfc bam-to-bed --input-dir ./bam_files --output-dir ./bed
blast-to-links	Convert BLAST tabular (outfmt 6) to link TSV	gfc blast-to-links --input-dir ./blast_results --output-dir ./links --min-length 100 --min-identity 30
delta-to-tab	Convert MUMmer .delta to tabular coordinates	gfc delta-to-tab --input-dir ./delta_files --output-dir ./tables
maf-to-xmfa	Convert MAF to XMFA (progressiveMauve format)	gfc maf-to-xmfa --input-dir ./maf_files --output-dir ./xmfa

## Variant Formats (VCF)
## Subcommand	Description	Example
vcf-to-bed	Convert VCF to BED intervals	gfc vcf-to-bed --input-dir ./vcf_files --output-dir ./bed
vcf-to-table	Convert VCF to tab‑separated table	gfc vcf-to-table --input-dir ./vcf_files --output-dir ./tables
vcf-to-consensus	Create consensus FASTA from VCF + reference	gfc vcf-to-consensus --input-dir ./data --output-dir ./consensus

## Phylogenetic Tree Formats
## Subcommand	Description	Example
tree-convert	Convert tree formats (newick, nexus, phyloxml)	gfc tree-convert --input-dir ./trees --output-dir ./converted --in-format newick --out-format nexus
annotate-tree	Add alignment sequences to tree (output NEXUS)	gfc annotate-tree --tree tree.nwk --aln alignment.fasta --output annotated.nex

## Scripts Overview
## The underlying Python scripts are located in src/genome_format_converters/converters/. Each script can also be run independently (though the gfc interface is recommended). Below is a quick reference of the scripts and their input/output formats.

## Script	Description	Input extensions	Output extension
gff3_to_gtf.py	GFF3 → GTF	.gff3, .gff	.gtf
gff3_to_bed.py	GFF3 → 6‑column BED	.gff3, .gff	.bed
genbank_to_gff3.py	GenBank → GFF3	.gbk, .gb	.gff3
gff3_to_table.py	GFF3 → tab‑separated feature table	.gff3, .gff	.tsv
gff3_to_protein.py	GFF3 + FASTA → protein FASTA	.gff3/.gff + .fasta/.fa	.faa
fasta_to_fastq.py	FASTA → FASTQ (with default quality)	.fasta, .fa, .fna, .fas	.fastq
fastq_to_fasta.py	FASTQ → FASTA	.fastq, .fq	.fasta
fasta_qual_to_fastq.py	Combine FASTA + QUAL → FASTQ	.fasta/.fa + .qual	.fastq
fastq_to_fasta_qual.py	Split FASTQ → FASTA + QUAL	.fastq, .fq	.fasta, .qual
convert_alignment.py	Alignment format converter (FASTA, PHYLIP, NEXUS, CLUSTAL)	any alignment file	user‑specified
fasta_to_table.py	FASTA → two‑column TSV (ID, sequence)	.fasta, .fa, .fna, .fas	.tsv
bam_to_bed.py	BAM/SAM → BED6	.bam, .sam	.bed
blast_tab_to_links.py	BLAST tabular (outfmt 6) → simplified link TSV	.tab	.links.tsv
delta_to_tab.py	MUMmer .delta → tabular alignment coordinates	.delta	.tsv
maf_to_xmfa.py	MAF → XMFA (progressiveMauve format)	.maf	.xmfa
vcf_to_bed.py	VCF/BCF → 1‑bp BED intervals	.vcf, .vcf.gz, .bcf	.bed
vcf_to_table.py	VCF/BCF → tab‑separated table (TSV)	.vcf, .vcf.gz, .bcf	.tsv
vcf_to_consensus.py	VCF + reference FASTA → consensus FASTA per sample	.vcf/.vcf.gz + .fasta	.fa
tree_convert.py	Newick ↔ NEXUS ↔ PhyloXML	.nwk, .nex, .xml	user‑specified
annotate_tree.py	Add alignment sequences to tree (NEXUS output)	.nwk + .fasta (aligned)	.nex
convert_all_gff_fasta_to_gbk.py	FASTA+GFF → GenBank	.fasta/.fa + .gff3/.gff	.gbk


### Testing
All scripts have been tested on small example datasets located in the tests/test_data/ directory. These test files cover the basic functionality of each converter. To run the tests yourself, install the package in development mode (pip install -e .) and execute the example commands from the Command Reference using the provided test data. For instance:

##bash
gfc gff3-to-gtf --input-dir tests/test_data --output-dir test_output
gfc fasta-to-fastq --input-dir tests/test_data --output-dir test_output
# ... etc.

## License
## This project is licensed under the MIT License – see the LICENSE file for details.

## Contributing
## Contributions are welcome! If you have a new converter or an improvement, please open an issue or submit a pull request.