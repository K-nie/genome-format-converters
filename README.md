Genome Format Converters
https://img.shields.io/badge/python-3.6+-blue.svg
https://img.shields.io/badge/License-MIT-yellow.svg

A collection of Python scripts for converting common bioinformatics file formats.
Each script follows a simple, uniform interface: you point it to an input directory, and it writes converted files to an output directory.

Table of Contents
Features

Installation

Usage

Scripts Overview

Annotation format conversions

Sequence format conversions

Alignment / mapping results

Variant formats (VCF)

Phylogenetic tree formats

License

Contributing

Features
Uniform interface: all scripts accept --input-dir and --output-dir arguments.

Batch processing: convert all files of a given type in a directory at once.

Lightweight: only requires a few well‑maintained Python libraries.

Well tested: each script has been tested on small example datasets.

Installation
Clone the repository:

bash
git clone https://github.com/K-nie/genome-format-converters.git
cd genome-format-converters
Install the required dependencies:

bash
pip install -r requirements.txt
Note: For scripts that work with BAM or VCF files, you also need pysam (included in requirements.txt).
For BLAST tabular conversion, you need BLAST+ installed separately (optional – only if you generate the input files).

Usage
All scripts are used in the same way:

bash
python script_name.py --input-dir /path/to/input/files --output-dir /path/to/output
The input directory should contain the files you want to convert.

The output directory will be created if it doesn’t exist.

Each script processes all files with recognised extensions in the input directory.

Scripts Overview
Annotation format conversions
Script	Description	Input extensions	Output extension
gff3_to_gtf.py	GFF3 → GTF	.gff3, .gff	.gtf
gff3_to_bed.py	GFF3 → 6‑column BED	.gff3, .gff	.bed
genbank_to_gff3.py	GenBank → GFF3	.gbk, .gb	.gff3
gff3_to_table.py	GFF3 → tab‑separated feature table	.gff3, .gff	.tsv
gff3_to_protein.py	GFF3 + FASTA → protein FASTA	.gff3, .gff + .fasta/.fa	.faa
Sequence format conversions
Script	Description	Input extensions	Output extension
fasta_to_fastq.py	FASTA → FASTQ (with default quality)	.fasta, .fa, .fna, .fas	.fastq
fastq_to_fasta.py	FASTQ → FASTA	.fastq, .fq	.fasta
fasta_qual_to_fastq.py	Combine FASTA + QUAL → FASTQ	.fasta/.fa + .qual	.fastq
fastq_to_fasta_qual.py	Split FASTQ → FASTA + QUAL	.fastq, .fq	.fasta, .qual
convert_alignment.py	Aligned FASTA → PHYLIP / NEXUS / CLUSTAL	.fasta (aligned)	.phylip / .nexus / .clustal
convert_any_alignment.py	Universal alignment converter (any format supported by Bio.AlignIO)	any alignment file	user‑specified
fasta_to_table.py	FASTA → two‑column TSV (ID, sequence)	.fasta, .fa, .fna, .fas	.tsv
Alignment / mapping results
Script	Description	Input extensions	Output extension
bam_to_bed.py	BAM/SAM → BED6	.bam, .sam	.bed
blast_tab_to_links.py	BLAST tabular (outfmt 6) → simplified link TSV	.tab	.links.tsv
delta_to_tab.py	MUMmer .delta → tabular alignment coordinates	.delta	.tsv
maf_to_xmfa.py	MAF → XMFA (progressiveMauve format)	.maf	.xmfa
Variant formats (VCF)
Script	Description	Input extensions	Output extension
vcf_to_bed.py	VCF/BCF → 1‑bp BED intervals	.vcf, .vcf.gz, .bcf	.bed
vcf_to_table.py	VCF/BCF → tab‑separated table (TSV)	.vcf, .vcf.gz, .bcf	.tsv
vcf_to_consensus.py	VCF + reference FASTA → consensus FASTA per sample	.vcf/.vcf.gz + .fasta	.fa
Phylogenetic tree formats
Script	Description	Input extensions	Output extension
tree_convert.py	Newick ↔ NEXUS ↔ PhyloXML	.nwk, .nex, .xml	user‑specified
annotate_tree.py	Add alignment sequences to tree (NEXUS output)	.nwk + .fasta (aligned)	.nex
License
This project is licensed under the MIT License – see the LICENSE file for details.

Contributing
Contributions are welcome! If you have a new converter or an improvement, please open an issue or submit a pull request.


