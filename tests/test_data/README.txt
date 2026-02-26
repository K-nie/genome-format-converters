General Usage
All scripts are used in the same way:

bash
python script_name.py --input-dir /path/to/input/files --output-dir /path/to/output
The input directory should contain the files you want to convert.

The output directory will be created if it doesn’t exist.

Each script processes all files with recognised extensions in the input directory.

Output files keep the same base name but receive a new extension (see individual script descriptions).

Scripts Overview
Annotation format conversions
Script	Description	Input extensions	Output extension
gff3_to_gtf.py	GFF3 → GTF	.gff3, .gff	.gtf
gff3_to_bed.py	GFF3 → 6‑column BED	.gff3, .gff	.bed
genbank_to_gff3.py	GenBank → GFF3	.gbk, .gb, .gbff	.gff3
gff3_to_table.py	GFF3 → tab‑separated feature table	.gff3, .gff	.tsv
gff3_to_protein.py	GFF3 + FASTA → protein FASTA	.gff3/.gff + .fasta/.fa	.faa
Examples:

bash
python gff3_to_gtf.py --input-dir ./my_gff_files --output-dir ./gtf_output
python gff3_to_protein.py --input-dir ./my_data --output-dir ./proteins
Sequence format conversions
Script	Description	Input extensions	Output extension
fasta_to_fastq.py	FASTA → FASTQ (with default quality)	.fasta, .fa, .fna, .fas	.fastq
fastq_to_fasta.py	FASTQ → FASTA	.fastq, .fq	.fasta
fasta_qual_to_fastq.py	Combine FASTA + QUAL → FASTQ	.fasta/.fa + .qual	.fastq
fastq_to_fasta_qual.py	Split FASTQ → FASTA + QUAL	.fastq, .fq	.fasta, .qual
convert_alignment.py	Alignment format converter (FASTA, PHYLIP, NEXUS, CLUSTAL)	any alignment file	user‑specified
fasta_to_table.py	FASTA → two‑column TSV (ID, sequence)	.fasta, .fa, .fna, .fas	.tsv
Examples:

bash
python fasta_to_fastq.py --input-dir ./fasta_files --output-dir ./fastq_output
python convert_alignment.py --input-dir ./aln --output-dir ./phylip --in-format fasta --out-format phylip
Alignment / mapping results
Script	Description	Input extensions	Output extension
bam_to_bed.py	BAM/SAM → BED6	.bam, .sam	.bed
blast_tab_to_links.py	BLAST tabular (outfmt 6) → simplified link TSV	.tab	.links.tsv
delta_to_tab.py	MUMmer .delta → tabular alignment coordinates	.delta	.tsv
maf_to_xmfa.py	MAF → XMFA (progressiveMauve format)	.maf	.xmfa
Examples:

bash
python bam_to_bed.py --input-dir ./bam_files --output-dir ./bed_output
python blast_tab_to_links.py --input-dir ./blast_results --output-dir ./links --min-length 100 --min-identity 30
Variant formats (VCF)
Script	Description	Input extensions	Output extension
vcf_to_bed.py	VCF/BCF → BED intervals	.vcf, .vcf.gz, .bcf	.bed
vcf_to_table.py	VCF/BCF → tab‑separated table (TSV)	.vcf, .vcf.gz, .bcf	.tsv
vcf_to_consensus.py	VCF + reference FASTA → consensus FASTA per sample	.vcf/.vcf.gz + .fasta	.fa
Examples:

bash
python vcf_to_bed.py --input-dir ./vcf_files --output-dir ./bed_output
python vcf_to_consensus.py --input-dir ./vcf_with_ref --output-dir ./consensus
Phylogenetic tree formats
Script	Description	Input extensions	Output extension
tree_convert.py	Newick ↔ NEXUS ↔ PhyloXML	any tree file	user‑specified
annotate_tree.py	Add alignment sequences to tree (NEXUS output)	.nwk + .fasta (aligned)	.nex
Examples:

bash
python tree_convert.py --input-dir ./trees --output-dir ./converted --in-format newick --out-format nexus
python annotate_tree.py --tree tree.nwk --aln alignment.fasta --output annotated.nex
Testing with Provided Data
The repository includes a set of small example files in the test_data/ directory. You can use them to verify that the scripts work on your system.

Navigate to the repository root:

bash
cd genome-format-converters
Create a test output directory (e.g., test_output/):

bash
mkdir test_output
Run any script on the test data. For example:

bash
python gff3_to_gtf.py --input-dir test_data --output-dir test_output
Inspect the output in test_output/. The generated files should contain correctly converted data.

All test files are described below:

File	Purpose
tiny.gff3	Minimal GFF3 with a two‑exon gene.
tiny.fasta	Corresponding genome sequence for tiny.gff3.
tiny.gbk	Minimal GenBank record for genbank_to_gff3.py.
test.fasta, test.fastq, test.qual	Small FASTA, FASTQ, and QUAL files for sequence conversions.
aln.fasta	Small aligned FASTA for alignment converters.
test.sam, test.bam	Single‑read SAM/BAM for bam_to_bed.py.
blast.tab	Two‑line BLAST tabular output.
test.delta	Simplified MUMmer .delta file.
test.maf	Two‑block MAF alignment.
test.vcf	Three‑variant VCF file for VCF scripts.
ref.fasta (symlink to tiny.fasta)	Reference for vcf_to_consensus.py.
test.nwk	Simple Newick tree.