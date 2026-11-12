"""Converter modules for genome format conversions.

Submodules are not imported here, so that importing one converter (or the
CLI) does not pull in every converter's dependencies. Import them directly,
for example ``from genome_format_converters.converters import gff3_to_gtf``.
"""

__all__ = [
    "gff3_to_gtf",
    "gff3_to_bed",
    "gff3_to_bed12",
    "genbank_to_gff3",
    "gff3_to_table",
    "gff3_to_protein",
    "gtf_to_gff3",
    "fasta_to_fastq",
    "fastq_to_fasta",
    "fasta_qual_to_fastq",
    "fastq_to_fasta_qual",
    "convert_alignment",
    "fasta_to_table",
    "stockholm_to_fasta",
    "fasta_to_stockholm",
    "bam_to_bed",
    "blast_tab_to_links",
    "delta_to_tab",
    "maf_to_xmfa",
    "hmmer_tblout_to_tsv",
    "orthogroups_to_fasta",
    "vcf_to_bed",
    "vcf_to_table",
    "vcf_to_consensus",
    "vcf_to_eigenstrat",
    "vcf_to_pseudohaploid",
    "vcf_to_plink",
    "vcf_to_snp_alignment",
    "eigenstrat_to_vcf",
    "tree_convert",
    "annotate_tree",
    "convert_all_gff_fasta_to_gbk",
]
