#!/usr/bin/env python3
"""
Unified command‑line interface for genome-format-converters.

Author: Benjamin Narh-Madey
Affiliation: Laboratory of Genetics (Hittinger Lab),
             Wisconsin Energy Institute,
             University of Wisconsin-Madison.
Email: narhmadey@swisc.edu

This package provides a collection of tools to convert common
bioinformatics file formats (GFF3, GTF, BED, GenBank, FASTA,
FASTQ, QUAL, alignment formats, BAM, VCF, MAF, MUMmer delta,
phylogenetic trees, etc.) with a simple, uniform interface.
"""

import argparse
import sys

# Import batch conversion functions from each module
from genome_format_converters.converters import (
    gff3_to_gtf,
    gff3_to_bed,
    genbank_to_gff3,
    gff3_to_table,
    gff3_to_protein,
    fasta_to_fastq,
    fastq_to_fasta,
    fasta_qual_to_fastq,
    fastq_to_fasta_qual,
    convert_alignment,
    fasta_to_table,
    bam_to_bed,
    blast_tab_to_links,
    delta_to_tab,
    maf_to_xmfa,
    vcf_to_bed,
    vcf_to_table,
    vcf_to_consensus,
    tree_convert,
    annotate_tree,
    convert_all_gff_fasta_to_gbk,
)

__version__ = "0.1.0"

def main():
    prog_description = (
        "Genome Format Converters – a suite of tools for converting "
        "bioinformatics file formats (GFF3, GTF, BED, GenBank, FASTA, "
        "FASTQ, QUAL, alignment formats, BAM, VCF, MAF, MUMmer delta, "
        "phylogenetic trees, etc.) with a uniform interface."
    )
    epilog_text = (
        "Author: Benjamin Narh-Madey\n"
        "Affiliation: Laboratory of Genetics (Hittinger Lab),\n"
        "             Wisconsin Energy Institute,\n"
        "             University of Wisconsin-Madison.\n"
        "Email: narhmadey@swisc.edu\n"
        f"\nVersion: {__version__}\n"
        "For bug reports or feature requests, please visit the GitHub repository."
    )

    parser = argparse.ArgumentParser(
        description=prog_description,
        epilog=epilog_text,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command", help="Subcommands")

    def add_io_args(subp):
        subp.add_argument("--input-dir", required=True,
                          help="Directory containing input files")
        subp.add_argument("--output-dir", required=True,
                          help="Output directory (created if it does not exist)")

    # GFF3 to GTF
    p = subparsers.add_parser("gff3-to-gtf", help="Convert GFF3 to GTF")
    add_io_args(p)

    # GFF3 to BED
    p = subparsers.add_parser("gff3-to-bed", help="Convert GFF3 to 6‑column BED")
    add_io_args(p)

    # GenBank to GFF3
    p = subparsers.add_parser("genbank-to-gff3", help="Convert GenBank to GFF3")
    add_io_args(p)

    # GFF3 to table
    p = subparsers.add_parser("gff3-to-table", help="Convert GFF3 to tab‑separated feature table")
    add_io_args(p)

    # GFF3+FASTA to protein
    p = subparsers.add_parser("gff3-to-protein", help="Extract protein sequences from GFF3 and FASTA")
    add_io_args(p)

    # FASTA to FASTQ
    p = subparsers.add_parser("fasta-to-fastq", help="Convert FASTA to FASTQ with default quality")
    add_io_args(p)
    p.add_argument("--qual-char", default="I",
                   help="Quality character (default I = 40, Illumina 1.8+ encoding)")

    # FASTQ to FASTA
    p = subparsers.add_parser("fastq-to-fasta", help="Convert FASTQ to FASTA (drop qualities)")
    add_io_args(p)

    # FASTA+QUAL to FASTQ
    p = subparsers.add_parser("fasta-qual-to-fastq", help="Combine FASTA and QUAL into FASTQ")
    add_io_args(p)

    # FASTQ to FASTA+QUAL
    p = subparsers.add_parser("fastq-to-fasta-qual", help="Split FASTQ into FASTA and QUAL files")
    add_io_args(p)

    # Alignment converter
    p = subparsers.add_parser("convert-alignment", help="Convert alignment between formats")
    add_io_args(p)
    p.add_argument("--in-format", required=True,
                   choices=["fasta", "phylip", "nexus", "clustal"],
                   help="Input alignment format")
    p.add_argument("--out-format", required=True,
                   choices=["fasta", "phylip", "nexus", "clustal"],
                   help="Output alignment format")

    # FASTA to table
    p = subparsers.add_parser("fasta-to-table", help="Convert FASTA to two‑column TSV (id, sequence)")
    add_io_args(p)

    # BAM to BED
    p = subparsers.add_parser("bam-to-bed", help="Convert BAM/SAM to BED6")
    add_io_args(p)

    # BLAST tabular to links
    p = subparsers.add_parser("blast-to-links", help="Convert BLAST tabular (outfmt 6) to simplified link TSV")
    add_io_args(p)
    p.add_argument("--min-length", type=int, default=0,
                   help="Minimum alignment length (default: 0 = no filter)")
    p.add_argument("--min-identity", type=float, default=0,
                   help="Minimum percent identity (default: 0 = no filter)")

    # Delta to tab
    p = subparsers.add_parser("delta-to-tab", help="Convert MUMmer .delta file to tabular coordinates")
    add_io_args(p)

    # MAF to XMFA
    p = subparsers.add_parser("maf-to-xmfa", help="Convert MAF to XMFA (progressiveMauve format)")
    add_io_args(p)

    # VCF to BED
    p = subparsers.add_parser("vcf-to-bed", help="Convert VCF to BED intervals")
    add_io_args(p)

    # VCF to table
    p = subparsers.add_parser("vcf-to-table", help="Convert VCF to tab‑separated table")
    add_io_args(p)

    # VCF to consensus
    p = subparsers.add_parser("vcf-to-consensus", help="Create consensus FASTA from VCF + reference")
    add_io_args(p)

    # Tree convert
    p = subparsers.add_parser("tree-convert", help="Convert tree formats (newick, nexus, phyloxml)")
    add_io_args(p)
    p.add_argument("--in-format", required=True,
                   choices=["newick", "nexus", "phyloxml"],
                   help="Input tree format")
    p.add_argument("--out-format", required=True,
                   choices=["newick", "nexus", "phyloxml"],
                   help="Output tree format")

    # Tree annotation
    p = subparsers.add_parser("annotate-tree", help="Add alignment sequences to tree (output NEXUS)")
    p.add_argument("--tree", required=True, help="Input tree file (Newick)")
    p.add_argument("--aln", required=True, help="Input alignment file (FASTA)")
    p.add_argument("--output", required=True, help="Output NEXUS file")

    # Convert FASTA+GFF to GenBank (original converter)
    p = subparsers.add_parser("fasta-gff-to-gbk",
                               help="Convert paired FASTA and GFF3 files to GenBank format")
    add_io_args(p)

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    # Dispatch to the appropriate batch function
    if args.command == "gff3-to-gtf":
        gff3_to_gtf.batch_convert(args.input_dir, args.output_dir)
    elif args.command == "gff3-to-bed":
        gff3_to_bed.batch_convert(args.input_dir, args.output_dir)
    elif args.command == "genbank-to-gff3":
        genbank_to_gff3.batch_convert(args.input_dir, args.output_dir)
    elif args.command == "gff3-to-table":
        gff3_to_table.batch_convert(args.input_dir, args.output_dir)
    elif args.command == "gff3-to-protein":
        gff3_to_protein.batch_convert(args.input_dir, args.output_dir)
    elif args.command == "fasta-to-fastq":
        fasta_to_fastq.batch_convert(args.input_dir, args.output_dir,
                                      qual_char=args.qual_char)
    elif args.command == "fastq-to-fasta":
        fastq_to_fasta.batch_convert(args.input_dir, args.output_dir)
    elif args.command == "fasta-qual-to-fastq":
        fasta_qual_to_fastq.batch_convert(args.input_dir, args.output_dir)
    elif args.command == "fastq-to-fasta-qual":
        fastq_to_fasta_qual.batch_convert(args.input_dir, args.output_dir)
    elif args.command == "convert-alignment":
        convert_alignment.batch_convert(args.input_dir, args.output_dir,
                                        args.in_format, args.out_format)
    elif args.command == "fasta-to-table":
        fasta_to_table.batch_convert(args.input_dir, args.output_dir)
    elif args.command == "bam-to-bed":
        bam_to_bed.batch_convert(args.input_dir, args.output_dir)
    elif args.command == "blast-to-links":
        blast_tab_to_links.batch_convert(args.input_dir, args.output_dir,
                                          min_length=args.min_length,
                                          min_identity=args.min_identity)
    elif args.command == "delta-to-tab":
        delta_to_tab.batch_convert(args.input_dir, args.output_dir)
    elif args.command == "maf-to-xmfa":
        maf_to_xmfa.batch_convert(args.input_dir, args.output_dir)
    elif args.command == "vcf-to-bed":
        vcf_to_bed.batch_convert(args.input_dir, args.output_dir)
    elif args.command == "vcf-to-table":
        vcf_to_table.batch_convert(args.input_dir, args.output_dir)
    elif args.command == "vcf-to-consensus":
        vcf_to_consensus.batch_convert(args.input_dir, args.output_dir)
    elif args.command == "tree-convert":
        tree_convert.batch_convert(args.input_dir, args.output_dir,
                                   args.in_format, args.out_format)
    elif args.command == "annotate-tree":
        annotate_tree.annotate_tree(args.tree, args.aln, args.output)
    elif args.command == "fasta-gff-to-gbk":
        convert_all_gff_fasta_to_gbk.batch_convert(args.input_dir, args.output_dir)
    else:
        print(f"Unknown command: {args.command}")
        sys.exit(1)

if __name__ == "__main__":
    main()