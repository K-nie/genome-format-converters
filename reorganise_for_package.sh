#!/bin/bash
# reorganise_for_package.sh – corrected version
# Run this script inside the 00_genome-format-converters directory.

set -e  # exit on error

echo "Reorganising files for Python package..."

# 1. Remove .sh extensions from test data files
if [ -d test_data ]; then
    echo "Renaming test files (removing .sh)..."
    for f in test_data/*.sh; do
        if [ -f "$f" ]; then
            mv "$f" "${f%.sh}"
        fi
    done
fi

# 2. Create the target directory structure
mkdir -p src/genome_format_converters/converters
mkdir -p tests/test_data

# 3. Move all Python scripts (except those already in the subdir) into converters/
echo "Moving Python scripts to src/genome_format_converters/converters/..."
for py in *.py; do
    if [ -f "$py" ]; then
        # Skip the reorganise script itself
        if [[ "$py" != "reorganise_for_package.sh" ]]; then
            mv "$py" src/genome_format_converters/converters/
        fi
    fi
done

# 4. Handle the existing genome-format-converters/ directory
if [ -d genome-format-converters ]; then
    echo "Moving contents of genome-format-converters/ to top level..."
    # Move everything (including hidden files) from that directory up
    mv genome-format-converters/* . 2>/dev/null || true
    mv genome-format-converters/.[!.]* . 2>/dev/null || true
    # Remove the now‑empty directory
    rmdir genome-format-converters 2>/dev/null || true
    # If still not empty, warn but continue
    if [ -d genome-format-converters ]; then
        echo "Warning: genome-format-converters/ not empty after move – you may need to clean it manually."
    fi
fi

# 5. Move test_data into tests/
if [ -d test_data ]; then
    echo "Moving test_data into tests/..."
    # Move all files (including hidden) from test_data to tests/test_data/
    shopt -s dotglob  # include hidden files
    mv test_data/* tests/test_data/ 2>/dev/null || true
    shopt -u dotglob
    # Remove the old test_data directory
    rm -rf test_data
fi

# 6. Create empty __init__.py files
touch src/genome_format_converters/__init__.py
touch src/genome_format_converters/converters/__init__.py

# 7. Create cli.py with the unified command-line interface (shortened version)
cat > src/genome_format_converters/cli.py << 'EOF'
#!/usr/bin/env python3
"""
Unified command‑line interface for genome-format-converters.
"""

import argparse
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Genome Format Converters – a collection of bioinformatics file conversion tools",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", help="Subcommands")

    def add_io_args(subp):
        subp.add_argument("--input-dir", required=True, help="Directory containing input files")
        subp.add_argument("--output-dir", required=True, help="Output directory")

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
    p = subparsers.add_parser("gff3-to-table", help="Convert GFF3 to tab‑separated table")
    add_io_args(p)

    # GFF3+FASTA to protein
    p = subparsers.add_parser("gff3-to-protein", help="Extract protein sequences from GFF3 and FASTA")
    add_io_args(p)

    # FASTA to FASTQ
    p = subparsers.add_parser("fasta-to-fastq", help="Convert FASTA to FASTQ with default quality")
    add_io_args(p)
    p.add_argument("--qual-char", default="I", help="Quality character (default I = 40)")

    # FASTQ to FASTA
    p = subparsers.add_parser("fastq-to-fasta", help="Convert FASTQ to FASTA")
    add_io_args(p)

    # FASTA+QUAL to FASTQ
    p = subparsers.add_parser("fasta-qual-to-fastq", help="Combine FASTA and QUAL into FASTQ")
    add_io_args(p)

    # FASTQ to FASTA+QUAL
    p = subparsers.add_parser("fastq-to-fasta-qual", help="Split FASTQ into FASTA and QUAL")
    add_io_args(p)

    # Alignment converter
    p = subparsers.add_parser("convert-alignment", help="Convert alignment between formats")
    add_io_args(p)
    p.add_argument("--in-format", required=True, choices=["fasta", "phylip", "nexus", "clustal"])
    p.add_argument("--out-format", required=True, choices=["fasta", "phylip", "nexus", "clustal"])

    # FASTA to table
    p = subparsers.add_parser("fasta-to-table", help="Convert FASTA to two‑column TSV")
    add_io_args(p)

    # BAM to BED
    p = subparsers.add_parser("bam-to-bed", help="Convert BAM/SAM to BED6")
    add_io_args(p)

    # BLAST tab to links
    p = subparsers.add_parser("blast-to-links", help="Convert BLAST tabular to link TSV")
    add_io_args(p)
    p.add_argument("--min-length", type=int, default=0, help="Minimum alignment length")
    p.add_argument("--min-identity", type=float, default=0, help="Minimum percent identity")

    # Delta to tab
    p = subparsers.add_parser("delta-to-tab", help="Convert MUMmer .delta to tabular")
    add_io_args(p)

    # MAF to XMFA
    p = subparsers.add_parser("maf-to-xmfa", help="Convert MAF to XMFA")
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
    p.add_argument("--in-format", required=True, choices=["newick", "nexus", "phyloxml"])
    p.add_argument("--out-format", required=True, choices=["newick", "nexus", "phyloxml"])

    # Tree annotation
    p = subparsers.add_parser("annotate-tree", help="Add alignment sequences to tree")
    p.add_argument("--tree", required=True, help="Input tree file (Newick)")
    p.add_argument("--aln", required=True, help="Input alignment file (FASTA)")
    p.add_argument("--output", required=True, help="Output NEXUS file")

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    # TODO: Import the appropriate batch functions and call them.
    print(f"Command '{args.command}' not yet implemented.")
    print("You need to refactor the original scripts into functions and import them here.")

if __name__ == "__main__":
    main()
EOF

echo "Done. New structure created."
echo "Next steps:"
echo " 1. Refactor each script to provide a batch function (e.g., batch_convert_*)."
echo " 2. Import those functions in cli.py and uncomment the calls."
echo " 3. Install the package in development mode: pip install -e ."
echo " 4. Test with the data in tests/test_data/."