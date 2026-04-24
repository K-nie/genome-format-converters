"""Categorised cheat sheet and per-subcommand examples.

Surfaced through the ``gfc guide`` and ``gfc examples <subcommand>`` CLI
subcommands.  Kept in-module (rather than reading README.md at runtime) so
the package is self-contained and the guide is available offline.
"""

from __future__ import annotations

from typing import Dict, List


# One-line example invocations per subcommand. The first entry is the
# canonical / happy-path example; later entries show common flags.
EXAMPLES: Dict[str, List[str]] = {
    # Annotation --------------------------------------------------------
    "gff3-to-gtf": [
        "gfc gff3-to-gtf --input-dir ./gff --output-dir ./gtf",
        "gfc gff3-to-gtf --input tiny.gff3 --output tiny.gtf",
    ],
    "gtf-to-gff3": [
        "gfc gtf-to-gff3 --input-dir ./gtf --output-dir ./gff3",
        "# Ensembl Fungi ships GTFs with no gene rows — the converter synthesises them.",
    ],
    "gff3-to-bed": [
        "gfc gff3-to-bed --input-dir ./gff --output-dir ./bed",
    ],
    "gff3-to-bed12": [
        "gfc gff3-to-bed12 --input-dir ./gff --output-dir ./bed12",
        "# BED12 is what IGV / UCSC expect for multi-exon transcripts.",
    ],
    "genbank-to-gff3": [
        "gfc genbank-to-gff3 --input-dir ./gbk --output-dir ./gff3",
    ],
    "gff3-to-table": [
        "gfc gff3-to-table --input-dir ./gff --output-dir ./tsv",
    ],
    "gff3-to-protein": [
        "gfc gff3-to-protein --input-dir ./genome --output-dir ./proteins",
        "# --input-dir must contain paired *.fasta + *.gff3 files.",
    ],
    "fasta-gff-to-gbk": [
        "gfc fasta-gff-to-gbk --input-dir ./genome --output-dir ./gbk",
    ],
    # Sequence ----------------------------------------------------------
    "fasta-to-fastq": [
        "gfc fasta-to-fastq --input-dir ./fa --output-dir ./fq --qual-char I",
    ],
    "fastq-to-fasta": [
        "gfc fastq-to-fasta --input-dir ./fq --output-dir ./fa",
        "gfc fastq-to-fasta --input reads.fastq --output reads.fasta",
    ],
    "fasta-qual-to-fastq": [
        "gfc fasta-qual-to-fastq --input-dir ./sanger --output-dir ./fq",
    ],
    "fastq-to-fasta-qual": [
        "gfc fastq-to-fasta-qual --input-dir ./fq --output-dir ./split",
    ],
    "convert-alignment": [
        "gfc convert-alignment --input-dir ./aln --output-dir ./phy "
        "--in-format fasta --out-format phylip",
    ],
    "fasta-to-table": [
        "gfc fasta-to-table --input-dir ./fa --output-dir ./tsv",
    ],
    "stockholm-to-fasta": [
        "gfc stockholm-to-fasta --input-dir ./rfam --output-dir ./fa",
    ],
    "fasta-to-stockholm": [
        "gfc fasta-to-stockholm --input-dir ./fa --output-dir ./sto",
    ],
    # Alignment / mapping -----------------------------------------------
    "bam-to-bed": [
        "gfc bam-to-bed --input-dir ./bam --output-dir ./bed",
    ],
    "blast-to-links": [
        "gfc blast-to-links --input-dir ./blast --output-dir ./links "
        "--min-length 100 --min-identity 30",
    ],
    "delta-to-tab": [
        "gfc delta-to-tab --input-dir ./delta --output-dir ./tsv",
    ],
    "maf-to-xmfa": [
        "gfc maf-to-xmfa --input-dir ./maf --output-dir ./xmfa",
    ],
    "hmmer-tblout-to-tsv": [
        "gfc hmmer-tblout-to-tsv --input-dir ./hmmer --output-dir ./tsv",
        "gfc hmmer-tblout-to-tsv --input-dir ./dom --output-dir ./tsv "
        "--hmmer-format domtblout",
    ],
    "orthogroups-to-fasta": [
        "gfc orthogroups-to-fasta --orthogroups Orthogroups.tsv "
        "--fasta-dir ./proteomes --output-dir ./per_og",
    ],
    # Variant -----------------------------------------------------------
    "vcf-to-bed": [
        "gfc vcf-to-bed --input-dir ./vcf --output-dir ./bed",
    ],
    "vcf-to-table": [
        "gfc vcf-to-table --input-dir ./vcf --output-dir ./tsv",
    ],
    "vcf-to-consensus": [
        "gfc vcf-to-consensus --input-dir ./data --output-dir ./consensus",
        "gfc vcf-to-consensus --input-dir ./data --output-dir ./out "
        "--het-policy iupac --mask-missing",
    ],
    "vcf-to-eigenstrat": [
        "gfc vcf-to-eigenstrat --input-dir ./vcf --output-dir ./eig",
        "gfc vcf-to-eigenstrat --input-dir ./vcf --output-dir ./eig "
        "--pop-map pops.tsv --chrom-map chroms.tsv --min-maf 0.05",
        "gfc vcf-to-eigenstrat --input-dir ./vcf --output-dir ./eig --also-plink",
    ],
    "vcf-to-pseudohaploid": [
        "gfc vcf-to-pseudohaploid --input-dir ./vcf --output-dir ./eig --seed 42",
        "# Per-file RNG seeds so order of files doesn't change the output.",
    ],
    "vcf-to-plink": [
        "gfc vcf-to-plink --input-dir ./vcf --output-dir ./plink --min-maf 0.05",
        "# Produces binary .bed/.bim/.fam; use plink --bfile <stem> downstream.",
    ],
    # Phylogenetic trees ------------------------------------------------
    "tree-convert": [
        "gfc tree-convert --input-dir ./trees --output-dir ./nex "
        "--in-format newick --out-format nexus",
    ],
    "annotate-tree": [
        "gfc annotate-tree --tree species.nwk --aln concat.fasta "
        "--output species.annotated.nex",
    ],
}


# Category grouping for `gfc guide`.
CATEGORIES: List[tuple] = [
    ("Annotation formats",
     ["gff3-to-gtf", "gtf-to-gff3", "gff3-to-bed", "gff3-to-bed12",
      "genbank-to-gff3", "gff3-to-table", "gff3-to-protein",
      "fasta-gff-to-gbk"]),
    ("Sequence formats",
     ["fasta-to-fastq", "fastq-to-fasta", "fasta-qual-to-fastq",
      "fastq-to-fasta-qual", "convert-alignment", "fasta-to-table",
      "stockholm-to-fasta", "fasta-to-stockholm"]),
    ("Alignment / mapping",
     ["bam-to-bed", "blast-to-links", "delta-to-tab", "maf-to-xmfa",
      "hmmer-tblout-to-tsv", "orthogroups-to-fasta"]),
    ("Variant formats (VCF / EIGENSTRAT / PLINK)",
     ["vcf-to-bed", "vcf-to-table", "vcf-to-consensus",
      "vcf-to-eigenstrat", "vcf-to-pseudohaploid", "vcf-to-plink"]),
    ("Phylogenetic trees",
     ["tree-convert", "annotate-tree"]),
]


# One-sentence summaries used by `gfc guide`.
SUMMARIES: Dict[str, str] = {
    "gff3-to-gtf": "GFF3 → GTF",
    "gtf-to-gff3": "GTF → GFF3 (synthesises gene rows)",
    "gff3-to-bed": "GFF3 → 6-column BED",
    "gff3-to-bed12": "GFF3 → BED12 with exon blocks per transcript",
    "genbank-to-gff3": "GenBank → GFF3",
    "gff3-to-table": "GFF3 → tab-separated feature table",
    "gff3-to-protein": "Paired GFF3 + FASTA → protein FASTA (phase-aware)",
    "fasta-gff-to-gbk": "Paired FASTA + GFF3 → GenBank",
    "fasta-to-fastq": "FASTA → FASTQ with constant quality",
    "fastq-to-fasta": "FASTQ → FASTA (drop qualities)",
    "fasta-qual-to-fastq": "Paired FASTA + QUAL → FASTQ",
    "fastq-to-fasta-qual": "FASTQ → paired FASTA + QUAL",
    "convert-alignment": "FASTA / PHYLIP / NEXUS / CLUSTAL interconversion",
    "fasta-to-table": "FASTA → two-column TSV (id, sequence)",
    "stockholm-to-fasta": "Stockholm (Rfam/Pfam) → FASTA alignment",
    "fasta-to-stockholm": "FASTA alignment → Stockholm",
    "bam-to-bed": "BAM / SAM / CRAM → BED6",
    "blast-to-links": "BLAST tabular (outfmt 6) → link TSV",
    "delta-to-tab": "MUMmer .delta → tabular coordinates",
    "maf-to-xmfa": "MAF → XMFA (progressiveMauve)",
    "hmmer-tblout-to-tsv": "HMMER --tblout / --domtblout → structured TSV",
    "orthogroups-to-fasta": "OrthoFinder Orthogroups.tsv → one FASTA per OG",
    "vcf-to-bed": "VCF / BCF → BED intervals",
    "vcf-to-table": "VCF / BCF → tab-separated table",
    "vcf-to-consensus": "VCF + reference FASTA → per-sample consensus FASTA",
    "vcf-to-eigenstrat": "VCF → EIGENSTRAT triplet (.geno/.snp/.ind)",
    "vcf-to-pseudohaploid": "VCF → pseudohaploid EIGENSTRAT (random allele per het)",
    "vcf-to-plink": "VCF → PLINK binary (.bed/.bim/.fam)",
    "tree-convert": "Newick ↔ NEXUS ↔ PhyloXML",
    "annotate-tree": "Tree + alignment → annotated NEXUS",
}


def render_guide() -> str:
    """Build the categorised cheat sheet shown by `gfc guide`."""
    lines: List[str] = []
    for cat_name, cmds in CATEGORIES:
        lines.append("")
        lines.append(f"{cat_name}")
        lines.append("=" * len(cat_name))
        width = max(len(c) for c in cmds)
        for cmd in cmds:
            summary = SUMMARIES.get(cmd, "")
            lines.append(f"  {cmd:<{width}}   {summary}")
    lines.append("")
    lines.append("Run `gfc examples <subcommand>` to see concrete invocations.")
    lines.append("Run `gfc <subcommand> --help` for full flag reference.")
    lines.append("")
    return "\n".join(lines)


def render_examples(subcommand: str) -> str:
    examples = EXAMPLES.get(subcommand)
    if examples is None:
        known = ", ".join(sorted(EXAMPLES.keys()))
        return (
            f"No examples for {subcommand!r}.\n"
            f"Known subcommands: {known}\n"
        )
    lines = [f"Examples for `{subcommand}`:", ""]
    for ex in examples:
        lines.append(f"    {ex}")
    lines.append("")
    return "\n".join(lines)
