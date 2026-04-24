#!/usr/bin/env python3
"""Unified command-line interface for genome-format-converters.

Author: Benjamin Narh-Madey
Affiliation: Laboratory of Genetics (Hittinger Lab),
             Wisconsin Energy Institute,
             University of Wisconsin-Madison.
Email: narhmadey@wisc.edu

This package provides a collection of tools to convert common
bioinformatics file formats (GFF3, GTF, BED, GenBank, FASTA,
FASTQ, QUAL, alignment formats, BAM, VCF, MAF, MUMmer delta,
EIGENSTRAT (.geno/.snp/.ind, diploid and pseudohaploid),
phylogenetic trees, etc.) with a simple, uniform interface.
"""

import argparse
import sys
from pathlib import Path

from genome_format_converters import __version__
from genome_format_converters._banner import get_banner, print_banner
from genome_format_converters._guide import (
    EXAMPLES as _EXAMPLES,
    render_examples,
    render_guide,
)
from genome_format_converters.converters import (
    gff3_to_gtf,
    gff3_to_bed,
    gff3_to_bed12,
    genbank_to_gff3,
    gff3_to_table,
    gff3_to_protein,
    gtf_to_gff3,
    fasta_to_fastq,
    fastq_to_fasta,
    fasta_qual_to_fastq,
    fastq_to_fasta_qual,
    convert_alignment,
    fasta_to_table,
    stockholm_to_fasta,
    fasta_to_stockholm,
    bam_to_bed,
    blast_tab_to_links,
    delta_to_tab,
    maf_to_xmfa,
    hmmer_tblout_to_tsv,
    orthogroups_to_fasta,
    vcf_to_bed,
    vcf_to_table,
    vcf_to_consensus,
    vcf_to_eigenstrat,
    vcf_to_pseudohaploid,
    vcf_to_plink,
    tree_convert,
    annotate_tree,
    convert_all_gff_fasta_to_gbk,
)
from genome_format_converters.converters._common import set_verbosity


# ---------- argparse helpers -------------------------------------------------

def _subparser_epilog(cmd: str) -> str:
    """Format a subparser epilog with the canonical examples for `cmd`."""
    examples = _EXAMPLES.get(cmd, [])
    if not examples:
        return ""
    lines = ["examples:"]
    for ex in examples:
        lines.append(f"  {ex}")
    return "\n".join(lines) + "\n"


def _add_io_args(sub, *, input_dir: bool = True, single_file: bool = True) -> None:
    """Add the standard I/O flag group.

    `--input-dir` / `--output-dir` is the batch mode. `--input` / `--output`
    is the single-file shortcut (only meaningful for 1:1 converters; not
    emitted when single_file=False).
    """
    g = sub.add_argument_group("input / output")
    if input_dir:
        g.add_argument("--input-dir", help="Directory containing input files")
        g.add_argument("--output-dir",
                       help="Output directory (created if it does not exist)")
    if single_file:
        g.add_argument("--input", "-i",
                       help="Single input file (alternative to --input-dir). "
                            "Use '-' for stdin where the converter supports it.")
        g.add_argument("--output", "-o",
                       help="Single output file (alternative to --output-dir). "
                            "Use '-' for stdout where the converter supports it.")


def _add_batch_flags(sub) -> None:
    g = sub.add_argument_group("batch controls")
    g.add_argument("--pattern", default=None,
                   help="Restrict batch mode to files matching this glob "
                        "(e.g. '*_clean.vcf.gz').")
    g.add_argument("--force", action="store_true",
                   help="Overwrite the output directory if it already exists "
                        "and contains files.")
    g.add_argument("--dry-run", action="store_true",
                   help="List what would be converted without writing output.")
    g.add_argument("--threads", type=int, default=1,
                   help="Number of parallel worker processes for batch mode "
                        "(default 1). Only affects 1:1 converters.")


def _add_eigenstrat_flags(sub) -> None:
    g = sub.add_argument_group("EIGENSTRAT options")
    g.add_argument("--pop-map", default=None,
                   help="Two-column TSV (sample <TAB> population) for column 3 of .ind")
    g.add_argument("--sex-map", default=None,
                   help="Two-column TSV (sample <TAB> sex M/F/U) for column 2 of .ind")
    g.add_argument("--chrom-map", default=None,
                   help="Two-column TSV (contig <TAB> integer). Required by legacy smartpca.")
    g.add_argument("--default-chrom", type=int, default=None,
                   help="Integer to use for any contig not covered by --chrom-map.")
    g.add_argument("--genetic-map", default=None,
                   help="PLINK-style genetic map (chrom snp_id cM bp) for the "
                        "Morgans column of .snp")
    g.add_argument("--ancestral-fasta", default=None,
                   help="FASTA providing the ancestral allele per site; used to "
                        "polarise ref/alt so EIGENSTRAT 2 == ancestral.")
    g.add_argument("--info-aa", action="store_true",
                   help="Polarise using the AA INFO field of each VCF record.")
    g.add_argument("--transversions-only", action="store_true",
                   help="Drop transitions (A<->G, C<->T). Useful against aDNA deamination.")
    g.add_argument("--min-maf", type=float, default=0.0,
                   help="Minimum minor-allele frequency (default 0 = no filter).")
    g.add_argument("--max-missing", type=float, default=1.0,
                   help="Maximum fraction of missing genotypes per site (default 1 = no filter).")
    g.add_argument("--strict-maps", action="store_true",
                   help="Fail with non-zero exit when a pop/sex/chrom map mentions "
                        "a sample or contig that isn't present in the VCF.")
    g.add_argument("--also-plink", action="store_true",
                   help="Also emit PLINK-compatible .ped / .map alongside the "
                        "EIGENSTRAT triplet (useful for plink --make-bed).")


# ---------- I/O resolution ---------------------------------------------------

def _resolve_io(args, default_exts=None):
    """Map the argparse namespace onto either batch-mode (input_dir, output_dir)
    or single-file-mode (input, output). Enforces mutual exclusion.

    Returns a dict with either `{"mode": "batch", "input_dir": ..., "output_dir": ...}`
    or `{"mode": "single", "input": ..., "output": ...}`.
    """
    has_dir = bool(getattr(args, "input_dir", None) or getattr(args, "output_dir", None))
    has_single = bool(getattr(args, "input", None) or getattr(args, "output", None))

    if has_dir and has_single:
        print("[error] --input-dir/--output-dir cannot be mixed with --input/--output.",
              file=sys.stderr)
        sys.exit(2)
    if not has_dir and not has_single:
        print("[error] Provide either --input-dir + --output-dir or --input + --output.",
              file=sys.stderr)
        sys.exit(2)

    if has_dir:
        if not (args.input_dir and args.output_dir):
            print("[error] Batch mode needs both --input-dir and --output-dir.",
                  file=sys.stderr)
            sys.exit(2)
        return {"mode": "batch", "input_dir": args.input_dir, "output_dir": args.output_dir}
    if not (getattr(args, "input", None) and getattr(args, "output", None)):
        print("[error] Single-file mode needs both --input and --output.",
              file=sys.stderr)
        sys.exit(2)
    return {"mode": "single", "input": args.input, "output": args.output}


def _single_file_batch_wrapper(input_path: str, output_path: str):
    """Package an `--input FILE --output FILE` pair into an ad-hoc
    `--input-dir DIR --output-dir DIR` view so we can reuse `batch_convert`
    for single-file runs without per-converter duplication."""
    in_p = Path(input_path)
    out_p = Path(output_path)
    # We fake a "batch" by pointing batch_convert at the parent dir of the
    # input file and passing --pattern equal to the exact filename. The
    # output is written to the parent dir of the requested output path, then
    # renamed to the exact requested path if the produced filename differs.
    in_parent = in_p.parent
    out_parent = out_p.parent
    pattern = in_p.name
    return str(in_parent), str(out_parent), pattern, out_p.name


# ---------- main -------------------------------------------------------------

def main():
    prog_description = (
        f"{get_banner()}\n"
        "Genome Format Converters – a suite of tools for converting "
        "bioinformatics file formats (GFF3, GTF, BED, GenBank, FASTA, "
        "FASTQ, QUAL, alignment formats, BAM, VCF, MAF, MUMmer delta, "
        "EIGENSTRAT (diploid and pseudohaploid), PLINK binary, Stockholm, "
        "HMMER tblout, OrthoFinder orthogroups, phylogenetic trees, etc.) "
        "with a uniform interface.\n\n"
        "Run `gfc guide` for a categorised cheat sheet or "
        "`gfc examples <subcommand>` for concrete invocations."
    )
    epilog_text = (
        "Author: Benjamin Narh-Madey\n"
        "Affiliation: Laboratory of Genetics (Hittinger Lab),\n"
        "             Wisconsin Energy Institute,\n"
        "             University of Wisconsin-Madison.\n"
        "Email: narhmadey@wisc.edu\n"
        f"\nVersion: {__version__}\n"
        "For bug reports or feature requests, please visit "
        "https://github.com/K-nie/genome-format-converters"
    )

    parser = argparse.ArgumentParser(
        description=prog_description,
        epilog=epilog_text,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("-q", "--quiet", action="store_true",
                        help="Suppress progress messages on stderr.")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Show per-record debug messages on stderr.")

    subparsers = parser.add_subparsers(dest="command", help="Subcommands")

    # ---- annotation conversions -----------------------------------------

    p = subparsers.add_parser("gff3-to-gtf", help="Convert GFF3 to GTF")
    _add_io_args(p); _add_batch_flags(p)

    p = subparsers.add_parser("gff3-to-bed", help="Convert GFF3 to 6-column BED")
    _add_io_args(p); _add_batch_flags(p)

    p = subparsers.add_parser("gff3-to-bed12",
                              help="Convert GFF3 to BED12 (exon blocks per transcript; IGV/UCSC-ready)")
    _add_io_args(p); _add_batch_flags(p)

    p = subparsers.add_parser("gtf-to-gff3",
                              help="Convert GTF (Ensembl-style) to GFF3 (synthesises gene rows)")
    _add_io_args(p); _add_batch_flags(p)

    p = subparsers.add_parser("genbank-to-gff3", help="Convert GenBank to GFF3")
    _add_io_args(p); _add_batch_flags(p)

    p = subparsers.add_parser("gff3-to-table", help="Convert GFF3 to tab-separated feature table")
    _add_io_args(p); _add_batch_flags(p)

    p = subparsers.add_parser("gff3-to-protein",
                              help="Extract protein sequences from paired GFF3 and FASTA")
    _add_io_args(p, single_file=False); _add_batch_flags(p)

    p = subparsers.add_parser("fasta-gff-to-gbk",
                              help="Convert paired FASTA and GFF3 files to GenBank")
    _add_io_args(p, single_file=False); _add_batch_flags(p)

    # ---- sequence conversions ------------------------------------------

    p = subparsers.add_parser("fasta-to-fastq", help="Convert FASTA to FASTQ with default quality")
    _add_io_args(p); _add_batch_flags(p)
    p.add_argument("--qual-char", default="I",
                   help="Quality character (default 'I' = phred 40, Illumina 1.8+)")

    p = subparsers.add_parser("fastq-to-fasta", help="Convert FASTQ to FASTA (drop qualities)")
    _add_io_args(p); _add_batch_flags(p)

    p = subparsers.add_parser("fasta-qual-to-fastq", help="Combine FASTA and QUAL into FASTQ")
    _add_io_args(p, single_file=False); _add_batch_flags(p)

    p = subparsers.add_parser("fastq-to-fasta-qual", help="Split FASTQ into FASTA and QUAL files")
    _add_io_args(p, single_file=False); _add_batch_flags(p)

    p = subparsers.add_parser("convert-alignment", help="Convert alignment between formats")
    _add_io_args(p, single_file=False); _add_batch_flags(p)
    p.add_argument("--in-format", required=True,
                   choices=["fasta", "phylip", "nexus", "clustal"],
                   help="Input alignment format")
    p.add_argument("--out-format", required=True,
                   choices=["fasta", "phylip", "nexus", "clustal"],
                   help="Output alignment format")

    p = subparsers.add_parser("fasta-to-table",
                              help="Convert FASTA to two-column TSV (id, sequence)")
    _add_io_args(p); _add_batch_flags(p)

    p = subparsers.add_parser("stockholm-to-fasta",
                              help="Convert Stockholm (Rfam/Pfam) alignments to FASTA")
    _add_io_args(p, single_file=False); _add_batch_flags(p)

    p = subparsers.add_parser("fasta-to-stockholm",
                              help="Convert FASTA alignment to Stockholm")
    _add_io_args(p); _add_batch_flags(p)

    # ---- alignment / mapping -------------------------------------------

    p = subparsers.add_parser("bam-to-bed", help="Convert BAM/SAM/CRAM to BED6")
    _add_io_args(p); _add_batch_flags(p)

    p = subparsers.add_parser("blast-to-links",
                              help="Convert BLAST tabular (outfmt 6) to link TSV")
    _add_io_args(p); _add_batch_flags(p)
    p.add_argument("--min-length", type=int, default=0,
                   help="Minimum alignment length (default 0 = no filter)")
    p.add_argument("--min-identity", type=float, default=0,
                   help="Minimum percent identity (default 0 = no filter)")

    p = subparsers.add_parser("delta-to-tab",
                              help="Convert MUMmer .delta to tabular coordinates")
    _add_io_args(p); _add_batch_flags(p)

    p = subparsers.add_parser("maf-to-xmfa",
                              help="Convert MAF to XMFA (progressiveMauve format)")
    _add_io_args(p); _add_batch_flags(p)

    p = subparsers.add_parser("hmmer-tblout-to-tsv",
                              help="Convert HMMER --tblout / --domtblout to structured TSV")
    _add_io_args(p); _add_batch_flags(p)
    p.add_argument("--hmmer-format", choices=["tblout", "domtblout"], default="tblout",
                   help="Layout of the HMMER output (default: tblout).")

    p = subparsers.add_parser("orthogroups-to-fasta",
                              help="OrthoFinder Orthogroups.tsv + per-species FASTAs -> one FASTA per orthogroup")
    p.add_argument("--orthogroups", required=True,
                   help="Path to the OrthoFinder Orthogroups.tsv file.")
    p.add_argument("--fasta-dir", required=True,
                   help="Directory containing one FASTA per species, matching the TSV column names.")
    p.add_argument("--output-dir", required=True,
                   help="Output directory for per-orthogroup FASTA files.")
    p.add_argument("--force", action="store_true",
                   help="Overwrite the output directory if non-empty.")

    # ---- variant formats ------------------------------------------------

    p = subparsers.add_parser("vcf-to-bed", help="Convert VCF to BED intervals")
    _add_io_args(p); _add_batch_flags(p)

    p = subparsers.add_parser("vcf-to-table", help="Convert VCF to tab-separated table")
    _add_io_args(p); _add_batch_flags(p)

    p = subparsers.add_parser("vcf-to-consensus",
                              help="Create per-sample consensus FASTA from VCF + reference")
    _add_io_args(p, single_file=False); _add_batch_flags(p)
    p.add_argument("--het-policy", choices=["iupac", "ref", "alt"], default="iupac",
                   help="How to resolve heterozygous SNPs in the consensus (default iupac).")
    p.add_argument("--mask-missing", action="store_true",
                   help="Replace missing-genotype positions with 'N'.")

    p = subparsers.add_parser("vcf-to-eigenstrat",
                              help="Convert VCF to EIGENSTRAT (.geno/.snp/.ind); biallelic SNPs only")
    _add_io_args(p, single_file=False); _add_batch_flags(p); _add_eigenstrat_flags(p)

    p = subparsers.add_parser("vcf-to-pseudohaploid",
                              help="Convert VCF to pseudohaploid EIGENSTRAT (one random allele per het)")
    _add_io_args(p, single_file=False); _add_batch_flags(p); _add_eigenstrat_flags(p)
    p.add_argument("--seed", type=int, default=None,
                   help="Base seed; per-file RNG is seed XOR stable hash(filename) "
                        "so file order doesn't affect output.")

    p = subparsers.add_parser("vcf-to-plink",
                              help="Convert VCF to PLINK binary (.bed/.bim/.fam)")
    _add_io_args(p, single_file=False); _add_batch_flags(p); _add_eigenstrat_flags(p)

    # ---- phylogenetic trees --------------------------------------------

    p = subparsers.add_parser("tree-convert",
                              help="Convert tree formats (newick, nexus, phyloxml)")
    _add_io_args(p, single_file=False); _add_batch_flags(p)
    p.add_argument("--in-format", required=True,
                   choices=["newick", "nexus", "phyloxml"],
                   help="Input tree format")
    p.add_argument("--out-format", required=True,
                   choices=["newick", "nexus", "phyloxml"],
                   help="Output tree format")

    p = subparsers.add_parser("annotate-tree",
                              help="Add alignment sequences to a tree (output NEXUS)")
    p.add_argument("--tree", required=True, help="Input tree file (Newick)")
    p.add_argument("--aln", required=True, help="Input alignment file (FASTA)")
    p.add_argument("--output", required=True, help="Output NEXUS file")

    # ---- help / discovery subcommands ----------------------------------

    p = subparsers.add_parser("guide",
                              help="Print a categorised cheat sheet of every subcommand")

    p = subparsers.add_parser("examples",
                              help="Show concrete example invocations for a subcommand")
    p.add_argument("target", nargs="?",
                   help="Subcommand name (omit to list all subcommands with examples)")

    # Attach canonical examples to every subparser's --help epilog in one
    # pass, rather than touching all 30 add_parser call sites above.
    for _cmd_name, _sub in subparsers.choices.items():
        _epi = _subparser_epilog(_cmd_name)
        if _epi:
            _sub.epilog = _epi
            _sub.formatter_class = argparse.RawDescriptionHelpFormatter

    # ---------- parse + dispatch ---------------------------------------------

    args = parser.parse_args()
    set_verbosity(verbose=getattr(args, "verbose", False),
                  quiet=getattr(args, "quiet", False))

    if args.command is None:
        print_banner(force=True)
        print(
            "Run `gfc --help` for full options, `gfc guide` for a cheat sheet,\n"
            "or `gfc examples <subcommand>` for worked invocations.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Discovery / help subcommands — emit on stdout so the output is pipe-safe.
    if args.command == "guide":
        sys.stdout.write(render_guide())
        return
    if args.command == "examples":
        if args.target is None:
            known = sorted(_EXAMPLES.keys())
            sys.stdout.write("Subcommands with worked examples:\n\n")
            for name in known:
                sys.stdout.write(f"  gfc examples {name}\n")
            return
        sys.stdout.write(render_examples(args.target))
        return

    # `annotate-tree` takes file paths directly; nothing to resolve.
    if args.command == "annotate-tree":
        annotate_tree.annotate_tree(args.tree, args.aln, args.output)
        return

    # `orthogroups-to-fasta` has its own I/O shape (a TSV + a FASTA dir).
    if args.command == "orthogroups-to-fasta":
        orthogroups_to_fasta.convert(
            args.orthogroups, args.fasta_dir, args.output_dir, force=args.force
        )
        return

    io = _resolve_io(args)

    # In single-file mode we reduce to a 1-item "batch" run. This works for
    # any converter whose batch_convert accepts --pattern (all of them do).
    if io["mode"] == "single":
        in_dir, out_dir, pattern, forced_output_name = _single_file_batch_wrapper(
            io["input"], io["output"]
        )
        # force=True so we don't refuse to write into a pre-existing out_dir.
        args.force = True
        # Record a rename step for after the conversion runs (see dispatch).
        rename_after = (out_dir, forced_output_name)
    else:
        in_dir = io["input_dir"]
        out_dir = io["output_dir"]
        pattern = getattr(args, "pattern", None)
        rename_after = None

    if getattr(args, "dry_run", False):
        from genome_format_converters.converters._common import iter_input_files
        # Best-effort listing; exact extension list depends on the subcommand.
        print(f"[dry-run] command: {args.command}", file=sys.stderr)
        print(f"[dry-run] input: {in_dir}  (pattern={pattern!r})", file=sys.stderr)
        print(f"[dry-run] output: {out_dir}", file=sys.stderr)
        # Skip actual conversion.
        return

    dispatch_table = {
        "gff3-to-gtf":        lambda: gff3_to_gtf.batch_convert(
            in_dir, out_dir, force=args.force, pattern=pattern),
        "gff3-to-bed":        lambda: gff3_to_bed.batch_convert(
            in_dir, out_dir, force=args.force, pattern=pattern),
        "gff3-to-bed12":      lambda: gff3_to_bed12.batch_convert(
            in_dir, out_dir, force=args.force, pattern=pattern),
        "gtf-to-gff3":        lambda: gtf_to_gff3.batch_convert(
            in_dir, out_dir, force=args.force, pattern=pattern),
        "genbank-to-gff3":    lambda: genbank_to_gff3.batch_convert(
            in_dir, out_dir, force=args.force, pattern=pattern),
        "gff3-to-table":      lambda: gff3_to_table.batch_convert(
            in_dir, out_dir, force=args.force, pattern=pattern),
        "gff3-to-protein":    lambda: gff3_to_protein.batch_convert(
            in_dir, out_dir, force=args.force, pattern=pattern),
        "fasta-to-fastq":     lambda: fasta_to_fastq.batch_convert(
            in_dir, out_dir, qual_char=args.qual_char,
            force=args.force, pattern=pattern),
        "fastq-to-fasta":     lambda: fastq_to_fasta.batch_convert(
            in_dir, out_dir, force=args.force, pattern=pattern),
        "fasta-qual-to-fastq": lambda: fasta_qual_to_fastq.batch_convert(
            in_dir, out_dir, force=args.force, pattern=pattern),
        "fastq-to-fasta-qual": lambda: fastq_to_fasta_qual.batch_convert(
            in_dir, out_dir, force=args.force, pattern=pattern),
        "convert-alignment":  lambda: convert_alignment.batch_convert(
            in_dir, out_dir, args.in_format, args.out_format,
            force=args.force, pattern=pattern),
        "fasta-to-table":     lambda: fasta_to_table.batch_convert(
            in_dir, out_dir, force=args.force, pattern=pattern),
        "stockholm-to-fasta": lambda: stockholm_to_fasta.batch_convert(
            in_dir, out_dir, force=args.force, pattern=pattern),
        "fasta-to-stockholm": lambda: fasta_to_stockholm.batch_convert(
            in_dir, out_dir, force=args.force, pattern=pattern),
        "bam-to-bed":         lambda: bam_to_bed.batch_convert(
            in_dir, out_dir, force=args.force, pattern=pattern),
        "blast-to-links":     lambda: blast_tab_to_links.batch_convert(
            in_dir, out_dir,
            min_length=args.min_length, min_identity=args.min_identity,
            force=args.force, pattern=pattern),
        "delta-to-tab":       lambda: delta_to_tab.batch_convert(
            in_dir, out_dir, force=args.force, pattern=pattern),
        "maf-to-xmfa":        lambda: maf_to_xmfa.batch_convert(
            in_dir, out_dir, force=args.force, pattern=pattern),
        "hmmer-tblout-to-tsv": lambda: hmmer_tblout_to_tsv.batch_convert(
            in_dir, out_dir, hmmer_format=args.hmmer_format,
            force=args.force, pattern=pattern),
        "vcf-to-bed":         lambda: vcf_to_bed.batch_convert(
            in_dir, out_dir, force=args.force, pattern=pattern),
        "vcf-to-table":       lambda: vcf_to_table.batch_convert(
            in_dir, out_dir, force=args.force, pattern=pattern),
        "vcf-to-consensus":   lambda: vcf_to_consensus.batch_convert(
            in_dir, out_dir,
            het_policy=args.het_policy, mask_missing=args.mask_missing,
            force=args.force, pattern=pattern),
        "vcf-to-eigenstrat":  lambda: vcf_to_eigenstrat.batch_convert(
            in_dir, out_dir,
            pop_map=args.pop_map, sex_map=args.sex_map,
            chrom_map=args.chrom_map, default_chrom=args.default_chrom,
            genetic_map=args.genetic_map,
            ancestral_fasta=args.ancestral_fasta, info_aa=args.info_aa,
            transversions_only=args.transversions_only,
            min_maf=args.min_maf, max_missing=args.max_missing,
            strict_maps=args.strict_maps, also_plink=args.also_plink,
            force=args.force, pattern=pattern),
        "vcf-to-pseudohaploid": lambda: vcf_to_pseudohaploid.batch_convert(
            in_dir, out_dir,
            pop_map=args.pop_map, sex_map=args.sex_map,
            chrom_map=args.chrom_map, default_chrom=args.default_chrom,
            genetic_map=args.genetic_map,
            ancestral_fasta=args.ancestral_fasta, info_aa=args.info_aa,
            transversions_only=args.transversions_only,
            min_maf=args.min_maf, max_missing=args.max_missing,
            strict_maps=args.strict_maps, also_plink=args.also_plink,
            seed=args.seed,
            force=args.force, pattern=pattern),
        "vcf-to-plink":       lambda: vcf_to_plink.batch_convert(
            in_dir, out_dir,
            pop_map=args.pop_map, sex_map=args.sex_map,
            chrom_map=args.chrom_map, default_chrom=args.default_chrom,
            genetic_map=args.genetic_map,
            ancestral_fasta=args.ancestral_fasta, info_aa=args.info_aa,
            transversions_only=args.transversions_only,
            min_maf=args.min_maf, max_missing=args.max_missing,
            strict_maps=args.strict_maps,
            force=args.force, pattern=pattern),
        "tree-convert":       lambda: tree_convert.batch_convert(
            in_dir, out_dir, args.in_format, args.out_format,
            force=args.force, pattern=pattern),
        "fasta-gff-to-gbk":   lambda: convert_all_gff_fasta_to_gbk.batch_convert(
            in_dir, out_dir, force=args.force, pattern=pattern),
    }

    if args.command not in dispatch_table:
        print(f"[error] Unknown command: {args.command}", file=sys.stderr)
        sys.exit(2)

    # --threads: where supported, batch_convert itself parallelises. For now
    # we accept the flag and warn when the converter hasn't been parallelised.
    if getattr(args, "threads", 1) > 1:
        print(
            f"[warn] --threads {args.threads} requested, but per-file parallelism "
            "is not yet wired through every converter. Running sequentially.",
            file=sys.stderr,
        )

    dispatch_table[args.command]()

    # Single-file mode rename — the batch_convert wrote into out_dir using
    # its own derived filename; rename to the user's requested --output path
    # if the names differ.
    if rename_after is not None:
        out_dir_p, forced_name = rename_after
        out_dir_p = Path(out_dir_p)
        candidates = [p for p in out_dir_p.iterdir()
                      if p.is_file() and p.name != forced_name]
        # If exactly one new file was produced and its name differs, rename it.
        if len(candidates) == 1 and candidates[0].name != forced_name:
            candidates[0].rename(out_dir_p / forced_name)


if __name__ == "__main__":
    main()
