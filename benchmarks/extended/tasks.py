"""Registry of the extended benchmark: the 24 subcommands not covered by T1-T8.

Each task names one gfc subcommand, the input directory that
``prepare_inputs.sh`` builds under ``benchmarks/data/extended/<id>/``, the
gfc command, the other tools that perform the same conversion, the tool
whose output is the reference, and the comparator in ``compare_ext.py``.

Command templates are formatted with:
    {IN}    the task's input directory
    {OUT}   a fresh output directory for this tool and replicate
    {BIOCONVERT}  path to the bioconvert executable (env GFC_BENCH_BIOCONVERT)

``reference`` is a tool name, or ``"input"`` for a reference-free check of
gfc's output against the input files (used where no other tool performs the
conversion, as for T8).

A tool whose executables are not on PATH is skipped and recorded as such.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List


@dataclass
class Tool:
    name: str
    cmd: str
    version: str = ""          # shell command printing the version
    needs: List[str] = field(default_factory=list)  # executables required
    route: str = "full"        # "full" or "partial" (Supplementary Table S7)


@dataclass
class Task:
    id: str
    subcommand: str
    conversion: str
    domain: str
    input_note: str
    gfc: str
    tools: List[Tool] = field(default_factory=list)
    reference: str = "input"
    comparator: str = ""
    options: Dict[str, str] = field(default_factory=dict)


def loop(ext: str, body: str, src: str = "{IN}") -> str:
    """Shell loop over every input file with extension ``ext``; ``$f`` is the
    file and ``$s`` its stem."""
    return (f'for f in {src}/*{ext}; do s=$(basename "$f" {ext}); {body}; done')


BC = "{BIOCONVERT}"
SEQRET_V = "seqret -version 2>&1 | tail -1"
SEQKIT_V = "seqkit version"
GFFREAD_V = "gffread --version 2>&1"
AGAT_V = ("agat_convert_sp_gff2gtf.pl --help 2>&1 | grep -i 'version' | head -1 "
          "| tr -s ' ' | cut -d':' -f2- | tr -d ' |'")
BC_V = "{BIOCONVERT} --version 2>&1 | tail -1"
ESL_V = "esl-reformat -h | sed -n 2p"
BCF_V = "bcftools --version | head -1"
PLINK2_V = "plink2 --version"


TASKS: List[Task] = [
    # ------------------------------------------------------------ sequences
    Task("E01", "fasta-to-fastq", "FASTA -> FASTQ", "Sequences & reads",
         "500,000 simulated 150-bp reads from S. cerevisiae R64, as FASTA",
         "gfc fasta-to-fastq --input-dir {IN} --output-dir {OUT} --force",
         [Tool("EMBOSS-seqret",
               loop(".fasta", 'seqret -sequence "$f" -outseq {OUT}/$s.fastq -osformat2 fastq-sanger -auto'),
               SEQRET_V, ["seqret"]),
          Tool("BioConvert", loop(".fasta", BC + ' fasta2fastq "$f" {OUT}/$s.fastq --force'),
               BC_V, [BC])],
         reference="EMBOSS-seqret", comparator="fastq", options={"fixed_quality_ok": "1"}),
    Task("E02", "fastq-to-fasta", "FASTQ -> FASTA", "Sequences & reads",
         "500,000 simulated 150-bp reads from S. cerevisiae R64 (ART, HiSeq 2500 profile)",
         "gfc fastq-to-fasta --input-dir {IN} --output-dir {OUT} --force",
         [Tool("seqkit", loop(".fastq", 'seqkit fq2fa "$f" -o {OUT}/$s.fasta'), SEQKIT_V, ["seqkit"]),
          Tool("EMBOSS-seqret",
               loop(".fastq", 'seqret -sequence "$f" -sformat1 fastq-sanger -outseq {OUT}/$s.fasta -osformat2 fasta -auto'),
               SEQRET_V, ["seqret"]),
          Tool("samtools", loop(".fastq", 'samtools import -0 "$f" | samtools fasta - > {OUT}/$s.fasta'),
               "samtools --version | head -1", ["samtools"]),
          Tool("BioConvert", loop(".fastq", BC + ' fastq2fasta "$f" {OUT}/$s.fasta --force'), BC_V, [BC])],
         reference="seqkit", comparator="fasta"),
    Task("E03", "fastq-to-fasta-qual", "FASTQ -> FASTA + QUAL", "Sequences & reads",
         "as E02",
         "gfc fastq-to-fasta-qual --input-dir {IN} --output-dir {OUT} --force",
         [Tool("py-ref", loop(".fastq", 'python {REFS}/fastq_to_fasta_qual.py "$f" {OUT}/$s'),
               "python -c 'import Bio; print(\"biopython\", Bio.__version__)'"),
          Tool("BioConvert",
               loop(".fastq", BC + ' fastq2fasta_qual "$f" {OUT}/$s.fasta {OUT}/$s.qual --force'),
               BC_V, [BC])],
         reference="py-ref", comparator="fasta_qual"),
    Task("E04", "fasta-qual-to-fastq", "FASTA + QUAL -> FASTQ", "Sequences & reads",
         "the E03 reads as a FASTA and a QUAL file",
         "gfc fasta-qual-to-fastq --input-dir {IN} --output-dir {OUT} --force",
         [Tool("py-ref", loop(".fasta", 'python {REFS}/fasta_qual_to_fastq.py "$f" {IN}/$s.qual {OUT}/$s.fastq'),
               "python -c 'import Bio; print(\"biopython\", Bio.__version__)'"),
          Tool("BioConvert",
               loop(".fasta", BC + ' fasta_qual2fastq "$f" {IN}/$s.qual {OUT}/$s.fastq --force'),
               BC_V, [BC])],
         reference="py-ref", comparator="fastq"),
    Task("E05", "fasta-to-table", "FASTA -> TSV", "Sequences & reads",
         "as E01",
         "gfc fasta-to-table --input-dir {IN} --output-dir {OUT} --force",
         [Tool("seqkit", loop(".fasta", 'seqkit fx2tab "$f" -o {OUT}/$s.tsv'), SEQKIT_V, ["seqkit"])],
         reference="seqkit", comparator="seq_table"),
    Task("E06", "bam-to-bed", "BAM -> BED6", "Sequences & reads",
         "the E02 reads aligned to R64 with minimap2, coordinate-sorted BAM",
         "gfc bam-to-bed --input-dir {IN} --output-dir {OUT} --pattern '*.bam' --force",
         [Tool("bedtools", loop(".bam", 'bedtools bamtobed -i "$f" > {OUT}/$s.bed'),
               "bedtools --version", ["bedtools"])],
         reference="bedtools", comparator="bed6"),
    # ------------------------------------------------------------ annotation
    Task("E07", "gtf-to-gff3", "GTF -> GFF3", "Annotation",
         "Ensembl 112 GTF of S. cerevisiae R64 and of GRCh38 chromosome 22",
         "gfc gtf-to-gff3 --input-dir {IN} --output-dir {OUT} --force",
         [Tool("gffread", loop(".gtf", 'gffread "$f" -o {OUT}/$s.gff3'), GFFREAD_V, ["gffread"]),
          Tool("AGAT", "cd {OUT} && " + loop(".gtf", 'agat_convert_sp_gxf2gxf.pl -g "$f" -o $s.gff3 > /dev/null'),
               AGAT_V, ["agat_convert_sp_gxf2gxf.pl"])],
         reference="gffread", comparator="gff3_models"),
    Task("E08", "gff3-to-bed", "GFF3 -> BED6", "Annotation",
         "Ensembl 112 GFF3 of S. cerevisiae R64 and of GRCh38 chromosome 22",
         "gfc gff3-to-bed --input-dir {IN} --output-dir {OUT} --pattern '*.gff3' --force",
         [Tool("gffread", loop(".gff3", 'gffread "$f" --bed | cut -f1-6 > {OUT}/$s.bed'),
               GFFREAD_V, ["gffread"], route="partial")],
         reference="input", comparator="gff3_bed6"),
    Task("E09", "gff3-to-table", "GFF3 -> TSV", "Annotation",
         "as E08",
         "gfc gff3-to-table --input-dir {IN} --output-dir {OUT} --pattern '*.gff3' --force",
         [Tool("gffread", loop(".gff3", 'gffread "$f" --table @chr,@start,@end,@strand,@id > {OUT}/$s.tsv'),
               GFFREAD_V, ["gffread"]),
          Tool("AGAT", "cd {OUT} && " + loop(".gff3", 'agat_convert_sp_gff2tsv.pl -gff "$f" -o $s.tsv > /dev/null'),
               AGAT_V, ["agat_convert_sp_gff2tsv.pl"])],
         reference="input", comparator="gff3_table"),
    Task("E10", "gff3-to-protein", "GFF3 + FASTA -> protein FASTA", "Annotation",
         "Ensembl 112 GFF3 and genome FASTA of S. cerevisiae R64 (nuclear chromosomes)",
         "gfc gff3-to-protein --input-dir {IN} --output-dir {OUT} --force",
         [Tool("gffread", loop(".gff3", 'gffread "$f" -g {IN}/$s.fasta -y {OUT}/$s.faa'),
               GFFREAD_V, ["gffread"]),
          Tool("AGAT", "cd {OUT} && " + loop(".gff3",
               'agat_sp_extract_sequences.pl -g "$f" -f {IN}/$s.fasta -p -o $s.faa > /dev/null'),
               AGAT_V, ["agat_sp_extract_sequences.pl"])],
         reference="gffread", comparator="protein"),
    Task("E11", "genbank-to-gff3", "GenBank -> GFF3", "Annotation",
         "NCBI RefSeq GenBank flat file of S. cerevisiae R64 (GCF_000146045.2)",
         "gfc genbank-to-gff3 --input-dir {IN} --output-dir {OUT} --force",
         [Tool("BioConvert", loop(".gbk", BC + ' genbank2gff3 "$f" {OUT}/$s.gff3 --force'), BC_V, [BC])],
         reference="input", comparator="genbank_vs_input", options={"compare_others": "1"}),
    # ------------------------------------------------------------ alignments
    Task("E12", "stockholm-to-fasta", "Stockholm -> aligned FASTA", "Alignments",
         "500 Pfam 37.0 seed alignments, one Stockholm file per family",
         "gfc stockholm-to-fasta --input-dir {IN} --output-dir {OUT} --force",
         [Tool("easel", loop(".sto", 'esl-reformat afa "$f" > {OUT}/$s.fasta'), ESL_V, ["esl-reformat"]),
          Tool("EMBOSS-seqret",
               loop(".sto", 'seqret -sequence "$f" -sformat1 stockholm -outseq {OUT}/$s.fasta -osformat2 fasta -auto'),
               SEQRET_V, ["seqret"]),
          # BioConvert has no direct route: Stockholm -> Clustal -> FASTA.
          Tool("BioConvert", loop(".sto", BC + ' stockholm2clustal "$f" {OUT}/$s.clw --force && '
                                  + BC + ' clustal2fasta {OUT}/$s.clw {OUT}/$s.fasta --force && rm {OUT}/$s.clw'),
               BC_V, [BC])],
         reference="easel", comparator="alignment", options={"cand_fmt": "fasta", "ref_fmt": "fasta"}),
    Task("E13", "fasta-to-stockholm", "aligned FASTA -> Stockholm", "Alignments",
         "the E12 families as aligned FASTA",
         "gfc fasta-to-stockholm --input-dir {IN} --output-dir {OUT} --force",
         [Tool("easel", loop(".fasta", 'esl-reformat stockholm "$f" > {OUT}/$s.sto'), ESL_V, ["esl-reformat"]),
          # BioConvert has no direct route: FASTA -> Clustal -> Stockholm.
          Tool("BioConvert", loop(".fasta", BC + ' fasta2clustal "$f" {OUT}/$s.clw --force && '
                                  + BC + ' clustal2stockholm {OUT}/$s.clw {OUT}/$s.sto --force && rm {OUT}/$s.clw'),
               BC_V, [BC])],
         reference="easel", comparator="alignment", options={"cand_fmt": "stockholm", "ref_fmt": "stockholm"}),
    Task("E14", "convert-alignment", "aligned FASTA -> PHYLIP", "Alignments",
         "as E13",
         "gfc convert-alignment --input-dir {IN} --output-dir {OUT} --in-format fasta --out-format phylip-relaxed --force",
         [Tool("EMBOSS-seqret",
               loop(".fasta", 'seqret -sequence "$f" -outseq {OUT}/$s.phylip -osformat2 phylip -auto'),
               SEQRET_V, ["seqret"]),
          Tool("easel", loop(".fasta", 'esl-reformat phylip "$f" > {OUT}/$s.phylip'), ESL_V,
               ["esl-reformat"]),
          Tool("BioConvert", loop(".fasta", BC + ' fasta2phylip "$f" {OUT}/$s.phylip --force'), BC_V, [BC])],
         reference="EMBOSS-seqret", comparator="alignment",
         options={"cand_fmt": "phylip-relaxed", "ref_fmt": "phylip"}),
    # ------------------------------------------------------------ phylogenetics
    Task("E15", "tree-convert", "Newick -> NEXUS", "Phylogenetics",
         "FastTree trees of 200 of the E12 Pfam families",
         "gfc tree-convert --input-dir {IN} --output-dir {OUT} --in-format newick --out-format nexus --force",
         [Tool("BioConvert", loop(".nwk", BC + ' newick2nexus "$f" {OUT}/$s.nexus --force'), BC_V, [BC])],
         reference="input", comparator="tree_vs_input", options={"compare_others": "1"}),
    Task("E16", "annotate-tree", "Newick + alignment -> NEXUS", "Phylogenetics",
         "the E15 trees with their aligned FASTA",
         loop(".nwk", 'gfc annotate-tree --tree "$f" --aln {IN}/$s.fasta --output {OUT}/$s.nex'),
         [], reference="input", comparator="annotated_tree"),
    # ------------------------------------------------------------ comparative genomics
    Task("E17", "blast-to-links", "BLAST tabular -> links", "Comparative genomics & search",
         "DIAMOND all-against-all blastp of the Ensembl 112 R64 proteome (outfmt 6)",
         "gfc blast-to-links --input-dir {IN} --output-dir {OUT} --force",
         [], reference="input", comparator="blast_links"),
    Task("E18", "delta-to-tab", "MUMmer delta -> TSV", "Comparative genomics & search",
         "nucmer alignment of GRCh37 against GRCh38 chromosome 22 (Ensembl 75 and 112)",
         "gfc delta-to-tab --input-dir {IN} --output-dir {OUT} --force",
         [Tool("show-coords", loop(".delta", 'show-coords -T -H "$f" > {OUT}/$s.coords'),
               "nucmer --version 2>&1 | tail -1", ["show-coords"])],
         reference="show-coords", comparator="delta"),
    Task("E19", "maf-to-xmfa", "MAF -> XMFA", "Comparative genomics & search",
         "LAST alignment of GRCh37 against a 5 Mb window of GRCh38 chromosome 22 (25-30 Mb), as MAF",
         "gfc maf-to-xmfa --input-dir {IN} --output-dir {OUT} --force",
         [], reference="input", comparator="xmfa"),
    # ------------------------------------------------------------ variants
    Task("E20", "vcf-to-bed", "VCF -> BED", "Variants",
         "1000 Genomes phase 3 chromosome 22, biallelic SNPs (as T1)",
         "gfc vcf-to-bed --input-dir {IN} --output-dir {OUT} --pattern '*.vcf.gz' --force",
         [Tool("bcftools", loop(".vcf.gz", "bcftools query -f '%CHROM\\t%POS0\\t%END\\t%ID\\n' \"$f\" > {OUT}/$s.bed"),
               BCF_V, ["bcftools"]),
          Tool("BioConvert", loop(".vcf.gz", BC + ' vcf2bed "$f" {OUT}/$s.bed --force'), BC_V, [BC])],
         reference="bcftools", comparator="vcf_bed"),
    Task("E21", "vcf-to-table", "VCF -> TSV", "Variants",
         "as E20",
         "gfc vcf-to-table --input-dir {IN} --output-dir {OUT} --pattern '*.vcf.gz' --force",
         [Tool("bcftools", loop(".vcf.gz",
               "bcftools query -H -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\n' \"$f\" > {OUT}/$s.tsv"),
               BCF_V, ["bcftools"])],
         reference="bcftools", comparator="vcf_table"),
    Task("E22", "vcf-to-consensus", "VCF + FASTA -> consensus FASTA", "Variants",
         "the first 10 samples of the E20 VCF with the GRCh37 chromosome 22 FASTA",
         "gfc vcf-to-consensus --input-dir {IN} --output-dir {OUT} --force",
         [Tool("bcftools",
               loop(".vcf.gz", "for x in $(bcftools query -l \"$f\"); do bcftools consensus -I -s $x "
                    "-f {IN}/$s.fasta \"$f\" 2>/dev/null | sed \"s/^>/>$x|/\"; done > {OUT}/$s.fa"),
               BCF_V, ["bcftools"])],
         reference="bcftools", comparator="consensus"),
    Task("E23", "vcf-to-snp-alignment", "VCF -> SNP alignment (PHYLIP)", "Variants",
         "the first 100,000 sites of the E20 VCF, all 2,504 samples",
         "gfc vcf-to-snp-alignment --input-dir {IN} --output-dir {OUT} --pattern '*.vcf.gz' --format phylip --het iupac --force",
         [Tool("plink2", loop(".vcf.gz", 'plink2 --vcf "$f" --export phylip --out {OUT}/$s --threads 1 --memory 640'),
               PLINK2_V, ["plink2"], route="partial")],
         reference="plink2", comparator="snp_alignment"),
    Task("E24", "eigenstrat-to-vcf", "EIGENSTRAT -> VCF", "Variants",
         "the E23 sites as EIGENSTRAT (text .geno for gfc, packed .geno for plink2)",
         "gfc eigenstrat-to-vcf --input-dir {IN}/text --output-dir {OUT} --force",
         [Tool("plink2", loop(".snp", '"${GFC_BENCH_PLINK2_EIG:-plink2}" --eigfile {IN}/packed/$s --export vcf --out {OUT}/$s '
                              '--threads 1 --memory 640', src="{IN}/packed"),
               PLINK2_V, ["plink2"])],
         reference="plink2", comparator="vcf_genotypes"),
]

BY_ID = {t.id: t for t in TASKS}
