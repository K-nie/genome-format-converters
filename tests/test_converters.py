"""Smoke + golden-file tests for every converter.

Run: ``pytest`` from the package root. Fixtures live under
``tests/test_data/``; each test copies the relevant fixture into a per-test
tempdir so the canonical fixtures stay pristine.
"""

from pathlib import Path
import shutil

import pytest

from genome_format_converters.converters import (
    annotate_tree,
    bam_to_bed,
    blast_tab_to_links,
    convert_alignment,
    convert_all_gff_fasta_to_gbk,
    delta_to_tab,
    fasta_qual_to_fastq,
    fasta_to_fastq,
    fasta_to_stockholm,
    fasta_to_table,
    fastq_to_fasta,
    fastq_to_fasta_qual,
    genbank_to_gff3,
    gff3_to_bed,
    gff3_to_bed12,
    gff3_to_gtf,
    gff3_to_protein,
    gff3_to_table,
    gtf_to_gff3,
    hmmer_tblout_to_tsv,
    maf_to_xmfa,
    orthogroups_to_fasta,
    stockholm_to_fasta,
    tree_convert,
    vcf_to_bed,
    vcf_to_consensus,
    vcf_to_eigenstrat,
    vcf_to_plink,
    vcf_to_pseudohaploid,
    vcf_to_table,
)

FIXTURES = Path(__file__).parent / "test_data"


# ---------- helpers ---------------------------------------------------------

def _stage(tmp_path: Path, *names: str) -> Path:
    in_dir = tmp_path / "in"
    in_dir.mkdir(exist_ok=True)
    for n in names:
        shutil.copy(FIXTURES / n, in_dir / n)
    return in_dir


# ---------- GFF3 family -----------------------------------------------------


def test_gff3_to_bed(tmp_path):
    in_dir = _stage(tmp_path, "tiny.gff3")
    out_dir = tmp_path / "out"
    gff3_to_bed.batch_convert(str(in_dir), str(out_dir))
    rows = (out_dir / "tiny.bed").read_text().strip().splitlines()
    assert any("\tgene1\t" in r for r in rows)
    # Strand column is '+' on this fixture; must be the 6th tab-delimited column.
    for r in rows:
        fields = r.split("\t")
        assert len(fields) == 6
        assert fields[5] in {"+", "-", "."}


def test_gff3_to_gtf_hierarchy(tmp_path):
    in_dir = _stage(tmp_path, "tiny.gff3")
    out_dir = tmp_path / "out"
    gff3_to_gtf.batch_convert(str(in_dir), str(out_dir))
    content = (out_dir / "tiny.gtf").read_text()
    # Gene rows carry gene_id only; transcript / exon / CDS rows carry both.
    gene_lines = [l for l in content.splitlines() if l.split("\t")[2] == "gene"]
    assert gene_lines
    assert all("transcript_id" not in l for l in gene_lines)
    transcript_lines = [l for l in content.splitlines()
                        if l.split("\t")[2] in {"mRNA", "exon", "CDS"}]
    assert transcript_lines
    assert all("gene_id" in l and "transcript_id" in l for l in transcript_lines)


def test_gff3_to_table(tmp_path):
    in_dir = _stage(tmp_path, "tiny.gff3")
    out_dir = tmp_path / "out"
    gff3_to_table.batch_convert(str(in_dir), str(out_dir))
    rows = (out_dir / "tiny.tsv").read_text().strip().splitlines()
    assert rows[0].startswith("seqid\tsource\ttype\t")
    assert len(rows) > 1


def test_genbank_to_gff3(tmp_path):
    in_dir = _stage(tmp_path, "tiny.gbk")
    out_dir = tmp_path / "out"
    genbank_to_gff3.batch_convert(str(in_dir), str(out_dir))
    assert (out_dir / "tiny.gff3").exists()


def test_fasta_gff_to_gbk(tmp_path):
    in_dir = _stage(tmp_path, "tiny.fasta", "tiny.gff3")
    out_dir = tmp_path / "out"
    convert_all_gff_fasta_to_gbk.batch_convert(str(in_dir), str(out_dir))
    content = (out_dir / "tiny.gbk").read_text()
    assert "LOCUS" in content


def test_gff3_to_protein(tmp_path):
    in_dir = _stage(tmp_path, "tiny.fasta", "tiny.gff3")
    out_dir = tmp_path / "out"
    gff3_to_protein.batch_convert(str(in_dir), str(out_dir))
    content = (out_dir / "tiny.faa").read_text()
    # Exactly one protein per mRNA (the concatenation fix) rather than one per CDS.
    # tiny.gff3 has a single mRNA with two CDS children, so expect exactly one `>` line.
    assert content.count(">") == 1


# ---------- sequence family -------------------------------------------------


def test_fasta_to_fastq_and_back(tmp_path):
    in_dir = _stage(tmp_path, "tiny.fasta")
    fastq_out = tmp_path / "fastq"
    fasta_to_fastq.batch_convert(str(in_dir), str(fastq_out))
    fastq_file = fastq_out / "tiny.fastq"
    assert fastq_file.exists()
    # Round-trip back to FASTA and compare record counts.
    fasta_again = tmp_path / "fa2"
    fastq_to_fasta.batch_convert(str(fastq_out), str(fasta_again))
    orig_count = (FIXTURES / "tiny.fasta").read_text().count(">")
    new_count = (fasta_again / "tiny.fasta").read_text().count(">")
    assert orig_count == new_count


def test_fastq_to_fasta_qual_roundtrip(tmp_path):
    in_dir = _stage(tmp_path, "test.fastq")
    split_out = tmp_path / "split"
    fastq_to_fasta_qual.batch_convert(str(in_dir), str(split_out))
    assert (split_out / "test.fasta").exists()
    assert (split_out / "test.qual").exists()
    # Recombine and check the FASTQ comes back.
    # Stage the split outputs under a unified name so the pairing works.
    combine_in = tmp_path / "combine"
    combine_in.mkdir()
    shutil.copy(split_out / "test.fasta", combine_in / "test.fasta")
    shutil.copy(split_out / "test.qual", combine_in / "test.qual")
    combine_out = tmp_path / "fq_round"
    fasta_qual_to_fastq.batch_convert(str(combine_in), str(combine_out))
    assert (combine_out / "test.fastq").exists()


def test_fasta_to_table(tmp_path):
    in_dir = _stage(tmp_path, "tiny.fasta")
    out_dir = tmp_path / "out"
    fasta_to_table.batch_convert(str(in_dir), str(out_dir))
    content = (out_dir / "tiny.tsv").read_text()
    assert content.startswith("id\tsequence\n")


def test_convert_alignment(tmp_path):
    in_dir = _stage(tmp_path, "aln.fasta")
    out_dir = tmp_path / "out"
    convert_alignment.batch_convert(str(in_dir), str(out_dir),
                                    in_format="fasta", out_format="phylip")
    assert (out_dir / "aln.phylip").exists()


# ---------- alignment / mapping ---------------------------------------------


def test_bam_to_bed(tmp_path):
    in_dir = _stage(tmp_path, "test.bam", "test.bam.bai")
    out_dir = tmp_path / "out"
    bam_to_bed.batch_convert(str(in_dir), str(out_dir))
    rows = (out_dir / "test.bed").read_text().strip().splitlines()
    assert rows
    for r in rows:
        assert len(r.split("\t")) == 6


def test_blast_tab_to_links(tmp_path):
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    # Synthesise a minimal BLAST outfmt-6 file.
    (in_dir / "hits.tab").write_text(
        "q1\ts1\t95.0\t100\t5\t0\t1\t100\t201\t300\t1e-50\t180\n"
        "q2\ts2\t80.0\t50\t2\t0\t1\t50\t301\t350\t1e-20\t90\n"
    )
    out_dir = tmp_path / "out"
    blast_tab_to_links.batch_convert(str(in_dir), str(out_dir),
                                      min_length=0, min_identity=0)
    rows = (out_dir / "hits.links.tsv").read_text().strip().splitlines()
    assert rows[0].startswith("query_name\t")
    assert len(rows) == 3  # header + 2 data rows


def test_delta_to_tab(tmp_path):
    in_dir = _stage(tmp_path, "test.delta")
    out_dir = tmp_path / "out"
    delta_to_tab.batch_convert(str(in_dir), str(out_dir))
    assert (out_dir / "test.tsv").exists()


def test_maf_to_xmfa(tmp_path):
    in_dir = _stage(tmp_path, "test.maf")
    out_dir = tmp_path / "out"
    maf_to_xmfa.batch_convert(str(in_dir), str(out_dir))
    content = (out_dir / "test.xmfa").read_text()
    # XMFA must preserve gaps; the fixture contains gapped sequence.
    assert "=" in content  # block terminator


# ---------- VCF family ------------------------------------------------------


def test_vcf_to_bed(tmp_path):
    in_dir = _stage(tmp_path, "tiny.vcf")
    out_dir = tmp_path / "out"
    vcf_to_bed.batch_convert(str(in_dir), str(out_dir))
    rows = (out_dir / "tiny.bed").read_text().strip().splitlines()
    assert rows
    # Fixture has 3 records (2 SNPs + 1 indel); all should emit exactly one BED row each.
    assert len(rows) == 3


def test_vcf_to_table(tmp_path):
    in_dir = _stage(tmp_path, "tiny.vcf")
    out_dir = tmp_path / "out"
    vcf_to_table.batch_convert(str(in_dir), str(out_dir))
    content = (out_dir / "tiny.tsv").read_text()
    assert content.startswith("CHROM\tPOS\tID\tREF\tALT\t")


def test_vcf_to_consensus(tmp_path):
    in_dir = _stage(tmp_path, "tiny.fasta", "tiny.vcf")
    out_dir = tmp_path / "out"
    vcf_to_consensus.batch_convert(str(in_dir), str(out_dir))
    content = (out_dir / "tiny_consensus.fa").read_text()
    # Two samples in tiny.vcf, one contig → 2 records.
    assert content.count(">") == 2


def test_vcf_to_eigenstrat_golden(tmp_path):
    in_dir = _stage(tmp_path, "tiny.vcf")
    out_dir = tmp_path / "out"
    vcf_to_eigenstrat.batch_convert(str(in_dir), str(out_dir))
    geno = (out_dir / "tiny.geno").read_text()
    snp = (out_dir / "tiny.snp").read_text()
    ind = (out_dir / "tiny.ind").read_text()
    assert geno == "10\n01\n"
    assert snp.startswith("ctg1_5\tctg1\t0.0\t5\tA\tT\n")
    assert "sample1\tU\tsample1" in ind
    assert "sample2\tU\tsample2" in ind


def test_vcf_to_pseudohaploid_reproducible(tmp_path):
    in_dir = _stage(tmp_path, "tiny.vcf")
    out_a = tmp_path / "a"
    out_b = tmp_path / "b"
    vcf_to_pseudohaploid.batch_convert(str(in_dir), str(out_a), seed=42)
    vcf_to_pseudohaploid.batch_convert(str(in_dir), str(out_b), seed=42)
    assert (out_a / "tiny.geno").read_text() == (out_b / "tiny.geno").read_text()
    # No heterozygous call may survive pseudohaploid: only 0, 2, 9.
    chars = set((out_a / "tiny.geno").read_text().replace("\n", ""))
    assert chars.issubset({"0", "2", "9"})


def test_vcf_to_eigenstrat_also_plink(tmp_path):
    in_dir = _stage(tmp_path, "tiny.vcf")
    out_dir = tmp_path / "out"
    vcf_to_eigenstrat.batch_convert(str(in_dir), str(out_dir), also_plink=True)
    assert (out_dir / "tiny.ped").exists()
    assert (out_dir / "tiny.map").exists()
    ped = (out_dir / "tiny.ped").read_text().strip().splitlines()
    assert len(ped) == 2  # 2 samples
    # Each PED row has 6 metadata cols + 2 alleles per SNP (2 SNPs kept → 4).
    assert len(ped[0].split()) == 6 + 2 * 2


def test_vcf_to_eigenstrat_maf_filter(tmp_path):
    in_dir = _stage(tmp_path, "tiny.vcf")
    out_dir = tmp_path / "out"
    # With --min-maf 0.5 the fixture's MAF=0.25 sites should be filtered out.
    vcf_to_eigenstrat.batch_convert(str(in_dir), str(out_dir), min_maf=0.5)
    geno = (out_dir / "tiny.geno").read_text()
    assert geno == ""  # every site filtered


# ---------- phylogenetic trees ----------------------------------------------


def test_tree_convert_roundtrip(tmp_path):
    in_dir = _stage(tmp_path, "test.nwk")
    nexus_dir = tmp_path / "nexus"
    tree_convert.batch_convert(str(in_dir), str(nexus_dir),
                               in_format="newick", out_format="nexus")
    assert any(nexus_dir.iterdir())


def test_annotate_tree(tmp_path):
    tree_src = FIXTURES / "test_aln.nwk"
    aln_src = FIXTURES / "aln.fasta"
    out_file = tmp_path / "annotated.nex"
    annotate_tree.annotate_tree(str(tree_src), str(aln_src), str(out_file))
    content = out_file.read_text()
    assert "#NEXUS" in content
    assert "begin trees;" in content
    assert "begin characters;" in content


# ---------- 0.1.4 additions -------------------------------------------------


def test_gff3_to_bed12(tmp_path):
    in_dir = _stage(tmp_path, "tiny.gff3")
    out_dir = tmp_path / "out"
    gff3_to_bed12.batch_convert(str(in_dir), str(out_dir))
    rows = (out_dir / "tiny.bed12").read_text().strip().splitlines()
    assert rows
    # tiny.gff3's mRNA has 2 exons at 100-150 and 150-200 (0-based).
    fields = rows[0].split("\t")
    assert len(fields) == 12
    assert fields[0] == "ctg1"                 # chrom
    assert fields[1] == "100"                  # chromStart (0-based)
    assert fields[2] == "200"                  # chromEnd
    assert fields[9] == "2"                    # blockCount
    assert fields[10].rstrip(",") == "50,50"   # blockSizes
    assert fields[11].rstrip(",") == "0,50"    # blockStarts


def test_stockholm_roundtrip(tmp_path):
    # FASTA -> Stockholm -> FASTA, record count preserved.
    in_dir = _stage(tmp_path, "aln.fasta")
    sto_dir = tmp_path / "sto"
    fasta_to_stockholm.batch_convert(str(in_dir), str(sto_dir))
    assert (sto_dir / "aln.sto").exists()

    round_dir = tmp_path / "round"
    stockholm_to_fasta.batch_convert(str(sto_dir), str(round_dir))
    fasta_out = round_dir / "aln.fasta"
    assert fasta_out.exists()
    assert fasta_out.read_text().count(">") == (FIXTURES / "aln.fasta").read_text().count(">")


def test_gtf_to_gff3(tmp_path):
    # Synthesise a tiny GTF-only fixture (no explicit gene row — converter
    # should synthesise one by grouping exon rows by gene_id).
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    (in_dir / "tiny.gtf").write_text(
        'ctg1\ttest\texon\t101\t150\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
        'ctg1\ttest\texon\t151\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
        'ctg1\ttest\tCDS\t101\t150\t.\t+\t0\tgene_id "g1"; transcript_id "t1";\n'
    )
    out_dir = tmp_path / "out"
    gtf_to_gff3.batch_convert(str(in_dir), str(out_dir))
    content = (out_dir / "tiny.gff3").read_text()
    assert content.startswith("##gff-version 3")
    assert "\tgene\t" in content
    assert "\tmRNA\t" in content
    # Child rows must reference the synthesised parent transcript.
    assert "Parent=t1" in content


def test_hmmer_tblout_to_tsv(tmp_path):
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    # Three columns of real tblout header + a representative data row.
    sample = (
        "#                                                                 --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----\n"
        "# target name        accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target\n"
        "#------------------- ---------- -------------------- ---------- --------- ------ -----   --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------\n"
        "YAL001C               -          PF12831.10           -          1.2e-45  154.2   0.3   2.1e-45  153.4   0.3   1.1   1   0   0   1   1   1   1 Uncharacterised protein\n"
    )
    (in_dir / "scan.tblout").write_text(sample)
    out_dir = tmp_path / "out"
    hmmer_tblout_to_tsv.batch_convert(str(in_dir), str(out_dir),
                                       hmmer_format="tblout")
    rows = (out_dir / "scan.tsv").read_text().strip().splitlines()
    assert rows[0].startswith("target_name\t")
    assert "YAL001C" in rows[1]
    # Description column must be preserved whole.
    assert rows[1].endswith("Uncharacterised protein")


def test_orthogroups_to_fasta(tmp_path):
    # 2 species, 2 orthogroups, checking >species|gene headers.
    fasta_dir = tmp_path / "fastas"
    fasta_dir.mkdir()
    (fasta_dir / "spA.fasta").write_text(">a1\nMKT\n>a2\nMRK\n")
    (fasta_dir / "spB.fasta").write_text(">b1\nMTK\n>b2\nMQQ\n")
    tsv = tmp_path / "Orthogroups.tsv"
    tsv.write_text(
        "Orthogroup\tspA\tspB\n"
        "OG0000000\ta1, a2\tb1\n"
        "OG0000001\t\tb2\n"
    )
    out_dir = tmp_path / "out"
    orthogroups_to_fasta.convert(str(tsv), str(fasta_dir), str(out_dir))
    og0 = (out_dir / "OG0000000.fa").read_text()
    assert ">spA|a1" in og0
    assert ">spA|a2" in og0
    assert ">spB|b1" in og0
    og1 = (out_dir / "OG0000001.fa").read_text()
    assert ">spB|b2" in og1


def test_vcf_to_plink(tmp_path):
    in_dir = _stage(tmp_path, "tiny.vcf")
    out_dir = tmp_path / "out"
    vcf_to_plink.batch_convert(str(in_dir), str(out_dir))
    bed = (out_dir / "tiny.bed").read_bytes()
    # 3 magic bytes + 2 SNPs * ceil(2/4)=1 byte = 5 bytes total.
    assert bed[:3] == b"\x6c\x1b\x01"
    assert len(bed) == 5
    # Per-SNP bytes computed by hand (see CHANGELOG / docstring).
    # SNP 1 row = "10" -> sample0=het(10), sample1=hom-alt(00); last two samples padded missing(01 each).
    # Byte = 01_01_00_10 = 0x52
    # SNP 2 row = "01" -> sample0=hom-alt(00), sample1=het(10); pad 01 01.
    # Byte = 01_01_10_00 = 0x58
    assert bed[3:] == b"\x52\x58"

    bim = (out_dir / "tiny.bim").read_text().splitlines()
    assert len(bim) == 2
    # .bim columns: chrom snp_id cM bp a1(minor) a2(major)
    assert bim[0].split("\t")[0] == "ctg1"
    assert bim[0].split("\t")[4] == "T"  # alt = minor
    assert bim[0].split("\t")[5] == "A"  # ref = major

    fam = (out_dir / "tiny.fam").read_text().splitlines()
    assert len(fam) == 2
    assert fam[0].split()[1] == "sample1"
    assert fam[1].split()[1] == "sample2"
