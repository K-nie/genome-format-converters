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


def test_genbank_to_gff3_joined_cds(tmp_path):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    rec = SeqRecord(Seq("ATGAAACCCGGGTTTAAATAG" + "C" * 9), id="chrT", name="chrT",
                    annotations={"molecule_type": "DNA"})
    # join(1..7,11..21): a 7-base exon, an intron, then 11 bases.
    loc = CompoundLocation([FeatureLocation(0, 7, 1), FeatureLocation(10, 21, 1)])
    rec.features.append(SeqFeature(loc, type="CDS",
                                   qualifiers={"locus_tag": ["T1"], "pseudo": [""]}))
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    SeqIO.write(rec, str(in_dir / "t.gbk"), "genbank")
    out_dir = tmp_path / "out"
    genbank_to_gff3.batch_convert(str(in_dir), str(out_dir))
    rows = [l.split("\t") for l in (out_dir / "t.gff3").read_text().splitlines()
            if not l.startswith("#")]
    assert rows[0][2:5] == ["region", "1", "30"]
    cds = [r for r in rows if r[2] == "CDS"]
    # One line per part, sharing an ID; the second part starts 2 bases
    # into a codon (7 % 3 == 1), so its phase is 2.
    assert [(r[3], r[4], r[7]) for r in cds] == [("1", "7", "0"), ("11", "21", "2")]
    assert all(r[8].startswith("ID=CDS:T1;") for r in cds)
    assert "pseudo=true" in cds[0][8]


def test_gtf_to_gff3_cds_includes_stop_codon(tmp_path):
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    a = 'gene_id "G{g}"; transcript_id "T{g}";'
    rows = [
        ("transcript", 100, 1000, "+", ".", 1), ("exon", 100, 400, "+", ".", 1),
        ("exon", 600, 1000, "+", ".", 1), ("CDS", 200, 400, "+", "0", 1),
        ("CDS", 600, 797, "+", "0", 1), ("stop_codon", 798, 800, "+", "0", 1),
        # Minus strand, stop codon split across the intron: one base after
        # the last CDS row, two in the next exon.
        ("transcript", 2000, 3000, "-", ".", 2), ("exon", 2000, 2101, "-", ".", 2),
        ("exon", 2500, 3000, "-", ".", 2), ("CDS", 2501, 2899, "-", "0", 2),
        ("stop_codon", 2500, 2500, "-", "0", 2), ("stop_codon", 2100, 2101, "-", "2", 2),
    ]
    (in_dir / "x.gtf").write_text("".join(
        f"22\tens\t{t}\t{s}\t{e}\t.\t{st}\t{ph}\t{a.format(g=g)}\n"
        for t, s, e, st, ph, g in rows))
    out_dir = tmp_path / "out"
    gtf_to_gff3.batch_convert(str(in_dir), str(out_dir))
    cds = sorted((int(r[3]), int(r[4]), r[7]) for r in
                 (l.split("\t") for l in (out_dir / "x.gff3").read_text().splitlines())
                 if len(r) > 2 and r[2] == "CDS")
    # GFF3 CDS includes the stop codon. 399 + 1 coding bases come before
    # the split piece, so both of its bases finish that codon: phase 2.
    assert cds == [(200, 400, "0"), (600, 800, "0"), (2100, 2101, "2"), (2500, 2899, "0")]


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


def _pfam_style_fasta(path: Path, names) -> None:
    rows = ["AC.GT-AC", "ACAGT-A.", "AC.GTTAC"]
    path.write_text("".join(f">{n}\n{r}\n" for n, r in zip(names, rows)))


def test_convert_alignment_phylip_dots_and_relaxed(tmp_path):
    from Bio import AlignIO
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    names = ["Q9XYZ1_HUMAN/12-150", "Q9XYZ1_HUMAN/160-300", "P12345_MOUSE/1-90"]
    _pfam_style_fasta(in_dir / "fam.fasta", names)
    out_dir = tmp_path / "relaxed"
    convert_alignment.batch_convert(str(in_dir), str(out_dir),
                                    in_format="fasta", out_format="phylip-relaxed")
    aln = AlignIO.read(str(out_dir / "fam.phylip"), "phylip-relaxed")
    assert [r.id for r in aln] == names
    # PHYLIP has no '.' gap symbol.
    assert [str(r.seq) for r in aln] == ["AC-GT-AC", "ACAGT-A-", "AC-GTTAC"]


def test_convert_alignment_strict_phylip_refuses_name_clash(tmp_path):
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    _pfam_style_fasta(in_dir / "fam.fasta",
                      ["Q9XYZ1_HUMAN/12-150", "Q9XYZ1_HUMAN/160-300", "P12345_MOUSE/1-90"])
    out_dir = tmp_path / "strict"
    with pytest.raises(SystemExit):
        convert_alignment.batch_convert(str(in_dir), str(out_dir),
                                        in_format="fasta", out_format="phylip")
    assert not any(out_dir.glob("*.phylip*"))


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
    rows = [r.split("\t") for r in (out_dir / "test.tsv").read_text().splitlines()]
    assert rows[0] == ["query_name", "query_start", "query_end",
                       "ref_name", "ref_start", "ref_end", "strand"]
    # S1 E1 are reference coordinates, S2 E2 query; a high-to-low query
    # interval is a reverse-strand alignment.
    assert rows[1] == ["qry1", "1", "100", "ref1", "1", "100", "+"]
    assert rows[2] == ["qry1", "301", "400", "ref1", "201", "300", "-"]
    assert len(rows) == 3


def test_maf_to_xmfa(tmp_path):
    in_dir = _stage(tmp_path, "test.maf")
    out_dir = tmp_path / "out"
    maf_to_xmfa.batch_convert(str(in_dir), str(out_dir))
    lines = (out_dir / "test.xmfa").read_text().splitlines()
    assert lines[0] == "#FormatVersion Mauve1"
    assert "#Sequence3File\tseq3" in lines
    heads = [l for l in lines if l.startswith(">")]
    # Mauve headers: index, 1-based forward-strand interval, strand, name.
    # seq3 is on the minus strand at MAF start 2, size 4, srcSize 10, so
    # it covers forward positions 5-8.
    assert heads[:3] == ["> 1:1-4 + seq1", "> 2:1-4 + seq2", "> 3:5-8 - seq3"]
    assert heads[3:] == ["> 1:5-7 + seq1", "> 2:7-9 + seq2"]
    assert lines.count("=") == 2
    assert "AC-TG" in lines  # gaps kept


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


def test_tree_convert_nexus_quotes_and_precision(tmp_path):
    from Bio import Phylo
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    (in_dir / "t.nwk").write_text("((A/1-2:0.1,B_x:0.2)0.926:0.3,'C D':0.412970504);\n")
    nexus_dir = tmp_path / "nexus"
    tree_convert.batch_convert(str(in_dir), str(nexus_dir),
                               in_format="newick", out_format="nexus")
    text = (nexus_dir / "t.nexus").read_text()
    assert "TaxLabels 'A/1-2' B_x 'C D';" in text
    assert "0.926" in text and "0.412970504" in text
    back = tmp_path / "back"
    tree_convert.batch_convert(str(nexus_dir), str(back),
                               in_format="nexus", out_format="newick")
    tree = Phylo.read(str(back / "t.newick"), "newick")
    assert sorted(t.name for t in tree.get_terminals()) == ["A/1-2", "B_x", "C D"]


def test_annotate_tree(tmp_path):
    tree_src = FIXTURES / "test_aln.nwk"
    aln_src = FIXTURES / "aln.fasta"
    out_file = tmp_path / "annotated.nex"
    annotate_tree.annotate_tree(str(tree_src), str(aln_src), str(out_file))
    content = out_file.read_text()
    assert "#NEXUS" in content
    assert "begin trees;" in content
    assert "begin characters;" in content


def test_annotate_tree_quotes_pfam_names(tmp_path):
    from Bio.Nexus import Nexus
    names = ["Q9XYZ1_HUMAN/12-150", "Q9XYZ1_HUMAN/160-300", "P12345_MOUSE/1-90"]
    aln = tmp_path / "fam.fasta"
    aln.write_text("".join(f">{n}\n{r}\n" for n, r in
                           zip(names, ["MK.LV-W", "MKALV-W", "MK.LVQW"])))
    tree = tmp_path / "fam.nwk"
    tree.write_text("('Q9XYZ1_HUMAN/12-150':0.1,'Q9XYZ1_HUMAN/160-300':0.2,"
                    "'P12345_MOUSE/1-90':0.3);\n")
    out_file = tmp_path / "fam.nex"
    annotate_tree.annotate_tree(str(tree), str(aln), str(out_file))
    text = out_file.read_text()
    assert "datatype=protein" in text
    nex = Nexus.Nexus(str(out_file))
    assert sorted(nex.taxlabels) == sorted(names)
    assert str(nex.matrix[names[0]]) == "MK-LV-W"
    from Bio import Phylo
    # Bio.Phylo keeps the NEXUS quotes on translated names.
    leaves = sorted(t.name.strip("'") for t in
                    Phylo.read(str(out_file), "nexus").get_terminals())
    assert leaves == sorted(names)


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


# Ensembl and tRNAscan-SE hang exons on many feature types besides mRNA
# (pseudogenic_transcript, miRNA, V_gene_segment, unconfirmed_transcript,
# tRNA pseudogene rows). Each must come out as its own transcript.
_ENS_TX = {f"transcript:T{i}" for i in range(1, 7)}


def _ens_dir(tmp_path):
    d = tmp_path / "ens"
    d.mkdir()
    shutil.copy(FIXTURES / "ensembl_biotypes.gff3", d / "e.gff3")
    return d


def test_gff3_to_gtf_every_exon_parent_is_a_transcript(tmp_path):
    out = tmp_path / "gtf"
    gff3_to_gtf.batch_convert(str(_ens_dir(tmp_path)), str(out))
    exon_tx = set()
    for line in (out / "e.gtf").read_text().splitlines():
        c = line.split("\t")
        if c[2] == "exon":
            assert 'transcript_id "' in c[8], line
            exon_tx.add(c[8].split('transcript_id "')[1].split('"')[0])
    assert exon_tx == _ENS_TX


def test_gff3_to_bed12_every_exon_parent_is_a_transcript(tmp_path):
    out = tmp_path / "bed"
    gff3_to_bed12.batch_convert(str(_ens_dir(tmp_path)), str(out))
    names = [r.split("\t")[3] for r in (out / "e.bed12").read_text().splitlines()]
    assert sorted(names) == sorted(_ENS_TX)


def test_single_file_output_is_renamed_when_the_folder_has_other_files(tmp_path):
    import subprocess
    import sys as _sys
    src = Path(__file__).parent / "test_data" / "tiny.gff3"
    (tmp_path / "in.gff3").write_text(src.read_text())
    (tmp_path / "unrelated.txt").write_text("x\n")
    subprocess.run([_sys.executable, "-m", "genome_format_converters.cli", "-q",
                    "gff3-to-gtf", "--input", str(tmp_path / "in.gff3"),
                    "--output", str(tmp_path / "out.gtf")], check=True)
    assert (tmp_path / "out.gtf").exists()
    assert not (tmp_path / "in.gtf").exists()
