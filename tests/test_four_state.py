"""The four-state comparators: conventions -> structural, defects -> disagreement."""

import csv
import shutil
import sys
from pathlib import Path

import pytest

pytest.importorskip("cyvcf2")

BENCH = Path(__file__).resolve().parents[1] / "benchmarks"
sys.path.insert(0, str(BENCH))

import four_state as fs  # noqa: E402
import mutation_test  # noqa: E402
import check_correctness  # noqa: E402


@pytest.fixture(scope="module")
def cases(tmp_path_factory):
    return mutation_test.run(tmp_path_factory.mktemp("mut"))


def test_every_case_gets_expected_state(cases):
    bad = [(c.task, c.name, c.expected, c.got, c.reason) for c in cases if not c.ok]
    assert not bad


def test_every_task_has_mutations(cases):
    tasks = {c.task for c in cases if c.kind == "mutation"}
    assert tasks == {f"T{i}" for i in range(1, 9)}


def test_populate_writes_four_state_sidecar(tmp_path):
    work = tmp_path / "w"
    work.mkdir()
    vcf_dir = work / "in"
    vcf_dir.mkdir()
    mutation_test.make_vcf(vcf_dir / "toy.vcf")
    bench = tmp_path / "bench"
    mutation_test.gfc("vcf-to-eigenstrat", "--input-dir", str(vcf_dir),
                      "--output-dir", str(bench / "gfc_rep1"), "--force")
    shutil.copytree(bench / "gfc_rep1", bench / "convertf_rep1")
    mutation_test.eig_chrom_label(bench / "convertf_rep1")
    tsv = tmp_path / "T1.tsv"
    with tsv.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["task", "tool", "version", "replicate", "wall_s",
                    "peak_rss_mb", "exit_code", "correct", "notes"])
        for tool in ("gfc", "convertf"):
            w.writerow(["T1", tool, "x", "1", "1", "1", "0", "", ""])
    check_correctness.populate("T1", bench, tsv)
    rows = list(csv.DictReader(tsv.open(), delimiter="\t"))
    assert {r["tool"]: r["correct"] for r in rows} == {"gfc": "1", "convertf": "1"}
    side = list(csv.DictReader(tsv.with_suffix(".verdicts.tsv").open(), delimiter="\t"))
    gfc_row = next(r for r in side if r["tool"] == "gfc")
    assert gfc_row["state"] == fs.STRUCTURAL
    assert gfc_row["reason"] == "chrom-label"


def test_populate_sidecar_names_conventions_and_reference(tmp_path):
    test_populate_writes_four_state_sidecar(tmp_path)
    side = list(csv.DictReader((tmp_path / "T1.verdicts.tsv").open(), delimiter="\t"))
    by_tool = {r["tool"]: r for r in side}
    assert by_tool["gfc"]["conventions"] == "chrom-label"
    assert by_tool["convertf"]["state"] == "reference"


# ---- split multi-allelic sites: two biallelic records at one position ----

_DUP_VCF = """##fileformat=VCFv4.2
##contig=<ID=22>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc\td\te\tf
22\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1\t0|1\t0|0\t1|0
22\t200\t.\tG\tT\t.\tPASS\t.\tGT\t0|1\t0|0\t0|0\t1|1\t0|1\t0|0
22\t200\t.\tG\tA\t.\tPASS\t.\tGT\t0|0\t1|1\t0|1\t0|0\t0|0\t0|1
22\t300\t.\tC\tT\t.\tPASS\t.\tGT\t1|1\t0|1\t0|0\t0|0\t1|0\t0|1
"""


@pytest.fixture()
def dup_vcf(tmp_path):
    d = tmp_path / "in"
    d.mkdir()
    (d / "d.vcf").write_text(_DUP_VCF)
    return d


def test_t1_split_multiallelic_is_not_a_disagreement(dup_vcf, tmp_path):
    out = tmp_path / "eig"
    mutation_test.gfc("vcf-to-eigenstrat", "--input-dir", str(dup_vcf),
                      "--output-dir", str(out), "--force")
    ref = tmp_path / "ref"
    shutil.copytree(out, ref)
    mutation_test.eig_allele_swap(ref)
    v = fs.compare_t1(out, ref)
    assert v.state == fs.STRUCTURAL, v.reason
    assert "allele-order" in v.conventions


def test_t2_split_multiallelic_is_not_a_disagreement(dup_vcf, tmp_path):
    out = tmp_path / "pl"
    mutation_test.gfc("vcf-to-plink", "--input-dir", str(dup_vcf),
                      "--output-dir", str(out), "--force")
    ref = tmp_path / "ref"
    shutil.copytree(out, ref)
    mutation_test.plink_allele_swap(ref)
    v = fs.compare_t2(out, ref)
    assert v.state == fs.STRUCTURAL, v.reason


def test_t2_still_catches_a_flipped_genotype_at_a_split_site(dup_vcf, tmp_path):
    out = tmp_path / "pl"
    mutation_test.gfc("vcf-to-plink", "--input-dir", str(dup_vcf),
                      "--output-dir", str(out), "--force")
    ref = tmp_path / "ref"
    shutil.copytree(out, ref)
    # Flip sample a at the third record (22:200 G>A): 00 (hom A1) <-> 11.
    bed = next(ref.glob("*.bed"))
    raw = bytearray(bed.read_bytes())
    i = 3 + 2 * 2  # 6 samples -> 2 bytes per variant
    raw[i] ^= 0b11
    bed.write_bytes(bytes(raw))
    v = fs.compare_t2(out, ref)
    assert v.state == fs.DISAGREEMENT and "22:200" in v.reason, v.reason


def test_pseudohaploid_check_handles_split_sites(dup_vcf, tmp_path):
    out = tmp_path / "ph"
    mutation_test.gfc("vcf-to-pseudohaploid", "--input-dir", str(dup_vcf),
                      "--output-dir", str(out), "--seed", "1", "--force")
    v = fs.check_pseudohaploid(out, dup_vcf / "d.vcf")
    assert v.state == fs.STRUCTURAL, v.reason
    # A hom-ref call turned into the ALT allele must be caught.
    geno = next(out.glob("*.geno"))
    rows = geno.read_text().splitlines()
    rows[2] = ("0" if rows[2][0] == "2" else "2") + rows[2][1:]  # sample a, 22:200 G>A, hom-ref
    geno.write_text("\n".join(rows) + "\n")
    assert fs.check_pseudohaploid(out, dup_vcf / "d.vcf").state == fs.DISAGREEMENT


def test_t8_reproducibility_with_threads_over_symlinks(dup_vcf, tmp_path):
    import subprocess
    res = subprocess.run(
        [sys.executable, str(BENCH / "t8_reproducibility.py"), "--vcf",
         str(dup_vcf / "d.vcf"), "--work", str(tmp_path / "w"), "--copies", "4",
         "--threads", "4", "--out", str(tmp_path / "r.tsv")],
        capture_output=True, text=True, timeout=300)
    assert res.returncode == 0, res.stderr
    rows = list(csv.DictReader((tmp_path / "r.tsv").open(), delimiter="\t"))
    assert all(r["pass"] == "yes" for r in rows), rows


def test_t1_reference_that_drops_a_colocated_record(dup_vcf, tmp_path):
    """convertf keeps one record per position (seen at chr22:43236421)."""
    out = tmp_path / "eig"
    mutation_test.gfc("vcf-to-eigenstrat", "--input-dir", str(dup_vcf),
                      "--output-dir", str(out), "--force")
    ref = tmp_path / "ref"
    shutil.copytree(out, ref)
    for ext in ("snp", "geno"):
        f = next(ref.glob("*." + ext))
        rows = f.read_text().splitlines()
        del rows[2]  # the second record at 22:200
        f.write_text("\n".join(rows) + "\n")
    v = fs.compare_t1(out, ref)
    assert v.state == fs.STRUCTURAL, v.reason
    assert "reference-drops-colocated" in v.conventions
    # A candidate site with no reference record at that position is still wrong.
    for ext in ("snp", "geno"):
        f = next(ref.glob("*." + ext))
        rows = f.read_text().splitlines()
        del rows[0]
        f.write_text("\n".join(rows) + "\n")
    assert fs.compare_t1(out, ref).state == fs.DISAGREEMENT


def test_t1_convertf_default_genetic_position(dup_vcf, tmp_path):
    """convertf turns a 0 genetic position into pos x 1e-8 Morgans."""
    out = tmp_path / "eig"
    mutation_test.gfc("vcf-to-eigenstrat", "--input-dir", str(dup_vcf),
                      "--output-dir", str(out), "--force")
    ref = tmp_path / "ref"
    shutil.copytree(out, ref)
    snp = next(ref.glob("*.snp"))

    def set_gpos(fn):
        rows = [l.split() for l in snp.read_text().splitlines()]
        snp.write_text("".join("\t".join(r[:2] + [fn(r)] + r[3:]) + "\n" for r in rows))

    set_gpos(lambda r: f"{int(r[3]) * 1e-8:.6f}")  # convertf prints 6 decimals
    v = fs.compare_t1(out, ref)
    assert v.state == fs.STRUCTURAL and "genetic-position-default" in v.conventions, v
    set_gpos(lambda r: "0.5")
    assert fs.compare_t1(out, ref).state == fs.DISAGREEMENT


# ---- transposable elements (yeast R64): references that skip or truncate them

_TE_GFF = (
    "##gff-version 3\n"
    "I\tsgd\tgene\t100\t400\t.\t+\t.\tID=gene:YAL001C\n"
    "I\tsgd\tmRNA\t100\t400\t.\t+\t.\tID=transcript:YAL001C_mRNA;Parent=gene:YAL001C\n"
    "I\tsgd\texon\t100\t400\t.\t+\t.\tParent=transcript:YAL001C_mRNA\n"
    "I\tsgd\ttransposable_element_gene\t1000\t2000\t.\t-\t.\tID=gene:YAR009C\n"
    "I\tsgd\ttransposable_element\t1000\t2000\t.\t-\t.\tID=transcript:YAR009C;Parent=gene:YAR009C\n"
    "I\tsgd\texon\t1000\t2000\t.\t-\t.\tParent=transcript:YAR009C\n"
)


def _gtf(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(
        f'I\tsgd\texon\t{s}\t{e}\t.\t{st}\t.\tgene_id "{g}"; transcript_id "{t}";\n'
        for t, g, s, e, st in rows))


def test_t4_reference_that_drops_te_transcripts(tmp_path):
    (tmp_path / "in").mkdir()
    (tmp_path / "in" / "te.gff3").write_text(_TE_GFF)
    mrna = ("transcript:YAL001C_mRNA", "gene:YAL001C", 100, 400, "+")
    te = ("transcript:YAR009C", "gene:YAR009C", 1000, 2000, "-")
    _gtf(tmp_path / "c" / "te.gtf", [mrna, te])
    _gtf(tmp_path / "r" / "te.gtf", [mrna])
    v = fs.compare_t4(tmp_path / "c", tmp_path / "r", input_dir=tmp_path / "in")
    assert v.state == fs.STRUCTURAL and "reference-drops-transposons" in v.conventions
    # Without the input, or when the extra transcript is not a TE, it is a defect.
    assert fs.compare_t4(tmp_path / "c", tmp_path / "r").state == fs.DISAGREEMENT
    _gtf(tmp_path / "c" / "te.gtf", [mrna, te, ("transcript:X", "gene:X", 5, 9, "+")])
    assert fs.compare_t4(tmp_path / "c", tmp_path / "r",
                         input_dir=tmp_path / "in").state == fs.DISAGREEMENT


def test_t3_reference_with_truncated_keys_and_gfc_insdc_keys(tmp_path):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    from genome_format_converters.converters import convert_all_gff_fasta_to_gbk as gbk

    ind = tmp_path / "in"
    ind.mkdir()
    (ind / "te.gff3").write_text(_TE_GFF)
    (ind / "te.fasta").write_text(">I\n" + "A" * 2100 + "\n")
    gbk.batch_convert(str(ind), str(tmp_path / "c"))
    out = next((tmp_path / "c").iterdir())
    text = out.read_text()
    assert "transposable_ele" not in text.split("ORIGIN")[0].replace(
        "GFF3 type: transposable_element", "")
    types = [f.type for f in SeqIO.read(str(out), "genbank").features]
    assert types.count("mobile_element") == 1 and all(len(t) <= 15 for t in types)

    # A Biopython-written reference keeps only top-level features and cuts
    # the TE gene's key to "transposable_ele".
    rec = SeqRecord(Seq("A" * 2100), id="I", annotations={"molecule_type": "DNA"})
    rec.features = [SeqFeature(FeatureLocation(99, 400, 1), type="gene"),
                    SeqFeature(FeatureLocation(999, 2000, -1),
                               type="transposable_element_gene")]
    (tmp_path / "r").mkdir()
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        SeqIO.write(rec, str(tmp_path / "r" / out.name), "genbank")
    v = fs.compare_t3(tmp_path / "c", tmp_path / "r")
    assert v.state == fs.STRUCTURAL, v.reason
    assert "reference-truncates-keys" in v.conventions


def test_t5_reference_that_drops_te_transcripts(tmp_path):
    (tmp_path / "in").mkdir()
    (tmp_path / "in" / "te.gff3").write_text(_TE_GFF)
    mrna = "I\t99\t400\ttranscript:YAL001C_mRNA\t0\t+\t99\t400\t0\t1\t301,\t0,\n"
    te = "I\t999\t2000\ttranscript:YAR009C\t0\t-\t999\t999\t0\t1\t1001,\t0,\n"
    for d, rows in (("c", mrna + te), ("r", mrna)):
        (tmp_path / d).mkdir()
        (tmp_path / d / "te.bed").write_text(rows)
    v = fs.compare_t5(tmp_path / "c", tmp_path / "r", input_dir=tmp_path / "in")
    assert v.state == fs.STRUCTURAL and "reference-drops-transposons" in v.conventions
    assert "row-order" not in v.conventions
    assert fs.compare_t5(tmp_path / "c", tmp_path / "r").state == fs.DISAGREEMENT


def test_t3_insdc_writer_files_so_types_as_misc_feature(tmp_path):
    import warnings
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord

    def write(d, feats):
        rec = SeqRecord(Seq("A" * 500), id="I", annotations={"molecule_type": "DNA"})
        rec.features = [SeqFeature(FeatureLocation(s, e, st), type=t) for t, s, e, st in feats]
        (tmp_path / d).mkdir()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            SeqIO.write(rec, str(tmp_path / d / "x.gb"), "genbank")

    write("r", [("gene", 10, 50, 1), ("ncRNA_gene", 60, 90, 1), ("pseudogene", 100, 200, -1)])
    write("c", [("gene", 10, 50, 1), ("misc_feature", 60, 90, 1), ("misc_feature", 100, 200, -1)])
    v = fs.compare_t3(tmp_path / "c", tmp_path / "r")
    assert v.state == fs.STRUCTURAL and "insdc-feature-keys" in v.conventions
    # A misc_feature in the wrong place is still a defect.
    write("c2", [("gene", 10, 50, 1), ("misc_feature", 61, 90, 1), ("misc_feature", 100, 200, -1)])
    assert fs.compare_t3(tmp_path / "c2", tmp_path / "r").state == fs.DISAGREEMENT


def _write_packed_geno(path, rows):
    """Packed EIGENSTRAT as plink2 --export eig writes it."""
    nind = len(rows[0])
    rlen = max(48, -(-nind // 4))
    head = f"GENO {nind} {len(rows)} 0 0".encode().ljust(rlen, b"\0")
    code = {"0": 0, "1": 1, "2": 2, "9": 3}
    body = b""
    for r in rows:
        vals = [code[x] for x in r] + [0] * (-len(r) % 4)
        rec = bytes((vals[k] << 6) | (vals[k + 1] << 4) | (vals[k + 2] << 2) | vals[k + 3]
                    for k in range(0, len(vals), 4))
        body += rec.ljust(rlen, b"\0")
    path.write_bytes(head + body)


def test_t1_packed_geno_from_plink2(tmp_path):
    rows = ["2109012", "0129221", "9999990"]
    snp = "".join(f"rs{i}\t22\t0.0\t{1000 + i}\tC\tG\n" for i in range(3))
    ind = "".join(f"S{i}\tU\tPOP\n" for i in range(7))
    for d in ("c", "r"):
        (tmp_path / d).mkdir()
        (tmp_path / d / "x.snp").write_text(snp)
        (tmp_path / d / "x.ind").write_text(ind)
    (tmp_path / "r" / "x.geno").write_text("".join(r + "\n" for r in rows))
    _write_packed_geno(tmp_path / "c" / "x.geno", rows)
    v = fs.compare_t1(tmp_path / "c", tmp_path / "r")
    assert v.state == fs.STRUCTURAL and "packed-geno" in v.conventions
    _write_packed_geno(tmp_path / "c" / "x.geno", ["2109012", "0129222", "9999990"])
    assert fs.compare_t1(tmp_path / "c", tmp_path / "r").state == fs.DISAGREEMENT
