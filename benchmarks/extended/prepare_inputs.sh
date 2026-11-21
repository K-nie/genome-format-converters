#!/usr/bin/env bash
# Build the inputs of the extended benchmark (E01-E24) under
# benchmarks/data/extended/<task>/. Idempotent: a task whose directory
# already holds its inputs is skipped.
#
# Needs the T1-T8 data first (download_1kg_chr22.sh,
# download_gene_annotations.sh) and these programs on PATH:
#   art_illumina minimap2 samtools bcftools tabix bgzip diamond nucmer
#   lastdb lastal last-split FastTree esl-reformat plink2 python curl
# All downloads are public and pinned to a release:
#   Ensembl 112 (GTF, proteome, GRCh38 chr22), Ensembl 75 (GRCh37 chr22),
#   NCBI RefSeq GCF_000146045.2 (R64 GenBank), Pfam 37.0 (seed alignments).
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
bench="$(cd "$here/.." && pwd)"
data="${GFC_BENCH_DATA:-$bench/data}"
ext="$data/extended"
dl="$ext/_downloads"
mkdir -p "$ext" "$dl"
threads="${GFC_BENCH_PREP_THREADS:-4}"
seed=11

r64_fa="$data/annotations/sc_r64/Saccharomyces_cerevisiae.R64-1-1.112.fasta"
r64_gff="$data/annotations/sc_r64/Saccharomyces_cerevisiae.R64-1-1.112.gff3"
hs_gff="$data/annotations/hs_chr22/Homo_sapiens.GRCh38.112.chromosome.22.gff3"
vcf="$data/1kg/chr22.biallelic.snps.filtered.vcf.gz"
for f in "$r64_fa" "$r64_gff" "$hs_gff" "$vcf"; do
    [[ -s "$f" ]] || { echo "[err ] missing $f: run the T1-T8 downloads first" >&2; exit 1; }
done

fetch() {  # fetch <url> <dest>; gunzips when the URL ends in .gz
    local url="$1" dest="$2"
    [[ -s "$dest" ]] && return
    echo "[get ] $url" >&2
    if [[ "$url" == *.gz ]]; then
        curl -fSL --retry 3 "$url" | gunzip -c > "$dest.tmp"
    else
        curl -fSL --retry 3 -o "$dest.tmp" "$url"
    fi
    mv "$dest.tmp" "$dest"
}
done_flag() { [[ -f "$ext/$1/.done" ]]; }
mark() { touch "$ext/$1/.done"; echo "[ok  ] $1: $(ls "$ext/$1" | grep -v '^\.' | tr '\n' ' ')" >&2; }
fresh() { rm -rf "$ext/$1"; mkdir -p "$ext/$1"; }

ens=https://ftp.ensembl.org/pub
fetch "$ens/release-112/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.112.gtf.gz" "$dl/sc_r64.gtf"
fetch "$ens/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.chr.gtf.gz" "$dl/hs_all.gtf"
fetch "$ens/release-112/fasta/saccharomyces_cerevisiae/pep/Saccharomyces_cerevisiae.R64-1-1.pep.all.fa.gz" "$dl/sc_r64.pep.fa"
fetch "$ens/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz" "$dl/GRCh38.chr22.fa"
fetch "$ens/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.22.fa.gz" "$dl/GRCh37.chr22.fa"
fetch "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gbff.gz" "$dl/R64.gbk"
fetch "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.0/Pfam-A.seed.gz" "$dl/Pfam-A.seed"

# ---- reads (E01-E06) --------------------------------------------------
if ! done_flag E02; then
    fresh E02
    # ART HiSeq 2500 error and quality profile, single-end 150 bp, 6x
    # coverage of R64, fixed seed.
    art_illumina -ss HS25 -i "$r64_fa" -l 150 -f 6 -rs "$seed" -na -q -o "$dl/reads" >/dev/null
    mv "$dl/reads.fq" "$ext/E02/reads.fastq"
    mark E02
fi
if ! done_flag E03; then fresh E03; ln -s "$ext/E02/reads.fastq" "$ext/E03/reads.fastq"; mark E03; fi
if ! done_flag E04; then
    fresh E04
    python "$here/refs/fastq_to_fasta_qual.py" "$ext/E02/reads.fastq" "$ext/E04/reads"
    mark E04
fi
for t in E01 E05; do
    if ! done_flag $t; then fresh $t; ln -s "$ext/E04/reads.fasta" "$ext/$t/reads.fasta"; mark $t; fi
done
if ! done_flag E06; then
    fresh E06
    minimap2 -t "$threads" -ax sr "$r64_fa" "$ext/E02/reads.fastq" 2>/dev/null \
        | samtools sort -@ "$threads" -o "$ext/E06/reads.bam" -
    samtools index "$ext/E06/reads.bam"
    mark E06
fi

# ---- annotation (E07-E11) ---------------------------------------------
if ! done_flag E07; then
    fresh E07
    cp "$dl/sc_r64.gtf" "$ext/E07/sc_r64.gtf"
    awk -F'\t' '/^#/ || $1=="22"' "$dl/hs_all.gtf" > "$ext/E07/hs_chr22.gtf"
    mark E07
fi
for t in E08 E09; do
    if ! done_flag $t; then
        fresh $t
        ln -s "$r64_gff" "$ext/$t/sc_r64.gff3"
        ln -s "$hs_gff" "$ext/$t/hs_chr22.gff3"
        mark $t
    fi
done
if ! done_flag E10; then
    fresh E10
    ln -s "$r64_gff" "$ext/E10/sc_r64.gff3"
    ln -s "$r64_fa" "$ext/E10/sc_r64.fasta"
    mark E10
fi
if ! done_flag E11; then fresh E11; cp "$dl/R64.gbk" "$ext/E11/R64.gbk"; mark E11; fi

# ---- alignments and trees (E12-E16) -----------------------------------
if ! done_flag E12; then
    fresh E12
    # The first 500 seed families with at least 10 sequences, one
    # Stockholm file each (Pfam files are Latin-1).
    python - "$dl/Pfam-A.seed" "$ext/E12" <<'PY'
import os, sys, pathlib
src, out = sys.argv[1], pathlib.Path(sys.argv[2])
n, block, acc, nseq = 0, [], None, 0
with open(src, encoding="latin-1") as fh:
    for line in fh:
        block.append(line)
        if line.startswith("#=GF AC"):
            acc = line.split()[2].split(".")[0]
        elif line.startswith("#=GF ID") and acc is None:
            acc = line.split()[2]
        elif line.strip() and not line.startswith(("#", "//")):
            nseq += 1
        if line.startswith("//"):
            if acc and nseq >= int(os.environ.get("GFC_EXT_MIN_SEQS", "10")):
                (out / f"{acc}.sto").write_text("".join(block), encoding="latin-1")
                n += 1
            block, acc, nseq = [], None, 0
            if n == 500:
                break
print(f"{n} families", file=sys.stderr)
PY
    mark E12
fi
for t in E13 E14; do
    if ! done_flag $t; then
        fresh $t
        for f in "$ext/E12"/*.sto; do esl-reformat afa "$f" > "$ext/$t/$(basename "$f" .sto).fasta"; done
        mark $t
    fi
done
if ! done_flag E15; then
    fresh E15
    for f in $(ls "$ext/E13"/*.fasta | head -200); do
        FastTree -quiet -nopr "$f" > "$ext/E15/$(basename "$f" .fasta).nwk" 2>/dev/null
    done
    mark E15
fi
if ! done_flag E16; then
    fresh E16
    for f in "$ext/E15"/*.nwk; do
        s=$(basename "$f" .nwk)
        ln -s "$f" "$ext/E16/$s.nwk"
        ln -s "$ext/E13/$s.fasta" "$ext/E16/$s.fasta"
    done
    mark E16
fi

# ---- comparative genomics (E17-E19) ------------------------------------
# Both assemblies call the chromosome "22"; give them distinct names so the
# alignments say which side is which.
[[ -s "$dl/GRCh38_22.fa" ]] || sed '1s/^>.*/>GRCh38_22/' "$dl/GRCh38.chr22.fa" > "$dl/GRCh38_22.fa"
[[ -s "$dl/GRCh37_22.fa" ]] || sed '1s/^>.*/>GRCh37_22/' "$dl/GRCh37.chr22.fa" > "$dl/GRCh37_22.fa"
if ! done_flag E17; then
    fresh E17
    diamond makedb --in "$dl/sc_r64.pep.fa" -d "$dl/sc_r64" --quiet
    diamond blastp -q "$dl/sc_r64.pep.fa" -d "$dl/sc_r64" -o "$ext/E17/sc_r64_allvsall.tab" \
        --outfmt 6 --max-target-seqs 25 --evalue 1e-5 --threads "$threads" --quiet
    mark E17
fi
if ! done_flag E18; then
    fresh E18
    nucmer -t "$threads" -l 100 -c 1000 --prefix "$ext/E18/chr22_GRCh38_vs_GRCh37" \
        "$dl/GRCh38_22.fa" "$dl/GRCh37_22.fa"
    mark E18
fi
if ! done_flag E19; then
    fresh E19
    # A 5 Mb window of GRCh38 chr22 (25-30 Mb by default) against the
    # GRCh37 sequence around it (500 kb margin each side, since coordinates
    # shift between the assemblies). The whole chromosome with -m 100 took
    # too long on the cluster. The GRCh37 window is cut into 1 Mb queries so
    # that lastal -P actually runs in parallel (it threads over queries),
    # and tantan masks simple repeats (-R11 -c) so they do not flood the
    # seeds. Each sequence name records its window, so MAF coordinates map
    # back to the chromosome by adding the window start.
    w0="${GFC_EXT_E19_START:-25000001}" w1="${GFC_EXT_E19_END:-30000000}"
    pad="${GFC_EXT_E19_PAD:-500000}" chunk="${GFC_EXT_E19_CHUNK:-1000000}"
    python - "$dl/GRCh38_22.fa" "$dl/GRCh37_22.fa" "$dl/E19_target.fa" "$dl/E19_query.fa" \
        "$w0" "$w1" "$pad" "$chunk" <<'PY'
import sys
t_in, q_in, t_out, q_out = sys.argv[1:5]
w0, w1, pad, chunk = (int(x) for x in sys.argv[5:9])
def seq(path):
    with open(path) as fh:
        return "".join(l.strip() for l in fh if not l.startswith(">"))
def write(fh, name, s):
    fh.write(f">{name}\n")
    for i in range(0, len(s), 60):
        fh.write(s[i:i + 60] + "\n")
t = seq(t_in)
w1 = min(w1, len(t))
with open(t_out, "w") as fh:
    write(fh, f"GRCh38_22_{w0}_{w1}", t[w0 - 1:w1])
q = seq(q_in)
q0, q1 = max(1, w0 - pad), min(len(q), w1 + pad)
with open(q_out, "w") as fh:
    for a in range(q0, q1 + 1, chunk):
        b = min(a + chunk - 1, q1)
        write(fh, f"GRCh37_22_{a}_{b}", q[a - 1:b])
print(f"E19 target GRCh38 {w0}-{w1}, query GRCh37 {q0}-{q1}", file=sys.stderr)
PY
    lastdb -P "$threads" -uNEAR -R11 -c "$dl/E19_target.lastdb" "$dl/E19_target.fa"
    lastal -P "$threads" -E 1e-20 "$dl/E19_target.lastdb" "$dl/E19_query.fa" \
        | last-split -m 1e-5 > "$ext/E19/chr22_GRCh38_vs_GRCh37.maf"
    mark E19
fi

# ---- variants (E20-E24) --------------------------------------------------
for t in E20 E21; do
    if ! done_flag $t; then
        fresh $t
        ln -s "$vcf" "$ext/$t/chr22.vcf.gz"
        [[ -f "$vcf.tbi" ]] && ln -s "$vcf.tbi" "$ext/$t/chr22.vcf.gz.tbi"
        mark $t
    fi
done
if ! done_flag E22; then
    fresh E22
    first10=$(bcftools query -l "$vcf" | head -10 | paste -sd,)
    bcftools view -s "$first10" -Oz -o "$ext/E22/chr22.vcf.gz" "$vcf"
    tabix -p vcf "$ext/E22/chr22.vcf.gz"
    cp "$dl/GRCh37.chr22.fa" "$ext/E22/chr22.fasta"
    samtools faidx "$ext/E22/chr22.fasta"
    mark E22
fi
if ! done_flag E23; then
    fresh E23
    bcftools view "$vcf" | awk '/^#/ {print; next} n < 100000 {print; n++}' \
        | bgzip > "$ext/E23/chr22_100k.vcf.gz"
    tabix -p vcf "$ext/E23/chr22_100k.vcf.gz"
    mark E23
fi
if ! done_flag E24; then
    fresh E24
    mkdir -p "$ext/E24/packed" "$ext/E24/text"
    "${GFC_BENCH_PLINK2_EIG:-plink2}" --vcf "$ext/E23/chr22_100k.vcf.gz" --export eig --out "$ext/E24/packed/chr22_100k" \
        --threads "$threads" --memory 4000 >/dev/null
    rm -f "$ext/E24/packed"/*.log
    # Text EIGENSTRAT with the same content: unpack the 2-bit .geno rows.
    python "$here/refs/unpack_geno.py" "$ext/E24/packed/chr22_100k" "$ext/E24/text/chr22_100k"
    mark E24
fi
echo "[done] extended inputs in $ext" >&2
