#!/usr/bin/env bash
# Pfam-A.hmm + tblout provisioning for the gfc benchmark suite (T6 input).
#
# Stage 2 path: pull the latest Pfam-A.hmm.gz from EBI, run hmmsearch on a
# small reference proteome (S. cerevisiae S288C, ~6,000 proteins) to
# generate a real ~10^4-row tblout, deposit it as
# benchmarks/data/hmmer/sample.tblout. ~1.4 GB Pfam download + ~5-30 min
# single-core hmmsearch.
#
# Stage 3 path: replace the S288C proteome with the 20-species Y1000+
# proteome slice for a ~10^5-row tblout. Adds ~200 MB of proteome FASTAs
# (still TODO; see benchmarks/data/download_y1000plus.sh Stage 2).
#
# Why not Stage 1: the Pfam-A.hmm.gz download is large enough that running
# it as part of every benchmark bootstrap is wasteful. Run this script
# once, manually, on the cluster submit host before the next paper-grade
# rerun.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
dest="$here/hmmer"
mkdir -p "$dest"

tblout="$dest/sample.tblout"
if [[ -s "$tblout" && "$(wc -l < "$tblout")" -gt 100 ]]; then
    echo "[skip] $tblout already populated ($(wc -l < "$tblout") lines)"
    exit 0
fi

pfam_url="https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
proteome_url="https://sgd-archive.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans_all.fasta.gz"

pfam_hmm="$dest/Pfam-A.hmm"
proteome_fa="$dest/S288C_orf_trans_all.fasta"

if [[ ! -s "$pfam_hmm" ]]; then
    echo "[get ] $pfam_url" >&2
    curl -fSL --retry 3 -o "${pfam_hmm}.gz" "$pfam_url"
    gunzip -f "${pfam_hmm}.gz"
    if command -v hmmpress >/dev/null 2>&1; then
        hmmpress -f "$pfam_hmm" >/dev/null
    fi
fi

if [[ ! -s "$proteome_fa" ]]; then
    echo "[get ] $proteome_url" >&2
    curl -fSL --retry 3 -o "${proteome_fa}.gz" "$proteome_url"
    gunzip -f "${proteome_fa}.gz"
fi

if command -v hmmsearch >/dev/null 2>&1; then
    threads="${GFC_BENCH_HMMER_THREADS:-1}"
    echo "[run ] hmmsearch --cpu $threads --tblout $tblout Pfam-A.hmm S288C_orf_trans_all.fasta" >&2
    hmmsearch --cpu "$threads" --tblout "$tblout" "$pfam_hmm" "$proteome_fa" > /dev/null
    echo "[done] $tblout populated ($(wc -l < "$tblout") lines)"
else
    echo "[fail] hmmsearch not on PATH; install HMMER (bioconda: hmmer) and rerun" >&2
    exit 1
fi
