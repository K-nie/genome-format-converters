#!/usr/bin/env bash
# Download / regenerate a HMMER --tblout fixture (Pfam-A × Y1000+ proteomes).
# Target: benchmarks/data/hmmer/pfam_vs_y1000plus.tblout
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
dest="$here/hmmer"
mkdir -p "$dest"

# TODO: pin a Zenodo mirror of the lab's pre-computed tblout so reviewers
# don't have to re-run hmmsearch. Size should be ~2 GB.
echo "[TODO] populate $dest/pfam_vs_y1000plus.tblout" >&2
echo "       Either deposit the lab's existing scan on Zenodo and curl it" >&2
echo "       here, or regenerate via:"                                    >&2
echo "         hmmsearch --tblout pfam_vs_y1000plus.tblout Pfam-A.hmm proteomes.fa" >&2
