#!/usr/bin/env python3
"""Convert BLAST tabular output (outfmt 6) to a simplified link TSV.

BLAST outfmt 6 has 12 columns: qseqid, sseqid, pident, length, mismatch,
gapopen, qstart, qend, sstart, send, evalue, bitscore.
"""

from pathlib import Path
from typing import Optional

from ._common import iter_input_files, log_info, log_warn, prepare_output_dir


def _convert_file(in_file: Path, out_file: Path, min_len: int, min_id: float) -> None:
    malformed = 0
    kept = 0
    with open(in_file) as f_in, open(out_file, "w") as f_out:
        f_out.write("query_name\tquery_start\tquery_end\tref_name\tref_start\tref_end\tidentity\n")
        for lineno, raw in enumerate(f_in, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 12:
                malformed += 1
                continue
            try:
                qname, rname = parts[0], parts[1]
                ident = float(parts[2])
                qstart = int(parts[6])
                qend = int(parts[7])
                rstart = int(parts[8])
                rend = int(parts[9])
            except ValueError:
                malformed += 1
                continue
            if (qend - qstart + 1) < min_len or ident < min_id:
                continue
            f_out.write(f"{qname}\t{qstart}\t{qend}\t{rname}\t{rstart}\t{rend}\t{ident}\n")
            kept += 1
    if malformed:
        log_warn(f"{in_file.name}: skipped {malformed} malformed or short rows.")
    log_info(f"{in_file.name}: wrote {kept} link rows.")


def batch_convert(input_dir: str, output_dir: str,
                  min_length: int = 0,
                  min_identity: float = 0.0,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)
    # Historically the package accepted `.tab` and any `.tab*` extension;
    # keep that behaviour, and also accept `.tsv` since many pipelines emit
    # BLAST tabular with that extension.
    exts = [".tab", ".tsv", ".tabular"]
    for tab in iter_input_files(in_path, exts, pattern=pattern):
        out_file = out_path / (tab.stem + ".links.tsv")
        log_info(f"Converting {tab.name} -> {out_file.name}")
        _convert_file(tab, out_file, min_length, min_identity)
