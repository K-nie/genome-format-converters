#!/usr/bin/env python3
"""
Convert BLAST tabular output (outfmt 6) to a simplified link TSV.
"""

from pathlib import Path

def _convert_file(in_file: Path, out_file: Path, min_len: int, min_id: float) -> None:
    """Convert one BLAST tabular file to link TSV, applying filters."""
    with open(in_file) as f_in, open(out_file, 'w') as f_out:
        f_out.write("query_name\tquery_start\tquery_end\tref_name\tref_start\tref_end\tidentity\n")
        for line in f_in:
            parts = line.strip().split("\t")
            if len(parts) < 12:
                continue
            qname = parts[0]
            rname = parts[1]
            ident = float(parts[2])
            qlen = int(parts[3])
            qstart = int(parts[6])
            qend   = int(parts[7])
            rstart = int(parts[8])
            rend   = int(parts[9])
            # Filter by length and identity
            if (qend - qstart + 1) < min_len or ident < min_id:
                continue
            f_out.write(f"{qname}\t{qstart}\t{qend}\t{rname}\t{rstart}\t{rend}\t{ident}\n")

def batch_convert(input_dir: str, output_dir: str, min_length: int = 0, min_identity: float = 0) -> None:
    """Convert all BLAST tabular files in input_dir to link TSV files."""
    in_path = Path(input_dir)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    for tab in in_path.glob("*.tab*"):
        out_file = out_path / tab.with_suffix(".links.tsv").name
        print(f"Converting {tab.name} -> {out_file.name}")
        _convert_file(tab, out_file, min_length, min_identity)