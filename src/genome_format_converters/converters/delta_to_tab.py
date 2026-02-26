#!/usr/bin/env python3
"""
Convert MUMmer .delta file to tabular format (TSV).
"""

from pathlib import Path

def _convert_file(in_file: Path, out_file: Path) -> None:
    """Convert a single .delta file to TSV."""
    with open(in_file) as f_in, open(out_file, 'w') as f_out:
        f_out.write("query_name\tquery_start\tquery_end\tref_name\tref_start\tref_end\n")
        current_ref = None
        current_query = None
        for line in f_in:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if not current_ref:
                parts = line.split()
                if len(parts) >= 2:
                    current_ref, current_query = parts[0], parts[1]
                continue
            if line.startswith('>'):
                continue
            try:
                parts = list(map(int, line.split()))
                if len(parts) >= 7:
                    qstart, qend, rstart, rend = parts[0], parts[1], parts[2], parts[3]
                    # Ensure forward orientation
                    if qstart > qend:
                        qstart, qend = qend, qstart
                    if rstart > rend:
                        rstart, rend = rend, rstart
                    f_out.write(f"{current_query}\t{qstart}\t{qend}\t{current_ref}\t{rstart}\t{rend}\n")
            except ValueError:
                # Reset for next alignment
                current_ref = current_query = None

def batch_convert(input_dir: str, output_dir: str) -> None:
    """Convert all .delta files in input_dir to TSV in output_dir."""
    in_path = Path(input_dir)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    for delta in in_path.glob("*.delta"):
        out_file = out_path / delta.with_suffix(".tsv").name
        print(f"Converting {delta.name} -> {out_file.name}")
        _convert_file(delta, out_file)