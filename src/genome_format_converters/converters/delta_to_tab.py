#!/usr/bin/env python3
"""Convert MUMmer .delta to TSV (per-alignment coordinates).

Each alignment block in a .delta file starts with a sequence-pair header
(`>ref_name query_name ref_len query_len`), followed by coordinate lines
(`S1 E1 S2 E2 errors similarity_errors stop_codons`), each followed by a
run of indel positions ending in `0`. S1 E1 are reference coordinates and
S2 E2 query coordinates; a reverse-strand alignment has its query (or, for
promer, either) interval written high to low. We emit one TSV row per
alignment, with both intervals as 1-based start <= end and the relative
orientation in a `strand` column, as `show-coords` reports them.
"""

from pathlib import Path
from typing import Optional

from ._common import iter_input_files, log_info, log_warn, prepare_output_dir


def _convert_file(in_file: Path, out_file: Path) -> None:
    with open(in_file) as f_in, open(out_file, "w") as f_out:
        f_out.write("query_name\tquery_start\tquery_end\tref_name\tref_start\tref_end\tstrand\n")
        current_ref = current_query = None
        for line in f_in:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith(">"):
                # '>ref query reflen querylen': sequence-pair header
                parts = line[1:].split()
                if len(parts) >= 2:
                    current_ref, current_query = parts[0], parts[1]
                else:
                    current_ref = current_query = None
                continue
            if current_ref is None:
                # First non-comment line of the file is the REF / QUERY paths.
                continue
            # Coordinate lines are exactly 7 integers; indel lines are a
            # single integer. Non-numeric lines indicate corruption: warn
            # and skip.
            tokens = line.split()
            try:
                parts = [int(t) for t in tokens]
            except ValueError:
                log_warn(f"{in_file.name}: skipping non-numeric line: {line!r}")
                continue
            if len(parts) != 7:
                continue
            rstart, rend, qstart, qend = parts[0], parts[1], parts[2], parts[3]
            strand = "+" if (rstart <= rend) == (qstart <= qend) else "-"
            if qstart > qend:
                qstart, qend = qend, qstart
            if rstart > rend:
                rstart, rend = rend, rstart
            f_out.write(f"{current_query}\t{qstart}\t{qend}\t{current_ref}\t{rstart}\t{rend}\t{strand}\n")


def batch_convert(input_dir: str, output_dir: str,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)
    for delta in iter_input_files(in_path, [".delta"], pattern=pattern):
        out_file = out_path / (delta.stem + ".tsv")
        log_info(f"Converting {delta.name} -> {out_file.name}")
        _convert_file(delta, out_file)
