#!/usr/bin/env python3
"""Convert HMMER ``--tblout`` / ``--domtblout`` files to structured TSV.

HMMER writes these as whitespace-aligned columns with an arbitrary-width
free-text ``description_of_target`` trailing column. The file also has
``#``-prefixed banner / footer lines that should be stripped.

Use ``--hmmer-format tblout`` (default) for per-sequence output and
``--hmmer-format domtblout`` for per-domain output.
"""

from pathlib import Path
from typing import List, Optional

from ._common import iter_input_files, log_info, log_warn, prepare_output_dir

# Column schema per the HMMER 3 user's guide.
_TBLOUT_COLUMNS: List[str] = [
    "target_name", "target_accession",
    "query_name", "query_accession",
    "full_evalue", "full_score", "full_bias",
    "best_domain_evalue", "best_domain_score", "best_domain_bias",
    "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc",
    "description_of_target",
]

_DOMTBLOUT_COLUMNS: List[str] = [
    "target_name", "target_accession", "tlen",
    "query_name", "query_accession", "qlen",
    "full_evalue", "full_score", "full_bias",
    "domain_num", "domain_of",
    "c_evalue", "i_evalue", "domain_score", "domain_bias",
    "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to",
    "acc", "description_of_target",
]


def _schema(fmt: str) -> List[str]:
    if fmt == "tblout":
        return _TBLOUT_COLUMNS
    if fmt == "domtblout":
        return _DOMTBLOUT_COLUMNS
    raise ValueError(f"Unknown HMMER format: {fmt!r}")


def _split_row(line: str, n_fixed: int) -> List[str]:
    # Fixed-width columns are whitespace-separated, and the final column is
    # free text that may itself contain whitespace. `split(maxsplit=n_fixed)`
    # pins the last column correctly.
    return line.split(None, n_fixed)


def _convert_file(in_file: Path, out_file: Path, fmt: str) -> None:
    columns = _schema(fmt)
    n_fixed = len(columns) - 1
    malformed = 0
    rows = 0
    with open(in_file) as f_in, open(out_file, "w") as f_out:
        f_out.write("\t".join(columns) + "\n")
        for raw in f_in:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = _split_row(line, n_fixed)
            if len(parts) < n_fixed:
                malformed += 1
                continue
            # Pad a missing description with empty string.
            while len(parts) < len(columns):
                parts.append("")
            f_out.write("\t".join(parts) + "\n")
            rows += 1
    if malformed:
        log_warn(f"{in_file.name}: skipped {malformed} malformed rows.")
    log_info(f"{in_file.name}: wrote {rows} {fmt} rows.")


def batch_convert(input_dir: str, output_dir: str,
                  hmmer_format: str = "tblout",
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)
    # Common extensions people use when saving these files.
    exts = [".tblout", ".domtblout", ".hmmer", ".hmmsearch", ".tbl", ".txt"]
    for in_file in iter_input_files(in_path, exts, pattern=pattern):
        out_file = out_path / (in_file.stem + ".tsv")
        log_info(f"Converting {in_file.name} -> {out_file.name}")
        _convert_file(in_file, out_file, hmmer_format)
