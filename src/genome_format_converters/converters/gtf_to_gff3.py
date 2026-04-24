#!/usr/bin/env python3
"""GTF → GFF3.

GTF (Ensembl-style) is the common output of RNA-seq pipelines and the
annotation format for Ensembl Fungi. It doesn't always emit explicit
``gene`` / ``transcript`` rows, so this converter synthesises them by
aggregating ``exon`` / ``CDS`` rows by ``gene_id`` / ``transcript_id``.
GFF3 output uses proper ``ID`` / ``Parent`` attributes so downstream tools
(IGV, apollo, bcbio-gff, etc.) can reconstruct the hierarchy.
"""

import re
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

from ._common import iter_input_files, log_info, log_warn, prepare_output_dir

# GTF row types that map 1:1 to GFF3 types. Anything not here passes through
# unchanged and under a parent ID derived from transcript_id if present.
_TYPE_MAP = {
    "gene": "gene",
    "transcript": "mRNA",
    "mRNA": "mRNA",
    "exon": "exon",
    "CDS": "CDS",
    "five_prime_utr": "five_prime_UTR",
    "3UTR": "three_prime_UTR",
    "5UTR": "five_prime_UTR",
    "three_prime_utr": "three_prime_UTR",
    "start_codon": "start_codon",
    "stop_codon": "stop_codon",
    "Selenocysteine": "Selenocysteine",
}

_ATTR_RX = re.compile(r'(\w+)\s+"([^"]*)"')


def _parse_gtf_attrs(raw: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for key, val in _ATTR_RX.findall(raw or ""):
        out.setdefault(key, val)
    return out


def _iter_gtf_rows(in_file: Path) -> Iterator[Tuple[str, str, str, int, int, str, str, str, Dict[str, str]]]:
    with open(in_file) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            chrom, source, ftype, start, end, score, strand, phase, attrs_raw = parts[:9]
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue
            attrs = _parse_gtf_attrs(attrs_raw)
            yield chrom, source, ftype, start_i, end_i, score, strand, phase, attrs


def _convert_file(in_file: Path, out_file: Path) -> None:
    # Pass 1: aggregate per-gene and per-transcript spans.
    gene_info: Dict[str, Dict] = {}
    tx_info: Dict[str, Dict] = {}
    child_rows: List[Tuple] = []  # non-gene/transcript rows

    for (chrom, source, ftype, start_i, end_i,
         score, strand, phase, attrs) in _iter_gtf_rows(in_file):
        gene_id = attrs.get("gene_id")
        tx_id = attrs.get("transcript_id")

        if ftype in {"gene"}:
            if gene_id:
                g = gene_info.setdefault(
                    gene_id,
                    {"chrom": chrom, "start": start_i, "end": end_i,
                     "strand": strand, "source": source,
                     "name": attrs.get("gene_name", gene_id)},
                )
                g["start"] = min(g["start"], start_i)
                g["end"] = max(g["end"], end_i)
            continue

        if ftype in {"transcript", "mRNA"}:
            if tx_id:
                t = tx_info.setdefault(
                    tx_id,
                    {"gene_id": gene_id or "", "chrom": chrom,
                     "start": start_i, "end": end_i,
                     "strand": strand, "source": source,
                     "name": attrs.get("transcript_name", tx_id)},
                )
                t["start"] = min(t["start"], start_i)
                t["end"] = max(t["end"], end_i)
            # Ensure the gene entry exists even if the GTF lacked a gene row.
            if gene_id:
                g = gene_info.setdefault(
                    gene_id,
                    {"chrom": chrom, "start": start_i, "end": end_i,
                     "strand": strand, "source": source,
                     "name": attrs.get("gene_name", gene_id)},
                )
                g["start"] = min(g["start"], start_i)
                g["end"] = max(g["end"], end_i)
            continue

        # Child feature (exon/CDS/etc.). Back-fill gene + transcript spans
        # from its coordinates so a GTF without explicit gene/transcript
        # rows still emits a valid GFF3 hierarchy.
        if gene_id:
            g = gene_info.setdefault(
                gene_id,
                {"chrom": chrom, "start": start_i, "end": end_i,
                 "strand": strand, "source": source,
                 "name": attrs.get("gene_name", gene_id)},
            )
            g["start"] = min(g["start"], start_i)
            g["end"] = max(g["end"], end_i)

        if tx_id:
            t = tx_info.setdefault(
                tx_id,
                {"gene_id": gene_id or "", "chrom": chrom,
                 "start": start_i, "end": end_i,
                 "strand": strand, "source": source,
                 "name": attrs.get("transcript_name", tx_id)},
            )
            t["start"] = min(t["start"], start_i)
            t["end"] = max(t["end"], end_i)

        child_rows.append(
            (chrom, source, ftype, start_i, end_i, score, strand, phase, attrs)
        )

    # Pass 2: write GFF3.
    child_id_counter: Dict[str, int] = {}  # per-type unique ID generator

    def _next_child_id(base: str, ftype: str) -> str:
        key = (base, ftype)
        n = child_id_counter.get(key, 0) + 1
        child_id_counter[key] = n
        return f"{base}.{ftype}.{n}"

    with open(out_file, "w") as fh:
        fh.write("##gff-version 3\n")
        for gene_id, g in gene_info.items():
            attrs_out = f"ID={gene_id}"
            if g.get("name") and g["name"] != gene_id:
                attrs_out += f";Name={g['name']}"
            fh.write(
                f"{g['chrom']}\t{g['source']}\tgene\t{g['start']}\t{g['end']}\t"
                f".\t{g['strand']}\t.\t{attrs_out}\n"
            )
            for tx_id, t in tx_info.items():
                if t["gene_id"] != gene_id:
                    continue
                tx_attrs = f"ID={tx_id};Parent={gene_id}"
                if t.get("name") and t["name"] != tx_id:
                    tx_attrs += f";Name={t['name']}"
                fh.write(
                    f"{t['chrom']}\t{t['source']}\tmRNA\t{t['start']}\t{t['end']}\t"
                    f".\t{t['strand']}\t.\t{tx_attrs}\n"
                )

        # Transcripts that claim no gene_id — emit as orphan mRNAs.
        orphan_tx = [tx_id for tx_id, t in tx_info.items() if not t["gene_id"]]
        for tx_id in orphan_tx:
            t = tx_info[tx_id]
            fh.write(
                f"{t['chrom']}\t{t['source']}\tmRNA\t{t['start']}\t{t['end']}\t"
                f".\t{t['strand']}\t.\tID={tx_id}\n"
            )

        for (chrom, source, ftype, start_i, end_i,
             score, strand, phase, attrs) in child_rows:
            gff_type = _TYPE_MAP.get(ftype, ftype)
            tx_id = attrs.get("transcript_id")
            parent = tx_id or attrs.get("gene_id", "")
            base = tx_id or attrs.get("gene_id", "feature")
            child_id = _next_child_id(base, gff_type)
            attrs_out = f"ID={child_id}"
            if parent:
                attrs_out += f";Parent={parent}"
            fh.write(
                f"{chrom}\t{source}\t{gff_type}\t{start_i}\t{end_i}\t"
                f"{score or '.'}\t{strand}\t{phase or '.'}\t{attrs_out}\n"
            )

    log_info(
        f"{in_file.name}: wrote {len(gene_info)} gene(s), {len(tx_info)} transcript(s), "
        f"{len(child_rows)} child row(s)."
    )


def batch_convert(input_dir: str, output_dir: str,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)
    for gtf in iter_input_files(in_path, [".gtf", ".gff2"], pattern=pattern):
        out_file = out_path / (gtf.stem + ".gff3")
        log_info(f"Converting {gtf.name} -> {out_file.name}")
        _convert_file(gtf, out_file)
