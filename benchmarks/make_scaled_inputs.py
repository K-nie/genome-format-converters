#!/usr/bin/env python3
"""Build size-scaled benchmark inputs, so start-up cost and throughput separate.

Two modes:

``gff3``  Replicate every feature of the input GFF3 files ``--factor`` times.
          Copy *k* gets seqid ``<seqid>_c<k>`` and every ``ID``/``Parent``
          value prefixed ``c<k>_``, so copies are distinct genes on distinct
          contigs and no tool can merge them. ``##sequence-region`` lines are
          copied per contig; ``##FASTA`` sections are dropped.

``vcf``   Keep the header and the first ``--sites`` records of a VCF
          (plain or gzip). Output is gzip-compressed when it ends in ``.gz``.

Examples::

    python benchmarks/make_scaled_inputs.py gff3 --input-dir data/annotations/hs_chr22 \\
        --factor 10 --output-dir scaled/hs_chr22_x10
    python benchmarks/make_scaled_inputs.py vcf --input data/1kg/chr22.vcf.gz \\
        --sites 100000 --output scaled/chr22_100k.vcf.gz
"""
from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path

_ATTR_IDS = re.compile(r"(?<![^;])(ID|Parent|Derives_from)=([^;\n]*)")


def _open(path: Path, mode: str):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode)


def _prefix_ids(attrs: str, tag: str) -> str:
    def sub(m):
        vals = ",".join(f"{tag}_{v}" for v in m.group(2).split(","))
        return f"{m.group(1)}={vals}"
    return _ATTR_IDS.sub(sub, attrs)


def scale_gff3(src: Path, dest: Path, factor: int) -> int:
    header, regions, body = [], [], []
    with _open(src, "r") as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("##sequence-region"):
                regions.append(line.split())
            elif line.startswith("#"):
                if not body:
                    header.append(line)
            elif line.strip():
                body.append(line.rstrip("\n").split("\t"))
    n = 0
    with _open(dest, "w") as out:
        for h in header:
            if not h.startswith("##sequence-region"):
                out.write(h)
        for k in range(1, factor + 1):
            tag = f"c{k}"
            for r in regions:
                out.write(f"##sequence-region {r[1]}_{tag} {' '.join(r[2:])}\n")
            for cols in body:
                if len(cols) != 9:
                    continue
                cols = list(cols)
                cols[0] = f"{cols[0]}_{tag}"
                cols[8] = _prefix_ids(cols[8], tag)
                out.write("\t".join(cols) + "\n")
                n += 1
    return n


def head_vcf(src: Path, dest: Path, sites: int) -> int:
    n = 0
    with _open(src, "r") as fh, _open(dest, "w") as out:
        for line in fh:
            if line.startswith("#"):
                out.write(line)
                continue
            if n >= sites:
                break
            out.write(line)
            n += 1
    return n


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="mode", required=True)
    g = sub.add_parser("gff3")
    g.add_argument("--input-dir", required=True)
    g.add_argument("--output-dir", required=True)
    g.add_argument("--factor", type=int, required=True)
    v = sub.add_parser("vcf")
    v.add_argument("--input", required=True)
    v.add_argument("--output", required=True)
    v.add_argument("--sites", type=int, required=True)
    args = ap.parse_args()

    if args.mode == "gff3":
        out_dir = Path(args.output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        files = sorted(Path(args.input_dir).glob("*.gff3"))
        if not files:
            sys.exit(f"no .gff3 files in {args.input_dir}")
        total = 0
        for f in files:
            total += scale_gff3(f, out_dir / f.name, args.factor)
        print(f"{len(files)} files x{args.factor}: {total} feature rows -> {out_dir}",
              file=sys.stderr)
    else:
        out = Path(args.output)
        out.parent.mkdir(parents=True, exist_ok=True)
        n = head_vcf(Path(args.input), out, args.sites)
        if n < args.sites:
            print(f"[warn] input has only {n} records (asked for {args.sites})",
                  file=sys.stderr)
        print(f"{n} sites -> {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
