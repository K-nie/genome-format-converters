#!/usr/bin/env python3
"""Convert a VCF to a SNP alignment (relaxed PHYLIP, FASTA or NEXUS).

One sequence per sample, one column per biallelic SNP, for phylogenetic
inference (IQ-TREE, RAxML-NG) from variant calls:

* homozygous REF / ALT -> that base;
* heterozygous -> the IUPAC ambiguity code of the two bases
  (``--het random`` instead picks one allele at random, seeded per file;
  ``--het missing`` writes ``N``);
* missing -> ``N``.

Indels, multi-allelic and non-ACGT sites are skipped. ``--min-samples``
drops sites called in fewer samples. Output names: ``<stem>.phy``,
``<stem>.fasta`` or ``<stem>.nex``.

When the alignment will be analysed with an ascertainment-bias correction
(IQ-TREE ``+ASC``), note that the default IUPAC coding keeps sites that are
invariant once ambiguity is ignored; use ``--drop-invariant`` to remove
columns where all called samples share one base.
"""

from __future__ import annotations

import zlib
from pathlib import Path
from typing import Optional

import numpy as np

try:
    import cyvcf2
except ImportError:
    cyvcf2 = None

from ._common import (
    iter_input_files,
    log_info,
    prepare_output_dir,
    strip_compound_suffix,
)
from ._eigenstrat_common import BASES
from ._vcf_fast import iter_records

_VCF_EXTS = [".vcf", ".vcf.gz", ".bcf"]
_IUPAC = {frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("CG"): "S",
          frozenset("AT"): "W", frozenset("GT"): "K", frozenset("AC"): "M"}
_EXT = {"phylip": ".phy", "fasta": ".fasta", "nexus": ".nex"}


def _columns(in_file: Path, het: str, min_samples: int, drop_invariant: bool,
             seed: Optional[int]):
    rng = np.random.default_rng(
        None if seed is None else seed ^ zlib.crc32(in_file.name.encode()))
    n_ = ord("N")
    with cyvcf2.VCF(str(in_file)) as vcf:
        samples = list(vcf.samples)
        cols = []
        for var in iter_records(in_file, vcf):
            if not var.ALT or len(var.ALT) != 1:
                continue
            ref, alt = var.REF.upper(), var.ALT[0].upper()
            if len(ref) != 1 or len(alt) != 1 or ref not in BASES or alt not in BASES:
                continue
            gt = np.asarray(var.gt_types)
            called = int((gt != 2).sum())
            if called < min_samples:
                continue
            col = np.full(len(samples), n_, dtype=np.uint8)
            col[gt == 0] = ord(ref)
            col[gt == 3] = ord(alt)
            hets = gt == 1
            if het == "iupac":
                col[hets] = ord(_IUPAC[frozenset((ref, alt))])
            elif het == "random":
                pick = rng.integers(0, 2, size=int(hets.sum()))
                col[hets] = np.where(pick == 0, ord(ref), ord(alt))
            if drop_invariant:
                bases = set(col[(col != n_)].tolist())
                if len(bases) <= 1:
                    continue
            cols.append(col)
    if cols:
        mat = np.stack(cols, axis=1)
    else:
        mat = np.zeros((len(samples), 0), dtype=np.uint8)
    return samples, mat


def _write(out_file: Path, fmt: str, samples, mat) -> None:
    seqs = [row.tobytes().decode("ascii") for row in mat]
    n_sites = mat.shape[1]
    with open(out_file, "w") as fh:
        if fmt == "phylip":
            fh.write(f"{len(samples)} {n_sites}\n")
            width = max((len(s) for s in samples), default=1) + 1
            for name, seq in zip(samples, seqs):
                fh.write(f"{name.ljust(width)}{seq}\n")
        elif fmt == "fasta":
            for name, seq in zip(samples, seqs):
                fh.write(f">{name}\n")
                for i in range(0, len(seq), 60):
                    fh.write(seq[i:i + 60] + "\n")
                if not seq:
                    fh.write("\n")
        else:
            fh.write("#NEXUS\nBEGIN DATA;\n")
            fh.write(f"  DIMENSIONS NTAX={len(samples)} NCHAR={n_sites};\n")
            fh.write("  FORMAT DATATYPE=DNA MISSING=N GAP=-;\n  MATRIX\n")
            width = max((len(s) for s in samples), default=1) + 2
            for name, seq in zip(samples, seqs):
                fh.write(f"  {name.ljust(width)}{seq}\n")
            fh.write("  ;\nEND;\n")


def batch_convert(input_dir: str, output_dir: str,
                  out_format: str = "phylip",
                  het: str = "iupac",
                  min_samples: int = 0,
                  drop_invariant: bool = False,
                  seed: Optional[int] = None,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    if cyvcf2 is None:
        raise ImportError("cyvcf2 is required for vcf-to-snp-alignment. "
                          "Install via `pip install cyvcf2`.")
    if out_format not in _EXT:
        raise ValueError(f"unknown format {out_format!r}")
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)
    for vcf in iter_input_files(in_path, _VCF_EXTS, pattern=pattern):
        stem = strip_compound_suffix(vcf, _VCF_EXTS)
        out_file = out_path / (stem + _EXT[out_format])
        log_info(f"Converting {vcf.name} -> {out_file.name}")
        samples, mat = _columns(vcf, het, min_samples, drop_invariant, seed)
        _write(out_file, out_format, samples, mat)
        log_info(f"{vcf.name}: {len(samples)} samples x {mat.shape[1]} SNP columns")
