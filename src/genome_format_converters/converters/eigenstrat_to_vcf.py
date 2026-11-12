#!/usr/bin/env python3
"""Convert EIGENSTRAT (``.geno`` + ``.snp`` + ``.ind``) back to VCF.

The reverse of ``vcf-to-eigenstrat``. Each ``<stem>.geno`` in the input
directory is paired with ``<stem>.snp`` and ``<stem>.ind``.

Encoding: an EIGENSTRAT genotype counts copies of the *first* allele in the
``.snp`` row (column 5), so that allele is written as VCF REF and column 6
as ALT:

    2 -> 0/0     1 -> 0/1     0 -> 1/1     9 -> ./.

With ``--haploid`` (for pseudohaploid data, which only holds 0/2/9) calls
are written as haploid ``0`` / ``1`` / ``.``. The EIGENSTRAT REF allele is
whatever the producing tool chose (gfc's polarisation, convertf's choice,
...), so it may not match the reference genome; pass ``--ref-fasta`` to
swap alleles where the ``.snp`` first allele disagrees with the genome.

Only the text EIGENSTRAT format is read. Packed ANCESTRYMAP / PACKEDPED
``.geno`` files (binary, starting with ``GENO``/``TGENO``) are rejected with
an explanation.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional

from genome_format_converters import __version__

from ._common import die, iter_input_files, log_info, log_warn, prepare_output_dir

_DIPLOID = {"2": "0/0", "1": "0/1", "0": "1/1", "9": "./."}
_DIPLOID_SWAP = {"2": "1/1", "1": "0/1", "0": "0/0", "9": "./."}
_HAPLOID = {"2": "0", "0": "1", "9": ".", "1": "."}
_HAPLOID_SWAP = {"2": "1", "0": "0", "9": ".", "1": "."}


def _read_ind(path: Path) -> List[str]:
    out = []
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if parts:
                out.append(parts[0])
    return out


def _convert(stem: str, geno: Path, snp: Path, ind: Path, out_file: Path,
             haploid: bool, ref_fasta) -> None:
    with open(geno, "rb") as fh:
        head = fh.read(5)
    if head.startswith(b"GENO") or head.startswith(b"TGENO"):
        die(f"{geno.name} is packed ANCESTRYMAP/TGENO (binary). Convert it to "
            "text EIGENSTRAT first (convertf outputformat: EIGENSTRAT).")
    samples = _read_ind(ind)
    n = len(samples)
    straight = _HAPLOID if haploid else _DIPLOID
    swapped = _HAPLOID_SWAP if haploid else _DIPLOID_SWAP

    # Collect contig names in first-seen order for the header.
    contigs: List[str] = []
    seen = set()
    with open(snp) as fh:
        for line in fh:
            p = line.split()
            if len(p) >= 6 and p[1] not in seen:
                seen.add(p[1])
                contigs.append(p[1])

    n_sites = n_swapped = n_het_haploid = 0
    with open(snp) as sfh, open(geno) as gfh, open(out_file, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write(f"##source=gfc {__version__} eigenstrat-to-vcf\n")
        for c in contigs:
            out.write(f"##contig=<ID={c}>\n")
        out.write('##INFO=<ID=CM,Number=1,Type=Float,'
                  'Description="Genetic position (Morgans) from the .snp file">\n')
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                  + "\t".join(samples) + "\n")
        for line_no, (sline, gline) in enumerate(zip(sfh, gfh), start=1):
            p = sline.split()
            row = gline.strip()
            if len(p) < 6:
                die(f"{snp.name}:{line_no}: expected 6 columns, got {len(p)}")
            if len(row) != n:
                die(f"{geno.name}:{line_no}: {len(row)} genotypes but "
                    f"{n} individuals in {ind.name}")
            snp_id, chrom, morgans, pos, a1, a2 = p[:6]
            table = straight
            ref, alt = a1, a2
            if ref_fasta is not None:
                base = ref_fasta.fetch(chrom, int(pos) - 1, int(pos)).upper()
                if base == a2.upper() and base != a1.upper():
                    ref, alt, table = a2, a1, swapped
                    n_swapped += 1
            if haploid and "1" in row:
                n_het_haploid += row.count("1")
            calls = "\t".join(table[c] for c in row)
            info = f"CM={morgans}" if float(morgans) != 0.0 else "."
            vid = snp_id if snp_id else "."
            out.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\tPASS\t{info}\tGT\t{calls}\n")
            n_sites += 1
    log_info(f"{stem}: {n_sites} sites x {n} samples written"
             + (f"; {n_swapped} sites re-oriented to the reference" if n_swapped else ""))
    if n_het_haploid:
        log_warn(f"{stem}: --haploid set but {n_het_haploid} heterozygous (1) "
                 "calls found; written as missing.")


def batch_convert(input_dir: str, output_dir: str,
                  haploid: bool = False,
                  ref_fasta: Optional[str] = None,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)
    fasta = None
    if ref_fasta:
        from ._common import require_pysam
        try:
            import pysam
        except ImportError:
            require_pysam()
        fasta = pysam.FastaFile(ref_fasta)
    try:
        for geno in iter_input_files(in_path, [".geno"], pattern=pattern):
            stem = geno.name[: -len(".geno")]
            snp = geno.with_name(stem + ".snp")
            ind = geno.with_name(stem + ".ind")
            if not snp.is_file() or not ind.is_file():
                log_warn(f"{geno.name}: missing {snp.name} or {ind.name}; skipped.")
                continue
            out_file = out_path / (stem + ".vcf")
            log_info(f"Converting {stem}.geno/.snp/.ind -> {out_file.name}")
            _convert(stem, geno, snp, ind, out_file, haploid, fasta)
    finally:
        if fasta is not None:
            fasta.close()
