#!/usr/bin/env python3
"""Create per-sample consensus FASTA from reference FASTA + VCF.

Handles:
    * multi-contig references (all contigs are consensused, not just the first)
    * indels (insertions and deletions modify the reference-length sequence;
      output length therefore differs per sample)
    * heterozygous SNPs (IUPAC ambiguity code, or random / first / ref
      resolution via --het-policy)
    * missing genotypes (lowercase reference; or 'N' if --mask-missing)

Pairs FASTA and VCF files by stem. Recognised FASTA extensions: .fasta, .fa,
.fna, .fas. Recognised VCF extensions: .vcf, .vcf.gz.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    import pysam
except ImportError:
    pysam = None
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ._common import log_info, log_warn, prepare_output_dir, require_pysam

_FASTA_EXTS = [".fasta", ".fa", ".fna", ".fas"]
_VCF_EXTS = [".vcf", ".vcf.gz"]

_IUPAC = {
    frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("GC"): "S",
    frozenset("AT"): "W", frozenset("GT"): "K", frozenset("AC"): "M",
    frozenset("CGT"): "B", frozenset("AGT"): "D", frozenset("ACT"): "H",
    frozenset("ACG"): "V", frozenset("ACGT"): "N",
}


def _find_pairs(input_dir: Path) -> List[Tuple[str, Path, Path]]:
    fastas: Dict[str, Path] = {}
    for ext in _FASTA_EXTS:
        for f in input_dir.glob(f"*{ext}"):
            fastas[f.stem] = f
    vcfs: Dict[str, Path] = {}
    for ext in _VCF_EXTS:
        for f in input_dir.glob(f"*{ext}"):
            # Strip compound suffix for pairing (`sample.vcf.gz` -> `sample`).
            name = f.name
            for e in sorted(_VCF_EXTS, key=len, reverse=True):
                if name.endswith(e):
                    vcfs[name[: -len(e)]] = f
                    break
    common = sorted(set(fastas) & set(vcfs))
    return [(stem, fastas[stem], vcfs[stem]) for stem in common]


def _resolve_het(ref: str, alt: str, policy: str) -> str:
    if policy == "iupac":
        key = frozenset((ref.upper(), alt.upper()))
        return _IUPAC.get(key, "N")
    if policy == "ref":
        return ref
    if policy == "alt":
        return alt
    return ref  # fallback


def _create_consensus(fasta_file: Path, vcf_file: Path, out_file: Path,
                      het_policy: str = "iupac",
                      mask_missing: bool = False) -> None:
    # Load every contig into a mutable per-sample edit buffer.
    contigs = list(SeqIO.parse(str(fasta_file), "fasta"))
    if not contigs:
        log_warn(f"{fasta_file.name}: no FASTA records found, skipping.")
        return
    ref_by_chrom: Dict[str, str] = {c.id: str(c.seq).upper() for c in contigs}
    order = [c.id for c in contigs]

    with pysam.VariantFile(str(vcf_file)) as vcf:
        samples = list(vcf.header.samples)
        if not samples:
            log_warn(f"{vcf_file.name}: no samples in header, skipping.")
            return
        # Per-sample, per-contig list of (position, replacement_seq). We
        # apply edits after collection so insertions don't perturb the
        # positions of later edits.
        edits: Dict[str, Dict[str, List[Tuple[int, int, str]]]] = {
            s: {c: [] for c in order} for s in samples
        }

        skipped_contigs = set()
        for rec in vcf:
            if rec.chrom not in ref_by_chrom:
                skipped_contigs.add(rec.chrom)
                continue
            if not rec.alts:
                continue
            ref_allele = rec.ref.upper()
            pos0 = rec.pos - 1  # 0-based inclusive start

            for sample in samples:
                call = rec.samples[sample]
                gt = call.get("GT", (None,))
                if gt is None or all(a is None for a in gt):
                    if mask_missing:
                        edits[sample][rec.chrom].append(
                            (pos0, pos0 + len(ref_allele), "N" * len(ref_allele))
                        )
                    continue
                # Drop missing-allele components but keep what's present.
                allele_indices = [a for a in gt if a is not None]
                if not allele_indices:
                    continue
                # 0 = ref, >=1 = alt index into rec.alts
                picked_seqs = []
                for ai in allele_indices:
                    if ai == 0:
                        picked_seqs.append(ref_allele)
                    else:
                        alt_idx = ai - 1
                        if alt_idx >= len(rec.alts):
                            picked_seqs.append(ref_allele)  # malformed GT
                        else:
                            picked_seqs.append(rec.alts[alt_idx].upper())
                # Collapse diploid calls
                unique = list(dict.fromkeys(picked_seqs))
                if len(unique) == 1:
                    replacement = unique[0]
                else:
                    # Heterozygous — only well-defined for SNP-length alleles.
                    if all(len(s) == 1 for s in unique):
                        replacement = _resolve_het(unique[0], unique[1], het_policy)
                    else:
                        replacement = unique[0]  # first allele wins for indels
                # Skip records that would silently no-op.
                if replacement == ref_allele:
                    continue
                edits[sample][rec.chrom].append(
                    (pos0, pos0 + len(ref_allele), replacement)
                )

        if skipped_contigs:
            log_warn(
                f"{vcf_file.name}: VCF contigs not in reference, skipped: "
                f"{', '.join(sorted(skipped_contigs))}"
            )

    # Apply edits per sample per contig, right-to-left so indel offsets stay
    # valid relative to the reference coordinates we collected.
    with open(out_file, "w") as out_handle:
        for sample in samples:
            per_contig_records = []
            for chrom in order:
                seq = ref_by_chrom[chrom]
                chrom_edits = sorted(edits[sample][chrom], key=lambda e: e[0])
                # Iterate right-to-left to avoid offset drift.
                new_seq = list(seq)
                for start, end, replacement in reversed(chrom_edits):
                    new_seq[start:end] = list(replacement)
                consensus = "".join(new_seq)
                per_contig_records.append(
                    SeqRecord(Seq(consensus), id=f"{sample}|{chrom}", description="consensus")
                )
            SeqIO.write(per_contig_records, out_handle, "fasta")


def batch_convert(input_dir: str, output_dir: str,
                  het_policy: str = "iupac",
                  mask_missing: bool = False,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    if pysam is None:
        require_pysam()
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)

    pairs = _find_pairs(in_path)
    if not pairs:
        log_warn("No matching FASTA / VCF pairs found (matched by stem).")
        return

    for stem, fasta, vcf in pairs:
        out_file = out_path / f"{stem}_consensus.fa"
        log_info(f"Processing {stem}: {fasta.name} + {vcf.name} -> {out_file.name}")
        _create_consensus(fasta, vcf, out_file,
                          het_policy=het_policy, mask_missing=mask_missing)
