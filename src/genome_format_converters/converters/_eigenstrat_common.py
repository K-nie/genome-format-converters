"""Shared helpers for the EIGENSTRAT-family converters.

Implements the bits that are identical between `vcf_to_eigenstrat` and
`vcf_to_pseudohaploid`: site filtering, allele polarisation, map-file
loading, PLINK-compatible sidecar output, and the `.snp` / `.ind` writers.
Only the per-sample genotype-encoding step differs between the two
converters, and that's left to the caller.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from ._common import log_warn, die

BASES = {"A", "C", "G", "T"}
_TRANSITIONS = {frozenset("AG"), frozenset("CT")}


# ---------- map-file loaders ------------------------------------------------


def load_genetic_map(path: Optional[str]) -> Dict[Tuple[str, int], float]:
    """Load a PLINK-style genetic map (`chrom snp_id cM bp`).

    Returns a dict keyed by `(chrom, bp)` -> cM. Morgans are 100× cM.
    Callers should handle missing sites by falling back to 0.0.
    """
    out: Dict[Tuple[str, int], float] = {}
    if path is None:
        return out
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            chrom = parts[0]
            try:
                cM = float(parts[2])
                bp = int(parts[3])
            except ValueError:
                continue
            out[(chrom, bp)] = cM
    return out


def morgans_for(snp_chrom: str, snp_pos: int,
                gmap: Dict[Tuple[str, int], float]) -> float:
    """Look up a genetic position in Morgans; 0.0 when unknown."""
    if not gmap:
        return 0.0
    cM = gmap.get((snp_chrom, snp_pos))
    if cM is None:
        return 0.0
    return cM / 100.0


def load_chrom_map(path: Optional[str]) -> Dict[str, str]:
    """Load a two-column TSV mapping contig names to smartpca-safe
    integer codes (or whatever the user supplies)."""
    out: Dict[str, str] = {}
    if path is None:
        return out
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            out[parts[0]] = parts[1]
    return out


def resolve_chrom(chrom: str,
                  chrom_map: Dict[str, str],
                  default_chrom: Optional[int]) -> str:
    """Apply the user's chrom-map; fall through to default when absent.

    Returns the mapped value as a string (to allow any label), or the
    original chromosome when no mapping is given.
    """
    if chrom in chrom_map:
        return chrom_map[chrom]
    if default_chrom is not None:
        return str(default_chrom)
    return chrom


# ---------- polarisation -----------------------------------------------------


class AncestralProvider:
    """Lazily opens an ancestral FASTA (pysam.FastaFile) and yields the
    ancestral base at a given position. Returns None when no ancestral
    call can be made (contig missing, out of bounds, non-ACGT)."""

    def __init__(self, path: Optional[str]):
        self._path = path
        self._fa = None

    def __enter__(self):
        if self._path is not None:
            import pysam
            self._fa = pysam.FastaFile(self._path)
        return self

    def __exit__(self, exc_type, exc, tb):
        if self._fa is not None:
            self._fa.close()

    def get(self, chrom: str, pos: int) -> Optional[str]:
        if self._fa is None:
            return None
        try:
            base = self._fa.fetch(chrom, pos - 1, pos).upper()
        except (KeyError, ValueError):
            return None
        if base and base[0] in BASES:
            return base[0]
        return None


def polarise(ref: str, alt: str, chrom: str, pos: int,
             ancestral: Optional[AncestralProvider],
             info_aa: Optional[str]) -> Tuple[str, str, bool]:
    """Return (major, minor, flipped) where `major` is the allele to be
    encoded as reference in EIGENSTRAT (i.e. contributing to genotype
    code `2`). `flipped` is True when we swapped ref/alt so the caller
    can invert genotype counts.

    Priority: INFO/AA field if requested and present, else ancestral
    FASTA if provided, else the VCF REF allele.
    """
    ancestral_base: Optional[str] = None

    if info_aa:
        raw = info_aa.strip().split("|")[0].upper()  # 1000G-style
        if raw and raw[0] in BASES:
            ancestral_base = raw[0]

    if ancestral_base is None and ancestral is not None:
        ancestral_base = ancestral.get(chrom, pos)

    if ancestral_base is None:
        return ref, alt, False

    if ancestral_base == ref:
        return ref, alt, False
    if ancestral_base == alt:
        return alt, ref, True
    # Ancestral call disagrees with both ref and alt — can't polarise.
    return ref, alt, False


# ---------- site filters -----------------------------------------------------


def is_transition(ref: str, alt: str) -> bool:
    return frozenset((ref, alt)) in _TRANSITIONS


def passes_maf(allele_codes: Sequence[str], min_maf: float) -> bool:
    """Check that the minor-allele frequency exceeds the threshold.
    `allele_codes` is the per-sample EIGENSTRAT digits ('0'/'1'/'2'/'9')
    written for this site.
    """
    if min_maf <= 0:
        return True
    n_ref = 0
    n_alt = 0
    for c in allele_codes:
        if c == "2":
            n_ref += 2
        elif c == "1":
            n_ref += 1
            n_alt += 1
        elif c == "0":
            n_alt += 2
        # "9" is missing
    total = n_ref + n_alt
    if total == 0:
        return False
    p = n_ref / total
    maf = min(p, 1.0 - p)
    return maf >= min_maf


def passes_missing(allele_codes: Sequence[str], max_missing: float) -> bool:
    if max_missing >= 1.0:
        return True
    if not allele_codes:
        return False
    miss = sum(1 for c in allele_codes if c == "9")
    return (miss / len(allele_codes)) <= max_missing


# ---------- map validation --------------------------------------------------


def validate_map_coverage(name: str,
                          mapping: Dict[str, str],
                          universe: Iterable[str],
                          strict: bool) -> None:
    """Warn (or fail under --strict-maps) when the map mentions keys that
    don't exist in the VCF samples / contigs."""
    known = set(universe)
    unknown = [k for k in mapping if k not in known]
    if not unknown:
        return
    msg = (f"{name}: {len(unknown)} key(s) in the map are not present "
           f"in the VCF: {', '.join(unknown[:5])}"
           + (" ..." if len(unknown) > 5 else ""))
    if strict:
        die(msg)
    else:
        log_warn(msg)


# ---------- writers ----------------------------------------------------------


def write_ind(ind_path: Path,
              samples: Sequence[str],
              pop_map: Dict[str, str],
              sex_map: Dict[str, str]) -> None:
    with open(ind_path, "w") as ind:
        for s in samples:
            sex = sex_map.get(s, "U")
            pop = pop_map.get(s, s)
            ind.write(f"{s}\t{sex}\t{pop}\n")


def write_plink_sidecar(out_prefix: Path,
                        samples: Sequence[str],
                        pop_map: Dict[str, str],
                        sex_map: Dict[str, str],
                        snp_rows: Sequence[Tuple[str, str, float, int, str, str]],
                        geno_rows: Sequence[str]) -> None:
    """Emit PLINK-compatible `.ped` + `.map` sidecar files.

    `.ped` is space-separated: family_id, individual_id, paternal_id,
    maternal_id, sex (1=M, 2=F, other=unknown), phenotype (-9), then two
    allele columns per SNP. `.map` is: chrom, snp_id, cM, bp.
    """
    ped_path = out_prefix.with_suffix(".ped")
    map_path = out_prefix.with_suffix(".map")
    sex_code = {"M": "1", "F": "2"}
    with open(map_path, "w") as mfh:
        for snp_id, chrom, morgans, bp, _ref, _alt in snp_rows:
            mfh.write(f"{chrom}\t{snp_id}\t{morgans * 100.0:.6f}\t{bp}\n")
    with open(ped_path, "w") as pfh:
        for s_idx, sample in enumerate(samples):
            fam = pop_map.get(sample, sample)
            sex_raw = sex_map.get(sample, "U")
            sex_col = sex_code.get(sex_raw, "0")
            row = [fam, sample, "0", "0", sex_col, "-9"]
            for snp_idx, (_snp_id, _chrom, _m, _bp, ref, alt) in enumerate(snp_rows):
                code = geno_rows[snp_idx][s_idx]
                if code == "2":
                    row += [ref, ref]
                elif code == "1":
                    row += [ref, alt]
                elif code == "0":
                    row += [alt, alt]
                else:  # "9"
                    row += ["0", "0"]
            pfh.write(" ".join(row) + "\n")
