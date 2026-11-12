"""Fast VCF record iterator for the EIGENSTRAT / PLINK / pseudohaploid paths.

cyvcf2 parses every FORMAT field of every sample into htslib's BCF
representation before Python sees the record, and ``Variant.gt_types`` then
walks that structure again. On 1000 Genomes-style VCFs (thousands of samples,
``FORMAT=GT`` only) that work dominates wall time, even though the
converters only need a 0/1/2/3 genotype class per sample.

This module reads plain or gzip/bgzip-compressed VCF text directly and, when
a record's FORMAT is exactly ``GT`` and every call is a diploid ``a/b`` or
``a|b`` with alleles in {0, 1, .}, decodes the whole sample block with a few
numpy operations on the raw bytes. Any record outside that shape (other
FORMAT fields, haploid or half-missing calls, multi-digit allele indices) is
handed to cyvcf2's own parser via ``Writer.variant_from_string`` so its
``gt_types`` is exactly what cyvcf2 would have produced. BCF input always
goes through cyvcf2.

The records are duck-typed to the subset of ``cyvcf2.Variant`` the
converters use: ``CHROM``, ``POS``, ``ID``, ``REF``, ``ALT``, ``INFO.get``,
and ``gt_types``. ``gt_types`` is computed lazily, so sites rejected by the
REF/ALT filters never pay for genotype decoding.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Iterator, Optional

import numpy as np

# gt_types codes, matching cyvcf2 (gts012=False).
HOM_REF, HET, UNKNOWN, HOM_ALT = 0, 1, 2, 3

_TAB = ord("\t")
_SLASH = ord("/")
_PIPE = ord("|")
_ZERO = ord("0")
_ONE = ord("1")
_DOT = ord(".")

# Lookup for a 2-allele call encoded as (a1 code, a2 code) with
# code 0 = '0', 1 = '1', 2 = '.'. Index = a1 * 3 + a2.
# Only the combinations accepted by the fast path are used; half-missing
# calls are routed to cyvcf2 because its semantics for them are irregular.
_ALLELE_CODE = np.full(256, 255, dtype=np.uint8)
_ALLELE_CODE[_ZERO] = 0
_ALLELE_CODE[_ONE] = 1
_ALLELE_CODE[_DOT] = 2
_PAIR_TO_GT = np.array(
    [HOM_REF, HET, 255,
     HET, HOM_ALT, 255,
     255, 255, UNKNOWN],
    dtype=np.int32,
)


class _Info:
    __slots__ = ("_raw", "_parsed")

    def __init__(self, raw: str):
        self._raw = raw
        self._parsed: Optional[dict] = None

    def get(self, key, default=None):
        if self._parsed is None:
            parsed = {}
            if self._raw not in ("", "."):
                for item in self._raw.split(";"):
                    k, sep, v = item.partition("=")
                    parsed[k] = v if sep else True
            self._parsed = parsed
        return self._parsed.get(key, default)


class FastRecord:
    """Minimal stand-in for ``cyvcf2.Variant``."""

    __slots__ = ("CHROM", "POS", "ID", "REF", "ALT", "_info", "_fmt",
                 "_samples", "_line", "_reader", "_gt")

    def __init__(self, fields, line, reader):
        self.CHROM = fields[0].decode()
        self.POS = int(fields[1])
        vid = fields[2].decode()
        self.ID = None if vid == "." else vid
        self.REF = fields[3].decode()
        alt = fields[4].decode()
        self.ALT = [] if alt == "." else alt.split(",")
        self._info = fields[7]
        self._fmt = fields[8] if len(fields) > 8 else b""
        self._samples = fields[9] if len(fields) > 9 else b""
        self._line = line
        self._reader = reader
        self._gt = None

    @property
    def INFO(self) -> _Info:
        return _Info(self._info.decode())

    @property
    def gt_types(self) -> np.ndarray:
        if self._gt is None:
            gt = None
            if self._fmt == b"GT":
                gt = _decode_gt_block(self._samples, self._reader.n_samples)
            if gt is None:
                gt = self._reader.slow_gt_types(self._line)
            self._gt = gt
        return self._gt


def _decode_gt_block(block: bytes, n_samples: int) -> Optional[np.ndarray]:
    """Vectorised decode of a GT-only sample block, or None if irregular."""
    if n_samples == 0 or len(block) != 4 * n_samples - 1:
        return None
    arr = np.frombuffer(block, dtype=np.uint8)
    if n_samples > 1 and not np.all(arr[3::4] == _TAB):
        return None
    sep = arr[1::4]
    if not np.all((sep == _SLASH) | (sep == _PIPE)):
        return None
    a1 = _ALLELE_CODE[arr[0::4]]
    a2 = _ALLELE_CODE[arr[2::4]]
    if a1.max() == 255 or a2.max() == 255:
        return None
    gt = _PAIR_TO_GT[a1.astype(np.int32) * 3 + a2]
    if gt.max() == 255:  # half-missing call such as 0/. -> defer to cyvcf2
        return None
    return gt


class _Reader:
    def __init__(self, path: Path, raw_header: str, n_samples: int):
        self.path = path
        self.raw_header = raw_header
        self.n_samples = n_samples
        self._writer = None

    def slow_gt_types(self, line: bytes) -> np.ndarray:
        if self._writer is None:
            import cyvcf2
            self._writer = cyvcf2.Writer.from_string("/dev/null", self.raw_header)
        var = self._writer.variant_from_string(line.decode().rstrip("\r\n"))
        # Copy while `var` is alive: gt_types can be a view into memory that
        # htslib frees with the record.
        return np.array(var.gt_types, dtype=np.int32, copy=True)


def _open_text(path: Path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rb")
    return open(path, "rb")


def supports_fast_path(path: Path) -> bool:
    """True for VCF text (plain or gzip/bgzip). BCF stays on cyvcf2."""
    name = path.name.lower()
    return not name.endswith(".bcf")


def iter_records(path: Path, vcf) -> Iterator:
    """Yield records from ``path``.

    ``vcf`` is an already-open ``cyvcf2.VCF`` on the same file; it supplies
    the header (samples, raw header for the slow path). For BCF input the
    cyvcf2 iterator itself is returned unchanged.
    """
    path = Path(path)
    if not supports_fast_path(path):
        yield from vcf
        return
    reader = _Reader(path, vcf.raw_header, len(vcf.samples))
    with _open_text(path) as fh:
        for line in fh:
            if not line or line[0] == 0x23:  # '#'
                continue
            fields = line.rstrip(b"\r\n").split(b"\t", 9)
            if len(fields) < 8:
                continue
            yield FastRecord(fields, line, reader)
