"""GenBank to GFF3 conversion.

Each GenBank record gives a `##sequence-region` line and one `region` line
for the whole sequence (carrying the organism and molecule type). Each
feature gives one GFF3 line per location part: a joined (spliced or
trans-spliced) location such as `join(1..100,201..300)` becomes several
lines that share one `ID`, as GFF3 requires for discontinuous features,
rather than one line spanning the intron. CDS parts get their phase from
`/codon_start` and the lengths of the parts before them. Qualifiers are
kept as attributes, escaped as GFF3 requires. `source` features are
covered by the `region` line and are not repeated.
"""

from pathlib import Path
from typing import Dict, Optional
from urllib.parse import quote

from Bio import SeqIO

from ._common import iter_input_files, log_info, prepare_output_dir

_GBK_EXTS = [".gbk", ".gb", ".gbff"]
# GFF3 column 9 reserves ; = & , and control characters; % escapes itself.
_SAFE = " !\"#$'()*+-./:<>?@[\\]^_`{|}~"


def _esc(value) -> str:
    return quote(str(value), safe=_SAFE)


def _strand(strand) -> str:
    return {1: "+", -1: "-"}.get(strand, ".")


def _feature_id(feat, used: Dict[str, int]) -> str:
    q = feat.qualifiers
    base = None
    for key in ("locus_tag", "gene", "protein_id", "label"):
        if q.get(key):
            base = f"{feat.type}:{q[key][0]}"
            break
    base = base or feat.type
    used[base] = used.get(base, 0) + 1
    return base if used[base] == 1 else f"{base}.{used[base]}"


def _write_record(rec, out) -> None:
    length = len(rec.seq) if rec.seq is not None else 0
    seqid = rec.id
    out.write(f"##sequence-region {seqid} 1 {length}\n")
    attrs = [f"ID={_esc(seqid)}"]
    ann = rec.annotations
    for key in ("organism", "molecule_type", "topology"):
        if ann.get(key):
            attrs.append(f"{key}={_esc(ann[key])}")
    out.write(f"{seqid}\tGenBank\tregion\t1\t{length}\t.\t+\t.\t{';'.join(attrs)}\n")
    used: Dict[str, int] = {}
    for feat in rec.features:
        if feat.type == "source":
            continue
        parts = feat.location.parts if feat.location is not None else []
        if not parts:
            continue
        fid = _feature_id(feat, used)
        attrs = [f"ID={_esc(fid)}"]
        for key, values in feat.qualifiers.items():
            if key in ("ID", "Parent"):
                key = key.lower()
            vals = values if isinstance(values, list) else [values]
            # Flag qualifiers (/pseudo, /trans_splicing) have no value;
            # GFF3 attributes need one.
            vals = [v if v != "" else "true" for v in vals]
            attrs.append(f"{_esc(key)}=" + ",".join(_esc(v) for v in vals))
        attr_text = ";".join(attrs)
        codon_start = int(feat.qualifiers.get("codon_start", ["1"])[0])
        done = 0  # coding bases before this part, after the codon_start offset
        for i, part in enumerate(parts):
            start, end = int(part.start) + 1, int(part.end)
            if feat.type == "CDS":
                phase = str(codon_start - 1 if i == 0 else (3 - done % 3) % 3)
                done += (end - start + 1) - (codon_start - 1 if i == 0 else 0)
            else:
                phase = "."
            out.write(f"{seqid}\tGenBank\t{feat.type}\t{start}\t{end}\t.\t"
                      f"{_strand(part.strand)}\t{phase}\t{attr_text}\n")


def _convert_file(in_file: Path, out_file: Path) -> None:
    with open(in_file) as in_handle, open(out_file, "w") as out_handle:
        out_handle.write("##gff-version 3\n")
        for rec in SeqIO.parse(in_handle, "genbank"):
            _write_record(rec, out_handle)


def batch_convert(input_dir: str, output_dir: str,
                  force: bool = False,
                  pattern: Optional[str] = None) -> None:
    in_path = Path(input_dir)
    out_path = prepare_output_dir(output_dir, force=force)
    for gbk in iter_input_files(in_path, _GBK_EXTS, pattern=pattern):
        out_file = out_path / (gbk.stem + ".gff3")
        log_info(f"Converting {gbk.name} -> {out_file.name}")
        _convert_file(gbk, out_file)
