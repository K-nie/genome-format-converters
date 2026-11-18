#!/usr/bin/env python3
"""Write a text EIGENSTRAT copy of a packed (PACKEDANCESTRYMAP) dataset.

Used only to build the E24 input, so that gfc (text .geno) and plink2
(packed .geno) read the same genotypes. A packed .geno starts with a
header record "GENO <nind> <nsnp> ..." of max(48, ceil(nind/4)) bytes, then
one record of that length per SNP, 2 bits per individual from the high bits
down; 3 means missing (written as 9).

Usage: unpack_geno.py <packed prefix> <text prefix>
"""
import shutil
import sys

src, dst = sys.argv[1], sys.argv[2]
lut = [bytes(b"0129"[(v >> s) & 3] for s in (6, 4, 2, 0)) for v in range(256)]
with open(src + ".geno", "rb") as fh, open(dst + ".geno", "wb") as out:
    head = fh.read(48)
    assert head.startswith(b"GENO "), "not a packed PACKEDANCESTRYMAP .geno"
    nind = int(head.split()[1])
    nbytes = -(-nind // 4)
    rlen = max(48, nbytes)
    fh.read(rlen - 48)
    while True:
        rec = fh.read(rlen)
        if not rec:
            break
        out.write(b"".join(lut[b] for b in rec[:nbytes])[:nind] + b"\n")
for suffix in (".snp", ".ind"):
    shutil.copyfile(src + suffix, dst + suffix)
