#!/usr/bin/env python3
# Author: Benjamin Narh-Madey
"""T6 reference: Bio.SearchIO.HmmerIO HMMER tblout -> TSV parser.

Citable competitor for T6 (HMMER tblout parsing). Exercises Biopython's
HmmerIO 'hmmer3-tab' parser, the most-cited Python parser for HMMER
outputs (Cock 2009, Bioinformatics 25:1422). The pyhmmer paper
(Larralde 2023) explicitly benchmarks against this code path; including
both pyhmmer and biopython rows in T6 lets the same comparison appear
in our table without re-running their numbers second-hand.

Output schema matches what gfc hmmer-tblout-to-tsv produces: one row
per hit, with target_name, target_accession, query_name,
query_accession, full_evalue, full_score, full_bias, best_evalue,
best_score, best_bias, exp, reg, clu, ov, env, dom, rep, inc,
description.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def parse_one(tblout: Path, output: Path) -> int:
    from Bio import SearchIO

    cols = [
        "target_name", "target_accession",
        "query_name", "query_accession",
        "full_evalue", "full_score", "full_bias",
        "best_evalue", "best_score", "best_bias",
        "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc",
        "description",
    ]
    n = 0
    with output.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for qresult in SearchIO.parse(str(tblout), "hmmer3-tab"):
            for hit in qresult:
                # hmmer3-tab puts per-hit values in hit.* and per-query
                # values in qresult.*. Description is on the hit object.
                row = [
                    hit.id,
                    getattr(hit, "accession", "-") or "-",
                    qresult.id,
                    getattr(qresult, "accession", "-") or "-",
                    f"{hit.evalue:g}",
                    f"{hit.bitscore:g}",
                    f"{getattr(hit, 'bias', 0.0):g}",
                    f"{hit.hsps[0].evalue:g}" if hit.hsps else "-",
                    f"{hit.hsps[0].bitscore:g}" if hit.hsps else "-",
                    f"{getattr(hit.hsps[0], 'bias', 0.0):g}" if hit.hsps else "-",
                    # exp/reg/clu/ov/env/dom/rep/inc — Biopython exposes
                    # the domain-estimation block via hit.domain_*
                    # attributes when present; fall back to "-".
                    f"{getattr(hit, 'domain_exp_num', '-')}",
                    f"{getattr(hit, 'domain_reg_num', '-')}",
                    f"{getattr(hit, 'domain_clu_num', '-')}",
                    f"{getattr(hit, 'domain_ov_num', '-')}",
                    f"{getattr(hit, 'domain_env_num', '-')}",
                    f"{getattr(hit, 'domain_obs_num', '-')}",
                    f"{getattr(hit, 'domain_reported_num', '-')}",
                    f"{getattr(hit, 'domain_included_num', '-')}",
                    (hit.description or "").strip(),
                ]
                fh.write("\t".join(str(v) for v in row) + "\n")
                n += 1
    return n


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tblout", required=True, type=Path,
                    help="Input HMMER tblout file")
    ap.add_argument("--output", required=True, type=Path,
                    help="Output TSV path")
    args = ap.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    n = parse_one(args.tblout, args.output)
    print(f"[t6_biopython_searchio] wrote {n} hits to {args.output}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
