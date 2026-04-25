# Benchmark datasets

Each script skips its target if the expected output file already exists,
so these are safe to rerun.

| Script | Pulls | Used by | Approx size |
|---|---|---|---|
| `download_y1000plus.sh` | 20-species slice + optional 1154-species full set from Shen 2018 / Opulente 2024 | T3, T4, T5, T7 | 0.4 GB / 20 GB |
| `download_1kg_chr22.sh` | 1000 Genomes Phase 3, chr22, 2504 samples, biallelic SNPs | T1, T2, T8 | 1.5 GB |
| `download_peter2018_yeast_vcf.sh` | Peter et al. 2018 ~100 S. cerevisiae isolates | T1, T2 | 0.5 GB |
| `download_pfam_scan.sh` | Pre-computed HMMER `--tblout` from scanning Pfam-A against Y1000+ proteomes | T6 | 2 GB |

**No raw data is committed to the repo.** The scripts place everything
under `benchmarks/data/` which is gitignored inside `benchmarks/`.

If the upstream URLs change, update them here and in the commit message;
the numbers in the paper depend on these being stable.
