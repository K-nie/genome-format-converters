#!/usr/bin/env python3
# Author: Benjamin Narh-Madey
"""Phase 3 stats: paired Wilcoxon + BH-FDR + Cliff's delta + bootstrap CI.

For every (task, competitor) pair, pairs replicates by ``replicate``
number against the matching gfc rep, then runs:

  * Paired Wilcoxon signed-rank on log10(metric) for both wall_s and
    peak_rss_mb. log-transform makes the test work on the same scale
    as the ratio claims in the paper (a 2x slowdown is the same
    distance as a 0.5x speedup once you log-transform).
  * Cliff's delta as a non-parametric effect size on the raw paired
    deltas (signed magnitude, range [-1, 1]).
  * Median ratio (competitor / gfc) with a 10000-resample bootstrap
    95 % CI. Median chosen over mean because the wall-time
    distributions can be right-skewed by occasional cluster
    interference; the median is robust to it.

Multiple-testing correction: Benjamini-Hochberg FDR at q = 0.05
across all (task x competitor x metric) tests; reported alongside
the raw p-value in ``p_adj_bh``.

Output: ``benchmarks/results/figures/stats_pairs.tsv`` with one
row per (task, competitor, metric) test. Columns:
  task, competitor, metric, n_pairs, wilcoxon_W, wilcoxon_p,
  p_adj_bh, cliffs_delta, median_ratio, ci_lo_95, ci_hi_95,
  gfc_wins

  ``gfc_wins`` is True when median_ratio > 1 (competitor takes
  more wall / RAM than gfc) AND the test is significant after
  BH-FDR correction. Strict definition: a "clean win" requires
  both directionality and significance.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

_RAW = Path(__file__).resolve().parent.parent / "results" / "raw" / "all.tsv"
_OUT = Path(__file__).resolve().parent.parent / "results" / "figures" / "stats_pairs.tsv"

_METRICS = ("wall_s", "peak_rss_mb")
_BOOT_N = 10_000
_FDR_Q = 0.05
_RNG = np.random.default_rng(20260503)


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """Cliff's delta on paired arrays. Range [-1, 1]; positive means
    a tends to exceed b, negative the reverse. Computed on the signed
    paired deltas (a - b) by sign tally, which is the paired-sample
    formulation; equivalent to the rank-biserial correlation of the
    Wilcoxon test statistic on the same paired deltas."""
    diff = np.asarray(a) - np.asarray(b)
    n = len(diff)
    if n == 0:
        return float("nan")
    pos = int((diff > 0).sum())
    neg = int((diff < 0).sum())
    return (pos - neg) / n


def median_ratio_with_ci(
    comp: np.ndarray, gfc: np.ndarray, n_boot: int = _BOOT_N
) -> tuple[float, float, float]:
    """Point estimate = median(comp/gfc) on the paired ratios.
    95 % CI = bootstrap percentile interval (resample paired indices
    with replacement n_boot times). The pairing is preserved on each
    resample so the CI reflects the dependency between gfc and
    competitor reps within the same replicate slot."""
    comp = np.asarray(comp, dtype=float)
    gfc = np.asarray(gfc, dtype=float)
    ratios = comp / gfc
    point = float(np.median(ratios))
    n = len(ratios)
    if n < 2:
        return point, float("nan"), float("nan")
    boot_medians = np.empty(n_boot)
    for i in range(n_boot):
        idx = _RNG.integers(0, n, size=n)
        boot_medians[i] = np.median(ratios[idx])
    lo, hi = np.percentile(boot_medians, [2.5, 97.5])
    return point, float(lo), float(hi)


def bh_adjust(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg step-up adjusted p-values. Returns array
    of same shape, monotonically corrected so larger raw p doesn't
    yield smaller adjusted p (standard BH formulation)."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * n / (np.arange(n) + 1)
    # enforce monotonicity from the largest p downward
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0.0, 1.0)
    out = np.empty_like(adj)
    out[order] = adj
    return out


def build_pair_rows(df: pd.DataFrame) -> pd.DataFrame:
    """One row per (task, competitor, metric)."""
    df = df[df["exit_code"] == 0].copy()
    rows: list[dict] = []

    tasks = sorted(df["task"].unique())
    for task in tasks:
        sub = df[df["task"] == task]
        if "gfc" not in set(sub["tool"]):
            continue
        gfc_rows = sub[sub["tool"] == "gfc"].set_index("replicate")
        comps = sorted(t for t in sub["tool"].unique() if t != "gfc")
        for comp in comps:
            comp_rows = sub[sub["tool"] == comp].set_index("replicate")
            paired_idx = sorted(set(gfc_rows.index) & set(comp_rows.index))
            if len(paired_idx) < 2:
                continue
            for metric in _METRICS:
                gfc_v = gfc_rows.loc[paired_idx, metric].to_numpy(dtype=float)
                comp_v = comp_rows.loc[paired_idx, metric].to_numpy(dtype=float)
                # log-transform for the Wilcoxon (ratio-scale geometry)
                log_diff = np.log10(comp_v) - np.log10(gfc_v)
                # scipy.stats.wilcoxon: paired one-sample on log_diff,
                # tests H0: median(diff) == 0. Default zero_method="wilcox"
                # drops zero-diffs; with n=10 we use the default.
                try:
                    res = stats.wilcoxon(log_diff, alternative="two-sided",
                                         zero_method="wilcox")
                    w_stat = float(res.statistic)
                    w_p = float(res.pvalue)
                except ValueError:
                    # all-zero diffs (identical reps) -> no test
                    w_stat = float("nan")
                    w_p = float("nan")
                delta = cliffs_delta(comp_v, gfc_v)
                med_r, lo, hi = median_ratio_with_ci(comp_v, gfc_v)
                rows.append({
                    "task": task,
                    "competitor": comp,
                    "metric": metric,
                    "n_pairs": len(paired_idx),
                    "wilcoxon_W": w_stat,
                    "wilcoxon_p": w_p,
                    "cliffs_delta": delta,
                    "median_ratio": med_r,
                    "ci_lo_95": lo,
                    "ci_hi_95": hi,
                })
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    valid = out["wilcoxon_p"].notna()
    out["p_adj_bh"] = float("nan")
    out.loc[valid, "p_adj_bh"] = bh_adjust(out.loc[valid, "wilcoxon_p"].to_numpy())
    out["gfc_wins"] = (out["median_ratio"] > 1.0) & (out["p_adj_bh"] < _FDR_Q)
    return out


def main() -> int:
    if not _RAW.exists():
        print(f"[error] {_RAW} missing.", file=sys.stderr)
        return 1
    df = pd.read_csv(_RAW, sep="\t")
    pair_df = build_pair_rows(df)
    if pair_df.empty:
        print("[error] no pair rows produced.", file=sys.stderr)
        return 1
    _OUT.parent.mkdir(parents=True, exist_ok=True)
    cols = ["task", "competitor", "metric", "n_pairs",
            "wilcoxon_W", "wilcoxon_p", "p_adj_bh",
            "cliffs_delta", "median_ratio", "ci_lo_95", "ci_hi_95",
            "gfc_wins"]
    pair_df[cols].to_csv(_OUT, sep="\t", index=False,
                         float_format="%.6g")
    print(f"[done] wrote {_OUT}", file=sys.stderr)

    # Summary to stdout: per-metric significant counts + clean-wins table.
    print()
    print(f"=== Phase 3 stats summary (BH-FDR q={_FDR_Q}) ===")
    print(f"Total tests: {len(pair_df)} "
          f"({pair_df['task'].nunique()} tasks x "
          f"{pair_df['competitor'].nunique()} competitors x "
          f"{len(_METRICS)} metrics)")
    print()
    for metric in _METRICS:
        sub = pair_df[pair_df["metric"] == metric]
        sig = sub[sub["p_adj_bh"] < _FDR_Q]
        print(f"  {metric:14s}: {len(sig)}/{len(sub)} significant after BH-FDR")
    print()
    wins = pair_df[pair_df["gfc_wins"]].copy()
    print(f"--- gfc clean wins (median_ratio > 1 AND p_adj_bh < {_FDR_Q}): "
          f"{len(wins)} of {len(pair_df)} tests ---")
    if not wins.empty:
        wins_view = wins[["task", "competitor", "metric", "n_pairs",
                          "median_ratio", "ci_lo_95", "ci_hi_95",
                          "cliffs_delta", "wilcoxon_p", "p_adj_bh"]]
        with pd.option_context("display.max_rows", None,
                               "display.width", 200,
                               "display.float_format", "{:.4g}".format):
            print(wins_view.to_string(index=False))
    print()
    print("--- non-wins (gfc loses or no significant difference) ---")
    losses = pair_df[~pair_df["gfc_wins"]].copy()
    if not losses.empty:
        losses_view = losses[["task", "competitor", "metric", "n_pairs",
                              "median_ratio", "ci_lo_95", "ci_hi_95",
                              "cliffs_delta", "wilcoxon_p", "p_adj_bh"]]
        with pd.option_context("display.max_rows", None,
                               "display.width", 200,
                               "display.float_format", "{:.4g}".format):
            print(losses_view.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
