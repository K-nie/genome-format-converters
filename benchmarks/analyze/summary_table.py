#!/usr/bin/env python3
# Author: Benjamin Narh-Madey
"""Headline-metrics summary table — gfc vs. best competitor per task.

Reads ``benchmarks/results/raw/all.tsv`` (the merged output of
``collect_results.py``) and writes:
  benchmarks/results/figures/summary_table.md   (markdown for paper / README)
  benchmarks/results/figures/summary_table.csv  (machine-readable)

For every task it reports:
  - gfc wall_s (mean ± sd, n)
  - best competitor wall_s (the competitor with the lowest mean wall_s)
  - speedup multiplier — median of the paired (competitor / gfc)
    ratios, with a 10000-resample bootstrap 95 % CI. >1 means gfc
    is faster, <1 means the competitor is faster (Phase 3 audit
    2026-05-03 switched both ratio columns from mean-of-ratios to
    median-with-CI for honest uncertainty reporting)
  - gfc peak_rss_mb (mean ± sd)
  - lowest-RAM competitor peak_rss_mb
  - memory ratio — median of paired (competitor / gfc) RSS ratios
    with bootstrap 95 % CI; same convention
  - Cliff's delta on the paired raw values (competitor vs gfc) for
    each ratio column, as the non-parametric effect-size complement
  - correctness flag from the gfc rep1 row (`1` / `0` / `skip` / blank)

Tasks with no competitor (only gfc rows) emit a row with `—` placeholders
in the competitor columns so the paper can show "we ran but had no
reference" honestly.

Failed rows (exit_code != 0) are excluded from speed/memory aggregation
but counted separately in a `failed_tools` cell so the table tells the
full story without the failures inflating the means.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_RAW = Path(__file__).resolve().parent.parent / "results" / "raw" / "all.tsv"
_FIG = Path(__file__).resolve().parent.parent / "results" / "figures"

_BOOT_N = 10_000
_RNG = np.random.default_rng(20260503)


def _agg_tool(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    """Mean + sd + n per (task, tool) on `metric`. Successful reps only."""
    ok = df[df["exit_code"] == 0]
    grouped = ok.groupby(["task", "tool"])[metric].agg(["mean", "std", "count"])
    return grouped.rename(columns={"mean": f"{metric}_mean",
                                   "std": f"{metric}_sd",
                                   "count": f"{metric}_n"})


def _paired_values(df_ok: pd.DataFrame, task: str, tool: str,
                   metric: str) -> np.ndarray | None:
    """Return the per-replicate values for (task, tool, metric) on
    successful reps, indexed by replicate number. None if missing."""
    sub = df_ok[(df_ok["task"] == task) & (df_ok["tool"] == tool)]
    if sub.empty:
        return None
    return sub.set_index("replicate")[metric]


def _median_ratio_ci(comp: pd.Series, gfc: pd.Series,
                     n_boot: int = _BOOT_N) -> tuple[float, float, float, int]:
    """Median of paired comp/gfc ratios with bootstrap 95 % CI.

    Pairs by replicate index intersection so missing reps don't get
    an unpaired ratio. Returns (median, ci_lo, ci_hi, n_pairs)."""
    paired = sorted(set(comp.index) & set(gfc.index))
    if len(paired) < 2:
        return float("nan"), float("nan"), float("nan"), len(paired)
    c = comp.loc[paired].to_numpy(dtype=float)
    g = gfc.loc[paired].to_numpy(dtype=float)
    ratios = c / g
    point = float(np.median(ratios))
    n = len(ratios)
    boot = np.empty(n_boot)
    for i in range(n_boot):
        idx = _RNG.integers(0, n, size=n)
        boot[i] = np.median(ratios[idx])
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return point, float(lo), float(hi), n


def _cliffs_delta(comp: pd.Series, gfc: pd.Series) -> float:
    """Cliff's delta on the paired comp - gfc differences."""
    paired = sorted(set(comp.index) & set(gfc.index))
    if not paired:
        return float("nan")
    diff = (comp.loc[paired].to_numpy(dtype=float)
            - gfc.loc[paired].to_numpy(dtype=float))
    n = len(diff)
    pos = int((diff > 0).sum())
    neg = int((diff < 0).sum())
    return (pos - neg) / n


def _fmt_ratio_ci(point: float, lo: float, hi: float) -> str:
    if not np.isfinite(point):
        return "—"
    if not (np.isfinite(lo) and np.isfinite(hi)):
        return f"{point:.2f}x"
    return f"{point:.2f}x [95% CI {lo:.2f}-{hi:.2f}]"


def _fmt_delta(d: float) -> str:
    if not np.isfinite(d):
        return "—"
    return f"{d:+.2f}"


def build_summary(df: pd.DataFrame) -> pd.DataFrame:
    """One row per task. gfc baseline + best competitor by wall + by RAM."""
    df_ok = df[df["exit_code"] == 0].copy()
    wall = _agg_tool(df, "wall_s")
    rss = _agg_tool(df, "peak_rss_mb")

    # Track failed tools for the audit cell.
    failed = (df[df["exit_code"] != 0]
              .groupby("task")["tool"].apply(lambda s: ",".join(sorted(set(s)))))

    # gfc correctness from rep1 (the row that drove check_correctness.py).
    correctness = (df[(df["tool"] == "gfc") & (df["replicate"] == 1)]
                   .set_index("task")["correct"].astype(str))

    tasks_with_data = sorted({t for (t, _) in wall.index})
    tasks_in_raw = sorted(df["task"].unique())

    rows = []
    for task in tasks_in_raw:
        if task not in tasks_with_data:
            rows.append({
                "task": task,
                "gfc_wall_s": "—",
                "best_speed_tool": "—",
                "best_speed_wall_s": "—",
                "speedup_x": "—",
                "speedup_cliffs_delta": "—",
                "gfc_rss_mb": "—",
                "best_mem_tool": "—",
                "best_mem_rss_mb": "—",
                "mem_ratio_x": "—",
                "mem_cliffs_delta": "—",
                "correct": "",
                "failed_tools": "(no successful reps)",
            })
            continue

        gfc_wall = wall.loc[(task, "gfc")] if (task, "gfc") in wall.index else None
        gfc_rss = rss.loc[(task, "gfc")] if (task, "gfc") in rss.index else None

        comp_wall = wall.loc[task].drop("gfc", errors="ignore")
        comp_rss = rss.loc[task].drop("gfc", errors="ignore")

        # Pick the fastest competitor by mean (preserves prior behaviour),
        # then bootstrap the median ratio of THAT competitor against gfc.
        if comp_wall.empty:
            best_speed_tool = "—"
            best_speed_mean = best_speed_sd = float("nan")
            speedup_str = "—"
            speedup_delta_str = "—"
        else:
            best_speed_tool = comp_wall["wall_s_mean"].idxmin()
            best_speed_mean = comp_wall.loc[best_speed_tool, "wall_s_mean"]
            best_speed_sd = comp_wall.loc[best_speed_tool, "wall_s_sd"]
            gfc_wall_vals = _paired_values(df_ok, task, "gfc", "wall_s")
            comp_wall_vals = _paired_values(df_ok, task, best_speed_tool, "wall_s")
            point, lo, hi, _ = _median_ratio_ci(comp_wall_vals, gfc_wall_vals)
            speedup_str = _fmt_ratio_ci(point, lo, hi)
            speedup_delta_str = _fmt_delta(_cliffs_delta(comp_wall_vals,
                                                          gfc_wall_vals))

        if comp_rss.empty:
            best_mem_tool = "—"
            best_mem_mean = best_mem_sd = float("nan")
            mem_ratio_str = "—"
            mem_delta_str = "—"
        else:
            best_mem_tool = comp_rss["peak_rss_mb_mean"].idxmin()
            best_mem_mean = comp_rss.loc[best_mem_tool, "peak_rss_mb_mean"]
            best_mem_sd = comp_rss.loc[best_mem_tool, "peak_rss_mb_sd"]
            gfc_rss_vals = _paired_values(df_ok, task, "gfc", "peak_rss_mb")
            comp_rss_vals = _paired_values(df_ok, task, best_mem_tool, "peak_rss_mb")
            point, lo, hi, _ = _median_ratio_ci(comp_rss_vals, gfc_rss_vals)
            mem_ratio_str = _fmt_ratio_ci(point, lo, hi)
            mem_delta_str = _fmt_delta(_cliffs_delta(comp_rss_vals, gfc_rss_vals))

        rows.append({
            "task": task,
            "gfc_wall_s": _fmt_mean_sd(gfc_wall, "wall_s"),
            "best_speed_tool": best_speed_tool,
            "best_speed_wall_s": _fmt_value_sd(best_speed_mean, best_speed_sd),
            "speedup_x": speedup_str,
            "speedup_cliffs_delta": speedup_delta_str,
            "gfc_rss_mb": _fmt_mean_sd(gfc_rss, "peak_rss_mb"),
            "best_mem_tool": best_mem_tool,
            "best_mem_rss_mb": _fmt_value_sd(best_mem_mean, best_mem_sd),
            "mem_ratio_x": mem_ratio_str,
            "mem_cliffs_delta": mem_delta_str,
            "correct": correctness.get(task, ""),
            "failed_tools": failed.get(task, ""),
        })

    return pd.DataFrame(rows)


def _fmt_mean_sd(row, metric: str) -> str:
    if row is None:
        return "—"
    mean = row[f"{metric}_mean"]
    sd = row[f"{metric}_sd"]
    n = int(row[f"{metric}_n"])
    sd_str = f"±{sd:.2f}" if pd.notna(sd) and n > 1 else ""
    return f"{mean:.2f}{sd_str} (n={n})"


def _fmt_value_sd(mean, sd) -> str:
    if mean is None or pd.isna(mean):
        return "—"
    sd_str = f"±{sd:.2f}" if pd.notna(sd) else ""
    return f"{mean:.2f}{sd_str}"


def render_markdown(table: pd.DataFrame) -> str:
    """Pipe-delimited markdown the README and the paper can paste verbatim."""
    cols = ["task", "gfc_wall_s", "best_speed_tool", "best_speed_wall_s",
            "speedup_x", "speedup_cliffs_delta",
            "gfc_rss_mb", "best_mem_tool", "best_mem_rss_mb",
            "mem_ratio_x", "mem_cliffs_delta",
            "correct", "failed_tools"]
    headers = ["Task", "gfc wall (s)", "Fastest competitor",
               "Competitor wall (s)", "Speedup (median, 95% CI)",
               "Speedup Cliff's delta", "gfc RSS (MB)",
               "Lowest-RAM competitor", "Competitor RSS (MB)",
               "Memory ratio (median, 95% CI)", "Memory Cliff's delta",
               "Correct", "Failed tools"]
    lines = ["| " + " | ".join(headers) + " |",
             "|" + "|".join("---" for _ in headers) + "|"]
    for _, row in table.iterrows():
        lines.append("| " + " | ".join(str(row[c]) for c in cols) + " |")
    return "\n".join(lines) + "\n"


def build_per_pair(df: pd.DataFrame) -> pd.DataFrame:
    """One row per (task, competitor) so wins against every individual
    competitor are visible — not just the single fastest one. Hides the
    'gfc loses to plink2 by 14x but beats AGAT by 13x on the same task'
    structure that the headline summary collapses.
    """
    df_ok = df[df["exit_code"] == 0].copy()
    wall = _agg_tool(df, "wall_s")
    rss = _agg_tool(df, "peak_rss_mb")
    rows = []
    for task in sorted({t for (t, _) in wall.index}):
        if (task, "gfc") not in wall.index:
            continue
        gfc_wall_mean = wall.loc[(task, "gfc"), "wall_s_mean"]
        gfc_rss_mean = rss.loc[(task, "gfc"), "peak_rss_mb_mean"]
        gfc_wall_vals = _paired_values(df_ok, task, "gfc", "wall_s")
        gfc_rss_vals = _paired_values(df_ok, task, "gfc", "peak_rss_mb")
        for comp_tool in sorted({tool for (t, tool) in wall.index
                                 if t == task and tool != "gfc"}):
            cw_mean = wall.loc[(task, comp_tool), "wall_s_mean"]
            cr_mean = rss.loc[(task, comp_tool), "peak_rss_mb_mean"]
            cw_vals = _paired_values(df_ok, task, comp_tool, "wall_s")
            cr_vals = _paired_values(df_ok, task, comp_tool, "peak_rss_mb")
            sp_pt, sp_lo, sp_hi, _ = _median_ratio_ci(cw_vals, gfc_wall_vals)
            mr_pt, mr_lo, mr_hi, _ = _median_ratio_ci(cr_vals, gfc_rss_vals)
            sp_delta = _cliffs_delta(cw_vals, gfc_wall_vals)
            mr_delta = _cliffs_delta(cr_vals, gfc_rss_vals)
            wins_speed = np.isfinite(sp_pt) and sp_pt > 1.0
            wins_mem = np.isfinite(mr_pt) and mr_pt > 1.0
            verdict = ("gfc wins" if wins_speed and wins_mem
                       else "speed only" if wins_speed
                       else "memory only" if wins_mem
                       else "loses")
            rows.append({
                "task": task,
                "competitor": comp_tool,
                "gfc_wall_s": f"{gfc_wall_mean:.2f}",
                "comp_wall_s": f"{cw_mean:.2f}",
                "speedup_x": _fmt_ratio_ci(sp_pt, sp_lo, sp_hi),
                "speedup_cliffs_delta": _fmt_delta(sp_delta),
                "gfc_rss_mb": f"{gfc_rss_mean:.2f}",
                "comp_rss_mb": f"{cr_mean:.2f}",
                "mem_ratio_x": _fmt_ratio_ci(mr_pt, mr_lo, mr_hi),
                "mem_cliffs_delta": _fmt_delta(mr_delta),
                "verdict": verdict,
            })
    return pd.DataFrame(rows)


def render_per_pair_markdown(table: pd.DataFrame) -> str:
    """Markdown for the per-(task, competitor) table — one row per
    pairwise comparison so the paper can read 'gfc beats AGAT 13x on
    T4' off the table without doing arithmetic."""
    cols = ["task", "competitor", "gfc_wall_s", "comp_wall_s",
            "speedup_x", "speedup_cliffs_delta",
            "gfc_rss_mb", "comp_rss_mb", "mem_ratio_x",
            "mem_cliffs_delta", "verdict"]
    headers = ["Task", "Competitor", "gfc wall (s)", "Comp wall (s)",
               "Speedup (median, 95% CI; >1 = gfc faster)",
               "Speedup Cliff's delta", "gfc RSS (MB)",
               "Comp RSS (MB)", "Mem ratio (median, 95% CI; >1 = gfc smaller)",
               "Memory Cliff's delta", "Verdict"]
    lines = ["| " + " | ".join(headers) + " |",
             "|" + "|".join("---" for _ in headers) + "|"]
    for _, row in table.iterrows():
        lines.append("| " + " | ".join(str(row[c]) for c in cols) + " |")
    return "\n".join(lines) + "\n"


def main() -> int:
    if not _RAW.exists():
        print(f"[error] {_RAW} missing. Run run_bench.sh + collect_results.py first.",
              file=sys.stderr)
        return 1

    df = pd.read_csv(_RAW, sep="\t")
    if df.empty:
        print(f"[error] {_RAW} is empty.", file=sys.stderr)
        return 1

    # Normalise the optional `correct` and `notes` columns: missing -> empty
    # string so the markdown / CSV doesn't render NaN cells.
    for col in ("correct", "notes"):
        if col in df.columns:
            df[col] = df[col].fillna("").astype(str)

    table = build_summary(df)
    pair_table = build_per_pair(df)
    _FIG.mkdir(parents=True, exist_ok=True)

    md_path = _FIG / "summary_table.md"
    csv_path = _FIG / "summary_table.csv"
    md_path.write_text(render_markdown(table))
    table.to_csv(csv_path, index=False)

    pair_md = _FIG / "summary_per_pair.md"
    pair_csv = _FIG / "summary_per_pair.csv"
    pair_md.write_text(render_per_pair_markdown(pair_table))
    pair_table.to_csv(pair_csv, index=False)

    print(f"[done] wrote {md_path}", file=sys.stderr)
    print(f"[done] wrote {csv_path}", file=sys.stderr)
    print(f"[done] wrote {pair_md}", file=sys.stderr)
    print(f"[done] wrote {pair_csv}", file=sys.stderr)
    print("\n--- per-pair summary (one row per competitor; >1x = gfc wins) ---",
          file=sys.stderr)
    print(render_per_pair_markdown(pair_table))
    print("\n--- headline summary (one row per task, fastest competitor only) ---",
          file=sys.stderr)
    print(render_markdown(table))
    return 0


if __name__ == "__main__":
    sys.exit(main())
