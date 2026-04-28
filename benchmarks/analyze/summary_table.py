#!/usr/bin/env python3
"""Headline-metrics summary table — gfc vs. best competitor per task.

Reads ``benchmarks/results/raw/all.tsv`` (the merged output of
``collect_results.py``) and writes:
  benchmarks/results/figures/summary_table.md   (markdown for paper / README)
  benchmarks/results/figures/summary_table.csv  (machine-readable)

For every task it reports:
  - gfc wall_s (mean ± sd, n)
  - best competitor wall_s (the competitor with the lowest mean wall_s)
  - speedup multiplier (best_competitor_mean / gfc_mean) — >1 means gfc
    is faster, <1 means the competitor is faster
  - gfc peak_rss_mb (mean ± sd)
  - lowest-RAM competitor peak_rss_mb
  - memory ratio (lowest_competitor_rss_mean / gfc_rss_mean) — >1 means
    gfc is more memory-efficient
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

import pandas as pd

_RAW = Path(__file__).resolve().parent.parent / "results" / "raw" / "all.tsv"
_FIG = Path(__file__).resolve().parent.parent / "results" / "figures"


def _agg_tool(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    """Mean + sd + n per (task, tool) on `metric`. Successful reps only."""
    ok = df[df["exit_code"] == 0]
    grouped = ok.groupby(["task", "tool"])[metric].agg(["mean", "std", "count"])
    return grouped.rename(columns={"mean": f"{metric}_mean",
                                   "std": f"{metric}_sd",
                                   "count": f"{metric}_n"})


def build_summary(df: pd.DataFrame) -> pd.DataFrame:
    """One row per task. gfc baseline + best competitor by wall + by RAM."""
    wall = _agg_tool(df, "wall_s")
    rss = _agg_tool(df, "peak_rss_mb")

    # Track failed tools for the audit cell.
    failed = (df[df["exit_code"] != 0]
              .groupby("task")["tool"].apply(lambda s: ",".join(sorted(set(s)))))

    # gfc correctness from rep1 (the row that drove check_correctness.py).
    correctness = (df[(df["tool"] == "gfc") & (df["replicate"] == 1)]
                   .set_index("task")["correct"].astype(str))

    # Tasks present in successful rows (after the exit_code filter inside
    # `_agg_tool`). If the raw TSV had a task whose every row was failed
    # OR malformed (e.g. the legacy convertf row with embedded newline
    # that pandas split into multi-row garbage), it'll be missing from
    # the aggregate index and pivoting on the raw `df["task"].unique()`
    # would KeyError on `wall.loc[task]`. Iterate the aggregated index
    # instead — emit a per-task "no data" row for tasks that are in
    # df but not in wall.
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
                "gfc_rss_mb": "—",
                "best_mem_tool": "—",
                "best_mem_rss_mb": "—",
                "mem_ratio_x": "—",
                "correct": "",
                "failed_tools": "(no successful reps)",
            })
            continue

        gfc_wall = wall.loc[(task, "gfc")] if (task, "gfc") in wall.index else None
        gfc_rss = rss.loc[(task, "gfc")] if (task, "gfc") in rss.index else None

        # Competitor candidates = every successful tool other than gfc.
        # `wall.loc[task]` raises KeyError if `task` isn't in the
        # aggregate's level-0; we already guarded above so this is safe.
        comp_wall = wall.loc[task].drop("gfc", errors="ignore")
        comp_rss = rss.loc[task].drop("gfc", errors="ignore")

        if comp_wall.empty:
            best_speed_tool = best_speed_mean = best_speed_sd = speedup = None
        else:
            best_speed_tool = comp_wall["wall_s_mean"].idxmin()
            best_speed_mean = comp_wall.loc[best_speed_tool, "wall_s_mean"]
            best_speed_sd = comp_wall.loc[best_speed_tool, "wall_s_sd"]
            speedup = (best_speed_mean / gfc_wall["wall_s_mean"]
                       if gfc_wall is not None else None)

        if comp_rss.empty:
            best_mem_tool = best_mem_mean = best_mem_sd = mem_ratio = None
        else:
            best_mem_tool = comp_rss["peak_rss_mb_mean"].idxmin()
            best_mem_mean = comp_rss.loc[best_mem_tool, "peak_rss_mb_mean"]
            best_mem_sd = comp_rss.loc[best_mem_tool, "peak_rss_mb_sd"]
            mem_ratio = (best_mem_mean / gfc_rss["peak_rss_mb_mean"]
                         if gfc_rss is not None else None)

        rows.append({
            "task": task,
            "gfc_wall_s": _fmt_mean_sd(gfc_wall, "wall_s"),
            "best_speed_tool": best_speed_tool or "—",
            "best_speed_wall_s": _fmt_value_sd(best_speed_mean, best_speed_sd),
            "speedup_x": f"{speedup:.2f}x" if speedup else "—",
            "gfc_rss_mb": _fmt_mean_sd(gfc_rss, "peak_rss_mb"),
            "best_mem_tool": best_mem_tool or "—",
            "best_mem_rss_mb": _fmt_value_sd(best_mem_mean, best_mem_sd),
            "mem_ratio_x": f"{mem_ratio:.2f}x" if mem_ratio else "—",
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
            "speedup_x", "gfc_rss_mb", "best_mem_tool", "best_mem_rss_mb",
            "mem_ratio_x", "correct", "failed_tools"]
    headers = ["Task", "gfc wall (s)", "Fastest competitor",
               "Competitor wall (s)", "Speedup", "gfc RSS (MB)",
               "Lowest-RAM competitor", "Competitor RSS (MB)",
               "Memory ratio", "Correct", "Failed tools"]
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
    _FIG.mkdir(parents=True, exist_ok=True)

    md_path = _FIG / "summary_table.md"
    csv_path = _FIG / "summary_table.csv"
    md_path.write_text(render_markdown(table))
    table.to_csv(csv_path, index=False)

    print(f"[done] wrote {md_path}", file=sys.stderr)
    print(f"[done] wrote {csv_path}", file=sys.stderr)
    print("\n--- summary ---", file=sys.stderr)
    print(render_markdown(table))
    return 0


if __name__ == "__main__":
    sys.exit(main())
