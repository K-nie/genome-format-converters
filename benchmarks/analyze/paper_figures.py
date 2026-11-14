#!/usr/bin/env python3
"""Build every paper figure and table from the paper-rerun results.

Input: a directory holding the per-phase tarballs written by
``condor/run_paper_rerun.sh`` (``main.<stamp>.tar.gz``, ``scaling.<stamp>.tar.gz``,
...) and/or already-extracted phase directories (``main/``, ``scaling/``, ...).
The newest tarball per phase is extracted to ``<results>/extracted/``. Phases
that are missing are skipped with a note, so the script can run as results
arrive.

Output (``--out``, default ``<results>/paper``)::

    figures/fig_runtime_memory.{pdf,png}   per-task wall time and peak RSS
    figures/fig_scaling.{pdf,png}          wall time vs input size
    figures/fig_correctness.{pdf,png}      four-state verdict matrix
    figures/fig_annotations.{pdf,png}      T4/T5 on real gene annotations
    tables/ratios_<phase>.tsv              competitor/gfc median ratios, 95% CI
    tables/per_cell_<phase>.tsv            median, CI, CV, n per (task, tool)
    tables/verdicts.tsv                    every four-state verdict + reason
    tables/t8_reproducibility.tsv          T8 checks
    tables/import_times.tsv                start-up time / RSS per subcommand
    tables/install_footprint.tsv           packages and disk per environment
    tables/cold_vs_warm.tsv                cold/warm median wall ratio per cell
    summary.md                             all tables in Markdown, for drafting

Statistics: medians with percentile bootstrap 95% CIs (10,000 resamples,
seeded). Ratios are median(competitor)/median(gfc) with the two samples
resampled independently, since replicate index is not a natural pairing.
No p-values. Only replicates with exit_code 0 are used; failed runs are
counted in the per-cell table.

Usage::

    python benchmarks/analyze/paper_figures.py --results ~/rerun_results
"""
from __future__ import annotations

import argparse
import re
import sys
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import ticker  # noqa: E402

PHASES = ("main", "annotations", "scaling", "cold", "import", "footprint")
SEED = 20261112
N_BOOT = 10_000

# Print palette: gfc is the one highlighted series; comparators are neutral
# ink, told apart by their tick labels / direct labels and marker shape.
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e4e3df"
GFC = "#2a78d6"
OTHER = "#7a7974"
MARKERS = ["o", "s", "^", "D", "v", "P", "X"]
STATE_COLOR = {"byte": "#1baf7a", "structural": "#2a78d6",
               "limitation": "#b5b4ae", "disagreement": "#e34948",
               "reference": "#ffffff"}
STATE_LABEL = {"byte": "byte-identical", "structural": "structural",
               "limitation": "limitation", "disagreement": "disagreement",
               "reference": "reference"}

TASK_TITLE = {
    "T1": "T1 VCF → EIGENSTRAT", "T2": "T2 VCF → PLINK",
    "T3": "T3 FASTA+GFF3 → GenBank", "T4": "T4 GFF3 → GTF",
    "T5": "T5 GFF3 → BED12", "T6": "T6 HMMER tblout → TSV",
    "T7": "T7 orthogroups → FASTA", "T8": "T8 VCF → pseudohaploid",
}
# Two-line panel titles for the narrow small multiples.
TASK_SHORT = {k: v.replace(" ", "\n", 1) for k, v in TASK_TITLE.items()}

plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.edgecolor": INK_2, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK_2, "ytick.color": INK_2, "axes.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 600,
})


# ------------------------------------------------------------ loading

def extract(results: Path) -> dict:
    """Return {phase: directory}, extracting the newest tarball per phase."""
    found = {}
    ex_root = results / "extracted"
    for phase in PHASES:
        tars = sorted(results.glob(f"{phase}.*.tar.gz"))
        if tars:
            dest = ex_root / tars[-1].name.replace(".tar.gz", "")
            if not dest.exists():
                dest.mkdir(parents=True)
                with tarfile.open(tars[-1]) as tf:
                    tf.extractall(dest, filter="data") if sys.version_info >= (3, 12) \
                        else tf.extractall(dest)  # noqa: S202 - our own archives
            sub = dest / phase
            found[phase] = sub if sub.is_dir() else dest
        elif (results / phase).is_dir():
            found[phase] = results / phase
    return found


def read_runs(paths) -> pd.DataFrame:
    frames = []
    for p in paths:
        df = pd.read_csv(p, sep="\t", dtype={"correct": str, "notes": str})
        df["source"] = p.name
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    df = pd.concat(frames, ignore_index=True)
    df = df[df["task"].astype(str).str.match(r"T\d")]
    for c in ("wall_s", "peak_rss_mb", "exit_code"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def task_tsvs(d: Path):
    """Per-task timing TSVs in a phase dir (no sidecars, no all.tsv)."""
    return sorted(p for p in d.glob("T*.tsv")
                  if not re.search(r"\.(verdicts|reproducibility)\.tsv$", p.name))


# ------------------------------------------------------------ statistics

_rng = np.random.default_rng(SEED)


def boot_median(x: np.ndarray):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return np.nan, np.nan, np.nan
    if len(x) == 1:
        return x[0], np.nan, np.nan
    idx = _rng.integers(0, len(x), size=(N_BOOT, len(x)))
    meds = np.median(x[idx], axis=1)
    return float(np.median(x)), float(np.percentile(meds, 2.5)), float(np.percentile(meds, 97.5))


def boot_ratio(num: np.ndarray, den: np.ndarray):
    num = np.asarray(num, float)[np.isfinite(num)]
    den = np.asarray(den, float)[np.isfinite(den)]
    if len(num) == 0 or len(den) == 0:
        return np.nan, np.nan, np.nan
    point = float(np.median(num) / np.median(den))
    if len(num) < 2 or len(den) < 2:
        return point, np.nan, np.nan
    a = np.median(num[_rng.integers(0, len(num), size=(N_BOOT, len(num)))], axis=1)
    b = np.median(den[_rng.integers(0, len(den), size=(N_BOOT, len(den)))], axis=1)
    r = a / b
    return point, float(np.percentile(r, 2.5)), float(np.percentile(r, 97.5))


def per_cell(df: pd.DataFrame, by=("task", "tool")) -> pd.DataFrame:
    rows = []
    for key, g in df.groupby(list(by), sort=False):
        ok = g[g["exit_code"] == 0]
        row = dict(zip(by, key if isinstance(key, tuple) else (key,)))
        row.update(n_ok=len(ok), n_failed=len(g) - len(ok))
        for metric in ("wall_s", "peak_rss_mb"):
            m, lo, hi = boot_median(ok[metric].values)
            x = ok[metric].values
            row[f"{metric}_median"] = m
            row[f"{metric}_ci_lo"] = lo
            row[f"{metric}_ci_hi"] = hi
            row[f"{metric}_cv"] = float(np.std(x, ddof=1) / np.mean(x)) if len(x) > 1 else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def ratios(df: pd.DataFrame, by=("task",)) -> pd.DataFrame:
    """median(competitor)/median(gfc): >1 means gfc is faster / leaner."""
    rows = []
    ok = df[df["exit_code"] == 0]
    for key, g in ok.groupby(list(by), sort=False):
        base = g[g["tool"] == "gfc"]
        if base.empty:
            continue
        for tool, c in g[g["tool"] != "gfc"].groupby("tool", sort=False):
            row = dict(zip(by, key if isinstance(key, tuple) else (key,)))
            row["competitor"] = tool
            for metric in ("wall_s", "peak_rss_mb"):
                r, lo, hi = boot_ratio(c[metric].values, base[metric].values)
                row[f"{metric}_ratio"] = r
                row[f"{metric}_ci_lo"] = lo
                row[f"{metric}_ci_hi"] = hi
            rows.append(row)
    return pd.DataFrame(rows)


# ------------------------------------------------------------ figures

def _save(fig, out: Path, name: str) -> None:
    out.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(out / f"{name}.{ext}", bbox_inches="tight",
                    dpi=600 if ext == "png" else None)
    plt.close(fig)
    print(f"[fig ] {out / name}.pdf", file=sys.stderr)


def _plain(v, _pos=None) -> str:
    return f"{v:g}"


def _style(ax, values=None) -> None:
    """Log y axis with plain-number labels, spanning at least one decade.

    The one-decade floor keeps panels comparable: two tools within a few
    percent of each other look close, instead of being stretched apart.
    """
    ax.set_yscale("log")
    if values is not None:
        v = np.asarray([x for x in values if np.isfinite(x) and x > 0])
        if len(v):
            lo, hi = np.log10(v.min()), np.log10(v.max())
            if hi - lo < 1:
                mid = (hi + lo) / 2
                lo, hi = mid - 0.5, mid + 0.5
            pad = 0.08 * (hi - lo)
            ax.set_ylim(10 ** (lo - pad), 10 ** (hi + pad))
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(_plain))
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())
    ax.grid(axis="y", color=GRID, linewidth=0.5)
    ax.set_axisbelow(True)


def fig_runtime_memory(df: pd.DataFrame, cells: pd.DataFrame, out: Path,
                       name: str, tasks=None) -> None:
    """Small multiples: one column per task, rows = wall time, peak RSS.

    Each tool is a point at its median with a 95% bootstrap CI whisker and
    the individual replicates as faint dots; gfc is the blue series.
    """
    tasks = tasks or [t for t in TASK_TITLE if t in set(cells["task"])]
    if not tasks:
        return
    fig, axes = plt.subplots(2, len(tasks), figsize=(1.55 * len(tasks) + 0.6, 4.2),
                             squeeze=False)
    ok = df[df["exit_code"] == 0]
    for j, task in enumerate(tasks):
        c = cells[cells["task"] == task].copy()
        tools = ["gfc"] + sorted(t for t in c["tool"] if t != "gfc")
        tools = [t for t in tools if t in set(c["tool"])]
        for i, (metric, label) in enumerate((("wall_s", "Wall time (s)"),
                                             ("peak_rss_mb", "Peak RSS (MB)"))):
            ax = axes[i][j]
            for x, tool in enumerate(tools):
                col = GFC if tool == "gfc" else OTHER
                reps = ok[(ok["task"] == task) & (ok["tool"] == tool)][metric]
                ax.scatter(np.full(len(reps), x) + _rng.uniform(-0.12, 0.12, len(reps)),
                           reps, s=6, color=col, alpha=0.35, linewidths=0, zorder=2)
                r = c[c["tool"] == tool].iloc[0]
                m, lo, hi = r[f"{metric}_median"], r[f"{metric}_ci_lo"], r[f"{metric}_ci_hi"]
                if np.isfinite(lo):
                    ax.plot([x, x], [lo, hi], color=col, linewidth=1.6,
                            solid_capstyle="round", zorder=3)
                ax.scatter([x], [m], s=26, color=col, edgecolor="white",
                           linewidth=0.8, zorder=4)
            vals = list(ok[(ok["task"] == task)][metric]) + \
                list(c[f"{metric}_ci_lo"]) + list(c[f"{metric}_ci_hi"])
            _style(ax, vals)
            ax.set_xticks(range(len(tools)))
            ax.set_xticklabels(tools if i == 1 else [], rotation=45, ha="right")
            ax.set_xlim(-0.6, len(tools) - 0.4)
            if j == 0:
                ax.set_ylabel(label)
            if i == 0:
                ax.set_title(TASK_SHORT.get(task, task), loc="left")
    fig.tight_layout(w_pad=0.6, h_pad=0.8)
    _save(fig, out, name)


def _size_value(label: str) -> float:
    m = re.match(r"sites(\d+)", label)
    if m:
        return float(m.group(1))
    m = re.search(r"_x(\d+)$", label)
    return float(m.group(1)) if m else np.nan


def fig_scaling(scal: pd.DataFrame, out: Path) -> pd.DataFrame:
    cells = per_cell(scal, by=("task", "size", "tool"))
    cells["size_value"] = cells["size"].map(_size_value)
    tasks = [t for t in ("T1", "T2", "T4", "T5") if t in set(cells["task"])]
    if not tasks:
        return cells
    fig, axes = plt.subplots(1, len(tasks), figsize=(2.0 * len(tasks) + 0.4, 2.5),
                             squeeze=False)
    for j, task in enumerate(tasks):
        ax = axes[0][j]
        c = cells[cells["task"] == task]
        others = sorted(t for t in set(c["tool"]) if t != "gfc")
        ends = []
        for k, tool in enumerate(["gfc"] + others):
            t = c[c["tool"] == tool].sort_values("size_value")
            if t.empty:
                continue
            col = GFC if tool == "gfc" else OTHER
            mk = "o" if tool == "gfc" else MARKERS[(k % (len(MARKERS) - 1)) + 1]
            ax.plot(t["size_value"], t["wall_s_median"], color=col,
                    linewidth=2 if tool == "gfc" else 1.2, marker=mk,
                    markersize=4.5, markeredgecolor="white", markeredgewidth=0.6)
            lo, hi = t["wall_s_ci_lo"].values, t["wall_s_ci_hi"].values
            good = np.isfinite(lo)
            if good.any():
                ax.vlines(t["size_value"].values[good], lo[good], hi[good],
                          color=col, linewidth=1)
            last = t.iloc[-1]
            ends.append([tool, last["size_value"], np.log10(last["wall_s_median"])])
        # Direct labels at line ends, nudged apart in log space.
        ends.sort(key=lambda e: e[2])
        for a in range(1, len(ends)):
            if ends[a][2] - ends[a - 1][2] < 0.12:
                ends[a][2] = ends[a - 1][2] + 0.12
        for tool, x, ly in ends:
            ax.annotate(tool, (x, 10 ** ly), xytext=(5, 0), textcoords="offset points",
                        va="center", fontsize=6.5, color=INK if tool == "gfc" else INK_2)
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(
            lambda v, _p: f"{v:,.0f}" if v >= 1 else f"{v:g}"))
        _style(ax, list(c["wall_s_median"]) + list(c["wall_s_ci_lo"]) + list(c["wall_s_ci_hi"]))
        ax.set_title(TASK_TITLE[task], loc="left")
        ax.set_xlabel("Sites" if task in ("T1", "T2") else "Annotation copies (x)")
        if j == 0:
            ax.set_ylabel("Wall time (s), median and 95% CI")
        ax.margins(x=0.25)
    fig.tight_layout(w_pad=1.0)
    _save(fig, out, "fig_scaling")
    return cells


def fig_correctness(verd: pd.DataFrame, out: Path) -> None:
    if verd.empty:
        return
    verd = verd.copy()
    rows = list(dict.fromkeys(verd["row"]))
    tools = ["gfc"] + sorted(t for t in set(verd["tool"]) if t != "gfc")
    fig, ax = plt.subplots(figsize=(2.4 + 1.15 * len(tools), 0.36 * len(rows) + 0.9))
    for i, row in enumerate(rows):
        for j, tool in enumerate(tools):
            v = verd[(verd["row"] == row) & (verd["tool"] == tool)]
            if v.empty:
                continue
            state = v.iloc[0]["state"]
            ax.add_patch(plt.Rectangle((j + 0.04, i + 0.06), 0.92, 0.88,
                                       facecolor=STATE_COLOR.get(state, "#ffffff"),
                                       edgecolor=INK_2 if state == "reference" else "none",
                                       linewidth=0.6, linestyle=(0, (2, 2))))
            conv = str(v.iloc[0].get("conventions", "") or "")
            txt = STATE_LABEL.get(state, state) + (f"\n{conv}" if conv else "")
            ax.text(j + 0.5, i + 0.5, txt, ha="center", va="center", fontsize=6,
                    color="white" if state in ("structural", "disagreement", "byte") else INK)
    ax.set_xlim(0, len(tools))
    ax.set_ylim(len(rows), 0)
    ax.set_xticks(np.arange(len(tools)) + 0.5)
    ax.set_xticklabels(tools, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(rows)) + 0.5)
    ax.set_yticklabels(rows)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    shown = [s for s in STATE_COLOR if s != "reference"]
    handles = [plt.Rectangle((0, 0), 1, 1, color=STATE_COLOR[s]) for s in shown]
    ax.legend(handles, [STATE_LABEL[s] for s in shown], ncol=4,
              loc="lower left", bbox_to_anchor=(0, 1.01), frameon=False,
              handlelength=1, columnspacing=1)
    fig.tight_layout()
    _save(fig, out, "fig_correctness")


# ------------------------------------------------------------ verdicts

_CONV = re.compile(r"\b(chrom-label|allele-order|attribute-format|gene-rows|exon-fusion|"
                   r"itemRgb|row-order|score|number-format|record-order|"
                   r"stochastic-draw|rendering|reference-drops-colocated|indels-excluded|"
                   r"site-order|monomorphic-allele-code|id-prefix|thick-span-noncoding|name-attribute|reference-omits-children)\b")


def load_verdicts(phase_dirs: dict) -> pd.DataFrame:
    frames = []
    for phase in ("main", "annotations"):
        d = phase_dirs.get(phase)
        if not d:
            continue
        for p in sorted(d.glob("T*.verdicts.tsv")):
            v = pd.read_csv(p, sep="\t", dtype=str).fillna("")
            label = p.name.split(".")[0]
            m = re.match(r"(T\d)_[^_].*?(?:__(.+))?$", label)
            task = m.group(1) if m else label[:2]
            data = m.group(2) if m and m.group(2) else ("Y1000+ slice" if task in ("T4", "T5") and phase == "main" else "")
            v["row"] = TASK_TITLE.get(task, task) + (f" ({data})" if data else "")
            v["phase"] = phase
            if "conventions" not in v.columns:  # sidecars from older runs
                v["conventions"] = v["reason"].map(
                    lambda r: ", ".join(dict.fromkeys(_CONV.findall(r))))
            frames.append(v)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


# ------------------------------------------------------------ main

def _fmt_ratio(r, lo, hi) -> str:
    if not np.isfinite(r):
        return ""
    ci = f" [{lo:.2f}, {hi:.2f}]" if np.isfinite(lo) else ""
    return f"{r:.2f}{ci}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--results", required=True, type=Path)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()
    out = args.out or args.results / "paper"
    figs, tabs = out / "figures", out / "tables"
    tabs.mkdir(parents=True, exist_ok=True)
    phases = extract(args.results)
    md = ["# Paper rerun: tables", ""]
    for p in PHASES:
        md.append(f"- {p}: {'found at ' + str(phases[p]) if p in phases else 'MISSING'}")
    md.append("")

    def write(df: pd.DataFrame, name: str, title: str) -> None:
        df.to_csv(tabs / f"{name}.tsv", sep="\t", index=False, float_format="%.4g")
        md.extend([f"## {title}", "", df.to_markdown(index=False, floatfmt=".3g")
                   if hasattr(df, "to_markdown") else df.to_string(index=False), ""])

    main_df = pd.DataFrame()
    if "main" in phases:
        main_df = read_runs(task_tsvs(phases["main"]))
        cells = per_cell(main_df)
        write(cells, "per_cell_main", "Main panel: per-cell medians (95% CI)")
        r = ratios(main_df)
        write(r, "ratios_main", "Main panel: competitor/gfc median ratios (>1 = gfc faster or leaner)")
        fig_runtime_memory(main_df, cells, figs, "fig_runtime_memory")
        rep = phases["main"] / "T8_pseudohaploid.reproducibility.tsv"
        if rep.is_file():
            write(pd.read_csv(rep, sep="\t"), "t8_reproducibility", "T8 reproducibility checks")

    if "annotations" in phases:
        ann = read_runs(task_tsvs(phases["annotations"]))
        if not ann.empty:
            ann["dataset"] = ann["source"].str.extract(r"__(.+)\.tsv$")[0]
            cells = per_cell(ann, by=("task", "dataset", "tool"))
            write(cells, "per_cell_annotations", "Real gene annotations: per-cell medians")
            write(ratios(ann, by=("task", "dataset")), "ratios_annotations",
                  "Real gene annotations: competitor/gfc ratios")
            for ds in sorted(ann["dataset"].dropna().unique()):
                sub = ann[ann["dataset"] == ds]
                fig_runtime_memory(sub, per_cell(sub), figs, f"fig_annotations_{ds}")

    verd = load_verdicts(phases)
    if not verd.empty:
        write(verd[["row", "tool", "reference", "state", "conventions", "reason"]],
              "verdicts", "Four-state verdicts")
        fig_correctness(verd, figs)

    if "scaling" in phases:
        scal = read_runs(task_tsvs(phases["scaling"]))
        if not scal.empty:
            scal["size"] = scal["source"].str.extract(r"__(.+)\.tsv$")[0]
            cells = fig_scaling(scal, figs)
            write(cells.drop(columns=["size_value"]), "per_cell_scaling", "Scaling: per-cell medians")
            write(ratios(scal, by=("task", "size")), "ratios_scaling", "Scaling: competitor/gfc ratios")

    if "cold" in phases and not main_df.empty:
        cold = read_runs(task_tsvs(phases["cold"]))
        if not cold.empty:
            rows = []
            okw, okc = main_df[main_df["exit_code"] == 0], cold[cold["exit_code"] == 0]
            for (task, tool), g in okc.groupby(["task", "tool"]):
                w = okw[(okw["task"] == task) & (okw["tool"] == tool)]
                rr, lo, hi = boot_ratio(g["wall_s"].values, w["wall_s"].values)
                rows.append(dict(task=task, tool=tool, cold_over_warm_wall=rr,
                                 ci_lo=lo, ci_hi=hi, n_cold=len(g), n_warm=len(w)))
            write(pd.DataFrame(rows), "cold_vs_warm", "Cold vs warm cache (median wall ratio)")
            write(ratios(cold), "ratios_cold", "Cold cache: competitor/gfc ratios")

    if "import" in phases:
        p = phases["import"] / "import_times.tsv"
        if p.is_file():
            imp = pd.read_csv(p, sep="\t")
            t = (imp.groupby(["mode", "target"])
                 .agg(wall_s_median=("wall_s", "median"), peak_rss_mb_median=("peak_rss_mb", "median"),
                      trials=("trial", "count"), heavy_modules=("modules", "first"))
                 .reset_index().sort_values("peak_rss_mb_median"))
            write(t, "import_times", "Start-up cost per subcommand")

    if "footprint" in phases:
        p = phases["footprint"] / "install_footprint.tsv"
        if p.is_file():
            write(pd.read_csv(p, sep="\t"), "install_footprint", "Install footprint")

    (out / "summary.md").write_text("\n".join(md))
    print(f"[done] {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
