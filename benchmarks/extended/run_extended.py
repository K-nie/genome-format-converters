#!/usr/bin/env python3
"""Run the extended benchmark (E01-E24): wall time, peak memory and a
four-state verdict for every subcommand not covered by T1-T8.

Usage:
    python benchmarks/extended/run_extended.py [--tasks E01 E05 ...] \
        [--reps 10] [--data benchmarks/data/extended] \
        [--out benchmarks/results/extended] [--tmp /scratch/gfc_ext]

Writes, under --out:
    raw/<task>.tsv          one row per tool and replicate (bench_one.py schema)
    verdicts_extended.tsv   one row per task and tool
    tasks_extended.tsv      the task table (subcommand, inputs, tools, reference)

Replicates run back to back, as in the main benchmark. Only replicate 1's
outputs are kept for the verdicts; later replicates are deleted as soon as
they finish so that scratch space stays small.
"""
from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
BENCH = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(BENCH))

from tasks import TASKS, BY_ID, Task, Tool  # noqa: E402
import compare_ext  # noqa: E402

BIOCONVERT = os.environ.get("GFC_BENCH_BIOCONVERT", "")


def fmt(template: str, IN: Path, OUT: Path) -> str:
    """Fill a command template (plain replacement, so awk/printf braces in
    the template are left alone)."""
    return (template.replace("{BIOCONVERT}", BIOCONVERT or "bioconvert-not-set")
            .replace("{REFS}", shlex.quote(str(HERE / "refs")))
            .replace("{IN}", shlex.quote(str(IN)))
            .replace("{OUT}", shlex.quote(str(OUT))))


def available(tool: Tool) -> str:
    """Empty string when every executable the tool needs is present,
    otherwise the name of the first missing one."""
    for exe in tool.needs:
        if exe == "{BIOCONVERT}":
            if not (BIOCONVERT and os.access(BIOCONVERT, os.X_OK)):
                return "bioconvert (set GFC_BENCH_BIOCONVERT)"
        elif shutil.which(exe) is None:
            return exe
    return ""


def version(cmd: str) -> str:
    if not cmd:
        return "unknown"
    cmd = cmd.replace("{BIOCONVERT}", BIOCONVERT or "false")
    try:
        out = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=120)
        v = (out.stdout or out.stderr).strip().splitlines()
        return " ".join(v[-1].split())[:80] if v else "unknown"
    except Exception:
        return "unknown"


def gfc_version() -> str:
    out = subprocess.run(["gfc", "--version"], capture_output=True, text=True)
    return (out.stdout.strip().split() or ["", "unknown"])[-1]


def run_task(task: Task, args, raw_dir: Path, verdict_rows: list, task_rows: list) -> None:
    in_dir = Path(args.data) / task.id
    if not (in_dir / ".done").exists():
        print(f"[skip] {task.id} {task.subcommand}: no inputs at {in_dir}", file=sys.stderr)
        task_rows.append([task.id, task.subcommand, task.conversion, task.domain, task.input_note,
                          "", task.reference, "inputs missing"])
        return
    work = Path(args.tmp) / task.id
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    tools = [Tool("gfc", task.gfc, "")] + list(task.tools)
    raw = raw_dir / f"{task.id}.tsv"
    with open(raw, "w") as fh:
        fh.write(subprocess.run([sys.executable, str(BENCH / "bench_one.py"), "--header"],
                                capture_output=True, text=True).stdout)
    ran, skipped = [], []
    for tool in tools:
        missing = available(tool)
        if missing:
            skipped.append(f"{tool.name} ({missing} not found)")
            print(f"[skip] {task.id} {tool.name}: {missing} not found", file=sys.stderr)
            continue
        ver = gfc_version() if tool.name == "gfc" else version(tool.version)
        for rep in range(1, args.reps + 1):
            out = work / f"{tool.name}_rep{rep}"
            out.mkdir(parents=True, exist_ok=True)
            cmd = fmt(tool.cmd, IN=in_dir, OUT=out)
            row = subprocess.run(
                [sys.executable, str(BENCH / "bench_one.py"), "--task", task.id,
                 "--tool", tool.name, "--version", ver, "--replicate", str(rep), "--cmd", cmd],
                capture_output=True, text=True)
            sys.stderr.write(row.stderr)
            with open(raw, "a") as fh:
                fh.write(row.stdout)
            if rep > 1:
                shutil.rmtree(out, ignore_errors=True)
        ran.append(tool.name)
        print(f"[done] {task.id} {tool.name} x{args.reps}", file=sys.stderr)
    task_rows.append([task.id, task.subcommand, task.conversion, task.domain, task.input_note,
                      ", ".join(t for t in ran if t != "gfc"), task.reference,
                      "; ".join(skipped)])
    # ---- verdicts on replicate 1 -------------------------------------------
    cmp = compare_ext.COMPARATORS[task.comparator]
    ref_dir = in_dir if task.reference == "input" else work / f"{task.reference}_rep1"
    if task.reference != "input" and task.reference not in ran:
        verdict_rows.append([task.id, task.subcommand, "gfc", task.reference, "limitation", "",
                             f"reference tool {task.reference} did not run"])
        return
    for name in ran:
        if name == task.reference:
            verdict_rows.append([task.id, task.subcommand, name, task.reference, "reference", "",
                                 "reference tool"])
            continue
        if task.reference == "input" and name != "gfc" and not task.options.get("compare_others"):
            verdict_rows.append([task.id, task.subcommand, name, "input", "not compared", "",
                                 "partial route with different content; timed only"])
            continue
        t0 = time.time()
        try:
            opts = {k: v for k, v in task.options.items() if k != "compare_others"}
            v = cmp(work / f"{name}_rep1", ref_dir, input_dir=in_dir, tool=name, **opts)
        except Exception as exc:  # a comparator crash must not stop the run
            v = compare_ext.Verdict("limitation", f"comparator error: {type(exc).__name__}: {exc}")
        verdict_rows.append([task.id, task.subcommand, name, task.reference, v.state,
                             ", ".join(v.conventions), v.reason])
        print(f"[verd] {task.id} {name} vs {task.reference}: {v.state} "
              f"{','.join(v.conventions)} ({time.time() - t0:.0f} s) {v.reason[:200]}", file=sys.stderr)
    if not args.keep:
        shutil.rmtree(work, ignore_errors=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tasks", nargs="*", default=[t.id for t in TASKS])
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--data", default=str(BENCH / "data" / "extended"))
    ap.add_argument("--out", default=str(BENCH / "results" / "extended"))
    ap.add_argument("--tmp", default=os.environ.get("GFC_BENCH_TMP", "/tmp/gfc_ext"))
    ap.add_argument("--keep", action="store_true", help="keep replicate-1 outputs")
    args = ap.parse_args()
    out = Path(args.out)
    raw_dir = out / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    verdict_rows, task_rows = [], []
    for tid in args.tasks:
        run_task(BY_ID[tid], args, raw_dir, verdict_rows, task_rows)
        # Rewrite the summaries after every task so a killed job keeps them.
        with open(out / "verdicts_extended.tsv", "w") as fh:
            fh.write("task\tsubcommand\ttool\treference\tstate\tconventions\treason\n")
            for r in verdict_rows:
                fh.write("\t".join(str(x).replace("\t", " ").replace("\n", " ") for x in r) + "\n")
        with open(out / "tasks_extended.tsv", "w") as fh:
            fh.write("task\tsubcommand\tconversion\tdomain\tinput\tother_tools\treference\tskipped\n")
            for r in task_rows:
                fh.write("\t".join(str(x) for x in r) + "\n")


if __name__ == "__main__":
    main()
