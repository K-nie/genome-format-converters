#!/usr/bin/env python3
"""Run one benchmark replicate and emit a single TSV row on stdout.

Called by the per-task shell scripts under ``benchmarks/tasks/``. One row
per invocation keeps the per-task TSVs append-only and trivially
concatenatable.

Measurements:
  wall_s       wall-clock seconds via ``time.perf_counter()``.
  peak_rss_mb  peak resident set size of the child process *and all its
               descendants*, polled every 50 ms via ``psutil``. Falls back
               to 0 if ``psutil`` isn't available.

Schema (matches benchmarks/README.md):
  task tool version replicate wall_s peak_rss_mb exit_code correct notes
"""
from __future__ import annotations

import argparse
import shlex
import subprocess
import sys
import threading
import time
from typing import Optional

try:
    import psutil  # type: ignore
except ImportError:  # pragma: no cover
    psutil = None


def _poll_peak_rss(pid: int, stop_event: threading.Event,
                   peak: dict, interval_s: float = 0.05) -> None:
    """Background thread: poll the process tree's RSS every `interval_s`
    seconds, tracking the highest value seen."""
    if psutil is None:
        return
    try:
        proc = psutil.Process(pid)
    except psutil.NoSuchProcess:
        return
    while not stop_event.is_set():
        try:
            total = proc.memory_info().rss
            for child in proc.children(recursive=True):
                try:
                    total += child.memory_info().rss
                except psutil.NoSuchProcess:
                    continue
            if total > peak["rss"]:
                peak["rss"] = total
        except psutil.NoSuchProcess:
            break
        time.sleep(interval_s)


def run_once(cmd: str) -> tuple[float, int, int]:
    """Run `cmd` under a shell. Return (wall_seconds, peak_rss_bytes, exit_code)."""
    peak = {"rss": 0}
    stop_event = threading.Event()

    start = time.perf_counter()
    # shell=True so the command string can use pipes / redirects when
    # tasks need them (e.g. `bcftools view | awk | gzip > out`).
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL)
    poller = threading.Thread(
        target=_poll_peak_rss, args=(proc.pid, stop_event, peak)
    )
    poller.start()

    rc = proc.wait()
    stop_event.set()
    poller.join(timeout=1.0)
    elapsed = time.perf_counter() - start
    return elapsed, peak["rss"], rc


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    # All args are optional at the argparse level; we enforce the "you need
    # these to actually run a benchmark" invariant in code below, so that
    # `--header` stays usable on its own.
    ap.add_argument("--task", help="Task id (e.g. T1)")
    ap.add_argument("--tool", help="Tool name (e.g. gfc)")
    ap.add_argument("--version", help="Tool version string")
    ap.add_argument("--replicate", type=int)
    ap.add_argument("--cmd", help="Full shell command string to time. Quote it once.")
    ap.add_argument("--correct", default="",
                    help="Optional correctness flag: '1' / '0' / empty.")
    ap.add_argument("--notes", default="")
    ap.add_argument("--header", action="store_true",
                    help="Just print the TSV header line and exit.")
    args = ap.parse_args()

    cols = ["task", "tool", "version", "replicate",
            "wall_s", "peak_rss_mb", "exit_code", "correct", "notes"]
    if args.header:
        print("\t".join(cols))
        return

    missing = [n for n in ("task", "tool", "version", "replicate", "cmd")
               if getattr(args, n) is None]
    if missing:
        ap.error(f"missing required args for a measurement run: {', '.join(missing)}")

    wall_s, peak_rss_b, rc = run_once(args.cmd)
    peak_rss_mb = peak_rss_b / (1024 * 1024)
    row = [
        args.task, args.tool, args.version, str(args.replicate),
        f"{wall_s:.4f}", f"{peak_rss_mb:.2f}",
        str(rc), args.correct, args.notes,
    ]
    print("\t".join(row))


if __name__ == "__main__":
    main()
