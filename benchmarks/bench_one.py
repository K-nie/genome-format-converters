#!/usr/bin/env python3
"""Run one benchmark replicate and emit a single TSV row on stdout.

Called by the per-task shell scripts under ``benchmarks/tasks/``. One row
per invocation keeps the per-task TSVs append-only and trivially
concatenatable.

Measurements:
  wall_s       wall-clock seconds via ``time.perf_counter()``.
  peak_rss_mb  ``ru_maxrss`` from the child process via ``os.wait4``. This
               is the kernel-tracked peak RSS the OS records for the
               specific child, not a polled snapshot of the live RSS — so
               mmap-heavy tools (plink2) report what they actually
               resident-set, not what they mmap'd. Units: KB on Linux,
               bytes on macOS / *BSD; normalised to bytes here.

Schema (matches benchmarks/README.md):
  task tool version replicate wall_s peak_rss_mb exit_code correct notes
"""
from __future__ import annotations

import argparse
import os
import resource
import subprocess
import sys
import time


def run_once(cmd: str, *, task: str, tool: str,
             replicate: int) -> tuple[float, int, int]:
    """Run `cmd` under a shell. Return (wall_seconds, peak_rss_bytes, exit_code).

    Uses ``os.wait4`` so we get the rusage of the specific child, not the
    accumulated rusage of all children since process start (which is what
    ``resource.getrusage(RUSAGE_CHILDREN)`` would give).

    Subprocess stderr is streamed to a per-run log file under
    ``$GFC_BENCH_LOG_DIR`` (default /tmp/gfc_bench_cmd_logs/). On non-zero
    exit, the tail of that log is echoed to this process's own stderr so
    it shows up in the Condor ``bench.<id>.err`` capture and the failure
    is debuggable without having to ssh to the execute slot. Successful
    runs delete the log file to avoid clutter.

    Streaming to a file rather than capturing via subprocess.PIPE keeps
    the parent's RSS low (no PIPE buffering in this process), so the
    benchmarked child's RSS measurement isn't perturbed.
    """
    log_dir = os.environ.get("GFC_BENCH_LOG_DIR", "/tmp/gfc_bench_cmd_logs")
    os.makedirs(log_dir, exist_ok=True)
    log_path = os.path.join(log_dir, f"{task}_{tool}_rep{replicate}.err")

    start = time.perf_counter()
    with open(log_path, "wb") as err_log:
        # shell=True so the command string can use pipes / redirects when
        # tasks need them (e.g. `bcftools view | awk | gzip > out`).
        proc = subprocess.Popen(
            cmd, shell=True,
            stdout=subprocess.DEVNULL,
            stderr=err_log,
        )
        _pid, status, ru = os.wait4(proc.pid, 0)
        elapsed = time.perf_counter() - start
    rc = os.waitstatus_to_exitcode(status)
    # Tell Popen the child is reaped so its destructor doesn't complain.
    proc.returncode = rc

    if rc != 0:
        try:
            with open(log_path, errors="replace") as f:
                tail = f.read()[-2000:]
        except OSError:
            tail = "(could not read stderr log)"
        print(
            f"[FAIL] {task}/{tool}/rep{replicate} exit={rc}\n"
            f"---- stderr tail ({log_path}) ----\n"
            f"{tail}\n"
            f"---- end {task}/{tool}/rep{replicate} stderr ----",
            file=sys.stderr,
        )
    else:
        try:
            os.unlink(log_path)
        except OSError:
            pass

    # ru_maxrss units differ by platform: KB on Linux, bytes on macOS/BSD.
    if sys.platform == "darwin" or sys.platform.startswith("freebsd"):
        peak_rss_bytes = ru.ru_maxrss
    else:
        peak_rss_bytes = ru.ru_maxrss * 1024

    return elapsed, peak_rss_bytes, rc


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

    wall_s, peak_rss_b, rc = run_once(
        args.cmd, task=args.task, tool=args.tool, replicate=args.replicate,
    )
    peak_rss_mb = peak_rss_b / (1024 * 1024)
    row = [
        args.task, args.tool, args.version, str(args.replicate),
        f"{wall_s:.4f}", f"{peak_rss_mb:.2f}",
        str(rc), args.correct, args.notes,
    ]
    print("\t".join(row))


if __name__ == "__main__":
    main()
