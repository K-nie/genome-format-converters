#!/usr/bin/env python3
"""Supplementary table — hardware fingerprint for the Methods section.

Parses ``benchmarks/results/hardware.log`` (written by run_on_condor.sh
at the start of every cluster smoke) and emits a paper-ready
markdown table summarising the execute-host spec: kernel, CPU model,
core count, RAM, filesystem, key tool versions.

Paper utility: every benchmark paper needs a Methods section that
states exactly where the numbers came from. This pulls the relevant
fields out of the raw dump so the user doesn't paste 200 lines of
``conda list --export`` into the paper by accident.

If hardware.log is missing or unparseable, exits 1 with a clear
message rather than producing a half-empty table.
"""
from pathlib import Path
import re
import sys

LOG = Path(__file__).resolve().parent.parent / "results" / "hardware.log"
FIG = Path(__file__).resolve().parent.parent / "results" / "figures"


def _section(text: str, header: str) -> str:
    """Pull the body of a `=== <header> ===` section out of hardware.log."""
    pattern = rf"=== {re.escape(header)} ===\n(.*?)(?=\n=== |\Z)"
    m = re.search(pattern, text, flags=re.DOTALL)
    return m.group(1).strip() if m else ""


def _first_match(text: str, pat: str, group: int = 1) -> str:
    m = re.search(pat, text)
    return m.group(group).strip() if m else ""


def main() -> int:
    if not LOG.exists():
        print(f"[error] {LOG} missing — run_on_condor.sh hasn't run yet.",
              file=sys.stderr)
        return 1
    text = LOG.read_text()

    # Pull the human-readable fields. Each is best-effort; missing fields
    # render as "—" rather than crashing the table.
    host = _section(text, "host").splitlines()
    hostname = host[0] if host else "—"
    kernel = host[1] if len(host) > 1 else "—"

    cpu_block = _section(text, "cpu")
    cpu_model = _first_match(cpu_block, r"Model name:\s*(.+)") or "—"
    cpus = _first_match(cpu_block, r"CPU\(s\):\s*(\d+)") or "—"
    sockets = _first_match(cpu_block, r"Socket\(s\):\s*(\d+)") or "—"
    cores_per_socket = _first_match(cpu_block, r"Core\(s\) per socket:\s*(\d+)") or "—"
    threads_per_core = _first_match(cpu_block, r"Thread\(s\) per core:\s*(\d+)") or "—"

    mem_block = _section(text, "memory")
    total_ram = _first_match(mem_block, r"Mem:\s+(\S+)") or "—"

    disk_block = _section(text, "disk (benchmarks/)")
    fs = _first_match(disk_block, r"^(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+",
                     group=1) or "—"

    tools = _section(text, "tool versions").splitlines()
    tool_lines = [t for t in tools if t and not t.startswith("[")]

    lines = [
        "| Field | Value |",
        "|---|---|",
        f"| Host | `{hostname}` |",
        f"| Kernel | `{kernel}` |",
        f"| CPU | {cpu_model} |",
        f"| Sockets × cores × threads | {sockets} × {cores_per_socket} × {threads_per_core} (= {cpus} logical CPUs) |",
        f"| RAM (total) | {total_ram} |",
        f"| Filesystem (benchmarks/) | {fs} |",
        "| Key tool versions | <ul>",
    ]
    for tl in tool_lines[:10]:  # cap to keep the table reasonable
        lines.append(f"| | `{tl.strip()}` |")
    lines.append("")

    out_md = "\n".join(lines) + "\n"
    FIG.mkdir(parents=True, exist_ok=True)
    md_path = FIG / "hardware_table.md"
    md_path.write_text(out_md)
    print(f"[done] wrote {md_path}", file=sys.stderr)
    print("\n--- hardware table ---", file=sys.stderr)
    print(out_md)
    return 0


if __name__ == "__main__":
    sys.exit(main())
