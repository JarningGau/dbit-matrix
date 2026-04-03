#!/usr/bin/env python3
"""Thin wrapper: forwards to scripts/methscan_run.py prepare (same CLI as before)."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def main() -> int:
    run_script = Path(__file__).resolve().parent / "methscan_run.py"
    cmd = [sys.executable, str(run_script), "prepare", *sys.argv[1:]]
    return subprocess.run(cmd).returncode


if __name__ == "__main__":
    raise SystemExit(main())
