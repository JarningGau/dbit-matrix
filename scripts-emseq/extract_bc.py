#!/usr/bin/env python3
"""Deprecated: EMSeq demux is unified into `scripts/extract_bc.py`.

This wrapper keeps backward compatibility with older entrypoints/automation.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def main() -> int:
    taps_extract = Path(__file__).resolve().parents[1] / "scripts" / "extract_bc.py"
    cmd = [sys.executable, str(taps_extract), *sys.argv[1:]]
    return subprocess.call(cmd)


if __name__ == "__main__":
    raise SystemExit(main())

