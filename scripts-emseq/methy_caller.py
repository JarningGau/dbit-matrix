#!/usr/bin/env python3
"""EMSeq methylation caller (pooled host or aggregate spike/mito)."""

from __future__ import annotations

import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from scripts import methy_call_core as core
from scripts.methy_caller import parse_args


def main(argv=None) -> int:
    args = parse_args(argv, description="EMSeq methylation caller")
    return core.run_caller(args, converted_means_methylated=False, log_prefix="emseq.methy_caller")


if __name__ == "__main__":
    raise SystemExit(main())
