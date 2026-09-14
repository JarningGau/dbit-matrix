#!/usr/bin/env python3
"""TAPS methylation caller (pooled host or aggregate spike/mito)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional, Tuple

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from scripts import methy_call_core as core
from scripts.methy_call_core import (
    bam_index_exists,
    ch_neighbor_matches,
    classify_reference_position,
    has_unclear_conversion_strand,
    load_barcode_whitelist,
    select_sample_positions,
    shuffle_host_parts,
    spot_id_from_cb,
    spot_is_allowed,
)

CONVERTED_MEANS_METHYLATED = True


def parse_args(
    argv: Optional[List[str]] = None,
    description: str = "TAPS methylation caller",
) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--mode", choices=("host", "aggregate"), required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("-c", "--chromosomes", required=True)
    parser.add_argument(
        "--context-mode", choices=("cg", "ch", "both"), default="cg"
    )
    parser.add_argument("--min-base-quality", type=int, default=30)
    parser.add_argument("--min-mapping-quality", type=int, default=10)
    parser.add_argument("--max-depth", type=int, default=50_000)
    parser.add_argument("--batch-size", type=int, default=10_000_000)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--sample-size", type=int, default=None)
    parser.add_argument("--r1-left-trimming", type=int, default=0)
    parser.add_argument("--r1-right-trimming", type=int, default=0)
    parser.add_argument("--r2-left-trimming", type=int, default=0)
    parser.add_argument("--r2-right-trimming", type=int, default=0)
    parser.add_argument("--cb-tag", default="CB")
    parser.add_argument("--barcodes", help="Whitelist TSV for host mode.")
    parser.add_argument(
        "--spot-roster",
        help="Optional .cov write-set (counts TSV or one XX_YY per line).",
    )
    parser.add_argument(
        "--counts-tsv",
        help="Host counts output. Default: <out-dir>/../per_spot_read_counts.tsv",
    )
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--force-recount", action="store_true")
    parser.add_argument("--out-dir", help="Host output directory (coverage/host).")
    parser.add_argument("--output", help="Aggregate CG/CH stem path without suffix.")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args(argv)


def classify_cg_observation(dinucleotide: str) -> Optional[bool]:
    return core.classify_cg_observation(dinucleotide, CONVERTED_MEANS_METHYLATED)


def classify_ch_observation(base: str, strand: str) -> Optional[Tuple[int, int]]:
    return core.classify_ch_observation(base, strand, CONVERTED_MEANS_METHYLATED)


def main(argv: Optional[List[str]] = None) -> int:
    return core.run_caller(parse_args(argv), True, "methy_caller")


if __name__ == "__main__":
    raise SystemExit(main())
