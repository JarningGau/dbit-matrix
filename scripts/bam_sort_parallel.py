#!/usr/bin/env python3
"""Sort per-spot BAM files in parallel and index the sorted outputs."""

from __future__ import annotations

import argparse
import shlex
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Find unsorted BAM files under a split output directory, then run "
            "samtools sort + samtools index in parallel and remove the original BAM."
        )
    )
    parser.add_argument(
        "--work-path",
        help=(
            "Sample work directory. If provided, BAMs are discovered under "
            "<work-path>/split_bams by default."
        ),
    )
    parser.add_argument(
        "--bam-dir",
        help=(
            "Directory containing per-spot BAM files. Overrides --work-path "
            "derived location when provided."
        ),
    )
    parser.add_argument(
        "--samtools-bin",
        default="samtools",
        help="samtools executable path or command name. Default: samtools.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=8,
        help="Maximum number of BAMs to sort concurrently. Default: 8.",
    )
    parser.add_argument(
        "--sort-threads",
        type=int,
        default=1,
        help="Thread count passed to each samtools sort -@. Default: 1.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands only; do not execute sorting.",
    )
    return parser.parse_args()


def quoted(parts: list[str]) -> str:
    return " ".join(shlex.quote(part) for part in parts)


def resolve_bam_dir(args: argparse.Namespace) -> Path:
    if args.bam_dir:
        return Path(args.bam_dir)
    if args.work_path:
        return Path(args.work_path) / "split_bams"
    raise ValueError("provide either --bam-dir or --work-path")


def discover_bams(bam_dir: Path) -> list[Path]:
    return sorted(
        path
        for path in bam_dir.rglob("*.bam")
        if not path.name.endswith(".sorted.bam")
    )


def build_commands(
    samtools_bin: str,
    bam_path: Path,
    sort_threads: int,
) -> tuple[list[str], list[str]]:
    sorted_bam = bam_path.with_suffix(".sorted.bam")
    sort_cmd = [
        samtools_bin,
        "sort",
        "-@",
        str(sort_threads),
        "-o",
        str(sorted_bam),
        str(bam_path),
    ]
    index_cmd = [samtools_bin, "index", str(sorted_bam)]
    return sort_cmd, index_cmd


def run_one_bam(
    samtools_bin: str,
    bam_path: Path,
    sort_threads: int,
) -> tuple[Path, Path]:
    sort_cmd, index_cmd = build_commands(samtools_bin, bam_path, sort_threads)
    subprocess.run(sort_cmd, check=True)
    subprocess.run(index_cmd, check=True)
    bam_path.unlink()
    return bam_path, bam_path.with_suffix(".sorted.bam")


def main() -> int:
    args = parse_args()
    if args.jobs <= 0:
        raise ValueError("jobs must be > 0")
    if args.sort_threads <= 0:
        raise ValueError("sort_threads must be > 0")

    bam_dir = resolve_bam_dir(args)
    bam_paths = discover_bams(bam_dir)
    if not bam_paths:
        raise ValueError(f"no unsorted BAM files found under: {bam_dir}")

    print(f"[bam_sort_parallel] bam_dir={bam_dir}")
    print(f"[bam_sort_parallel] bam_count={len(bam_paths)}")
    print(f"[bam_sort_parallel] jobs={args.jobs}")
    print(f"[bam_sort_parallel] sort_threads={args.sort_threads}")

    for bam_path in bam_paths:
        sort_cmd, index_cmd = build_commands(args.samtools_bin, bam_path, args.sort_threads)
        print(f"[bam_sort_parallel] command={quoted(sort_cmd)}")
        print(f"[bam_sort_parallel] command={quoted(index_cmd)}")
        print(f"[bam_sort_parallel] command=rm {shlex.quote(str(bam_path))}")

    if args.dry_run:
        print("[bam_sort_parallel] dry_run=1")
        return 0

    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(
                run_one_bam,
                args.samtools_bin,
                bam_path,
                args.sort_threads,
            ): bam_path
            for bam_path in bam_paths
        }
        for future in as_completed(futures):
            bam_path, sorted_bam = future.result()
            print(f"[bam_sort_parallel] sorted={sorted_bam}")
            print(f"[bam_sort_parallel] removed={bam_path}")

    print("[bam_sort_parallel] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
