#!/usr/bin/env python3
"""Run fastp to split paired FASTQ files into chunks."""

from __future__ import annotations

import argparse
import shlex
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Split paired FASTQ files with fastp for downstream DBiT workflow."
    )
    parser.add_argument("--r1", required=True, help="Input R1 FASTQ(.gz) path.")
    parser.add_argument("--r2", required=True, help="Input R2 FASTQ(.gz) path.")
    parser.add_argument(
        "--work-path",
        required=True,
        help="Work directory. Outputs are written into <work-path>/shard_fastq.",
    )
    parser.add_argument(
        "--fastp-threads",
        type=int,
        default=8,
        help="Thread count passed to fastp -w. Default: 8.",
    )
    parser.add_argument(
        "--number-of-split-parts",
        type=int,
        required=True,
        help="Value passed to fastp --split.",
    )
    parser.add_argument(
        "--fastp-bin",
        default="fastp",
        help="fastp executable path or name in PATH. Default: fastp.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Only print the generated command without executing it.",
    )
    return parser.parse_args()


def build_fastp_cmd(args: argparse.Namespace) -> tuple[list[str], Path]:
    work_path = Path(args.work_path)
    shard_dir = work_path / "shard_fastq"
    shard_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        args.fastp_bin,
        "-i",
        args.r1,
        "-I",
        args.r2,
        "-o",
        str(shard_dir / "R1.fq.gz"),
        "-O",
        str(shard_dir / "R2.fq.gz"),
        "-w",
        str(args.fastp_threads),
        "--split",
        str(args.number_of_split_parts),
        "--disable_adapter_trimming",
        "-h",
        str(shard_dir / "fastp.html"),
        "-j",
        str(shard_dir / "fastp.json"),
    ]
    return cmd, shard_dir


def main() -> int:
    args = parse_args()
    cmd, shard_dir = build_fastp_cmd(args)
    cmd_string = " ".join(shlex.quote(part) for part in cmd)

    print(f"[fastp_split] output_dir={shard_dir}")
    print(f"[fastp_split] command={cmd_string}")

    if args.dry_run:
        return 0

    subprocess.run(cmd, check=True)
    print("[fastp_split] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
