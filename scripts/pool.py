#!/usr/bin/env python3
"""Pool host/spike-in BAM shards for DBiT workflow."""

from __future__ import annotations

import argparse
import shlex
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Pool BAM shards under <work-path>/align_shards into "
            "<work-path>/pooled outputs."
        )
    )
    parser.add_argument(
        "--work-path",
        required=True,
        help="Sample work directory containing align_shards/.",
    )
    parser.add_argument(
        "--samtools-bin",
        default="samtools",
        help="samtools executable path or command name. Default: samtools.",
    )
    parser.add_argument(
        "--samtools-threads",
        type=int,
        default=4,
        help="Thread count for samtools sort -@. Default: 4.",
    )
    parser.add_argument(
        "--host-sort-mem",
        default="16G",
        help="Memory per thread for host sort (-m). Default: 16G.",
    )
    parser.add_argument(
        "--spike-in-name",
        action="append",
        default=[],
        help=(
            "Spike-in name to pool (e.g. lambda). "
            "May be specified multiple times."
        ),
    )
    parser.add_argument(
        "--mode",
        choices=["all", "host", "spike"],
        default="all",
        help="Run host + spike pooling, or host/spike only. Default: all.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands only; do not execute pooling.",
    )
    return parser.parse_args()


def quoted(parts: list[str]) -> str:
    return " ".join(shlex.quote(part) for part in parts)


def run_command(command: list[str]) -> None:
    subprocess.run(command, check=True)


def discover_host_shards(align_dir: Path) -> list[Path]:
    return sorted(align_dir.glob("*.cb.bam"))


def discover_spike_shards(align_dir: Path, spike_name: str) -> list[Path]:
    return sorted(align_dir.glob(f"*.{spike_name}.bam"))


def main() -> int:
    args = parse_args()
    if args.samtools_threads <= 0:
        raise ValueError("samtools_threads must be > 0")
    spike_names = [name.strip() for name in args.spike_in_name if name.strip()]
    if args.mode in ("all", "spike") and not spike_names:
        raise ValueError(
            "spike mode requires at least one --spike-in-name (NAME from spike_in_index)"
        )

    work_path = Path(args.work_path)
    align_dir = work_path / "align_shards"
    pooled_dir = work_path / "pooled"
    pooled_dir.mkdir(parents=True, exist_ok=True)

    print(f"[pool] mode={args.mode}")
    print(f"[pool] work_path={work_path}")
    print(f"[pool] align_dir={align_dir}")
    print(f"[pool] pooled_dir={pooled_dir}")

    # Keep local behavior predictable: spike-in first, then host.
    if args.mode in ("all", "spike"):
        print(f"[pool] spike_in_count={len(spike_names)}")
        for spike_name in spike_names:
            spike_shards = discover_spike_shards(align_dir, spike_name)
            if not spike_shards:
                raise ValueError(
                    f"no spike-in shards found for '{spike_name}': "
                    f"{align_dir}/*.{spike_name}.bam"
                )
            spike_tmp = pooled_dir / f"pooled.{spike_name}.bam"
            spike_sorted = pooled_dir / f"pooled.{spike_name}.sorted.bam"
            cat_cmd = [
                args.samtools_bin,
                "cat",
                "-o",
                str(spike_tmp),
                *[str(path) for path in spike_shards],
            ]
            sort_cmd = [
                args.samtools_bin,
                "sort",
                "-@",
                str(args.samtools_threads),
                "-o",
                str(spike_sorted),
                str(spike_tmp),
            ]
            index_cmd = [args.samtools_bin, "index", str(spike_sorted)]
            rm_cmd = ["rm", str(spike_tmp)]
            print(f"[pool] spike_in={spike_name}")
            print(f"[pool] spike_shard_count={len(spike_shards)}")
            print(f"[pool] command={quoted(cat_cmd)}")
            print(f"[pool] command={quoted(sort_cmd)}")
            print(f"[pool] command={quoted(index_cmd)}")
            print(f"[pool] command={quoted(rm_cmd)}")
            if not args.dry_run:
                run_command(cat_cmd)
                run_command(sort_cmd)
                run_command(index_cmd)
                spike_tmp.unlink(missing_ok=True)

    if args.mode in ("all", "host"):
        host_shards = discover_host_shards(align_dir)
        if not host_shards:
            raise ValueError(f"no host shards found: {align_dir}/*.cb.bam")
        host_tmp = pooled_dir / "pooled.cb.bam"
        host_sorted = pooled_dir / "pooled.byCB.bam"
        cat_cmd = [
            args.samtools_bin,
            "cat",
            "-o",
            str(host_tmp),
            *[str(path) for path in host_shards],
        ]
        sort_cmd = [
            args.samtools_bin,
            "sort",
            "-m",
            args.host_sort_mem,
            "-@",
            str(args.samtools_threads),
            "-t",
            "CB",
            "-o",
            str(host_sorted),
            str(host_tmp),
        ]
        rm_cmd = ["rm", str(host_tmp)]
        print("[pool] target=host")
        print(f"[pool] host_shard_count={len(host_shards)}")
        print(f"[pool] command={quoted(cat_cmd)}")
        print(f"[pool] command={quoted(sort_cmd)}")
        print(f"[pool] command={quoted(rm_cmd)}")
        if not args.dry_run:
            run_command(cat_cmd)
            run_command(sort_cmd)
            host_tmp.unlink(missing_ok=True)

    print("[pool] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
