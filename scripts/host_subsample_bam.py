#!/usr/bin/env python3
"""Shared helpers for preparing host subsampled BAMs."""

from __future__ import annotations

import shlex
from dataclasses import dataclass
from pathlib import Path

HOST_SUBSAMPLE_SEED = 11
# Cap alignment records after fractional subsampling (not user-configurable).
HOST_SUBSAMPLE_MAX_READS = 10_000_000
HOST_SUBSAMPLED_BAM_BASENAME = "host.subsampled.bam"
HOST_SUBSAMPLED_SORTED_BAM_BASENAME = "host.subsampled.sorted.bam"


@dataclass(frozen=True)
class HostSubsamplePaths:
    subsampled_bam: Path
    sorted_bam: Path
    sorted_bam_index: Path


def get_host_subsample_paths(output_dir: Path) -> HostSubsamplePaths:
    subsampled_bam = output_dir / HOST_SUBSAMPLED_BAM_BASENAME
    sorted_bam = output_dir / HOST_SUBSAMPLED_SORTED_BAM_BASENAME
    sorted_bam_index = output_dir / f"{HOST_SUBSAMPLED_SORTED_BAM_BASENAME}.bai"
    return HostSubsamplePaths(
        subsampled_bam=subsampled_bam,
        sorted_bam=sorted_bam,
        sorted_bam_index=sorted_bam_index,
    )


def build_subsample_fraction(seed: int, fraction: float) -> str:
    if fraction <= 0 or fraction > 1:
        raise ValueError("host-subsample-fraction must be in (0, 1]")
    if fraction == 1:
        return ""
    decimal_part = f"{fraction:.8f}".split(".", 1)[1].rstrip("0")
    if not decimal_part:
        decimal_part = "0"
    return f"{seed}.{decimal_part}"


def _build_subsample_pipeline_shell(
    samtools_bin: str,
    samtools_threads: int,
    pooled_host_bam: Path,
    out_bam: Path,
    host_subsample_fraction: float,
    host_subsample_seed: int,
) -> str:
    """SAM: optional -s subsample, then cap alignment lines; convert to BAM with header."""
    subsample_arg = build_subsample_fraction(host_subsample_seed, host_subsample_fraction)
    header_cmd = shlex.join([samtools_bin, "view", "-H", str(pooled_host_bam)])
    view1: list[str] = [samtools_bin, "view", "-@", str(samtools_threads)]
    if subsample_arg:
        view1.extend(["-s", subsample_arg])
    view1.append(str(pooled_host_bam))
    tail_cmd = shlex.join(
        [
            samtools_bin,
            "view",
            "-@",
            str(samtools_threads),
            "-b",
            "-o",
            str(out_bam),
            "-",
        ]
    )
    inner = f"{shlex.join(view1)} | head -n {HOST_SUBSAMPLE_MAX_READS}"
    return f"{{ {header_cmd}; {inner}; }} | {tail_cmd}"


def build_prepare_host_subsample_commands(
    pooled_host_bam: Path,
    paths: HostSubsamplePaths,
    samtools_bin: str,
    samtools_threads: int,
    host_subsample_fraction: float,
    host_subsample_seed: int = HOST_SUBSAMPLE_SEED,
) -> list[tuple[list[str], Path]]:
    pipeline = _build_subsample_pipeline_shell(
        samtools_bin=samtools_bin,
        samtools_threads=samtools_threads,
        pooled_host_bam=pooled_host_bam,
        out_bam=paths.subsampled_bam,
        host_subsample_fraction=host_subsample_fraction,
        host_subsample_seed=host_subsample_seed,
    )
    subsample_cmd = ["/bin/sh", "-c", pipeline]
    sort_cmd = [
        samtools_bin,
        "sort",
        "-@",
        str(samtools_threads),
        "-o",
        str(paths.sorted_bam),
        str(paths.subsampled_bam),
    ]
    index_cmd = [samtools_bin, "index", str(paths.sorted_bam)]
    return [
        (subsample_cmd, paths.subsampled_bam),
        (sort_cmd, paths.sorted_bam),
        (index_cmd, paths.sorted_bam_index),
    ]
