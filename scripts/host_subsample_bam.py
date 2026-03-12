#!/usr/bin/env python3
"""Shared helpers for preparing host subsampled BAMs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

HOST_SUBSAMPLE_SEED = 11
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


def build_prepare_host_subsample_commands(
    pooled_host_bam: Path,
    paths: HostSubsamplePaths,
    samtools_bin: str,
    samtools_threads: int,
    host_subsample_fraction: float,
    host_subsample_seed: int = HOST_SUBSAMPLE_SEED,
) -> list[tuple[list[str], Path]]:
    subsample_arg = build_subsample_fraction(host_subsample_seed, host_subsample_fraction)
    if subsample_arg:
        subsample_cmd = [
            samtools_bin,
            "view",
            "-@",
            str(samtools_threads),
            "-s",
            subsample_arg,
            "-b",
            "-o",
            str(paths.subsampled_bam),
            str(pooled_host_bam),
        ]
    else:
        subsample_cmd = [
            samtools_bin,
            "view",
            "-@",
            str(samtools_threads),
            "-b",
            "-o",
            str(paths.subsampled_bam),
            str(pooled_host_bam),
        ]
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
