#!/usr/bin/env python3
"""Shared helpers for preparing host subsampled BAMs."""

from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path

HOST_SUBSAMPLE_SEED = 11
# Cap expected alignment records via idxstats-adjusted -s (not user-configurable).
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


def count_mapped_reads(bam: Path, samtools_bin: str) -> int:
    """Sum mapped reads from ``samtools idxstats`` column 3."""
    result = subprocess.run(
        [samtools_bin, "idxstats", str(bam)],
        check=True,
        capture_output=True,
        text=True,
    )
    mapped = 0
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        mapped += int(parts[2])
    return mapped


def resolve_effective_subsample_fraction(
    requested: float,
    mapped: int,
    max_reads: int = HOST_SUBSAMPLE_MAX_READS,
) -> float:
    """Return name-hash ``-s`` fraction, capped so expected reads stay near max_reads."""
    if requested <= 0 or requested > 1:
        raise ValueError("host-subsample-fraction must be in (0, 1]")
    if max_reads <= 0:
        raise ValueError("max_reads must be > 0")
    if mapped <= 0:
        return requested
    return min(requested, max_reads / mapped)


def resolve_host_subsample_fraction(
    pooled_host_bam: Path,
    requested_fraction: float,
    samtools_bin: str,
    max_reads: int = HOST_SUBSAMPLE_MAX_READS,
) -> float:
    """Resolve effective fraction from idxstats when the BAM is available."""
    if not pooled_host_bam.is_file():
        return requested_fraction
    try:
        mapped = count_mapped_reads(pooled_host_bam, samtools_bin)
    except (OSError, subprocess.CalledProcessError, ValueError):
        return requested_fraction
    return resolve_effective_subsample_fraction(
        requested_fraction, mapped, max_reads=max_reads
    )


def build_subsample_view_command(
    samtools_bin: str,
    samtools_threads: int,
    pooled_host_bam: Path,
    out_bam: Path,
    host_subsample_fraction: float,
    host_subsample_seed: int,
) -> list[str]:
    """Build ``samtools view [-s] -b -o`` (no coordinate-ordered head cap)."""
    subsample_arg = build_subsample_fraction(
        host_subsample_seed, host_subsample_fraction
    )
    command = [
        samtools_bin,
        "view",
        "-@",
        str(samtools_threads),
        "-b",
        "-o",
        str(out_bam),
    ]
    if subsample_arg:
        command.extend(["-s", subsample_arg])
    command.append(str(pooled_host_bam))
    return command


def build_prepare_host_subsample_commands(
    pooled_host_bam: Path,
    paths: HostSubsamplePaths,
    samtools_bin: str,
    samtools_threads: int,
    host_subsample_fraction: float,
    host_subsample_seed: int = HOST_SUBSAMPLE_SEED,
    effective_fraction: float | None = None,
) -> list[tuple[list[str], Path]]:
    fraction = (
        host_subsample_fraction
        if effective_fraction is None
        else effective_fraction
    )
    subsample_cmd = build_subsample_view_command(
        samtools_bin=samtools_bin,
        samtools_threads=samtools_threads,
        pooled_host_bam=pooled_host_bam,
        out_bam=paths.subsampled_bam,
        host_subsample_fraction=fraction,
        host_subsample_seed=host_subsample_seed,
    )
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
