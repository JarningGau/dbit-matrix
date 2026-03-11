#!/usr/bin/env python3
"""Run M-bias QC for host and spike-in BAMs."""

from __future__ import annotations

import argparse
import csv
import shlex
import subprocess
from collections import defaultdict
from pathlib import Path

import pysam


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate M-bias QC tables from pooled BAMs under a sample work directory. "
            "Host uses fixed-fraction subsampling; spike-ins use full BAMs."
        )
    )
    parser.add_argument(
        "--work-path",
        required=True,
        help="Sample work directory containing pooled/ outputs.",
    )
    parser.add_argument(
        "--mode",
        choices=["all", "host", "spike"],
        default="all",
        help="Run host + spike, or host/spike only. Default: all.",
    )
    parser.add_argument(
        "--spike-in-name",
        action="append",
        default=[],
        help=(
            "Spike-in name to process (e.g. lambda). "
            "May be specified multiple times; if omitted in spike mode, auto-discover."
        ),
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
        help="Threads for samtools view/sort. Default: 4.",
    )
    parser.add_argument(
        "--host-subsample-fraction",
        type=float,
        default=0.1,
        help="Host subsampling fraction in (0, 1]. Default: 0.1.",
    )
    parser.add_argument(
        "--host-subsample-seed",
        type=int,
        default=11,
        help="Host subsampling seed. Default: 11.",
    )
    parser.add_argument(
        "--max-cycle",
        type=int,
        default=150,
        help="Maximum read cycle to report. Default: 150.",
    )
    parser.add_argument(
        "--min-mapping-quality",
        type=int,
        default=1,
        help="Minimum mapping quality for counting. Default: 1.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory. Default: <work_path>/qc/mbias.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands only; do not execute.",
    )
    return parser.parse_args()


def quoted(parts: list[str]) -> str:
    return " ".join(shlex.quote(part) for part in parts)


def run_command(command: list[str], dry_run: bool) -> None:
    print(f"[mbias] command={quoted(command)}")
    if not dry_run:
        subprocess.run(command, check=True)


def run_or_skip(command: list[str], output_path: Path, dry_run: bool) -> bool:
    if output_path.exists():
        print(f"[mbias] skip_existing_output={output_path}")
        return False
    run_command(command, dry_run)
    return True


def discover_spike_bams(pooled_dir: Path) -> dict[str, Path]:
    spikes: dict[str, Path] = {}
    for bam_path in sorted(pooled_dir.glob("pooled.*.sorted.bam")):
        stem = bam_path.name[: -len(".sorted.bam")]
        if stem == "pooled.byCB":
            continue
        if not stem.startswith("pooled."):
            continue
        spike_name = stem.split(".", 1)[1].strip()
        if spike_name:
            spikes[spike_name] = bam_path
    return spikes


def build_subsample_fraction(seed: int, fraction: float) -> str:
    if fraction <= 0 or fraction > 1:
        raise ValueError("host-subsample-fraction must be in (0, 1]")
    if fraction == 1:
        return ""
    decimal_part = f"{fraction:.8f}".split(".", 1)[1].rstrip("0")
    if not decimal_part:
        decimal_part = "0"
    return f"{seed}.{decimal_part}"


def count_mbias_rows(
    bam_path: Path,
    max_cycle: int,
    min_mapping_quality: int,
) -> list[dict[str, object]]:
    counts: dict[tuple[str, int, str], list[int]] = defaultdict(lambda: [0, 0])
    with pysam.AlignmentFile(str(bam_path), "rb") as bam_file:
        for record in bam_file.fetch(until_eof=True):
            if record.is_unmapped:
                continue
            if record.is_secondary or record.is_supplementary:
                continue
            if record.mapping_quality < min_mapping_quality:
                continue
            read_seq = record.query_sequence
            if not read_seq:
                continue
            read_label = "R1" if record.is_read1 else "R2" if record.is_read2 else "R0"
            read_len = len(read_seq)
            if read_len < 2:
                continue
            for idx in range(read_len - 1):
                cycle = idx + 1 if not record.is_reverse else read_len - idx
                if cycle < 1 or cycle > max_cycle:
                    continue
                dinucleotide = (read_seq[idx] + read_seq[idx + 1]).upper()
                if dinucleotide in ("TG", "CA"):
                    key = (read_label, cycle, dinucleotide)
                    counts[key][0] += 1
                elif dinucleotide == "CG":
                    key = (read_label, cycle, dinucleotide)
                    counts[key][1] += 1

    rows: list[dict[str, object]] = []
    for (read_label, cycle, context), (methylated, unmethylated) in sorted(counts.items()):
        coverage = methylated + unmethylated
        methylation_rate = (methylated / coverage) if coverage > 0 else 0.0
        rows.append(
            {
                "read": read_label,
                "cycle": cycle,
                "context": context,
                "methylated_count": methylated,
                "unmethylated_count": unmethylated,
                "coverage": coverage,
                "methylation_rate": round(methylation_rate, 6),
            }
        )
    return rows


def write_mbias_tsv(output_path: Path, rows: list[dict[str, object]]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "read",
        "cycle",
        "context",
        "methylated_count",
        "unmethylated_count",
        "coverage",
        "methylation_rate",
    ]
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def run_host_mode(args: argparse.Namespace, work_path: Path, output_dir: Path) -> None:
    pooled_host = work_path / "pooled" / "pooled.byCB.bam"
    if not pooled_host.exists():
        raise ValueError(f"missing host pooled BAM: {pooled_host}")

    output_dir.mkdir(parents=True, exist_ok=True)
    host_subsampled = output_dir / "host.subsampled.bam"
    host_sorted = output_dir / "host.subsampled.sorted.bam"
    host_index = output_dir / "host.subsampled.sorted.bam.bai"
    host_tsv = output_dir / "host.mbias.tsv"

    subsample_arg = build_subsample_fraction(
        args.host_subsample_seed, args.host_subsample_fraction
    )
    if subsample_arg:
        subsample_cmd = [
            args.samtools_bin,
            "view",
            "-@",
            str(args.samtools_threads),
            "-s",
            subsample_arg,
            "-b",
            "-o",
            str(host_subsampled),
            str(pooled_host),
        ]
    else:
        subsample_cmd = [
            args.samtools_bin,
            "view",
            "-@",
            str(args.samtools_threads),
            "-b",
            "-o",
            str(host_subsampled),
            str(pooled_host),
        ]
    sort_cmd = [
        args.samtools_bin,
        "sort",
        "-@",
        str(args.samtools_threads),
        "-o",
        str(host_sorted),
        str(host_subsampled),
    ]
    index_cmd = [args.samtools_bin, "index", str(host_sorted)]

    run_or_skip(subsample_cmd, host_subsampled, args.dry_run)
    run_or_skip(sort_cmd, host_sorted, args.dry_run)
    run_or_skip(index_cmd, host_index, args.dry_run)

    if host_tsv.exists():
        print(f"[mbias] skip_existing_output={host_tsv}")
    else:
        print(f"[mbias] host_input={host_sorted}")
        print(f"[mbias] host_output={host_tsv}")
        if not args.dry_run:
            rows = count_mbias_rows(
                host_sorted,
                max_cycle=args.max_cycle,
                min_mapping_quality=args.min_mapping_quality,
            )
            write_mbias_tsv(host_tsv, rows)


def run_spike_mode(args: argparse.Namespace, work_path: Path, output_dir: Path) -> None:
    pooled_dir = work_path / "pooled"
    discovered_spikes = discover_spike_bams(pooled_dir)
    if not discovered_spikes:
        raise ValueError(f"no spike-in BAMs found under: {pooled_dir}/pooled.*.sorted.bam")

    requested = [name.strip() for name in args.spike_in_name if name.strip()]
    spike_names = requested if requested else sorted(discovered_spikes.keys())
    print(f"[mbias] spike_in_count={len(spike_names)}")
    for spike_name in spike_names:
        spike_bam = discovered_spikes.get(spike_name)
        if spike_bam is None:
            raise ValueError(
                f"requested spike-in BAM missing for '{spike_name}': "
                f"{pooled_dir}/pooled.{spike_name}.sorted.bam"
            )
        output_path = output_dir / f"{spike_name}.mbias.tsv"
        if output_path.exists():
            print(f"[mbias] skip_existing_output={output_path}")
            continue
        print(f"[mbias] spike_in={spike_name}")
        print(f"[mbias] spike_input={spike_bam}")
        print(f"[mbias] output={output_path}")
        if not args.dry_run:
            rows = count_mbias_rows(
                spike_bam,
                max_cycle=args.max_cycle,
                min_mapping_quality=args.min_mapping_quality,
            )
            write_mbias_tsv(output_path, rows)


def main() -> int:
    args = parse_args()
    if args.samtools_threads <= 0:
        raise ValueError("samtools-threads must be > 0")
    if args.host_subsample_seed < 0:
        raise ValueError("host-subsample-seed must be >= 0")
    if args.max_cycle <= 0:
        raise ValueError("max-cycle must be > 0")
    if args.min_mapping_quality < 0:
        raise ValueError("min-mapping-quality must be >= 0")

    work_path = Path(args.work_path)
    output_dir = (
        Path(args.output_dir) if args.output_dir else work_path / "qc" / "mbias"
    )
    print(f"[mbias] mode={args.mode}")
    print(f"[mbias] work_path={work_path}")
    print(f"[mbias] output_dir={output_dir}")
    print(f"[mbias] host_subsample_fraction={args.host_subsample_fraction}")
    print(f"[mbias] host_subsample_seed={args.host_subsample_seed}")
    print(f"[mbias] max_cycle={args.max_cycle}")
    print(f"[mbias] min_mapping_quality={args.min_mapping_quality}")
    if args.dry_run:
        print("[mbias] dry_run=1")

    if args.mode in ("all", "host"):
        run_host_mode(args, work_path, output_dir)
    if args.mode in ("all", "spike"):
        run_spike_mode(args, work_path, output_dir)

    print("[mbias] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
