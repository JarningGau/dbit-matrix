#!/usr/bin/env python3
"""Run methylation calling for host spots and spike-ins."""

from __future__ import annotations

import argparse
import shlex
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pysam


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Call methylation sites from split host BAMs and pooled spike-in BAMs "
            "under a sample work directory."
        )
    )
    parser.add_argument(
        "--work-path",
        required=True,
        help="Sample work directory containing split_bams/ and pooled/.",
    )
    parser.add_argument(
        "--mode",
        choices=["all", "host", "spike"],
        default="all",
        help="Run host + spike calling, or host/spike only. Default: all.",
    )
    parser.add_argument(
        "--spike-in-name",
        action="append",
        default=[],
        help=(
            "Spike-in name to call (e.g. lambda). "
            "May be specified multiple times; if omitted in spike mode, auto-discover."
        ),
    )
    parser.add_argument(
        "--reference-file",
        required=True,
        help="Host reference FASTA path passed to methy_caller.",
    )
    parser.add_argument(
        "--spike-reference",
        action="append",
        default=[],
        help=(
            "Spike-in reference in NAME=FASTA format. "
            "May be specified multiple times and is required for spike mode."
        ),
    )
    parser.add_argument(
        "--spike-chromosomes",
        action="append",
        default=[],
        help=(
            "Optional spike chromosome override in NAME=CHR1,CHR2 format. "
            "If omitted, spike mode uses all contigs from each spike reference."
        ),
    )
    parser.add_argument(
        "--chromosomes",
        required=True,
        help="Comma-separated chromosome list for host/spike main outputs.",
    )
    parser.add_argument(
        "--mito-chromosomes",
        default="chrM",
        help="Comma-separated chromosome list used for per-spot host_mito outputs.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=8,
        help="Maximum host spot jobs to run concurrently. Default: 8.",
    )
    parser.add_argument(
        "--min-base-quality",
        type=int,
        default=30,
        help="Forwarded to methy_caller. Default: 30.",
    )
    parser.add_argument(
        "--min-mapping-quality",
        type=int,
        default=1,
        help="Forwarded to methy_caller. Default: 1.",
    )
    parser.add_argument(
        "--sample-size",
        type=int,
        default=None,
        help="Forwarded to methy_caller sample_size.",
    )
    parser.add_argument(
        "--max-depth",
        type=int,
        default=250,
        help="Forwarded to methy_caller max_depth. Default: 250.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=10_000_000,
        help="Forwarded to methy_caller batch_size. Default: 10,000,000.",
    )
    parser.add_argument(
        "--r1-left-trimming",
        type=int,
        default=0,
        help="Forwarded to methy_caller R1 left trimming. Default: 0.",
    )
    parser.add_argument(
        "--r1-right-trimming",
        type=int,
        default=0,
        help="Forwarded to methy_caller R1 right trimming. Default: 0.",
    )
    parser.add_argument(
        "--r2-left-trimming",
        type=int,
        default=0,
        help="Forwarded to methy_caller R2 left trimming. Default: 0.",
    )
    parser.add_argument(
        "--r2-right-trimming",
        type=int,
        default=0,
        help="Forwarded to methy_caller R2 right trimming. Default: 0.",
    )
    parser.add_argument(
        "--caller-script",
        default="scripts/methy_caller.py",
        help="Path to methy_caller script. Default: scripts/methy_caller.py.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands only; do not execute calling.",
    )
    return parser.parse_args()


def quoted(parts: list[str]) -> str:
    return " ".join(shlex.quote(part) for part in parts)


def discover_host_spot_bams(split_dir: Path) -> list[Path]:
    return sorted(split_dir.rglob("*.sorted.bam"))


def discover_spike_bams(pooled_dir: Path) -> dict[str, Path]:
    spikes: dict[str, Path] = {}
    for bam in sorted(pooled_dir.glob("pooled.*.sorted.bam")):
        stem = bam.name[: -len(".sorted.bam")]
        if stem == "pooled.byCB":
            continue
        if not stem.startswith("pooled."):
            continue
        spike_name = stem.split(".", 1)[1]
        if spike_name:
            spikes[spike_name] = bam
    return spikes


def build_caller_command(
    args: argparse.Namespace,
    bam_file: Path,
    output_file: Path,
    chromosomes: str,
    reference_file: str,
) -> list[str]:
    command = [
        sys.executable,
        args.caller_script,
        "-f",
        reference_file,
        "-i",
        str(bam_file),
        "-c",
        chromosomes,
        "-o",
        str(output_file),
        "--min-base-quality",
        str(args.min_base_quality),
        "--min-mapping-quality",
        str(args.min_mapping_quality),
        "--max-depth",
        str(args.max_depth),
        "--batch-size",
        str(args.batch_size),
        "--r1-left-trimming",
        str(args.r1_left_trimming),
        "--r1-right-trimming",
        str(args.r1_right_trimming),
        "--r2-left-trimming",
        str(args.r2_left_trimming),
        "--r2-right-trimming",
        str(args.r2_right_trimming),
    ]
    if args.sample_size is not None:
        command.extend(["--sample-size", str(args.sample_size)])
    return command


def parse_chromosome_csv(chromosome_csv: str) -> list[str]:
    chromosomes = [item.strip() for item in chromosome_csv.split(",") if item.strip()]
    if not chromosomes:
        raise ValueError("chromosome list is empty after parsing")
    return chromosomes


def validate_chromosomes_in_reference(
    reference_file: str,
    chromosome_csv: str,
    label: str,
) -> list[str]:
    requested = parse_chromosome_csv(chromosome_csv)
    with pysam.FastaFile(reference_file) as fasta:
        available = set(fasta.references)
    missing = [chrom for chrom in requested if chrom not in available]
    if missing:
        missing_str = ",".join(missing)
        raise ValueError(
            f"{label} chromosomes not found in reference '{reference_file}': {missing_str}"
        )
    return requested


def load_reference_chromosomes(reference_file: str) -> list[str]:
    with pysam.FastaFile(reference_file) as fasta:
        chromosomes = list(fasta.references)
    if not chromosomes:
        raise ValueError(f"reference has no contigs: {reference_file}")
    return chromosomes


def run_command(command: list[str], dry_run: bool) -> None:
    print(f"[call] command={quoted(command)}")
    if not dry_run:
        subprocess.run(command, check=True)


def run_or_skip(command: list[str], output_path: Path, dry_run: bool) -> bool:
    if output_path.exists():
        print(f"[call] skip_existing_output={output_path}")
        return False
    run_command(command, dry_run)
    return True


def run_host_spot(
    args: argparse.Namespace,
    work_path: Path,
    host_bam: Path,
) -> tuple[Path, Path]:
    split_root = work_path / "split_bams"
    relative = host_bam.relative_to(split_root)
    base_name = relative.name[: -len(".sorted.bam")]
    host_out = work_path / "coverage" / "host" / relative.parent / f"{base_name}.CG.cov"
    mito_out = (
        work_path / "coverage" / "host_mito" / relative.parent / f"{base_name}.CG.cov"
    )
    host_out.parent.mkdir(parents=True, exist_ok=True)
    mito_out.parent.mkdir(parents=True, exist_ok=True)

    host_cmd = build_caller_command(
        args,
        host_bam,
        host_out,
        args.chromosomes,
        args.reference_file,
    )
    mito_cmd = build_caller_command(
        args,
        host_bam,
        mito_out,
        args.mito_chromosomes,
        args.reference_file,
    )
    run_or_skip(host_cmd, host_out, args.dry_run)
    run_or_skip(mito_cmd, mito_out, args.dry_run)
    return host_out, mito_out


def run_host_mode(args: argparse.Namespace, work_path: Path) -> None:
    split_dir = work_path / "split_bams"
    host_bams = discover_host_spot_bams(split_dir)
    if not host_bams:
        raise ValueError(f"no host spot BAMs found under: {split_dir}/**/*.sorted.bam")

    print(f"[call] host_spot_count={len(host_bams)}")
    print(f"[call] host_jobs={args.jobs}")
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(run_host_spot, args, work_path, host_bam): host_bam
            for host_bam in host_bams
        }
        for future in as_completed(futures):
            host_bam = futures[future]
            host_out, mito_out = future.result()
            print(f"[call] host_bam={host_bam}")
            print(f"[call] output_host={host_out}")
            print(f"[call] output_host_mito={mito_out}")


def run_spike_mode(args: argparse.Namespace, work_path: Path) -> None:
    pooled_dir = work_path / "pooled"
    discovered_spikes = discover_spike_bams(pooled_dir)
    if not discovered_spikes:
        raise ValueError(f"no spike-in BAMs found under: {pooled_dir}/pooled.*.sorted.bam")

    spike_reference_map = parse_spike_references(args.spike_reference)
    spike_chromosome_map = parse_spike_chromosomes(args.spike_chromosomes)
    requested_names = [name.strip() for name in args.spike_in_name if name.strip()]
    spike_names = requested_names if requested_names else sorted(discovered_spikes.keys())
    print(f"[call] spike_in_count={len(spike_names)}")
    for spike_name in spike_names:
        spike_bam = discovered_spikes.get(spike_name)
        if spike_bam is None:
            raise ValueError(
                f"requested spike-in BAM missing for '{spike_name}': "
                f"{pooled_dir}/pooled.{spike_name}.sorted.bam"
            )
        spike_out = work_path / "coverage" / f"{spike_name}.CG.cov"
        spike_out.parent.mkdir(parents=True, exist_ok=True)
        spike_reference = spike_reference_map.get(spike_name)
        if spike_reference is None:
            raise ValueError(
                f"missing spike reference for '{spike_name}'. "
                "Provide --spike-reference NAME=FASTA."
            )
        spike_chromosome_csv = spike_chromosome_map.get(spike_name)
        if spike_chromosome_csv is None:
            spike_chromosomes = load_reference_chromosomes(spike_reference)
            spike_chromosome_csv = ",".join(spike_chromosomes)
        else:
            validate_chromosomes_in_reference(
                spike_reference,
                spike_chromosome_csv,
                f"spike '{spike_name}'",
            )
        command = build_caller_command(
            args,
            spike_bam,
            spike_out,
            spike_chromosome_csv,
            spike_reference,
        )
        run_or_skip(command, spike_out, args.dry_run)
        print(f"[call] spike_in={spike_name}")
        print(f"[call] output={spike_out}")


def parse_spike_references(items: list[str]) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for item in items:
        name, sep, reference_path = item.partition("=")
        if not sep or not name.strip() or not reference_path.strip():
            raise ValueError(
                f"invalid --spike-reference value (expected NAME=FASTA): {item}"
            )
        parsed[name.strip()] = reference_path.strip()
    return parsed


def parse_spike_chromosomes(items: list[str]) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for item in items:
        name, sep, chromosome_csv = item.partition("=")
        if not sep or not name.strip() or not chromosome_csv.strip():
            raise ValueError(
                f"invalid --spike-chromosomes value (expected NAME=CHR1,CHR2): {item}"
            )
        parse_chromosome_csv(chromosome_csv.strip())
        parsed[name.strip()] = chromosome_csv.strip()
    return parsed


def main() -> int:
    args = parse_args()
    if args.jobs <= 0:
        raise ValueError("jobs must be > 0")
    if args.min_base_quality < 0:
        raise ValueError("min-base-quality must be >= 0")
    if args.min_mapping_quality < 0:
        raise ValueError("min-mapping-quality must be >= 0")
    if args.max_depth <= 0:
        raise ValueError("max-depth must be > 0")
    if args.batch_size <= 0:
        raise ValueError("batch-size must be > 0")
    if args.sample_size is not None and args.sample_size <= 0:
        raise ValueError("sample-size must be > 0 when provided")
    if args.r1_left_trimming < 0:
        raise ValueError("r1-left-trimming must be >= 0")
    if args.r1_right_trimming < 0:
        raise ValueError("r1-right-trimming must be >= 0")
    if args.r2_left_trimming < 0:
        raise ValueError("r2-left-trimming must be >= 0")
    if args.r2_right_trimming < 0:
        raise ValueError("r2-right-trimming must be >= 0")

    work_path = Path(args.work_path)
    print(f"[call] mode={args.mode}")
    print(f"[call] work_path={work_path}")
    print(f"[call] host_reference={args.reference_file}")
    print(f"[call] chromosomes={args.chromosomes}")
    print(f"[call] mito_chromosomes={args.mito_chromosomes}")
    validate_chromosomes_in_reference(args.reference_file, args.chromosomes, "host")
    validate_chromosomes_in_reference(
        args.reference_file,
        args.mito_chromosomes,
        "host_mito",
    )
    if args.dry_run:
        print("[call] dry_run=1")

    if args.mode in ("all", "host"):
        run_host_mode(args, work_path)
    if args.mode in ("all", "spike"):
        run_spike_mode(args, work_path)

    print("[call] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
