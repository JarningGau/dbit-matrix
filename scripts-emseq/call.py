#!/usr/bin/env python3
"""Run EMSeq methylation calling for host spots and spike-ins.

Uses the pooled pysam caller. Does not run biscuit pileup, vcf2bed, mergecg,
or tabix.
"""

from __future__ import annotations

import argparse
import shlex
import subprocess
import sys
from pathlib import Path

import pysam

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from scripts.methy_call_core import bam_index_exists


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Call EMSeq methylation sites from pooled/pooled.host.sorted.bam (host) and pooled "
            "spike-in BAMs under a sample work directory. Host calling writes "
            "coverage/per_spot_read_counts.tsv and supports --smoke. Does not use biscuit."
        )
    )
    parser.add_argument(
        "--work-path",
        required=True,
        help="Sample work directory containing pooled/.",
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
        help="Comma-separated chromosome list for host nuclear calling.",
    )
    parser.add_argument(
        "--mito-chromosomes",
        default="chrM",
        help="Comma-separated chromosome list used for aggregate host_mito output.",
    )
    parser.add_argument(
        "--context-mode",
        choices=["cg", "ch", "both"],
        default="cg",
        help="Emit CpG outputs, CH outputs, or both. Default: cg.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=8,
        help="Interval workers inside methy_caller. Default: 8.",
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
        default=10,
        help="Forwarded to methy_caller. Default: 10.",
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
        default=50000,
        help="Forwarded to methy_caller max_depth. Default: 50000.",
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
        default="scripts-emseq/methy_caller.py",
        help="Path to methy_caller script. Default: scripts-emseq/methy_caller.py.",
    )
    parser.add_argument(
        "--barcodes",
        help="Barcode whitelist TSV (required for host calling).",
    )
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="Forwarded to methy_caller: limit host calling to a small spot subset.",
    )
    parser.add_argument(
        "--force-recount",
        action="store_true",
        help="Forwarded to methy_caller: rewrite per-spot read counts.",
    )
    parser.add_argument(
        "--spot-roster",
        help="Optional spot roster forwarded to methy_caller (counts TSV or one CB per line).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands only; do not execute calling.",
    )
    return parser.parse_args()


def quoted(parts: list[str]) -> str:
    return " ".join(shlex.quote(part) for part in parts)


def discover_spike_bams(pooled_dir: Path) -> dict[str, Path]:
    spikes: dict[str, Path] = {}
    for bam in sorted(pooled_dir.glob("pooled.*.sorted.bam")):
        stem = bam.name[: -len(".sorted.bam")]
        if not stem.startswith("pooled."):
            continue
        spike_name = stem.split(".", 1)[1]
        if spike_name in {"host", "byCB"}:
            continue
        if spike_name:
            spikes[spike_name] = bam
    return spikes


def context_suffixes(context_mode: str) -> list[str]:
    suffixes = []
    if context_mode in {"cg", "both"}:
        suffixes.append("CG")
    if context_mode in {"ch", "both"}:
        suffixes.append("CH")
    return suffixes


def aggregate_outputs(stem: Path, context_mode: str) -> list[Path]:
    return [Path(str(stem) + f".{suffix}.cov") for suffix in context_suffixes(context_mode)]


def build_caller_command(
    args: argparse.Namespace,
    *,
    mode: str,
    bam: Path,
    reference: str,
    chromosomes: str,
    output: Path | None = None,
    out_dir: Path | None = None,
    barcodes: Path | None = None,
) -> list[str]:
    command = [
        sys.executable,
        args.caller_script,
        "--mode",
        mode,
        "--bam",
        str(bam),
        "--reference",
        reference,
        "--chromosomes",
        chromosomes,
        "--context-mode",
        args.context_mode,
        "--min-base-quality",
        str(args.min_base_quality),
        "--min-mapping-quality",
        str(args.min_mapping_quality),
        "--max-depth",
        str(args.max_depth),
        "--batch-size",
        str(args.batch_size),
        "--jobs",
        str(args.jobs),
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
    if mode == "host":
        assert out_dir is not None and barcodes is not None
        command.extend(
            [
                "--out-dir",
                str(out_dir),
                "--barcodes",
                str(barcodes),
                "--counts-tsv",
                str(out_dir.parent / "per_spot_read_counts.tsv"),
            ]
        )
        if args.smoke:
            command.append("--smoke")
        if args.force_recount:
            command.append("--force-recount")
        if args.spot_roster:
            command.extend(["--spot-roster", str(args.spot_roster)])
    else:
        assert output is not None
        command.extend(["--output", str(output)])
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


def run_or_skip_aggregate(
    command: list[str],
    stem: Path,
    context_mode: str,
    dry_run: bool,
) -> bool:
    outputs = aggregate_outputs(stem, context_mode)
    if outputs and all(path.exists() for path in outputs):
        for path in outputs:
            print(f"[call] skip_existing_output={path}")
        return False
    run_command(command, dry_run)
    return True


def ensure_host_mito_bam(work_path: Path) -> Path:
    pooled_host = work_path / "pooled" / "pooled.host.sorted.bam"
    if not pooled_host.is_file():
        raise ValueError(f"missing pooled host BAM for host_mito: {pooled_host}")
    if not bam_index_exists(pooled_host):
        raise ValueError(f"missing BAM index for host_mito input: {pooled_host}")
    print(f"[call] host_mito_bam={pooled_host}")
    return pooled_host


def run_host_mito_aggregate(args: argparse.Namespace, work_path: Path) -> list[Path]:
    host_mito_bam = ensure_host_mito_bam(work_path)
    stem = work_path / "coverage" / "host_mito"
    stem.parent.mkdir(parents=True, exist_ok=True)
    command = build_caller_command(
        args,
        mode="aggregate",
        bam=host_mito_bam,
        reference=args.reference_file,
        chromosomes=args.mito_chromosomes,
        output=stem,
    )
    run_or_skip_aggregate(command, stem, args.context_mode, args.dry_run)
    return aggregate_outputs(stem, args.context_mode)


def run_host_mode(args: argparse.Namespace, work_path: Path) -> None:
    if not args.barcodes:
        raise ValueError("host calling requires --barcodes")
    pooled = work_path / "pooled" / "pooled.host.sorted.bam"
    out_dir = work_path / "coverage" / "host"
    if not pooled.is_file():
        raise ValueError(f"missing pooled host BAM: {pooled}")
    command = build_caller_command(
        args,
        mode="host",
        bam=pooled,
        reference=args.reference_file,
        chromosomes=args.chromosomes,
        out_dir=out_dir,
        barcodes=Path(args.barcodes),
    )
    run_command(command, args.dry_run)
    mito_outputs = run_host_mito_aggregate(args, work_path)
    for mito_out in mito_outputs:
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
        print(f"[call] spike_in={spike_name}")
        stem = work_path / "coverage" / spike_name
        stem.parent.mkdir(parents=True, exist_ok=True)
        command = build_caller_command(
            args,
            mode="aggregate",
            bam=spike_bam,
            reference=spike_reference,
            chromosomes=spike_chromosome_csv,
            output=stem,
        )
        run_or_skip_aggregate(command, stem, args.context_mode, args.dry_run)
        for spike_out in aggregate_outputs(stem, args.context_mode):
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
    if args.context_mode not in {"cg", "ch", "both"}:
        raise ValueError("context-mode must be one of: cg, ch, both")
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
    print(f"[call] context_mode={args.context_mode}")
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
