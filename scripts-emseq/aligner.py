#!/usr/bin/env python3
"""Align EMSeq demux/spike-in FASTQ pairs using biscuit."""

from __future__ import annotations

import argparse
import shlex
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run biscuit align on EMSeq demux FASTQ pairs and pipe host alignments to "
            "sinto nametotag to produce CB-tagged BAM files."
        )
    )
    parser.add_argument(
        "--work-path",
        help=(
            "Sample work directory. If provided, scan <work-path>/demux/*.R1.demux.fq.gz "
            "and write outputs to <work-path>/align_shards/<chunk>.cb.bam and "
            "<chunk>.<spike_name>.bam."
        ),
    )
    parser.add_argument(
        "--chunk",
        help="Optional chunk name filter in work-path mode (e.g. 0001).",
    )
    parser.add_argument("--r1", help="Input R1 demux FASTQ(.gz) for single-chunk mode.")
    parser.add_argument("--r2", help="Input R2 demux FASTQ(.gz) for single-chunk mode.")
    parser.add_argument("--output-bam", help="Output BAM path for single-chunk mode.")
    parser.add_argument(
        "--biscuit-reference",
        required=True,
        help="Reference FASTA passed to biscuit align.",
    )
    parser.add_argument(
        "--biscuit-threads",
        type=int,
        default=2,
        help="Thread count passed to biscuit align -@. Default: 2.",
    )
    parser.add_argument(
        "--biscuit-batch-size",
        type=int,
        default=1,
        help="Batch size passed to biscuit align -b. Default: 1.",
    )
    parser.add_argument(
        "--biscuit-bin",
        default="biscuit",
        help="biscuit executable path or command name. Default: biscuit.",
    )
    parser.add_argument(
        "--sinto-bin",
        default="sinto",
        help="sinto executable path or command name. Default: sinto.",
    )
    parser.add_argument(
        "--samtools-bin",
        default="samtools",
        help="samtools executable path or command name. Default: samtools.",
    )
    parser.add_argument(
        "--spike-in-index",
        action="append",
        default=[],
        help=(
            "Spike-in reference in NAME=INDEX format. "
            "May be specified multiple times; execution order follows argument order."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands only; do not execute alignment.",
    )
    return parser.parse_args()


def quoted(parts: list[str]) -> str:
    return " ".join(shlex.quote(part) for part in parts)


def parse_spike_in_index_items(items: list[str]) -> list[tuple[str, str]]:
    parsed: list[tuple[str, str]] = []
    for item in items:
        name, sep, index = item.partition("=")
        if not sep:
            raise ValueError(f"invalid --spike-in-index value (expected NAME=INDEX): {item}")
        name = name.strip()
        index = index.strip()
        if not name or not index:
            raise ValueError(f"invalid --spike-in-index value (expected NAME=INDEX): {item}")
        parsed.append((name, index))
    return parsed


def discover_demux_inputs(work_path: Path) -> list[tuple[str, Path, Path]]:
    demux_dir = work_path / "demux"
    pairs: list[tuple[str, Path, Path]] = []
    for r1 in sorted(demux_dir.glob("*.R1.demux.fq.gz")):
        chunk = r1.name[: -len(".R1.demux.fq.gz")]
        r2 = demux_dir / f"{chunk}.R2.demux.fq.gz"
        if not r2.exists():
            raise ValueError(f"missing paired R2 for chunk '{chunk}': {r2}")
        pairs.append((chunk, r1, r2))
    return pairs


def resolve_jobs(args: argparse.Namespace) -> list[tuple[str, Path, Path]]:
    has_work_mode = args.work_path is not None
    has_single_mode = any(value is not None for value in (args.r1, args.r2, args.output_bam))
    if has_work_mode and has_single_mode:
        raise ValueError("use either --work-path mode or --r1/--r2/--output-bam mode")
    if not has_work_mode and not has_single_mode:
        raise ValueError("missing input mode: provide --work-path or --r1 --r2 --output-bam")

    if has_work_mode:
        work_path = Path(args.work_path)
        pairs = discover_demux_inputs(work_path)
        if not pairs:
            raise ValueError(f"no demux chunks found under: {work_path / 'demux'}/*.R1.demux.fq.gz")
        if args.chunk:
            pairs = [job for job in pairs if job[0] == args.chunk]
            if not pairs:
                raise ValueError(f"requested chunk not found in demux outputs: {args.chunk}")
        return pairs

    if not (args.r1 and args.r2 and args.output_bam):
        raise ValueError("single-chunk mode requires --r1 --r2 --output-bam")
    return [("single", Path(args.r1), Path(args.r2))]


def run_pipeline(first_cmd: list[str], second_cmd: list[str]) -> None:
    first_proc = subprocess.Popen(first_cmd, stdout=subprocess.PIPE)
    if first_proc.stdout is None:
        raise RuntimeError("failed to capture command stdout")
    try:
        second_proc = subprocess.Popen(second_cmd, stdin=first_proc.stdout)
    finally:
        first_proc.stdout.close()

    second_code = second_proc.wait()
    first_code = first_proc.wait()
    if second_code != 0:
        raise subprocess.CalledProcessError(second_code, second_cmd)
    if first_code != 0:
        raise subprocess.CalledProcessError(first_code, first_cmd)


def main() -> int:
    args = parse_args()
    if args.biscuit_threads <= 0:
        raise ValueError("biscuit_threads must be > 0")
    if args.biscuit_batch_size <= 0:
        raise ValueError("biscuit_batch_size must be > 0")

    spike_indexes = parse_spike_in_index_items(args.spike_in_index)
    jobs = resolve_jobs(args)
    print(f"[align_emseq] job_count={len(jobs)}")
    print(f"[align_emseq] spike_in_count={len(spike_indexes)}")

    for index, (chunk, r1, r2) in enumerate(jobs, start=1):
        if args.work_path:
            work_path = Path(args.work_path)
            demux_dir = work_path / "demux"
            align_dir = work_path / "align_shards"
            host_out_bam = align_dir / f"{chunk}.cb.bam"
            spike_r1 = demux_dir / f"{chunk}.R1.spike-in.fq.gz"
            spike_r2 = demux_dir / f"{chunk}.R2.spike-in.fq.gz"
        else:
            host_out_bam = Path(args.output_bam)
            align_dir = host_out_bam.parent
            spike_r1 = None
            spike_r2 = None

        align_dir.mkdir(parents=True, exist_ok=True)
        print(f"[align_emseq] ({index}/{len(jobs)}) chunk={chunk}")

        for spike_name, spike_index in spike_indexes:
            if spike_r1 is None or spike_r2 is None:
                raise ValueError("spike-in alignment requires --work-path mode")
            if not spike_r1.exists() or not spike_r2.exists():
                raise ValueError(
                    f"missing spike-in FASTQ pair for chunk '{chunk}': {spike_r1}, {spike_r2}"
                )
            spike_out_bam = align_dir / f"{chunk}.{spike_name}.bam"
            spike_align_cmd = [
                args.biscuit_bin,
                "align",
                "-@",
                str(args.biscuit_threads),
                "-b",
                str(args.biscuit_batch_size),
                spike_index,
                str(spike_r1),
                str(spike_r2),
            ]
            spike_samtools_cmd = [
                args.samtools_bin,
                "view",
                "-b",
                "-o",
                str(spike_out_bam),
                "-",
            ]
            spike_command = f"{quoted(spike_align_cmd)} | {quoted(spike_samtools_cmd)}"
            print(f"[align_emseq] spike_in={spike_name}")
            print(f"[align_emseq] output_bam={spike_out_bam}")
            print(f"[align_emseq] command={spike_command}")
            if not args.dry_run:
                run_pipeline(spike_align_cmd, spike_samtools_cmd)

        host_align_cmd = [
            args.biscuit_bin,
            "align",
            "-@",
            str(args.biscuit_threads),
            "-b",
            str(args.biscuit_batch_size),
            args.biscuit_reference,
            str(r1),
            str(r2),
        ]
        sinto_cmd = [
            args.sinto_bin,
            "nametotag",
            "-b",
            "-",
            "-O",
            "b",
            "-o",
            str(host_out_bam),
        ]
        host_command = f"{quoted(host_align_cmd)} | {quoted(sinto_cmd)}"
        print("[align_emseq] target=host")
        print(f"[align_emseq] output_bam={host_out_bam}")
        print(f"[align_emseq] command={host_command}")
        if not args.dry_run:
            run_pipeline(host_align_cmd, sinto_cmd)

    print("[align_emseq] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

