#!/usr/bin/env python3
"""Run M-bias QC for host and spike-in BAMs."""

from __future__ import annotations

import argparse
import csv
import shlex
import subprocess
from collections import defaultdict
from pathlib import Path

import matplotlib
import pysam

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HOST_SUBSAMPLE_SEED = 11


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
        default="spike",
        help="Run host + spike, or host/spike only. Default: spike.",
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
        "--reference-file",
        default=None,
        help="Host reference FASTA used to locate true CpG sites.",
    )
    parser.add_argument(
        "--spike-reference",
        action="append",
        default=[],
        help=(
            "Spike-in reference in NAME=FASTA format. "
            "May be specified multiple times."
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


TARGET_FLAGS = {99, 147, 83, 163}


def resolve_cycle(
    record: pysam.AlignedSegment,
    query_pos: int,
    read_len: int,
    max_cycle: int,
    right_align_r1: bool,
) -> int:
    if right_align_r1 and record.is_read1:
        # Right-align R1 to the requested cycle window while keeping
        # the first aligned base at cycle 1 for full-length reads.
        return max_cycle - read_len + query_pos + 1
    if record.is_reverse:
        return read_len - query_pos
    return query_pos + 1


def count_mbias_rows(
    bam_path: Path,
    reference_file: str,
    max_cycle: int,
    min_mapping_quality: int,
    right_align_r1: bool = False,
) -> list[dict[str, object]]:
    counts: dict[tuple[str, int], list[int]] = defaultdict(lambda: [0, 0])
    with (
        pysam.AlignmentFile(str(bam_path), "rb") as bam_file,
        pysam.FastaFile(reference_file) as fasta,
    ):
        current_contig = None
        current_sequence = ""
        for record in bam_file.fetch(until_eof=True):
            if record.is_unmapped:
                continue
            if record.is_secondary or record.is_supplementary:
                continue
            if record.flag not in TARGET_FLAGS:
                continue
            if record.mapping_quality < min_mapping_quality:
                continue
            if not record.is_read1 and not record.is_read2:
                continue
            if record.reference_name is None:
                continue
            read_seq = record.query_sequence
            if not read_seq:
                continue
            read_label = "R1" if record.is_read1 else "R2"
            read_len = len(read_seq)
            if read_len < 2:
                continue
            if record.reference_name != current_contig:
                current_contig = record.reference_name
                current_sequence = fasta.fetch(current_contig).upper()
            for query_pos, ref_pos in record.get_aligned_pairs(matches_only=False):
                if query_pos is None or ref_pos is None:
                    continue
                if query_pos + 1 >= read_len:
                    continue
                if ref_pos < 0 or ref_pos + 1 >= len(current_sequence):
                    continue
                if current_sequence[ref_pos : ref_pos + 2] != "CG":
                    continue
                cycle = resolve_cycle(
                    record=record,
                    query_pos=query_pos,
                    read_len=read_len,
                    max_cycle=max_cycle,
                    right_align_r1=right_align_r1,
                )
                if cycle < 1 or cycle > max_cycle:
                    continue
                dinucleotide = (read_seq[query_pos] + read_seq[query_pos + 1]).upper()
                if dinucleotide in ("TG", "CA"):
                    key = (read_label, cycle)
                    counts[key][0] += 1
                elif dinucleotide == "CG":
                    key = (read_label, cycle)
                    counts[key][1] += 1

    rows: list[dict[str, object]] = []
    for (read_label, cycle), (methylated, unmethylated) in sorted(counts.items()):
        coverage = methylated + unmethylated
        methylation_rate = (methylated / coverage) if coverage > 0 else 0.0
        rows.append(
            {
                "read": read_label,
                "cycle": cycle,
                "context": "CG",
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


def write_mbias_png(
    output_path: Path,
    rows: list[dict[str, object]],
    title: str,
    max_cycle: int,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    grouped: dict[tuple[str, str], list[tuple[int, float]]] = defaultdict(list)
    for row in rows:
        read_label = str(row["read"])
        context = str(row["context"])
        cycle = int(row["cycle"])
        rate_percent = float(row["methylation_rate"]) * 100.0
        grouped[(read_label, context)].append((cycle, rate_percent))
    for key in grouped:
        grouped[key].sort(key=lambda item: item[0])

    palette = {
        ("R1", "CG"): "#2ca02c",
        ("R2", "CG"): "#17becf",
    }
    fig, axis = plt.subplots(figsize=(12, 6.4), dpi=120)
    for read_label in ("R1", "R2"):
        for context in ("CG",):
            series = grouped.get((read_label, context))
            if not series:
                continue
            color = palette.get((read_label, context), "#333333")
            axis.plot(
                [cycle for cycle, _ in series],
                [rate for _, rate in series],
                label=f"{read_label}-{context}",
                color=color,
                linewidth=1.8,
            )
    axis.set_title(title)
    axis.set_xlabel("Cycle")
    axis.set_ylabel("Methylation Rate (%)")
    axis.set_xlim(1, max_cycle)
    axis.set_ylim(0, 100)
    axis.grid(alpha=0.25)
    if axis.lines:
        axis.legend(loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(output_path, format="png")
    plt.close(fig)


def run_host_mode(args: argparse.Namespace, work_path: Path, output_dir: Path) -> None:
    if not args.reference_file:
        raise ValueError("host mode requires --reference-file")
    pooled_host = work_path / "pooled" / "pooled.byCB.bam"
    if not pooled_host.exists():
        raise ValueError(f"missing host pooled BAM: {pooled_host}")

    output_dir.mkdir(parents=True, exist_ok=True)
    host_subsampled = output_dir / "host.subsampled.bam"
    host_sorted = output_dir / "host.subsampled.sorted.bam"
    host_index = output_dir / "host.subsampled.sorted.bam.bai"
    host_tsv = output_dir / "host.mbias.tsv"
    host_png = output_dir / "host.mbias.png"

    subsample_arg = build_subsample_fraction(
        HOST_SUBSAMPLE_SEED, args.host_subsample_fraction
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

    if host_tsv.exists() and host_png.exists():
        print(f"[mbias] skip_existing_output={host_tsv}")
        print(f"[mbias] skip_existing_output={host_png}")
    else:
        rows: list[dict[str, object]] = []
        print(f"[mbias] host_input={host_sorted}")
        print(f"[mbias] host_output={host_tsv}")
        print(f"[mbias] host_plot={host_png}")
        if not args.dry_run:
            rows = count_mbias_rows(
                host_sorted,
                reference_file=args.reference_file,
                max_cycle=args.max_cycle,
                min_mapping_quality=args.min_mapping_quality,
                right_align_r1=True,
            )
        if not host_tsv.exists() and not args.dry_run:
            write_mbias_tsv(host_tsv, rows)
        if not host_png.exists() and not args.dry_run:
            write_mbias_png(host_png, rows, "M-bias: host", args.max_cycle)


def run_spike_mode(args: argparse.Namespace, work_path: Path, output_dir: Path) -> None:
    spike_reference_map = parse_spike_references(args.spike_reference)
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
        spike_reference = spike_reference_map.get(spike_name)
        if spike_reference is None:
            raise ValueError(
                f"missing spike reference for '{spike_name}'. "
                "Provide --spike-reference NAME=FASTA."
            )
        output_path = output_dir / f"{spike_name}.mbias.tsv"
        plot_path = output_dir / f"{spike_name}.mbias.png"
        if output_path.exists() and plot_path.exists():
            print(f"[mbias] skip_existing_output={output_path}")
            print(f"[mbias] skip_existing_output={plot_path}")
            continue
        rows: list[dict[str, object]] = []
        print(f"[mbias] spike_in={spike_name}")
        print(f"[mbias] spike_input={spike_bam}")
        print(f"[mbias] output={output_path}")
        print(f"[mbias] plot={plot_path}")
        if not args.dry_run:
            rows = count_mbias_rows(
                spike_bam,
                reference_file=spike_reference,
                max_cycle=args.max_cycle,
                min_mapping_quality=args.min_mapping_quality,
            )
        if not output_path.exists() and not args.dry_run:
            write_mbias_tsv(output_path, rows)
        if not plot_path.exists() and not args.dry_run:
            write_mbias_png(plot_path, rows, f"M-bias: {spike_name}", args.max_cycle)


def main() -> int:
    args = parse_args()
    if args.samtools_threads <= 0:
        raise ValueError("samtools-threads must be > 0")
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
    print(f"[mbias] host_subsample_seed={HOST_SUBSAMPLE_SEED}")
    print(f"[mbias] max_cycle={args.max_cycle}")
    print(f"[mbias] min_mapping_quality={args.min_mapping_quality}")
    if args.reference_file:
        print(f"[mbias] host_reference={args.reference_file}")
    if args.spike_reference:
        print(f"[mbias] spike_reference_count={len(args.spike_reference)}")
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
