#!/usr/bin/env python3
"""Split pooled CB-tagged BAM into per-spot BAM files."""

from __future__ import annotations

import argparse
import random
import sys
from collections import defaultdict
from pathlib import Path

import pysam


SMOKE_SPOT_LIMIT = 16


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Split pooled.byCB.bam into per-spot BAM files using the CB tag. "
            "Supports both legacy 16bp concatenated CB and the current X+Y format."
        )
    )
    parser.add_argument(
        "--in-bam",
        required=True,
        help="Input pooled.byCB.bam, typically already sorted by CB.",
    )
    parser.add_argument(
        "--barcodes",
        required=True,
        help="Barcode whitelist TSV (first column is the barcode sequence).",
    )
    parser.add_argument(
        "--out-dir",
        required=True,
        help="Output directory, e.g. /work/split_bams.",
    )
    parser.add_argument(
        "--cb-tag",
        default="CB",
        help="Tag name containing the combined X/Y barcode. Default: CB.",
    )
    parser.add_argument(
        "--threads-read",
        type=int,
        default=1,
        help="Reader threads for pysam.AlignmentFile. Default: 1.",
    )
    parser.add_argument(
        "--threads-write",
        type=int,
        default=1,
        help="Writer threads for pysam.AlignmentFile. Default: 1.",
    )
    parser.add_argument(
        "--smoke",
        action="store_true",
        help=(
            f"Smoke-test mode: randomly select up to {SMOKE_SPOT_LIMIT} observed spots "
            "and emit BAMs only for those spots."
        ),
    )
    return parser.parse_args()


def read_barcodes(tsv_path: str) -> tuple[list[str], dict[str, int]]:
    barcodes: list[str] = []
    with open(tsv_path, encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            barcode = line.split("\t", 1)[0].strip().upper()
            if barcode:
                barcodes.append(barcode)
    return barcodes, {barcode: index for index, barcode in enumerate(barcodes)}


def parse_cb_value(cb_value: str, barcode_len: int) -> tuple[str, str] | None:
    cb_value = cb_value.strip().upper()
    if not cb_value:
        return None
    if "+" in cb_value:
        parts = cb_value.split("+")
        if len(parts) != 2:
            return None
        x_barcode, y_barcode = parts
        if len(x_barcode) != barcode_len or len(y_barcode) != barcode_len:
            return None
        return x_barcode, y_barcode
    if len(cb_value) != barcode_len * 2:
        return None
    return cb_value[:barcode_len], cb_value[barcode_len:]


def resolve_spot(
    record: pysam.AlignedSegment,
    cb_tag: str,
    barcode_to_index: dict[str, int],
    barcode_len: int,
) -> tuple[int, int] | None:
    try:
        cb_value = record.get_tag(cb_tag)
    except KeyError:
        return None
    parsed = parse_cb_value(cb_value, barcode_len)
    if parsed is None:
        return None
    x_barcode, y_barcode = parsed
    x_index = barcode_to_index.get(x_barcode)
    y_index = barcode_to_index.get(y_barcode)
    if x_index is None or y_index is None:
        return None
    return x_index, y_index


def collect_spot_counts(
    in_bam: str,
    cb_tag: str,
    barcode_to_index: dict[str, int],
    barcode_len: int,
    threads_read: int,
) -> tuple[dict[tuple[int, int], int], int, int, int]:
    counts: dict[tuple[int, int], int] = defaultdict(int)
    total_reads = 0
    assigned_reads = 0
    skipped_reads = 0
    with pysam.AlignmentFile(in_bam, "rb", threads=threads_read) as bam_in:
        for record in bam_in.fetch(until_eof=True):
            total_reads += 1
            spot = resolve_spot(record, cb_tag, barcode_to_index, barcode_len)
            if spot is None:
                skipped_reads += 1
                continue
            counts[spot] += 1
            assigned_reads += 1
    return counts, total_reads, assigned_reads, skipped_reads


def choose_smoke_spots(counts: dict[tuple[int, int], int]) -> set[tuple[int, int]]:
    observed_spots = sorted(counts)
    if len(observed_spots) <= SMOKE_SPOT_LIMIT:
        return set(observed_spots)
    return set(random.sample(observed_spots, k=SMOKE_SPOT_LIMIT))


def emit_bams(
    in_bam: str,
    out_dir: Path,
    cb_tag: str,
    barcode_to_index: dict[str, int],
    barcode_len: int,
    threads_read: int,
    threads_write: int,
    selected_spots: set[tuple[int, int]] | None,
) -> dict[tuple[int, int], int]:
    emitted_counts: dict[tuple[int, int], int] = defaultdict(int)
    current_path: Path | None = None
    current_out: pysam.AlignmentFile | None = None

    with pysam.AlignmentFile(in_bam, "rb", threads=threads_read) as bam_in:
        header = bam_in.header
        for record in bam_in.fetch(until_eof=True):
            spot = resolve_spot(record, cb_tag, barcode_to_index, barcode_len)
            if spot is None:
                continue
            if selected_spots is not None and spot not in selected_spots:
                continue

            x_index, y_index = spot
            subdir = out_dir / f"{x_index:02d}"
            out_path = subdir / f"{x_index:02d}_{y_index:02d}.bam"
            if out_path != current_path:
                if current_out is not None:
                    current_out.close()
                subdir.mkdir(parents=True, exist_ok=True)
                current_out = pysam.AlignmentFile(
                    str(out_path),
                    "wb",
                    header=header,
                    threads=threads_write,
                )
                current_path = out_path

            current_out.write(record)
            emitted_counts[spot] += 1

    if current_out is not None:
        current_out.close()
    return emitted_counts


def split_all_bams(
    in_bam: str,
    out_dir: Path,
    cb_tag: str,
    barcode_to_index: dict[str, int],
    barcode_len: int,
    threads_read: int,
    threads_write: int,
) -> tuple[dict[tuple[int, int], int], int, int, int]:
    emitted_counts: dict[tuple[int, int], int] = defaultdict(int)
    total_reads = 0
    assigned_reads = 0
    skipped_reads = 0
    current_path: Path | None = None
    current_out: pysam.AlignmentFile | None = None

    with pysam.AlignmentFile(in_bam, "rb", threads=threads_read) as bam_in:
        header = bam_in.header
        for record in bam_in.fetch(until_eof=True):
            total_reads += 1
            spot = resolve_spot(record, cb_tag, barcode_to_index, barcode_len)
            if spot is None:
                skipped_reads += 1
                continue

            x_index, y_index = spot
            subdir = out_dir / f"{x_index:02d}"
            out_path = subdir / f"{x_index:02d}_{y_index:02d}.bam"
            if out_path != current_path:
                if current_out is not None:
                    current_out.close()
                subdir.mkdir(parents=True, exist_ok=True)
                current_out = pysam.AlignmentFile(
                    str(out_path),
                    "wb",
                    header=header,
                    threads=threads_write,
                )
                current_path = out_path

            current_out.write(record)
            emitted_counts[spot] += 1
            assigned_reads += 1

    if current_out is not None:
        current_out.close()
    return emitted_counts, total_reads, assigned_reads, skipped_reads


def write_stats(out_dir: Path, counts: dict[tuple[int, int], int]) -> Path:
    stats_path = out_dir / "per_spot_read_counts.tsv"
    with stats_path.open("w", encoding="utf-8") as handle:
        handle.write("X_index\tY_index\tspot\treads\n")
        for (x_index, y_index), read_count in sorted(counts.items()):
            handle.write(f"{x_index}\t{y_index}\t{x_index:02d}_{y_index:02d}\t{read_count}\n")
    return stats_path


def main() -> int:
    args = parse_args()
    barcodes, barcode_to_index = read_barcodes(args.barcodes)
    if not barcodes:
        sys.exit("ERROR: empty barcode list.")

    barcode_len = len(barcodes[0])
    if any(len(barcode) != barcode_len for barcode in barcodes):
        sys.exit("ERROR: inconsistent barcode lengths in whitelist.")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pysam.set_verbosity(0)

    selected_spots: set[tuple[int, int]] | None = None
    if args.smoke:
        all_counts, total_reads, assigned_reads, skipped_reads = collect_spot_counts(
            in_bam=args.in_bam,
            cb_tag=args.cb_tag,
            barcode_to_index=barcode_to_index,
            barcode_len=barcode_len,
            threads_read=args.threads_read,
        )
        selected_spots = choose_smoke_spots(all_counts)
        emitted_counts = emit_bams(
            in_bam=args.in_bam,
            out_dir=out_dir,
            cb_tag=args.cb_tag,
            barcode_to_index=barcode_to_index,
            barcode_len=barcode_len,
            threads_read=args.threads_read,
            threads_write=args.threads_write,
            selected_spots=selected_spots,
        )
        print(
            f"[done] smoke_mode=1 observed_spots={len(all_counts)} "
            f"selected_spots={len(selected_spots)}"
        )
    else:
        emitted_counts, total_reads, assigned_reads, skipped_reads = split_all_bams(
            in_bam=args.in_bam,
            out_dir=out_dir,
            cb_tag=args.cb_tag,
            barcode_to_index=barcode_to_index,
            barcode_len=barcode_len,
            threads_read=args.threads_read,
            threads_write=args.threads_write,
        )

    stats_path = write_stats(out_dir, emitted_counts)
    print(f"[done] total={total_reads}, assigned={assigned_reads}, skipped={skipped_reads}")
    print(f"[done] emitted_spots={len(emitted_counts)}")
    print(f"[done] outputs under: {out_dir}")
    print(f"[done] counts: {stats_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
