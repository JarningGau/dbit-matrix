#!/usr/bin/env python3
"""Extract DBiT-RNA barcodes and UMI into clean FASTQ."""

from __future__ import annotations

import argparse
import gzip
import json
import sys
import time
from dataclasses import dataclass
from pathlib import Path


@dataclass
class FastqRecord:
    name: str
    seq: str
    plus: str
    qual: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract BC2-BC1-UMI from R1 and write clean FASTQ pairs."
    )
    parser.add_argument("--r1", required=True, help="Input read1 FASTQ(.gz).")
    parser.add_argument("--r2", required=True, help="Input read2 FASTQ(.gz).")
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="Output prefix, writes <prefix>.R1.clean.fq.gz / <prefix>.R2.clean.fq.gz / <prefix>.stats.json",
    )
    parser.add_argument("--barcode1-whitelist", required=True, help="Barcode1 whitelist path.")
    parser.add_argument("--barcode2-whitelist", required=True, help="Barcode2 whitelist path.")
    parser.add_argument("--linker-bc", required=True, help="Linker sequence between BC2 and BC1.")
    parser.add_argument("--umi-left", required=True, help="Left anchor before UMI.")
    parser.add_argument("--umi-len", required=True, type=int, help="UMI length.")
    parser.add_argument("--barcode-hamming-distance", type=int, default=1)
    parser.add_argument("--linker-edit-distance", type=int, default=1)
    parser.add_argument("--gzip-level", type=int, default=1)
    parser.add_argument(
        "--progress",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Show live progress while processing reads (default: true).",
    )
    parser.add_argument(
        "--progress-interval",
        type=int,
        default=20000,
        help="Refresh interval in reads for non-TTY progress updates.",
    )
    parser.add_argument("--dry-run", action="store_true", help="Print outputs and exit.")
    return parser.parse_args()


def read_whitelist(path: str) -> tuple[set[str], int]:
    barcodes: set[str] = set()
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            barcodes.add(line.split()[0])
    if not barcodes:
        raise ValueError(f"Empty whitelist: {path}")
    lengths = {len(x) for x in barcodes}
    if len(lengths) != 1:
        raise ValueError(f"Whitelist has inconsistent barcode lengths: {path}")
    return barcodes, next(iter(lengths))


def open_text_maybe_gzip(path: str, mode: str):
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t", encoding="utf-8")
    return open(path, mode, encoding="utf-8")


def fastq_iter(path: str):
    with open_text_maybe_gzip(path, "r") as handle:
        while True:
            name = handle.readline()
            if not name:
                return
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not seq or not plus or not qual:
                raise ValueError(f"Incomplete FASTQ record: {path}")
            yield FastqRecord(
                name=name.rstrip("\n"),
                seq=seq.rstrip("\n"),
                plus=plus.rstrip("\n"),
                qual=qual.rstrip("\n"),
            )


def hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(1 for x, y in zip(a, b) if x != y)


def resolve_barcode(raw: str, whitelist: set[str], max_dist: int) -> str | None:
    if raw in whitelist:
        return raw
    best = None
    best_dist = max_dist + 1
    ties = 0
    for candidate in whitelist:
        dist = hamming(raw, candidate)
        if dist < best_dist:
            best = candidate
            best_dist = dist
            ties = 1
        elif dist == best_dist:
            ties += 1
    if best is None or best_dist > max_dist or ties > 1:
        return None
    return best


def locate_linker(seq: str, linker: str, max_dist: int) -> int:
    m = len(linker)
    for idx in range(0, max(0, len(seq) - m + 1)):
        window = seq[idx : idx + m]
        if hamming(window, linker) <= max_dist:
            return idx
    return -1


class ProgressReporter:
    def __init__(self, enabled: bool, interval: int):
        self.enabled = enabled
        self.interval = max(1, interval)
        self.is_tty = sys.stderr.isatty()
        self.start = time.time()
        self.last_render = 0.0
        self.frames = ["-", "\\", "|", "/"]

    def _line(self, total: int, kept: int) -> str:
        elapsed = max(1e-9, time.time() - self.start)
        rate = total / elapsed
        keep_rate = (kept / total * 100.0) if total else 0.0
        frame = self.frames[total % len(self.frames)]
        return (
            f"\r[{frame}] reads={total:,} kept={kept:,} "
            f"keep_rate={keep_rate:5.1f}% rate={rate:,.0f} r/s"
        )

    def update(self, total: int, kept: int) -> None:
        if not self.enabled:
            return
        if self.is_tty:
            now = time.time()
            if now - self.last_render < 0.1:
                return
            self.last_render = now
            sys.stderr.write(self._line(total, kept))
            sys.stderr.flush()
            return
        if total % self.interval == 0:
            sys.stderr.write(self._line(total, kept).lstrip("\r") + "\n")
            sys.stderr.flush()

    def finish(self, total: int, kept: int) -> None:
        if not self.enabled:
            return
        final = self._line(total, kept).lstrip("\r")
        if self.is_tty:
            sys.stderr.write("\r" + final + "\n")
        else:
            if total % self.interval != 0:
                sys.stderr.write(final + "\n")
        sys.stderr.flush()


def main() -> int:
    args = parse_args()
    b1_whitelist, b1_len = read_whitelist(args.barcode1_whitelist)
    b2_whitelist, b2_len = read_whitelist(args.barcode2_whitelist)

    out_prefix = Path(args.output_prefix)
    out_r1 = out_prefix.with_suffix("")
    out_r1 = Path(str(out_prefix) + ".R1.clean.fq.gz")
    out_r2 = Path(str(out_prefix) + ".R2.clean.fq.gz")
    out_stats = Path(str(out_prefix) + ".stats.json")

    if args.dry_run:
        print(out_r1)
        print(out_r2)
        print(out_stats)
        return 0

    out_r1.parent.mkdir(parents=True, exist_ok=True)

    stats = {
        "total_reads": 0,
        "kept_reads": 0,
        "reject_linker_not_found": 0,
        "reject_bad_barcode2": 0,
        "reject_bad_barcode1": 0,
        "reject_umi_left_not_found": 0,
        "reject_umi_too_short": 0,
    }
    progress = ProgressReporter(enabled=args.progress, interval=args.progress_interval)

    with gzip.open(out_r1, "wt", compresslevel=args.gzip_level, encoding="utf-8") as w1, gzip.open(
        out_r2, "wt", compresslevel=args.gzip_level, encoding="utf-8"
    ) as w2:
        for r1, r2 in zip(fastq_iter(args.r1), fastq_iter(args.r2)):
            stats["total_reads"] += 1
            linker_pos = locate_linker(r1.seq, args.linker_bc, args.linker_edit_distance)
            if linker_pos < 0:
                stats["reject_linker_not_found"] += 1
                progress.update(stats["total_reads"], stats["kept_reads"])
                continue

            bc2_raw = r1.seq[max(0, linker_pos - b2_len) : linker_pos]
            if len(bc2_raw) != b2_len:
                stats["reject_bad_barcode2"] += 1
                progress.update(stats["total_reads"], stats["kept_reads"])
                continue
            bc2 = resolve_barcode(bc2_raw, b2_whitelist, args.barcode_hamming_distance)
            if bc2 is None:
                stats["reject_bad_barcode2"] += 1
                progress.update(stats["total_reads"], stats["kept_reads"])
                continue

            bc1_start = linker_pos + len(args.linker_bc)
            bc1_raw = r1.seq[bc1_start : bc1_start + b1_len]
            if len(bc1_raw) != b1_len:
                stats["reject_bad_barcode1"] += 1
                progress.update(stats["total_reads"], stats["kept_reads"])
                continue
            bc1 = resolve_barcode(bc1_raw, b1_whitelist, args.barcode_hamming_distance)
            if bc1 is None:
                stats["reject_bad_barcode1"] += 1
                progress.update(stats["total_reads"], stats["kept_reads"])
                continue

            umi_left_start = bc1_start + b1_len
            umi_left_end = umi_left_start + len(args.umi_left)
            if r1.seq[umi_left_start:umi_left_end] != args.umi_left:
                stats["reject_umi_left_not_found"] += 1
                progress.update(stats["total_reads"], stats["kept_reads"])
                continue

            umi_start = umi_left_end
            umi_end = umi_start + args.umi_len
            umi = r1.seq[umi_start:umi_end]
            if len(umi) != args.umi_len:
                stats["reject_umi_too_short"] += 1
                progress.update(stats["total_reads"], stats["kept_reads"])
                continue

            q_bc2 = "I" * len(bc2)
            q_bc1 = "I" * len(bc1)
            q_umi = r1.qual[umi_start:umi_end]
            if len(q_umi) != args.umi_len:
                q_umi = "I" * args.umi_len
            clean_seq = bc2 + bc1 + umi
            clean_qual = q_bc2 + q_bc1 + q_umi
            clean_name = r1.name

            w1.write(f"{clean_name}\n{clean_seq}\n+\n{clean_qual}\n")
            w2.write(f"{clean_name}\n{r2.seq}\n+\n{r2.qual}\n")
            stats["kept_reads"] += 1
            progress.update(stats["total_reads"], stats["kept_reads"])

    denom = max(1, stats["total_reads"])
    stats["keep_rate"] = stats["kept_reads"] / denom
    progress.finish(stats["total_reads"], stats["kept_reads"])
    out_stats.write_text(json.dumps(stats, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
