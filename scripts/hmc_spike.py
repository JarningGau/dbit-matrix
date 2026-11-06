#!/usr/bin/env python3
"""Score 5hmC oligo conversion (hmC→T) on demux spike-in R1 FASTQ.

Locate a workflow-configured oligo with exact backbone match and [CT] at every
C in the oligo sequence, then score configured hmC positions (C=retained/
protected, T=converted). TAPSbeta should protect 5hmC (low conversion);
classical TAPS converts both 5mC and 5hmC.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from glob import glob
from pathlib import Path
from typing import Iterator, TextIO

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_int_list(text: str) -> list[int]:
    parts = [p.strip() for p in text.split(",") if p.strip()]
    if not parts:
        raise argparse.ArgumentTypeError("expected a comma-separated list of integers")
    try:
        return [int(p) for p in parts]
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"invalid integer list: {text!r}"
        ) from exc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Locate a 5hmC oligo on demux spike-in R1 (exact non-C bases; "
            "[CT] at every C) and compute site-level hmC→T conversion rates."
        )
    )
    parser.add_argument(
        "--r1-glob",
        required=True,
        help=(
            "Glob for demux spike-in R1 FASTQ(.gz), e.g. "
            "'work/<sample>/demux/*.R1.spike-in.fq.gz'."
        ),
    )
    parser.add_argument(
        "--out-dir",
        required=True,
        help="Output directory for summary / site table / plot.",
    )
    parser.add_argument(
        "--name",
        required=True,
        help="Spike name used as summary column prefix (e.g. ClaI).",
    )
    parser.add_argument(
        "--sequence",
        required=True,
        help=(
            "Oligo sequence to match (every C must be C or T in the read; "
            "all other bases must match exactly)."
        ),
    )
    parser.add_argument(
        "--hmc-positions",
        type=parse_int_list,
        required=True,
        help="Comma-separated 0-based positions of 5hmC sites within --sequence.",
    )
    parser.add_argument(
        "--score-positions",
        type=parse_int_list,
        required=True,
        help=(
            "Comma-separated 0-based positions used for mean hmC→T "
            "(must be a subset of --hmc-positions)."
        ),
    )
    parser.add_argument(
        "--max-reads",
        type=int,
        default=0,
        help="Stop after N reads across all files (0 = all).",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Parallel workers over input FASTQ files. Default: 1.",
    )
    parser.add_argument(
        "--progress-interval-seconds",
        type=float,
        default=30.0,
        help="Progress report interval. Set 0 to disable.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned inputs/outputs and exit without writing.",
    )
    return parser.parse_args()


def open_text(path: str, mode: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t", encoding="utf-8")
    return open(path, mode, encoding="utf-8")


def fastq_iter(handle: TextIO) -> Iterator[tuple[str, str]]:
    while True:
        head = handle.readline()
        if not head:
            break
        seq = handle.readline()
        plus = handle.readline()
        qual = handle.readline()
        if not (seq and plus and qual):
            raise ValueError("incomplete FASTQ record detected")
        yield head.rstrip("\n"), seq.rstrip("\n").upper()


def compile_hmc_motif_regex(sequence: str) -> re.Pattern[str]:
    """Exact non-C bases; every C in sequence must be C or T."""
    parts: list[str] = []
    for base in sequence.upper():
        if base == "C":
            parts.append("[CT]")
        else:
            parts.append(re.escape(base))
    return re.compile("".join(parts))


def find_hmc_motif(
    text: str,
    sequence: str,
    start_pos: int,
    end_pos: int,
    regex: re.Pattern[str] | None = None,
) -> tuple[int | None, str]:
    """Find unique motif hit in [start_pos, end_pos).

    Returns (position, status) where status is 'hit', 'none', or 'ambiguous'.
    """
    start_pos = max(0, start_pos)
    end_pos = min(end_pos, len(text))
    pattern_len = len(sequence)
    if pattern_len == 0 or start_pos + pattern_len > end_pos:
        return None, "none"

    if regex is None:
        regex = compile_hmc_motif_regex(sequence)

    window = text[start_pos:end_pos]
    hits = [start_pos + m.start() for m in regex.finditer(window)]
    hits = [p for p in hits if p + pattern_len <= end_pos]
    if len(hits) == 1:
        return hits[0], "hit"
    if len(hits) > 1:
        return None, "ambiguous"
    return None, "none"


def write_site_plot(
    name: str,
    positions: list[int],
    conversion_rates: list[float | None],
    out_png: Path,
) -> None:
    labels = [str(p) for p in positions]
    values = [0.0 if r is None else r for r in conversion_rates]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar(labels, values, color="#2c7fb8", edgecolor="white", linewidth=0.4)
    ax.set_xlabel(f"{name} C position (0-based)")
    ax.set_ylabel("hmC→T conversion rate (%)")
    ax.set_title(f"{name} 5hmC site conversion")
    ax.set_ylim(0, 105)
    for label, rate in zip(labels, conversion_rates):
        if rate is None:
            continue
        ax.text(
            label,
            rate + 2.0,
            f"{rate:.1f}",
            ha="center",
            va="bottom",
            fontsize=8,
        )
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def validate_positions(
    sequence: str, hmc_positions: list[int], score_positions: list[int]
) -> None:
    if not sequence:
        raise ValueError("--sequence must be non-empty")
    if not hmc_positions:
        raise ValueError("--hmc-positions must be non-empty")
    if not score_positions:
        raise ValueError("--score-positions must be non-empty")
    seq_len = len(sequence)
    hmc_set = set(hmc_positions)
    if len(hmc_set) != len(hmc_positions):
        raise ValueError("--hmc-positions must not contain duplicates")
    if len(set(score_positions)) != len(score_positions):
        raise ValueError("--score-positions must not contain duplicates")
    for pos in hmc_positions:
        if pos < 0 or pos >= seq_len:
            raise ValueError(f"hmc position {pos} out of range for sequence")
        if sequence[pos] != "C":
            raise ValueError(
                f"hmc position {pos} is {sequence[pos]!r}, expected 'C'"
            )
    for pos in score_positions:
        if pos not in hmc_set:
            raise ValueError(
                f"score position {pos} is not in --hmc-positions"
            )


@dataclass
class ScanCounters:
    site_retained: list[int]
    site_converted: list[int]
    site_other: list[int]
    n_scanned: int = 0
    n_hit: int = 0
    n_none: int = 0
    n_ambiguous: int = 0
    n_too_short: int = 0

    @classmethod
    def zeros(cls, n_sites: int) -> ScanCounters:
        return cls(
            site_retained=[0] * n_sites,
            site_converted=[0] * n_sites,
            site_other=[0] * n_sites,
        )

    def add(self, other: ScanCounters) -> None:
        for index in range(len(self.site_retained)):
            self.site_retained[index] += other.site_retained[index]
            self.site_converted[index] += other.site_converted[index]
            self.site_other[index] += other.site_other[index]
        self.n_scanned += other.n_scanned
        self.n_hit += other.n_hit
        self.n_none += other.n_none
        self.n_ambiguous += other.n_ambiguous
        self.n_too_short += other.n_too_short


def scan_fastq_file(
    path: str,
    sequence: str,
    hmc_positions: list[int],
    max_reads: int,
) -> ScanCounters:
    regex = compile_hmc_motif_regex(sequence)
    counters = ScanCounters.zeros(len(hmc_positions))
    pattern_len = len(sequence)
    with open_text(path, "r") as handle:
        for _head, seq in fastq_iter(handle):
            if max_reads and counters.n_scanned >= max_reads:
                break
            counters.n_scanned += 1
            if len(seq) < pattern_len:
                counters.n_too_short += 1
                continue
            hit_pos, status = find_hmc_motif(
                seq,
                sequence,
                start_pos=0,
                end_pos=len(seq),
                regex=regex,
            )
            if status == "none":
                counters.n_none += 1
                continue
            if status == "ambiguous":
                counters.n_ambiguous += 1
                continue
            assert hit_pos is not None
            counters.n_hit += 1
            obs = seq[hit_pos : hit_pos + pattern_len]
            for index, cpos in enumerate(hmc_positions):
                base = obs[cpos]
                if base == "C":
                    counters.site_retained[index] += 1
                elif base == "T":
                    counters.site_converted[index] += 1
                else:
                    counters.site_other[index] += 1
    return counters


def _scan_fastq_file_task(
    args: tuple[str, str, list[int], int],
) -> ScanCounters:
    path, sequence, hmc_positions, max_reads = args
    return scan_fastq_file(path, sequence, hmc_positions, max_reads)


def main() -> int:
    args = parse_args()
    if not args.name.strip():
        raise ValueError("--name must be non-empty")
    if args.max_reads < 0:
        raise ValueError("--max-reads must be >= 0")
    if args.workers < 1:
        raise ValueError("--workers must be >= 1")

    sequence = args.sequence.upper()
    hmc_positions = list(args.hmc_positions)
    score_positions = list(args.score_positions)
    validate_positions(sequence, hmc_positions, score_positions)

    r1_paths = sorted(glob(args.r1_glob))
    out_dir = Path(args.out_dir)

    if args.dry_run:
        print(
            f"[hmc_spike] dry-run name={args.name} sequence_len={len(sequence)} "
            f"hmc_positions={hmc_positions} score_positions={score_positions} "
            f"workers={args.workers} r1_glob={args.r1_glob} "
            f"r1_matches={len(r1_paths)} out_dir={out_dir}",
            flush=True,
        )
        return 0

    if not r1_paths:
        raise FileNotFoundError(f"no files matched --r1-glob: {args.r1_glob}")

    out_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    totals = ScanCounters.zeros(len(hmc_positions))
    workers = min(args.workers, len(r1_paths))

    if args.max_reads and workers > 1:
        # Preserve global --max-reads semantics by scanning files in order.
        remaining = args.max_reads
        for path in r1_paths:
            if remaining <= 0:
                break
            part = scan_fastq_file(
                path,
                sequence,
                hmc_positions,
                remaining,
            )
            totals.add(part)
            remaining = args.max_reads - totals.n_scanned
    elif workers == 1:
        remaining = args.max_reads
        last_report = t0
        for path in r1_paths:
            if args.max_reads and totals.n_scanned >= args.max_reads:
                break
            part = scan_fastq_file(
                path,
                sequence,
                hmc_positions,
                remaining if args.max_reads else 0,
            )
            totals.add(part)
            if args.max_reads:
                remaining = args.max_reads - totals.n_scanned
            if args.progress_interval_seconds > 0:
                now = time.time()
                if now - last_report >= args.progress_interval_seconds:
                    elapsed = now - t0
                    rate = totals.n_scanned / elapsed if elapsed > 0 else 0.0
                    print(
                        f"[hmc_spike] scanned={totals.n_scanned} hit={totals.n_hit} "
                        f"none={totals.n_none} ambiguous={totals.n_ambiguous} "
                        f"rate={rate:.0f} reads/s file={Path(path).name}",
                        flush=True,
                    )
                    last_report = now
    else:
        task_args = [
            (path, sequence, hmc_positions, 0)
            for path in r1_paths
        ]
        with ProcessPoolExecutor(max_workers=workers) as pool:
            for part in pool.map(_scan_fastq_file_task, task_args):
                totals.add(part)
                if args.progress_interval_seconds > 0:
                    elapsed = time.time() - t0
                    rate = totals.n_scanned / elapsed if elapsed > 0 else 0.0
                    print(
                        f"[hmc_spike] scanned={totals.n_scanned} hit={totals.n_hit} "
                        f"none={totals.n_none} ambiguous={totals.n_ambiguous} "
                        f"rate={rate:.0f} reads/s",
                        flush=True,
                    )

    site_rows: list[dict[str, int | float | None]] = []
    pos_to_index = {pos: i for i, pos in enumerate(hmc_positions)}
    for i, cpos in enumerate(hmc_positions):
        retained = totals.site_retained[i]
        converted = totals.site_converted[i]
        other = totals.site_other[i]
        denom = retained + converted
        conversion_rate = (
            (converted / denom) * 100.0 if denom > 0 else None
        )
        site_rows.append(
            {
                "pos": cpos,
                "retained_C": retained,
                "converted_T": converted,
                "other": other,
                "conversion_rate": conversion_rate,
            }
        )

    score_rates: list[float] = []
    score_retained = 0
    score_converted = 0
    for pos in score_positions:
        idx = pos_to_index[pos]
        retained = totals.site_retained[idx]
        converted = totals.site_converted[idx]
        score_retained += retained
        score_converted += converted
        denom = retained + converted
        if denom > 0:
            score_rates.append((converted / denom) * 100.0)
    mean_conversion_rate = (
        sum(score_rates) / len(score_rates) if score_rates else None
    )
    pooled_denom = score_retained + score_converted
    pooled_conversion_rate = (
        (score_converted / pooled_denom) * 100.0
        if pooled_denom > 0
        else None
    )

    site_tsv = out_dir / "site_conversion.tsv"
    with site_tsv.open("w", encoding="utf-8") as handle:
        handle.write("pos\tretained_C\tconverted_T\tother\tconversion_rate\n")
        for row in site_rows:
            rate = row["conversion_rate"]
            rate_s = "" if rate is None else f"{rate:.6f}"
            handle.write(
                f"{row['pos']}\t{row['retained_C']}\t{row['converted_T']}\t"
                f"{row['other']}\t{rate_s}\n"
            )

    site_png = out_dir / "site_conversion.png"
    write_site_plot(
        args.name,
        hmc_positions,
        [row["conversion_rate"] for row in site_rows],  # type: ignore[misc]
        site_png,
    )

    summary = {
        "name": args.name,
        "r1_glob": args.r1_glob,
        "r1_files": r1_paths,
        "sequence": sequence,
        "hmc_positions": hmc_positions,
        "score_positions": score_positions,
        "max_reads": args.max_reads,
        "workers": workers,
        "reads_scanned": totals.n_scanned,
        "reads_matched": totals.n_hit,
        "reads_no_hit": totals.n_none,
        "reads_ambiguous": totals.n_ambiguous,
        "reads_too_short": totals.n_too_short,
        "matched_fraction": (totals.n_hit / totals.n_scanned) if totals.n_scanned else None,
        "site_counts": {
            str(row["pos"]): {
                "retained_C": row["retained_C"],
                "converted_T": row["converted_T"],
                "other": row["other"],
                "conversion_rate": row["conversion_rate"],
            }
            for row in site_rows
        },
        "mean_conversion_rate": mean_conversion_rate,
        "pooled_conversion_rate": pooled_conversion_rate,
        "pooled_retained_C": score_retained,
        "pooled_converted_T": score_converted,
        "elapsed_seconds": time.time() - t0,
    }
    summary_path = out_dir / "summary.json"
    with summary_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")

    print(
        f"[hmc_spike] done name={args.name} scanned={totals.n_scanned} matched={totals.n_hit} "
        f"mean_hmCtoT={mean_conversion_rate} pooled_hmCtoT={pooled_conversion_rate} "
        f"out={out_dir}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
