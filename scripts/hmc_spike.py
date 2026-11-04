#!/usr/bin/env python3
"""Score 5hmC oligo conversion (hmC→T) on demux spike-in R1 FASTQ.

Locate a workflow-configured oligo via C→N masking, then score configured
hmC positions (C=retained/protected, T=converted). TAPSbeta should protect
5hmC (low conversion); classical TAPS converts both 5mC and 5hmC.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import time
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
            "Locate a 5hmC oligo on demux spike-in R1 (C→N mask) and compute "
            "site-level hmC→T conversion rates."
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
        help="Oligo sequence to match (Cs are treated as candidate 5hmC sites).",
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
        "--edit-distance",
        type=int,
        default=1,
        help=(
            "Max mismatches on non-N positions of the C→N-masked sequence. "
            "Default: 1."
        ),
    )
    parser.add_argument(
        "--max-reads",
        type=int,
        default=0,
        help="Stop after N reads across all files (0 = all).",
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


def mask_c_to_n(seq: str) -> str:
    return seq.upper().replace("C", "N")


def masked_mismatch_count(obs: str, mask: str) -> int:
    if len(obs) != len(mask):
        raise ValueError("obs and mask length mismatch")
    dist = 0
    for ob, mb in zip(obs, mask):
        if mb == "N":
            continue
        if ob != mb:
            dist += 1
    return dist


def find_masked_spike(
    text: str,
    spike_mask: str,
    start_pos: int,
    end_pos: int,
    max_mismatch: int,
) -> tuple[int | None, str]:
    """Find unique best C→N-masked hit in [start_pos, end_pos).

    Returns (position, status) where status is 'hit', 'none', or 'ambiguous'.
    """
    start_pos = max(0, start_pos)
    end_pos = min(end_pos, len(text))
    pattern_len = len(spike_mask)
    if pattern_len == 0 or start_pos + pattern_len > end_pos:
        return None, "none"

    regex = re.compile("".join("." if b == "N" else re.escape(b) for b in spike_mask))
    window = text[start_pos:end_pos]
    exact_hits = [start_pos + m.start() for m in regex.finditer(window)]
    exact_hits = [p for p in exact_hits if p + pattern_len <= end_pos]
    if len(exact_hits) == 1:
        return exact_hits[0], "hit"
    if len(exact_hits) > 1:
        return None, "ambiguous"

    if max_mismatch <= 0:
        return None, "none"

    best_pos: list[int] = []
    best_dist = max_mismatch + 1
    max_start = end_pos - pattern_len
    for pos in range(start_pos, max_start + 1):
        dist = masked_mismatch_count(text[pos : pos + pattern_len], spike_mask)
        if dist < best_dist:
            best_dist = dist
            best_pos = [pos]
        elif dist == best_dist:
            best_pos.append(pos)
    if best_dist > max_mismatch:
        return None, "none"
    if len(best_pos) != 1:
        return None, "ambiguous"
    return best_pos[0], "hit"


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
    ax.set_ylabel("hmC→T conversion rate")
    ax.set_title(f"{name} 5hmC site conversion")
    ax.set_ylim(0, 1.05)
    for label, rate in zip(labels, conversion_rates):
        if rate is None:
            continue
        ax.text(
            label,
            rate + 0.02,
            f"{rate:.3f}",
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


def main() -> int:
    args = parse_args()
    if not args.name.strip():
        raise ValueError("--name must be non-empty")
    if args.edit_distance < 0:
        raise ValueError("--edit-distance must be >= 0")
    if args.max_reads < 0:
        raise ValueError("--max-reads must be >= 0")

    sequence = args.sequence.upper()
    hmc_positions = list(args.hmc_positions)
    score_positions = list(args.score_positions)
    validate_positions(sequence, hmc_positions, score_positions)

    sequence_mask = mask_c_to_n(sequence)
    r1_paths = sorted(glob(args.r1_glob))
    out_dir = Path(args.out_dir)

    if args.dry_run:
        print(
            f"[hmc_spike] dry-run name={args.name} sequence_len={len(sequence)} "
            f"hmc_positions={hmc_positions} score_positions={score_positions} "
            f"edit_distance={args.edit_distance} r1_glob={args.r1_glob} "
            f"r1_matches={len(r1_paths)} out_dir={out_dir}",
            flush=True,
        )
        return 0

    if not r1_paths:
        raise FileNotFoundError(f"no files matched --r1-glob: {args.r1_glob}")

    out_dir.mkdir(parents=True, exist_ok=True)

    site_retained = [0] * len(hmc_positions)
    site_converted = [0] * len(hmc_positions)
    site_other = [0] * len(hmc_positions)

    n_scanned = 0
    n_hit = 0
    n_none = 0
    n_ambiguous = 0
    n_too_short = 0

    t0 = time.time()
    last_report = t0
    stop = False

    for path in r1_paths:
        if stop:
            break
        with open_text(path, "r") as handle:
            for _head, seq in fastq_iter(handle):
                if args.max_reads and n_scanned >= args.max_reads:
                    stop = True
                    break
                n_scanned += 1

                if len(seq) < len(sequence):
                    n_too_short += 1
                    continue

                hit_pos, status = find_masked_spike(
                    seq,
                    sequence_mask,
                    start_pos=0,
                    end_pos=len(seq),
                    max_mismatch=args.edit_distance,
                )
                if status == "none":
                    n_none += 1
                    continue
                if status == "ambiguous":
                    n_ambiguous += 1
                    continue
                assert hit_pos is not None
                n_hit += 1

                obs = seq[hit_pos : hit_pos + len(sequence)]
                for i, cpos in enumerate(hmc_positions):
                    base = obs[cpos]
                    if base == "C":
                        site_retained[i] += 1
                    elif base == "T":
                        site_converted[i] += 1
                    else:
                        site_other[i] += 1

                if args.progress_interval_seconds > 0:
                    now = time.time()
                    if now - last_report >= args.progress_interval_seconds:
                        elapsed = now - t0
                        rate = n_scanned / elapsed if elapsed > 0 else 0.0
                        print(
                            f"[hmc_spike] scanned={n_scanned} hit={n_hit} "
                            f"none={n_none} ambiguous={n_ambiguous} "
                            f"rate={rate:.0f} reads/s file={Path(path).name}",
                            flush=True,
                        )
                        last_report = now

    site_rows: list[dict[str, int | float | None]] = []
    pos_to_index = {pos: i for i, pos in enumerate(hmc_positions)}
    for i, cpos in enumerate(hmc_positions):
        retained = site_retained[i]
        converted = site_converted[i]
        other = site_other[i]
        denom = retained + converted
        conversion_rate = (converted / denom) if denom > 0 else None
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
        retained = site_retained[idx]
        converted = site_converted[idx]
        score_retained += retained
        score_converted += converted
        denom = retained + converted
        if denom > 0:
            score_rates.append(converted / denom)
    mean_conversion_rate = (
        sum(score_rates) / len(score_rates) if score_rates else None
    )
    pooled_denom = score_retained + score_converted
    pooled_conversion_rate = (
        score_converted / pooled_denom if pooled_denom > 0 else None
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
        "sequence_mask": sequence_mask,
        "hmc_positions": hmc_positions,
        "score_positions": score_positions,
        "edit_distance": args.edit_distance,
        "max_reads": args.max_reads,
        "reads_scanned": n_scanned,
        "reads_matched": n_hit,
        "reads_no_hit": n_none,
        "reads_ambiguous": n_ambiguous,
        "reads_too_short": n_too_short,
        "matched_fraction": (n_hit / n_scanned) if n_scanned else None,
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
        f"[hmc_spike] done name={args.name} scanned={n_scanned} matched={n_hit} "
        f"mean_hmCtoT={mean_conversion_rate} pooled_hmCtoT={pooled_conversion_rate} "
        f"out={out_dir}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
