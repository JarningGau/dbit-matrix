#!/usr/bin/env python3
"""Extract DBiT barcodes from paired FASTQ and write demultiplexed reads."""

from __future__ import annotations

import argparse
import gzip
import json
import time
from itertools import zip_longest
from pathlib import Path
from typing import Iterator, TextIO

from fuzzysearch import find_near_matches


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Locate linker2/Tn5 on R1, extract barcodeA/barcodeB, "
            "filter by whitelist, and output paired FASTQ."
        )
    )
    parser.add_argument("r1", help="Input chunk R1 FASTQ(.gz).")
    parser.add_argument("r2", help="Input chunk R2 FASTQ(.gz).")
    parser.add_argument(
        "-b1",
        "--barcode1-whitelist",
        required=True,
        help="Whitelist file for barcodeA (one barcode per line).",
    )
    parser.add_argument(
        "-b2",
        "--barcode2-whitelist",
        required=True,
        help="Whitelist file for barcodeB (one barcode per line).",
    )
    parser.add_argument(
        "-o",
        "--output-prefix",
        required=True,
        help=(
            "Output prefix. Writes demux/spike-in FASTQ as "
            "<prefix>.R1.demux.fq.gz, <prefix>.R1.spike-in.fq.gz, "
            "<prefix>.R2.demux.fq.gz, <prefix>.R2.spike-in.fq.gz."
        ),
    )
    parser.add_argument(
        "--linker1",
        default="GTGGCCGATGTTTCG",
        help="Expected linker1 sequence.",
    )
    parser.add_argument(
        "--linker2",
        default="ATCCACGTGCTTGAGAGGCCAGAGCATTCG",
        help="Expected linker2 sequence.",
    )
    parser.add_argument(
        "--tn5",
        default="CATCGGCGTACGACTAGATGTGTATAAGAGACAG",
        help="Expected Tn5 mosaic end sequence on R1.",
    )
    parser.add_argument(
        "--linker-edit-distance",
        type=int,
        default=1,
        help="Max edit distance for linker/Tn5 matching fallback. Default: 1.",
    )
    parser.add_argument(
        "--barcode-hamming-distance",
        type=int,
        default=1,
        help="Max Hamming distance for whitelist fallback. Default: 1.",
    )
    parser.add_argument(
        "--progress-interval-seconds",
        type=float,
        default=1.0,
        help="Report processing speed every N seconds. Set 0 to disable.",
    )
    parser.add_argument(
        "--gzip-level",
        type=int,
        default=6,
        help="gzip compress level for output FASTQ files (0-9). Default: 6.",
    )
    return parser.parse_args()


def open_text(path: str, mode: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t", encoding="utf-8")
    return open(path, mode, encoding="utf-8")


def read_whitelist(path: str) -> set[str]:
    barcodes: set[str] = set()
    with open_text(path, "r") as handle:
        for raw_line in handle:
            line = raw_line.strip().upper()
            if line:
                barcodes.add(line)
    if not barcodes:
        raise ValueError(f"empty barcode whitelist: {path}")
    return barcodes


def fastq_iter(handle: TextIO) -> Iterator[tuple[str, str, str, str]]:
    while True:
        head = handle.readline()
        if not head:
            break
        seq = handle.readline()
        plus = handle.readline()
        qual = handle.readline()
        if not (seq and plus and qual):
            raise ValueError("incomplete FASTQ record detected")
        yield head.rstrip("\n"), seq.rstrip("\n"), plus.rstrip("\n"), qual.rstrip("\n")


def annotate_header(header: str, bc1: str, bc2: str) -> str:
    if not header.startswith("@"):
        raise ValueError(f"invalid FASTQ header: {header}")
    parts = header.split(" ", 1)
    original = parts[0][1:]
    parts[0] = f"@{bc1}+{bc2}:{original}"
    if len(parts) == 1:
        return parts[0]
    return f"{parts[0]} {parts[1]}"


def levenshtein_distance(a: str, b: str, max_dist: int | None = None) -> int:
    if a == b:
        return 0
    if not a:
        return len(b)
    if not b:
        return len(a)

    if max_dist is not None and abs(len(a) - len(b)) > max_dist:
        return max_dist + 1

    if len(a) > len(b):
        a, b = b, a

    prev = list(range(len(a) + 1))
    for i, bc in enumerate(b, start=1):
        curr = [i]
        row_min = i
        for j, ac in enumerate(a, start=1):
            cost = 0 if ac == bc else 1
            curr_val = min(
                prev[j] + 1,
                curr[j - 1] + 1,
                prev[j - 1] + cost,
            )
            curr.append(curr_val)
            row_min = min(row_min, curr_val)
        if max_dist is not None and row_min > max_dist:
            return max_dist + 1
        prev = curr
    return prev[-1]


def find_pattern_layered(
    text: str,
    pattern: str,
    start_pos: int,
    end_pos: int,
    max_edit_distance: int,
) -> int | None:
    if not pattern:
        return None
    if start_pos < 0:
        start_pos = 0
    if end_pos > len(text):
        end_pos = len(text)
    if start_pos >= end_pos:
        return None

    max_start_exact = min(end_pos - len(pattern), len(text) - len(pattern))
    if max_start_exact < start_pos:
        max_start_exact = start_pos - 1

    # 1) Exact fast path.
    if max_start_exact >= start_pos:
        exact = text.find(pattern, start_pos, max_start_exact + len(pattern))
        if exact >= 0:
            return exact

    if max_edit_distance <= 0:
        return None

    # 2) Mismatch-only fallback (substitution).
    if max_start_exact >= start_pos:
        for pos in range(start_pos, max_start_exact + 1):
            if (
                hamming_distance(
                    text[pos : pos + len(pattern)],
                    pattern,
                    stop_at=max_edit_distance,
                )
                <= max_edit_distance
            ):
                return pos

    # 3) Full edit-distance fallback (handles indel) via fuzzysearch.
    near_matches = find_near_matches(
        pattern,
        text[start_pos:end_pos],
        max_l_dist=max_edit_distance,
    )
    if near_matches:
        return start_pos + near_matches[0].start
    return None


def hamming_distance(a: str, b: str, stop_at: int | None = None) -> int:
    if len(a) != len(b):
        raise ValueError("hamming distance requires equal-length strings")
    dist = 0
    for ac, bc in zip(a, b):
        if ac != bc:
            dist += 1
            if stop_at is not None and dist > stop_at:
                return dist
    return dist


def match_whitelist(obs: str, allow_set: set[str], max_hamming: int) -> str | None:
    if obs in allow_set:
        return obs
    if max_hamming <= 0:
        return None

    best: str | None = None
    best_dist = max_hamming + 1
    tie = False
    for candidate in allow_set:
        if len(candidate) != len(obs):
            continue
        dist = hamming_distance(obs, candidate, stop_at=best_dist - 1)
        if dist < best_dist:
            best = candidate
            best_dist = dist
            tie = False
        elif dist == best_dist:
            tie = True
    if best is None or best_dist > max_hamming or tie:
        return None
    return best


def find_linker2(s1_u: str, linker2: str, max_edit_distance: int) -> tuple[int, int] | None:
    # Exact path searches globally to maximize recall.
    pos = s1_u.find(linker2)
    if pos >= 0:
        return pos, pos + len(linker2)
    if max_edit_distance <= 0:
        return None

    # Fuzzy path is limited to a small front window (DBiT linker2 is near read start).
    fuzzy_start = 0
    fuzzy_end = min(len(s1_u), len(linker2) + 24)
    pos = find_pattern_layered(
        s1_u,
        linker2,
        start_pos=fuzzy_start,
        end_pos=fuzzy_end,
        max_edit_distance=max_edit_distance,
    )
    if pos is None:
        return None
    return pos, pos + len(linker2)


def find_tn5(
    s1_u: str, tn5: str, after_pos: int, max_edit_distance: int
) -> tuple[int, int] | None:
    if len(tn5) > 15:
        seed = tn5[-15:]
        offset = len(tn5) - len(seed)
    else:
        seed = tn5
        offset = 0
    seed_start_min = after_pos + offset
    if seed_start_min >= len(s1_u):
        return None
    seed_pos = find_pattern_layered(
        s1_u,
        seed,
        start_pos=seed_start_min,
        end_pos=len(s1_u),
        max_edit_distance=max_edit_distance,
    )
    if seed_pos is None:
        return None
    tn5_start = seed_pos - offset
    tn5_end = tn5_start + len(tn5)
    if tn5_start < 0 or tn5_end > len(s1_u):
        return None
    return tn5_start, tn5_end


def main() -> int:
    args = parse_args()
    if args.gzip_level < 0 or args.gzip_level > 9:
        raise ValueError("--gzip-level must be between 0 and 9")
    if args.linker_edit_distance < 0:
        raise ValueError("--linker-edit-distance must be >= 0")
    if args.barcode_hamming_distance < 0:
        raise ValueError("--barcode-hamming-distance must be >= 0")
    linker1 = args.linker1.upper()
    linker2 = args.linker2.upper()
    tn5 = args.tn5.upper()

    bc1_allow = read_whitelist(args.barcode1_whitelist)
    bc2_allow = read_whitelist(args.barcode2_whitelist)
    bc1_len = len(next(iter(bc1_allow)))
    bc2_len = len(next(iter(bc2_allow)))

    if any(len(x) != bc1_len for x in bc1_allow):
        raise ValueError("barcode1 whitelist contains inconsistent barcode lengths")
    if any(len(x) != bc2_len for x in bc2_allow):
        raise ValueError("barcode2 whitelist contains inconsistent barcode lengths")

    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    out_r1_demux = f"{prefix}.R1.demux.fq.gz"
    out_r1_spike = f"{prefix}.R1.spike-in.fq.gz"
    out_r2_demux = f"{prefix}.R2.demux.fq.gz"
    out_r2_spike = f"{prefix}.R2.spike-in.fq.gz"
    out_stats = f"{prefix}.stats.json"

    total = 0
    kept = 0
    reject_counts: dict[str, int] = {
        "short_r1": 0,
        "structure_mismatch": 0,
        "linker2_not_found": 0,
        "linker1_mismatch": 0,
        "tn5_not_found": 0,
        "barcode1_not_in_whitelist": 0,
        "barcode2_not_in_whitelist": 0,
    }

    prefix_len = bc2_len + len(linker2) + bc1_len + len(linker1) + len(tn5)
    t0 = time.monotonic()
    t_last = t0
    reads_last = 0

    def report_progress(force: bool = False) -> None:
        nonlocal t_last, reads_last
        if args.progress_interval_seconds <= 0:
            return
        now = time.monotonic()
        if not force and (now - t_last) < args.progress_interval_seconds:
            return
        dt = now - t_last
        if dt <= 0:
            return
        reads_now = total - reads_last
        rate = reads_now / dt
        print(
            f"[extract_bc] progress reads={total} kept={kept} speed={rate:.0f} reads/s"
        )
        t_last = now
        reads_last = total

    with open_text(args.r1, "r") as r1_in, open_text(args.r2, "r") as r2_in:
        with gzip.open(
            out_r1_demux, "wt", encoding="utf-8", compresslevel=args.gzip_level
        ) as r1_demux_out, gzip.open(
            out_r1_spike, "wt", encoding="utf-8", compresslevel=args.gzip_level
        ) as r1_spike_out, gzip.open(
            out_r2_demux, "wt", encoding="utf-8", compresslevel=args.gzip_level
        ) as r2_demux_out, gzip.open(
            out_r2_spike, "wt", encoding="utf-8", compresslevel=args.gzip_level
        ) as r2_spike_out:
            for rec1, rec2 in zip_longest(fastq_iter(r1_in), fastq_iter(r2_in)):
                if rec1 is None or rec2 is None:
                    raise ValueError("R1 and R2 FASTQ record counts are inconsistent")
                total += 1
                h1, s1, p1, q1 = rec1
                h2, s2, p2, q2 = rec2
                s1_u = s1.upper()
                min_len = prefix_len
                if len(s1_u) < min_len or len(q1) < min_len:
                    reject_counts["short_r1"] += 1
                    r1_spike_out.write(f"{h1}\n{s1}\n{p1}\n{q1}\n")
                    r2_spike_out.write(f"{h2}\n{s2}\n{p2}\n{q2}\n")
                    report_progress()
                    continue

                linker2_pos = find_linker2(
                    s1_u, linker2, max_edit_distance=args.linker_edit_distance
                )
                if linker2_pos is None:
                    reject_counts["structure_mismatch"] += 1
                    reject_counts["linker2_not_found"] += 1
                    r1_spike_out.write(f"{h1}\n{s1}\n{p1}\n{q1}\n")
                    r2_spike_out.write(f"{h2}\n{s2}\n{p2}\n{q2}\n")
                    report_progress()
                    continue
                linker2_start, linker2_end = linker2_pos
                bc2_start = linker2_start - bc2_len
                bc1_end = linker2_end + bc1_len
                linker1_start = bc1_end
                linker1_end = linker1_start + len(linker1)
                if bc2_start < 0 or linker1_end > len(s1_u):
                    reject_counts["structure_mismatch"] += 1
                    reject_counts["linker1_mismatch"] += 1
                    r1_spike_out.write(f"{h1}\n{s1}\n{p1}\n{q1}\n")
                    r2_spike_out.write(f"{h2}\n{s2}\n{p2}\n{q2}\n")
                    report_progress()
                    continue

                barcode2_obs = s1_u[bc2_start:linker2_start]
                barcode1_obs = s1_u[linker2_end:bc1_end]
                linker1_obs = s1_u[linker1_start:linker1_end]
                # linker1 fast path: exact first, then mismatch-only, then full edit.
                linker1_ok = linker1_obs == linker1
                if not linker1_ok and args.linker_edit_distance > 0:
                    if (
                        hamming_distance(
                            linker1_obs, linker1, stop_at=args.linker_edit_distance
                        )
                        <= args.linker_edit_distance
                    ):
                        linker1_ok = True
                    elif (
                        levenshtein_distance(
                            linker1_obs,
                            linker1,
                            max_dist=args.linker_edit_distance,
                        )
                        <= args.linker_edit_distance
                    ):
                        linker1_ok = True
                if not linker1_ok:
                    reject_counts["structure_mismatch"] += 1
                    reject_counts["linker1_mismatch"] += 1
                    r1_spike_out.write(f"{h1}\n{s1}\n{p1}\n{q1}\n")
                    r2_spike_out.write(f"{h2}\n{s2}\n{p2}\n{q2}\n")
                    report_progress()
                    continue

                tn5_pos = find_tn5(
                    s1_u,
                    tn5,
                    after_pos=linker1_end,
                    max_edit_distance=args.linker_edit_distance,
                )
                if tn5_pos is None:
                    reject_counts["structure_mismatch"] += 1
                    reject_counts["tn5_not_found"] += 1
                    r1_spike_out.write(f"{h1}\n{s1}\n{p1}\n{q1}\n")
                    r2_spike_out.write(f"{h2}\n{s2}\n{p2}\n{q2}\n")
                    report_progress()
                    continue
                _tn5_start, tn5_end = tn5_pos

                barcode1 = match_whitelist(
                    barcode1_obs, bc1_allow, args.barcode_hamming_distance
                )
                if barcode1 is None:
                    reject_counts["barcode1_not_in_whitelist"] += 1
                    r1_spike_out.write(f"{h1}\n{s1}\n{p1}\n{q1}\n")
                    r2_spike_out.write(f"{h2}\n{s2}\n{p2}\n{q2}\n")
                    report_progress()
                    continue
                barcode2 = match_whitelist(
                    barcode2_obs, bc2_allow, args.barcode_hamming_distance
                )
                if barcode2 is None:
                    reject_counts["barcode2_not_in_whitelist"] += 1
                    r1_spike_out.write(f"{h1}\n{s1}\n{p1}\n{q1}\n")
                    r2_spike_out.write(f"{h2}\n{s2}\n{p2}\n{q2}\n")
                    report_progress()
                    continue

                kept += 1
                annotated_h1 = annotate_header(h1, barcode1, barcode2)
                annotated_h2 = annotate_header(h2, barcode1, barcode2)
                trimmed_s1 = s1[tn5_end:]
                trimmed_q1 = q1[tn5_end:]
                r1_demux_out.write(f"{annotated_h1}\n{trimmed_s1}\n{p1}\n{trimmed_q1}\n")
                r2_demux_out.write(f"{annotated_h2}\n{s2}\n{p2}\n{q2}\n")
                report_progress()

            report_progress(force=True)

    stats = {
        "input_r1": args.r1,
        "input_r2": args.r2,
        "output_r1_demux": out_r1_demux,
        "output_r1_spike_in": out_r1_spike,
        "output_r2_demux": out_r2_demux,
        "output_r2_spike_in": out_r2_spike,
        "total_reads": total,
        "kept_reads": kept,
        "spike_in_reads": total - kept,
        "kept_fraction": (kept / total) if total else 0.0,
        "reject_counts": reject_counts,
    }
    Path(out_stats).write_text(json.dumps(stats, indent=2), encoding="utf-8")

    print(f"[extract_bc] input_r1={args.r1}")
    print(f"[extract_bc] input_r2={args.r2}")
    print(f"[extract_bc] output_r1_demux={out_r1_demux}")
    print(f"[extract_bc] output_r1_spike_in={out_r1_spike}")
    print(f"[extract_bc] output_r2_demux={out_r2_demux}")
    print(f"[extract_bc] output_r2_spike_in={out_r2_spike}")
    print(f"[extract_bc] stats={out_stats}")
    elapsed = max(time.monotonic() - t0, 1e-9)
    print(
        f"[extract_bc] kept={kept}/{total} spike_in={total - kept} "
        f"avg_speed={total / elapsed:.1f} reads/s"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
