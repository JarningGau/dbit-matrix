#!/usr/bin/env python3
"""TAPS v2 demux: C→N-masked methylated-linker locate + optional all-T filter.

Reuses shared barcode/linker helpers from ``scripts/extract_bc.py``. The
methylated-linker is the 5mC-bearing ``insert_left`` anchor on R1 (conversion
QC), not a reference spike-in (lambda/pUC19). Rejected reads still go to the
existing ``*.spike-in.fq.gz`` contract paths.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
import time
from itertools import zip_longest
from pathlib import Path

_SCRIPTS_DIR = Path(__file__).resolve().parent.parent / "scripts"
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from extract_bc import (  # noqa: E402
    annotate_header,
    fastq_iter,
    find_linker2,
    match_whitelist,
    open_text,
    read_whitelist,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Locate barcode-linked linker and C→N-masked insert_left "
            "(methylated-linker) on R1, extract barcodes, optionally require "
            "selected C sites to be T, and output paired FASTQ."
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
        "--linker-bc",
        required=True,
        help="Expected linker sequence between barcode2 and barcode1.",
    )
    parser.add_argument(
        "--insert-left",
        required=True,
        help=(
            "Expected insert-left / methylated-linker anchor sequence on R1 "
            "(equivalent to previous --tn5)."
        ),
    )
    parser.add_argument(
        "--linker-edit-distance",
        type=int,
        default=1,
        help="Max edit distance for linker matching fallback. Default: 1.",
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
        default=60,
        help="Report processing speed every N seconds. Set 0 to disable.",
    )
    parser.add_argument(
        "--gzip-level",
        type=int,
        default=6,
        help="gzip compress level for output FASTQ files (0-9). Default: 6.",
    )
    parser.add_argument(
        "--spike-edit-distance",
        type=int,
        default=1,
        help=(
            "Max mismatches on non-N positions of the C→N-masked insert_left "
            "(methylated-linker). Default: 1."
        ),
    )
    parser.add_argument(
        "--require-c-all-t",
        default="",
        help=(
            "Comma-separated 0-based C positions in insert_left that must all "
            "be T to keep a read (e.g. 3,6,10). Empty = no conversion filter."
        ),
    )
    return parser.parse_args()


def mask_c_to_n(seq: str) -> str:
    return seq.upper().replace("C", "N")


def parse_require_c_all_t(raw: str, insert_left: str) -> list[int]:
    if not raw.strip():
        return []
    positions: list[int] = []
    seen: set[int] = set()
    for part in raw.split(","):
        token = part.strip()
        if not token:
            continue
        pos = int(token)
        if pos < 0 or pos >= len(insert_left):
            raise ValueError(
                f"require-c-all-t position {pos} out of range for "
                f"insert_left length {len(insert_left)}"
            )
        if insert_left[pos] != "C":
            raise ValueError(
                f"require-c-all-t position {pos} is '{insert_left[pos]}' "
                f"in insert_left, expected C"
            )
        if pos not in seen:
            positions.append(pos)
            seen.add(pos)
    if not positions:
        raise ValueError("--require-c-all-t is empty after parsing")
    return positions


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
) -> int | None:
    """Find unique best C→N-masked methylated-linker hit in [start_pos, end_pos)."""
    start_pos = max(0, start_pos)
    end_pos = min(end_pos, len(text))
    pattern_len = len(spike_mask)
    if pattern_len == 0 or start_pos + pattern_len > end_pos:
        return None

    regex = re.compile("".join("." if b == "N" else re.escape(b) for b in spike_mask))
    window = text[start_pos:end_pos]
    exact_hits = [start_pos + m.start() for m in regex.finditer(window)]
    exact_hits = [p for p in exact_hits if p + pattern_len <= end_pos]
    if len(exact_hits) == 1:
        return exact_hits[0]
    if len(exact_hits) > 1:
        return None

    if max_mismatch <= 0:
        return None

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
    if best_dist > max_mismatch or len(best_pos) != 1:
        return None
    return best_pos[0]


def find_insert_left_masked(
    s1_u: str,
    insert_left: str,
    after_pos: int,
    max_mismatch: int,
    window_extra: int = 40,
) -> tuple[int, int] | None:
    """Locate insert_left via C→N mask after barcode1 end (allows linker1 gap)."""
    spike_mask = mask_c_to_n(insert_left)
    window_end = min(len(s1_u), after_pos + len(insert_left) + window_extra)
    spike_start = find_masked_spike(
        s1_u,
        spike_mask,
        start_pos=after_pos,
        end_pos=window_end,
        max_mismatch=max_mismatch,
    )
    if spike_start is None:
        return None
    spike_end = spike_start + len(insert_left)
    if spike_end > len(s1_u):
        return None
    return spike_start, spike_end


def sites_all_t(obs_spike: str, positions: list[int]) -> bool:
    for pos in positions:
        if pos >= len(obs_spike) or obs_spike[pos] != "T":
            return False
    return True


def c_positions_in_insert(insert_left: str) -> list[int]:
    return [i for i, base in enumerate(insert_left) if base == "C"]


def assay_positions_3_6_10(insert_left: str) -> list[int]:
    return [
        p
        for p in (3, 6, 10)
        if p < len(insert_left) and insert_left[p] == "C"
    ]


def score_conversion(
    obs_spike: str, c_positions: list[int]
) -> tuple[int, int, list[str]]:
    """Score T=converted, C=retained; other bases ignored in counts."""
    converted = 0
    retained = 0
    site_bases: list[str] = []
    for pos in c_positions:
        if pos >= len(obs_spike):
            site_bases.append(".")
            continue
        base = obs_spike[pos]
        site_bases.append(base)
        if base == "T":
            converted += 1
        elif base == "C":
            retained += 1
    return converted, retained, site_bases


def main() -> int:
    args = parse_args()
    if args.gzip_level < 0 or args.gzip_level > 9:
        raise ValueError("--gzip-level must be between 0 and 9")
    if args.linker_edit_distance < 0:
        raise ValueError("--linker-edit-distance must be >= 0")
    if args.barcode_hamming_distance < 0:
        raise ValueError("--barcode-hamming-distance must be >= 0")
    if args.spike_edit_distance < 0:
        raise ValueError("--spike-edit-distance must be >= 0")
    linker_bc = args.linker_bc.upper()
    insert_left = args.insert_left.upper()
    require_c_all_t = parse_require_c_all_t(args.require_c_all_t, insert_left)
    all_c_positions = c_positions_in_insert(insert_left)
    assay_sites = assay_positions_3_6_10(insert_left)

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
    scored_fully_converted = 0
    spike_scored_reads = 0
    spike_all_t_3_6_10 = 0
    spike_rate_sum_3_6_10 = 0.0
    spike_rate_n_3_6_10 = 0
    site_converted = [0] * len(all_c_positions)
    site_retained = [0] * len(all_c_positions)
    site_other = [0] * len(all_c_positions)
    reject_counts: dict[str, int] = {
        "short_r1": 0,
        "structure_mismatch": 0,
        "linker2_not_found": 0,
        "tn5_not_found": 0,
        "spikein_c_not_all_t": 0,
        "barcode1_not_in_whitelist": 0,
        "barcode2_not_in_whitelist": 0,
    }

    prefix_len = bc2_len + len(linker_bc) + bc1_len + len(insert_left)
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
        kept_fraction = (kept / total) if total else 0.0
        print(
            f"[extract_bc] progress reads={total} kept={kept} "
            f"keep_rate={kept_fraction:.4f} speed={rate:.0f} reads/s"
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

                linker_bc_window_start = 0
                linker_bc_window_end = min(
                    len(s1_u), bc2_len + len(linker_bc) + 12 + 50
                )
                linker2_pos = find_linker2(
                    s1_u,
                    linker_bc,
                    max_edit_distance=args.linker_edit_distance,
                    window_start=linker_bc_window_start,
                    window_end=linker_bc_window_end,
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
                if bc2_start < 0 or bc1_end > len(s1_u):
                    reject_counts["structure_mismatch"] += 1
                    r1_spike_out.write(f"{h1}\n{s1}\n{p1}\n{q1}\n")
                    r2_spike_out.write(f"{h2}\n{s2}\n{p2}\n{q2}\n")
                    report_progress()
                    continue

                barcode2_obs = s1_u[bc2_start:linker2_start]
                barcode1_obs = s1_u[linker2_end:bc1_end]

                tn5_pos = find_insert_left_masked(
                    s1_u,
                    insert_left,
                    after_pos=bc1_end,
                    max_mismatch=args.spike_edit_distance,
                )
                if tn5_pos is None:
                    reject_counts["structure_mismatch"] += 1
                    reject_counts["tn5_not_found"] += 1
                    r1_spike_out.write(f"{h1}\n{s1}\n{p1}\n{q1}\n")
                    r2_spike_out.write(f"{h2}\n{s2}\n{p2}\n{q2}\n")
                    report_progress()
                    continue
                tn5_start, tn5_end = tn5_pos
                obs_spike = s1_u[tn5_start:tn5_end]

                # Conversion assay: all located methylated-linker reads.
                spike_scored_reads += 1
                _, _, site_bases = score_conversion(obs_spike, all_c_positions)
                for i, base in enumerate(site_bases):
                    if base == "T":
                        site_converted[i] += 1
                    elif base == "C":
                        site_retained[i] += 1
                    else:
                        site_other[i] += 1
                if assay_sites:
                    conv_a, ret_a, _ = score_conversion(obs_spike, assay_sites)
                    denom_a = conv_a + ret_a
                    if denom_a > 0:
                        spike_rate_sum_3_6_10 += conv_a / denom_a
                        spike_rate_n_3_6_10 += 1
                    if sites_all_t(obs_spike, assay_sites):
                        spike_all_t_3_6_10 += 1

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

                if require_c_all_t:
                    if not sites_all_t(obs_spike, require_c_all_t):
                        reject_counts["spikein_c_not_all_t"] += 1
                        r1_spike_out.write(f"{h1}\n{s1}\n{p1}\n{q1}\n")
                        r2_spike_out.write(f"{h2}\n{s2}\n{p2}\n{q2}\n")
                        report_progress()
                        continue
                    scored_fully_converted += 1
                else:
                    if assay_sites and sites_all_t(obs_spike, assay_sites):
                        scored_fully_converted += 1

                kept += 1
                annotated_h1 = annotate_header(h1, barcode1, barcode2)
                annotated_h2 = annotate_header(h2, barcode1, barcode2)
                trimmed_s1 = s1[tn5_end:]
                trimmed_q1 = q1[tn5_end:]
                r1_demux_out.write(f"{annotated_h1}\n{trimmed_s1}\n{p1}\n{trimmed_q1}\n")
                r2_demux_out.write(f"{annotated_h2}\n{s2}\n{p2}\n{q2}\n")
                report_progress()

            report_progress(force=True)

    spike_mCtoT_site_counts: dict[str, dict[str, int]] = {}
    spike_mCtoT_rate_per_c_pos: dict[str, float | None] = {}
    for i, pos in enumerate(all_c_positions):
        key = str(pos)
        converted_t = site_converted[i]
        retained_c = site_retained[i]
        other = site_other[i]
        spike_mCtoT_site_counts[key] = {
            "converted_T": converted_t,
            "retained_C": retained_c,
            "other": other,
        }
        informative = converted_t + retained_c
        spike_mCtoT_rate_per_c_pos[key] = (
            converted_t / informative if informative else None
        )
    mean_rate_3_6_10: float | None = (
        spike_rate_sum_3_6_10 / spike_rate_n_3_6_10 if spike_rate_n_3_6_10 else None
    )
    fraction_all_t_3_6_10: float | None = (
        spike_all_t_3_6_10 / spike_scored_reads if spike_scored_reads else None
    )

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
        "require_c_all_t": require_c_all_t,
        "insert_left_mask_c_to_n": mask_c_to_n(insert_left),
        "spike_edit_distance": args.spike_edit_distance,
        "scored_fully_converted": scored_fully_converted,
        "spike_scored_reads": spike_scored_reads,
        "spike_mCtoT_assay_positions": assay_sites,
        "spike_mCtoT_site_counts": spike_mCtoT_site_counts,
        "spike_mCtoT_rate_per_c_pos": spike_mCtoT_rate_per_c_pos,
        "spike_mCtoT_mean_rate_3_6_10": mean_rate_3_6_10,
        "spike_mCtoT_fraction_all_T_3_6_10": fraction_all_t_3_6_10,
        "reject_counts": reject_counts,
    }
    Path(out_stats).write_text(json.dumps(stats, indent=2), encoding="utf-8")

    print(f"[extract_bc] input_r1={args.r1}")
    print(f"[extract_bc] input_r2={args.r2}")
    print(f"[extract_bc] output_r1_demux={out_r1_demux}")
    print(f"[extract_bc] output_r1_spike_in={out_r1_spike}")
    print(f"[extract_bc] output_r2_demux={out_r2_demux}")
    print(f"[extract_bc] output_r2_spike_in={out_r2_spike}")
    print(f"[extract_bc] require_c_all_t={require_c_all_t}")
    print(f"[extract_bc] stats={out_stats}")
    elapsed = max(time.monotonic() - t0, 1e-9)
    keep_rate = (kept / total) if total else 0.0
    print(
        f"[extract_bc] kept={kept}/{total} keep_rate={keep_rate:.4f} "
        f"spike_in={total - kept} fully_converted={scored_fully_converted} "
        f"avg_speed={total / elapsed:.1f} reads/s"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
