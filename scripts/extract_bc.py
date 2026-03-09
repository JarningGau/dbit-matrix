#!/usr/bin/env python3
"""Extract DBiT barcodes from paired FASTQ and write demultiplexed reads."""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
from typing import Iterator, TextIO


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Parse R1 as barcodeB-linker2-barcodeA-linker1-Tn5-insert, "
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
        help="Output prefix. Writes <prefix>.R1.fq.gz and <prefix>.R2.fq.gz.",
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
    parts[0] = f"{parts[0]}|BC1={bc1}|BC2={bc2}"
    if len(parts) == 1:
        return parts[0]
    return f"{parts[0]} {parts[1]}"


def main() -> int:
    args = parse_args()
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
    out_r1 = f"{prefix}.R1.fq.gz"
    out_r2 = f"{prefix}.R2.fq.gz"
    out_stats = f"{prefix}.stats.json"

    total = 0
    kept = 0
    reject_counts: dict[str, int] = {
        "short_r1": 0,
        "structure_mismatch": 0,
        "barcode1_not_in_whitelist": 0,
        "barcode2_not_in_whitelist": 0,
    }

    prefix_len = bc2_len + len(linker2) + bc1_len + len(linker1) + len(tn5)

    with open_text(args.r1, "r") as r1_in, open_text(args.r2, "r") as r2_in:
        with gzip.open(out_r1, "wt", encoding="utf-8") as r1_out, gzip.open(
            out_r2, "wt", encoding="utf-8"
        ) as r2_out:
            for rec1, rec2 in zip(fastq_iter(r1_in), fastq_iter(r2_in)):
                total += 1
                h1, s1, p1, q1 = rec1
                h2, s2, p2, q2 = rec2
                s1_u = s1.upper()
                min_len = prefix_len
                if len(s1_u) < min_len or len(q1) < min_len:
                    reject_counts["short_r1"] += 1
                    continue

                barcode2 = s1_u[:bc2_len]
                linker2_obs = s1_u[bc2_len : bc2_len + len(linker2)]
                bc1_start = bc2_len + len(linker2)
                barcode1 = s1_u[bc1_start : bc1_start + bc1_len]
                linker1_start = bc1_start + bc1_len
                linker1_obs = s1_u[linker1_start : linker1_start + len(linker1)]
                tn5_start = linker1_start + len(linker1)
                tn5_obs = s1_u[tn5_start : tn5_start + len(tn5)]

                if linker1_obs != linker1 or linker2_obs != linker2 or tn5_obs != tn5:
                    reject_counts["structure_mismatch"] += 1
                    continue
                if barcode1 not in bc1_allow:
                    reject_counts["barcode1_not_in_whitelist"] += 1
                    continue
                if barcode2 not in bc2_allow:
                    reject_counts["barcode2_not_in_whitelist"] += 1
                    continue

                kept += 1
                annotated_h1 = annotate_header(h1, barcode1, barcode2)
                annotated_h2 = annotate_header(h2, barcode1, barcode2)
                trimmed_s1 = s1[prefix_len:]
                trimmed_q1 = q1[prefix_len:]
                r1_out.write(f"{annotated_h1}\n{trimmed_s1}\n{p1}\n{trimmed_q1}\n")
                r2_out.write(f"{annotated_h2}\n{s2}\n{p2}\n{q2}\n")

            try:
                next(fastq_iter(r1_in))
                leftover_r1 = True
            except StopIteration:
                leftover_r1 = False
            try:
                next(fastq_iter(r2_in))
                leftover_r2 = True
            except StopIteration:
                leftover_r2 = False
            if leftover_r1 or leftover_r2:
                raise ValueError("R1 and R2 FASTQ record counts are inconsistent")

    stats = {
        "input_r1": args.r1,
        "input_r2": args.r2,
        "output_r1": out_r1,
        "output_r2": out_r2,
        "total_reads": total,
        "kept_reads": kept,
        "kept_fraction": (kept / total) if total else 0.0,
        "reject_counts": reject_counts,
    }
    Path(out_stats).write_text(json.dumps(stats, indent=2), encoding="utf-8")

    print(f"[extract_bc] input_r1={args.r1}")
    print(f"[extract_bc] input_r2={args.r2}")
    print(f"[extract_bc] output_r1={out_r1}")
    print(f"[extract_bc] output_r2={out_r2}")
    print(f"[extract_bc] stats={out_stats}")
    print(f"[extract_bc] kept={kept}/{total}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
