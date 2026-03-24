#!/usr/bin/env python3
"""Run STARsolo for DBiT-RNA sample-level matrix output."""

from __future__ import annotations

import argparse
import shutil
import shlex
import subprocess
from pathlib import Path


def load_barcode_index_map(path: str) -> dict[str, int]:
    mapping: dict[str, int] = {}
    with open(path, "r", encoding="utf-8") as handle:
        for index, raw in enumerate(handle):
            line = raw.strip()
            if not line:
                continue
            barcode = line.split()[0]
            if barcode in mapping:
                raise ValueError(f"Duplicate barcode in whitelist: {barcode}")
            mapping[barcode] = index
    if not mapping:
        raise ValueError(f"Empty whitelist: {path}")
    return mapping


def write_barcodes_pos(barcodes_tsv: Path, index_map: dict[str, int]) -> None:
    out_path = barcodes_tsv.with_name("barcodes_pos.tsv")
    with open(barcodes_tsv, "r", encoding="utf-8") as src, open(
        out_path, "w", encoding="utf-8"
    ) as dst:
        dst.write("barcode\tx\ty\n")
        for raw in src:
            barcode = raw.strip()
            if not barcode:
                continue
            if len(barcode) != 16:
                raise ValueError(
                    f"Expected 16bp barcode in {barcodes_tsv}, got {barcode!r} (len={len(barcode)})"
                )
            x_barcode = barcode[:8]
            y_barcode = barcode[8:]
            if x_barcode not in index_map or y_barcode not in index_map:
                raise ValueError(
                    f"Barcode {barcode!r} not found in whitelist split mapping from {barcodes_tsv}"
                )
            dst.write(f"{barcode}\t{index_map[x_barcode]}\t{index_map[y_barcode]}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run STARsolo on demuxed DBiT-RNA reads and write matrix outputs."
    )
    parser.add_argument("--work-path", required=True, help="Sample work directory.")
    parser.add_argument(
        "--star-bin",
        default="STAR",
        help="STAR executable path or command name. Default: STAR.",
    )
    parser.add_argument(
        "--star-genome-dir",
        required=True,
        help="STAR genome index directory.",
    )
    parser.add_argument(
        "--solo-cb-whitelist",
        required=True,
        help="STARsolo --soloCBwhitelist value (e.g. None or whitelist path).",
    )
    parser.add_argument(
        "--barcode-whitelist",
        required=True,
        help="8bp barcode whitelist TSV for barcode->space coordinate mapping.",
    )
    parser.add_argument(
        "--solo-cell-filter",
        default="EmptyDrops_CR",
        help="STARsolo --soloCellFilter value. Default: EmptyDrops_CR.",
    )
    parser.add_argument(
        "--out-tmp-dir",
        help="Optional STAR --outTmpDir (useful on WSL).",
    )
    parser.add_argument("--threads", type=int, default=1, help="STAR runThreadN.")
    parser.add_argument("--dry-run", action="store_true", help="Print outputs and exit.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.threads <= 0:
        raise ValueError("--threads must be > 0")

    work_path = Path(args.work_path)
    demux_dir = work_path / "demux"
    r1_candidates = sorted(demux_dir.glob("*.R1.clean.fq.gz"))
    if not r1_candidates:
        raise FileNotFoundError(f"No demux R1 clean FASTQ found under: {demux_dir}")
    if len(r1_candidates) > 1:
        raise ValueError(
            f"Expected exactly one demux R1 clean FASTQ, found {len(r1_candidates)} under: {demux_dir}"
        )
    in_r1 = r1_candidates[0]
    in_r2 = demux_dir / in_r1.name.replace(".R1.clean.fq.gz", ".R2.clean.fq.gz")
    solo_dir = work_path / "solo"
    star_prefix = solo_dir / "star."
    gene_raw_barcodes = solo_dir / "star.Solo.out" / "Gene" / "raw" / "barcodes.tsv"
    gene_filtered_barcodes = solo_dir / "star.Solo.out" / "Gene" / "filtered" / "barcodes.tsv"

    star_cmd = [
        args.star_bin,
        "--genomeDir",
        args.star_genome_dir,
        "--runMode",
        "alignReads",
        "--runThreadN",
        str(args.threads),
        "--outSAMtype",
        "None",
        "--outFileNamePrefix",
        str(star_prefix),
        "--readFilesIn",
        str(in_r2),
        str(in_r1),
        "--readFilesCommand",
        "zcat",
        "--soloType",
        "CB_UMI_Simple",
        "--soloCBstart",
        "1",
        "--soloCBlen",
        "16",
        "--soloUMIstart",
        "17",
        "--soloUMIlen",
        "10",
        "--soloBarcodeReadLength",
        "0",
        "--soloCBwhitelist",
        args.solo_cb_whitelist,
        "--soloCellFilter",
        args.solo_cell_filter,
        "--soloFeatures",
        "Gene",
    ]
    if args.out_tmp_dir:
        star_cmd.extend(["--outTmpDir", args.out_tmp_dir])

    if args.dry_run:
        print(solo_dir)
        print(gene_raw_barcodes.with_name("barcodes_pos.tsv"))
        print(gene_filtered_barcodes.with_name("barcodes_pos.tsv"))
        print(" ".join(shlex.quote(part) for part in star_cmd))
        return 0

    if not in_r1.exists() or not in_r2.exists():
        raise FileNotFoundError(
            f"Missing clean FASTQ inputs: {in_r1} / {in_r2}"
        )
    if not Path(args.star_genome_dir).exists():
        raise FileNotFoundError(f"Missing STAR genomeDir: {args.star_genome_dir}")
    if args.out_tmp_dir:
        out_tmp_dir = Path(args.out_tmp_dir)
        if out_tmp_dir.exists():
            shutil.rmtree(out_tmp_dir)
    solo_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(star_cmd, check=True)
    if not gene_raw_barcodes.exists() or not gene_filtered_barcodes.exists():
        raise FileNotFoundError(
            "STARsolo outputs missing barcode tables under Gene/raw or Gene/filtered."
        )
    barcode_index_map = load_barcode_index_map(args.barcode_whitelist)
    write_barcodes_pos(gene_raw_barcodes, barcode_index_map)
    write_barcodes_pos(gene_filtered_barcodes, barcode_index_map)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
