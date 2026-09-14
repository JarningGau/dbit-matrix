#!/usr/bin/env python3
"""TAPS methylation caller (pooled host or aggregate spike/mito)."""

from __future__ import annotations

import argparse
import multiprocessing
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import BinaryIO, Dict, List, Optional, Set, Tuple

import pysam

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from scripts.split_bams import parse_cb_value, read_barcodes

TOP_FLAGS: Set[int] = {99, 147}
BOT_FLAGS: Set[int] = {83, 163}
PAIRED_FLAGS: Set[int] = TOP_FLAGS | BOT_FLAGS
CH_FORWARD_CONTEXTS: Set[str] = {"CA", "CC", "CT"}
CH_REVERSE_MAP: Dict[str, str] = {"T": "CA", "G": "CC", "A": "CT"}
CH_REVERSE_NEIGHBOR: Dict[str, str] = {"CA": "T", "CC": "G", "CT": "A"}

WORKER_CHROMOSOME_SEQUENCE: Optional[str] = None
BarcodeMap = Dict[str, int]
CgCounts = Dict[str, Dict[int, List[int]]]
ChCounts = Dict[str, Dict[int, List[object]]]


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="TAPS methylation caller")
    parser.add_argument("--mode", choices=("host", "aggregate"), required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("-c", "--chromosomes", required=True)
    parser.add_argument(
        "--context-mode", choices=("cg", "ch", "both"), default="cg"
    )
    parser.add_argument("--min-base-quality", type=int, default=30)
    parser.add_argument("--min-mapping-quality", type=int, default=10)
    parser.add_argument("--max-depth", type=int, default=50_000)
    parser.add_argument("--batch-size", type=int, default=10_000_000)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--sample-size", type=int, default=None)
    parser.add_argument("--r1-left-trimming", type=int, default=0)
    parser.add_argument("--r1-right-trimming", type=int, default=0)
    parser.add_argument("--r2-left-trimming", type=int, default=0)
    parser.add_argument("--r2-right-trimming", type=int, default=0)
    parser.add_argument("--cb-tag", default="CB")
    parser.add_argument("--barcodes", help="Whitelist TSV for host mode.")
    parser.add_argument(
        "--spot-roster",
        help="split_bams/per_spot_read_counts.tsv for host mode.",
    )
    parser.add_argument("--out-dir", help="Host output directory (coverage/host).")
    parser.add_argument("--output", help="Aggregate CG/CH stem path without suffix.")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args(argv)


def parse_chromosomes(csv: str) -> List[str]:
    chromosomes = [item.strip() for item in csv.split(",") if item.strip()]
    if not chromosomes:
        raise ValueError("chromosome list is empty")
    return chromosomes


def validate_chromosomes(reference: str, chromosomes: List[str]) -> None:
    with pysam.FastaFile(reference) as fasta:
        available = set(fasta.references)
    missing = [chrom for chrom in chromosomes if chrom not in available]
    if missing:
        raise ValueError(
            f"chromosomes not found in reference '{reference}': {','.join(missing)}"
        )


def bam_index_exists(path: Path) -> bool:
    return any(
        candidate.is_file()
        for candidate in (
            Path(str(path) + ".bai"),
            path.with_suffix(".bai"),
            Path(str(path) + ".csi"),
        )
    )


def classify_cg_observation(dinucleotide: str) -> Optional[bool]:
    if dinucleotide in {"TG", "CA"}:
        return True
    if dinucleotide == "CG":
        return False
    return None


def classify_ch_observation(base: str, strand: str) -> Optional[Tuple[int, int]]:
    if strand == "+":
        if base not in {"C", "T"}:
            return None
        methylated = base == "T"
    else:
        if base not in {"G", "A"}:
            return None
        methylated = base == "A"
    return (1, 0) if methylated else (0, 1)


def ch_neighbor_matches(context: str, strand: str, base: str) -> bool:
    observed = base.upper()
    if context == "CC":
        return observed in (("C", "T") if strand == "+" else ("G", "A"))
    expected = context[1] if strand == "+" else CH_REVERSE_NEIGHBOR[context]
    return observed == expected


def classify_reference_position(
    sequence: str, pos: int
) -> Tuple[bool, Optional[Tuple[str, str]]]:
    base = sequence[pos]
    next_base = sequence[pos + 1] if pos + 1 < len(sequence) else ""
    if base == "C":
        dinucleotide = base + next_base
        if dinucleotide == "CG":
            return True, None
        if dinucleotide in CH_FORWARD_CONTEXTS:
            return False, (dinucleotide, "+")
        return False, None
    if base == "G" and pos > 0:
        context = CH_REVERSE_MAP.get(sequence[pos - 1])
        if context is not None:
            return False, (context, "-")
    return False, None


def is_sample_site(is_cpg: bool, ch_target: Optional[Tuple[str, str]], context_mode: str) -> bool:
    if context_mode in {"cg", "both"} and is_cpg:
        return True
    if context_mode in {"ch", "both"} and ch_target is not None:
        return True
    return False


def select_sample_positions(
    sequence: str, context_mode: str, sample_size: Optional[int]
) -> Optional[Set[int]]:
    """First N CG/CH reference sites on this chromosome, or None to keep all."""
    if sample_size is None:
        return None
    positions: List[int] = []
    for pos in range(len(sequence)):
        is_cpg, ch_target = classify_reference_position(sequence, pos)
        if is_sample_site(is_cpg, ch_target, context_mode):
            positions.append(pos)
            if len(positions) >= sample_size:
                break
    return set(positions)


def has_unclear_conversion_strand(record: pysam.AlignedSegment) -> bool:
    try:
        return record.get_tag("YD") == "u"
    except KeyError:
        return False


def is_trimmed(
    record: pysam.AlignedSegment,
    query_pos: int,
    read_len: int,
    r1_left: int,
    r1_right: int,
    r2_left: int,
    r2_right: int,
) -> bool:
    if record.is_read1:
        left_trim, right_trim = r1_left, r1_right
    elif record.is_read2:
        left_trim, right_trim = r2_left, r2_right
    else:
        return False
    if left_trim == 0 and right_trim == 0:
        return False
    if record.is_reverse:
        left_cycle = read_len - query_pos
        right_cycle = query_pos + 1
    else:
        left_cycle = query_pos + 1
        right_cycle = read_len - query_pos
    return left_cycle <= left_trim or right_cycle <= right_trim


def spot_id_from_cb(
    cb: str, barcode_to_index: BarcodeMap, barcode_len: int
) -> Optional[str]:
    parsed = parse_cb_value(cb, barcode_len)
    if parsed is None:
        return None
    x_barcode, y_barcode = parsed
    x_index = barcode_to_index.get(x_barcode)
    y_index = barcode_to_index.get(y_barcode)
    if x_index is None or y_index is None:
        return None
    return f"{x_index:02d}_{y_index:02d}"


def spot_is_allowed(spot_id: Optional[str], allowed: Optional[Set[str]]) -> bool:
    if spot_id is None:
        return False
    if allowed is None:
        return True
    return spot_id in allowed


def load_barcode_whitelist(path: str) -> Tuple[List[str], BarcodeMap, int]:
    whitelist = Path(path)
    if not whitelist.is_file():
        raise ValueError(f"barcode whitelist not found: {whitelist}")
    barcodes, barcode_to_index = read_barcodes(path)
    if not barcodes:
        raise ValueError("barcode whitelist is empty")
    barcode_len = len(barcodes[0])
    if any(len(barcode) != barcode_len for barcode in barcodes):
        raise ValueError("barcode whitelist has unequal lengths")
    return barcodes, barcode_to_index, barcode_len


def load_spot_roster(path: str) -> List[str]:
    roster_path = Path(path)
    if not roster_path.is_file():
        raise ValueError(f"spot roster not found: {roster_path}")
    spots: List[str] = []
    with roster_path.open(encoding="utf-8") as handle:
        header = handle.readline()
        if not header:
            raise ValueError(f"spot roster is empty: {path}")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            spots.append(fields[2])
    if not spots:
        raise ValueError(f"spot roster has no spots: {path}")
    return spots


def host_cov_path(out_dir: Path, spot_id: str, suffix: str) -> Path:
    row = spot_id.split("_", 1)[0]
    return out_dir / row / f"{spot_id}.{suffix}.cov"


def format_cg_line(chrom: str, pos: int, meth: int, unmeth: int) -> str:
    total = meth + unmeth
    pct = round((meth / total) * 100, 2) if total > 0 else 0.0
    return f"{chrom}\t{pos}\t{pos}\t{pct:.2f}\t{meth}\t{unmeth}\n"


def format_ch_line(
    chrom: str, pos: int, meth: int, unmeth: int, context: str, strand: str
) -> str:
    total = meth + unmeth
    pct = round((meth / total) * 100, 2) if total > 0 else 0.0
    return (
        f"{chrom}\t{pos}\t{pos}\t{pct:.2f}\t{meth}\t{unmeth}\t{context}\t{strand}\n"
    )


def effective_batch_size(
    chromosome_length: int, configured_batch_size: int, jobs: int
) -> int:
    return min(configured_batch_size, max(1, chromosome_length // max(jobs, 1)))


def count_cg_column(
    pileup_col,
    dest: Dict[int, List[int]] | CgCounts,
    host: bool,
    cb_tag: str,
    barcode_to_index: Optional[BarcodeMap],
    barcode_len: int,
    allowed_spots: Optional[Set[str]],
    skip_spots: Set[str],
    r1_left: int,
    r1_right: int,
    r2_left: int,
    r2_right: int,
) -> None:
    pos = pileup_col.pos
    for pileup_read in pileup_col.pileups:
        record = pileup_read.alignment
        if record.flag not in PAIRED_FLAGS:
            continue
        if has_unclear_conversion_strand(record):
            continue
        if pileup_read.is_del or pileup_read.is_refskip:
            continue
        qpos = pileup_read.query_position
        if qpos is None:
            continue
        seq = record.query_sequence
        if not seq or qpos + 1 >= len(seq):
            continue
        if is_trimmed(record, qpos, len(seq), r1_left, r1_right, r2_left, r2_right):
            continue
        if is_trimmed(
            record, qpos + 1, len(seq), r1_left, r1_right, r2_left, r2_right
        ):
            continue
        methylated = classify_cg_observation((seq[qpos] + seq[qpos + 1]).upper())
        if methylated is None:
            continue
        key = 0 if methylated else 1
        if not host:
            dest.setdefault(pos, [0, 0])[key] += 1
            continue
        try:
            cb = record.get_tag(cb_tag)
        except KeyError:
            continue
        assert barcode_to_index is not None
        spot = spot_id_from_cb(str(cb), barcode_to_index, barcode_len)
        if not spot_is_allowed(spot, allowed_spots) or spot in skip_spots:
            continue
        assert isinstance(spot, str)
        slot = dest.setdefault(spot, {}).setdefault(pos, [0, 0])
        slot[key] += 1


def count_ch_column(
    pileup_col,
    context: str,
    strand: str,
    dest,
    host: bool,
    cb_tag: str,
    barcode_to_index: Optional[BarcodeMap],
    barcode_len: int,
    allowed_spots: Optional[Set[str]],
    skip_spots: Set[str],
    r1_left: int,
    r1_right: int,
    r2_left: int,
    r2_right: int,
) -> None:
    pos = pileup_col.pos
    strand_flags = TOP_FLAGS if strand == "+" else BOT_FLAGS
    for pileup_read in pileup_col.pileups:
        record = pileup_read.alignment
        if record.flag not in strand_flags:
            continue
        if has_unclear_conversion_strand(record):
            continue
        if pileup_read.is_del or pileup_read.is_refskip:
            continue
        qpos = pileup_read.query_position
        if qpos is None:
            continue
        seq = record.query_sequence
        if not seq or qpos >= len(seq):
            continue
        if is_trimmed(record, qpos, len(seq), r1_left, r1_right, r2_left, r2_right):
            continue
        if strand == "+":
            if qpos + 1 >= len(seq):
                continue
            neighbor_qpos = qpos + 1
            neighbor = seq[neighbor_qpos]
        else:
            if qpos == 0:
                continue
            neighbor_qpos = qpos - 1
            neighbor = seq[neighbor_qpos]
        if is_trimmed(
            record, neighbor_qpos, len(seq), r1_left, r1_right, r2_left, r2_right
        ):
            continue
        if not ch_neighbor_matches(context, strand, neighbor):
            continue
        observation = classify_ch_observation(seq[qpos].upper(), strand)
        if observation is None:
            continue
        meth, unmeth = observation
        if not host:
            values = dest.setdefault((pos, context, strand), [0, 0])
            values[0] += meth
            values[1] += unmeth
            continue
        try:
            cb = record.get_tag(cb_tag)
        except KeyError:
            continue
        assert barcode_to_index is not None
        spot = spot_id_from_cb(str(cb), barcode_to_index, barcode_len)
        if not spot_is_allowed(spot, allowed_spots) or spot in skip_spots:
            continue
        assert isinstance(spot, str)
        slot = dest.setdefault(spot, {}).setdefault(pos, [0, 0, context, strand])
        slot[0] += meth
        slot[1] += unmeth


def process_interval(
    bam_path: str,
    chrom: str,
    batch_index: int,
    start: int,
    stop: int,
    part_dir: str,
    context_mode: str,
    host: bool,
    cb_tag: str,
    barcode_to_index: Optional[BarcodeMap],
    barcode_len: int,
    allowed_spots: Optional[Set[str]],
    skip_spots: Set[str],
    min_base_quality: int,
    min_mapping_quality: int,
    max_depth: int,
    r1_left: int,
    r1_right: int,
    r2_left: int,
    r2_right: int,
    allowed_positions: Optional[Set[int]],
) -> Tuple[int, Dict[str, str]]:
    sequence = WORKER_CHROMOSOME_SEQUENCE
    if sequence is None:
        raise RuntimeError("worker chromosome sequence was not initialized")
    cg_host: CgCounts = defaultdict(dict)
    ch_host: ChCounts = defaultdict(dict)
    cg_agg: Dict[int, List[int]] = defaultdict(lambda: [0, 0])
    ch_agg: Dict[Tuple[int, str, str], List[int]] = defaultdict(lambda: [0, 0])
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for pileup_col in bam.pileup(
            chrom,
            start,
            stop,
            truncate=True,
            ignore_overlaps=True,
            min_base_quality=min_base_quality,
            stepper="samtools",
            redo_baq=False,
            max_depth=max_depth,
            ignore_orphans=True,
            compute_baq=True,
            min_mapping_quality=min_mapping_quality,
        ):
            if allowed_positions is not None and pileup_col.pos not in allowed_positions:
                continue
            is_cpg, ch_target = classify_reference_position(sequence, pileup_col.pos)
            if is_cpg and context_mode in {"cg", "both"}:
                count_cg_column(
                    pileup_col,
                    cg_host if host else cg_agg,
                    host,
                    cb_tag,
                    barcode_to_index,
                    barcode_len,
                    allowed_spots,
                    skip_spots,
                    r1_left,
                    r1_right,
                    r2_left,
                    r2_right,
                )
            if ch_target is not None and context_mode in {"ch", "both"}:
                ctx, strand = ch_target
                count_ch_column(
                    pileup_col,
                    ctx,
                    strand,
                    ch_host if host else ch_agg,
                    host,
                    cb_tag,
                    barcode_to_index,
                    barcode_len,
                    allowed_spots,
                    skip_spots,
                    r1_left,
                    r1_right,
                    r2_left,
                    r2_right,
                )
    part_root = Path(part_dir)
    paths: Dict[str, str] = {}
    if context_mode in {"cg", "both"}:
        cg_path = part_root / f"batch_{batch_index:06d}.CG.part"
        with cg_path.open("w", encoding="utf-8") as handle:
            if host:
                for spot_id, pos_map in cg_host.items():
                    for pos, (meth, unmeth) in pos_map.items():
                        if meth + unmeth <= 0:
                            continue
                        handle.write(
                            format_cg_line(chrom, pos, meth, unmeth).rstrip("\n")
                            + f"\t{spot_id}\n"
                        )
            else:
                for pos, (meth, unmeth) in sorted(cg_agg.items()):
                    if meth + unmeth <= 0:
                        continue
                    handle.write(format_cg_line(chrom, pos, meth, unmeth))
        paths["CG"] = str(cg_path)
    if context_mode in {"ch", "both"}:
        ch_path = part_root / f"batch_{batch_index:06d}.CH.part"
        with ch_path.open("w", encoding="utf-8") as handle:
            if host:
                for spot_id, pos_map in ch_host.items():
                    for pos, values in pos_map.items():
                        meth, unmeth, ctx, strand = values
                        if int(meth) + int(unmeth) <= 0:
                            continue
                        handle.write(
                            format_ch_line(
                                chrom, pos, int(meth), int(unmeth), str(ctx), str(strand)
                            ).rstrip("\n")
                            + f"\t{spot_id}\n"
                        )
            else:
                for (pos, ctx, strand), (meth, unmeth) in sorted(ch_agg.items()):
                    if meth + unmeth <= 0:
                        continue
                    handle.write(format_ch_line(chrom, pos, meth, unmeth, ctx, strand))
        paths["CH"] = str(ch_path)
    return batch_index, paths


def append_parts(
    part_paths: Dict[str, str],
    handles: Dict[str, BinaryIO],
) -> None:
    for context, path_str in part_paths.items():
        path = Path(path_str)
        if path.stat().st_size == 0:
            continue
        with path.open("rb") as part_handle:
            shutil.copyfileobj(part_handle, handles[context], length=16 * 1024 * 1024)


HOST_CG_FIELDS = 7
HOST_CH_FIELDS = 9


def sort_part_by_spot(part_path: Path, field_index: int, sorted_path: Path) -> None:
    """Sort a combined part file by the spot_id column into sorted_path."""
    try:
        subprocess.run(
            [
                "sort",
                "-s",
                "-t",
                "\t",
                "-k",
                f"{field_index},{field_index}",
                "-T",
                str(part_path.parent),
                "-o",
                str(sorted_path),
                str(part_path),
            ],
            check=True,
            env={**os.environ, "LC_ALL": "C"},
        )
    except FileNotFoundError as exc:
        raise ValueError("GNU sort is required for host coverage shuffle") from exc


def replace_atomic(temp_path: Path, final_path: Path) -> None:
    final_path.parent.mkdir(parents=True, exist_ok=True)
    os.replace(temp_path, final_path)


def write_empty_cov(path: Path) -> None:
    temp_path = path.with_name(path.name + f".tmp.{os.getpid()}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temp_path.write_text("", encoding="utf-8")
    os.replace(temp_path, path)


def stream_sorted_parts(
    part_path: Path,
    suffix: str,
    expected_fields: int,
    out_dir: Path,
    allowed_spots: Set[str],
    skip_spots: Set[str],
) -> None:
    if not part_path.is_file() or part_path.stat().st_size == 0:
        return
    sorted_path = part_path.with_name(part_path.name + f".sorted.{os.getpid()}")
    current_spot: Optional[str] = None
    current_temp: Optional[Path] = None
    current_handle = None
    try:
        sort_part_by_spot(part_path, expected_fields, sorted_path)

        def close_current() -> None:
            nonlocal current_spot, current_temp, current_handle
            if current_handle is not None:
                current_handle.close()
                current_handle = None
            if current_spot is not None and current_temp is not None:
                replace_atomic(current_temp, host_cov_path(out_dir, current_spot, suffix))
            current_spot = None
            current_temp = None

        with sorted_path.open(encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) != expected_fields:
                    continue
                spot_id = fields[-1]
                if spot_id not in allowed_spots or spot_id in skip_spots:
                    continue
                if spot_id != current_spot:
                    close_current()
                    current_spot = spot_id
                    final_path = host_cov_path(out_dir, spot_id, suffix)
                    current_temp = final_path.with_name(
                        final_path.name + f".tmp.{os.getpid()}"
                    )
                    final_path.parent.mkdir(parents=True, exist_ok=True)
                    current_handle = current_temp.open("w", encoding="utf-8")
                assert current_handle is not None
                current_handle.write("\t".join(fields[:-1]) + "\n")
        close_current()
    finally:
        if current_handle is not None:
            current_handle.close()
        if current_temp is not None:
            try:
                current_temp.unlink()
            except FileNotFoundError:
                pass
        try:
            sorted_path.unlink()
        except FileNotFoundError:
            pass


def shuffle_host_parts(
    combined_cg: Path,
    combined_ch: Optional[Path],
    out_dir: Path,
    roster: List[str],
    skip_spots: Set[str],
    context_mode: str,
) -> None:
    allowed_spots = set(roster)
    if context_mode in {"cg", "both"}:
        stream_sorted_parts(
            combined_cg, "CG", HOST_CG_FIELDS, out_dir, allowed_spots, skip_spots
        )
    if context_mode in {"ch", "both"} and combined_ch is not None:
        stream_sorted_parts(
            combined_ch, "CH", HOST_CH_FIELDS, out_dir, allowed_spots, skip_spots
        )
    suffixes = []
    if context_mode in {"cg", "both"}:
        suffixes.append("CG")
    if context_mode in {"ch", "both"}:
        suffixes.append("CH")
    for spot_id in roster:
        if spot_id in skip_spots:
            continue
        for suffix in suffixes:
            path = host_cov_path(out_dir, spot_id, suffix)
            if not path.exists():
                write_empty_cov(path)


def process_chromosome(
    args: argparse.Namespace,
    chrom: str,
    part_root: Path,
    host: bool,
    barcode_to_index: Optional[BarcodeMap],
    barcode_len: int,
    allowed_spots: Optional[Set[str]],
    skip_spots: Set[str],
    combined_handles: Dict[str, BinaryIO],
) -> None:
    global WORKER_CHROMOSOME_SEQUENCE
    with pysam.FastaFile(args.reference) as fasta:
        sequence = fasta.fetch(chrom).upper()
    WORKER_CHROMOSOME_SEQUENCE = sequence
    chrom_len = len(sequence)
    allowed_positions = select_sample_positions(
        sequence, args.context_mode, args.sample_size
    )
    batch = effective_batch_size(chrom_len, args.batch_size, args.jobs)
    intervals = [
        (index, start, min(start + batch, chrom_len))
        for index, start in enumerate(range(0, chrom_len, batch))
    ]
    chrom_dir = part_root / chrom
    chrom_dir.mkdir(parents=True, exist_ok=True)
    interval_args = [
        (
            args.bam,
            chrom,
            batch_index,
            start,
            stop,
            str(chrom_dir),
            args.context_mode,
            host,
            args.cb_tag,
            barcode_to_index,
            barcode_len,
            allowed_spots,
            skip_spots,
            args.min_base_quality,
            args.min_mapping_quality,
            args.max_depth,
            args.r1_left_trimming,
            args.r1_right_trimming,
            args.r2_left_trimming,
            args.r2_right_trimming,
            allowed_positions,
        )
        for batch_index, start, stop in intervals
    ]
    results = []
    if args.jobs == 1 or len(intervals) == 1:
        for worker_args in interval_args:
            results.append(process_interval(*worker_args))
    else:
        fork_context = multiprocessing.get_context("fork")
        with ProcessPoolExecutor(
            max_workers=min(args.jobs, len(intervals)),
            mp_context=fork_context,
        ) as executor:
            futures = [
                executor.submit(process_interval, *worker_args)
                for worker_args in interval_args
            ]
            for future in as_completed(futures):
                results.append(future.result())
    for _, part_paths in sorted(results):
        append_parts(part_paths, combined_handles)
    shutil.rmtree(chrom_dir)
    WORKER_CHROMOSOME_SEQUENCE = None


def existing_skip_spots(
    out_dir: Path, roster: List[str], context_mode: str
) -> Set[str]:
    skip: Set[str] = set()
    suffixes = []
    if context_mode in {"cg", "both"}:
        suffixes.append("CG")
    if context_mode in {"ch", "both"}:
        suffixes.append("CH")
    for spot_id in roster:
        if all(host_cov_path(out_dir, spot_id, suffix).is_file() for suffix in suffixes):
            skip.add(spot_id)
    return skip


def validate_runtime(args: argparse.Namespace) -> None:
    if args.jobs <= 0:
        raise ValueError("--jobs must be > 0")
    if args.batch_size <= 0 or args.max_depth <= 0:
        raise ValueError("--batch-size and --max-depth must be > 0")
    if args.min_base_quality < 0 or args.min_mapping_quality < 0:
        raise ValueError("quality filters must be >= 0")
    for name in (
        "r1_left_trimming",
        "r1_right_trimming",
        "r2_left_trimming",
        "r2_right_trimming",
    ):
        if getattr(args, name) < 0:
            raise ValueError(f"--{name.replace('_', '-')} must be >= 0")
    if args.sample_size is not None and args.sample_size <= 0:
        raise ValueError("--sample-size must be > 0 when provided")
    if args.mode == "host":
        if not args.barcodes or not args.spot_roster or not args.out_dir:
            raise ValueError("host mode requires --barcodes, --spot-roster, and --out-dir")
    elif not args.output:
        raise ValueError("aggregate mode requires --output")


def print_resolved_parameters(args: argparse.Namespace) -> None:
    print(f"[methy_caller] mode={args.mode}")
    print(f"[methy_caller] bam={args.bam}")
    print(f"[methy_caller] reference={args.reference}")
    print(f"[methy_caller] chromosomes={args.chromosomes}")
    print(f"[methy_caller] context_mode={args.context_mode}")
    print(f"[methy_caller] jobs={args.jobs}")
    print(f"[methy_caller] r1_left_trimming={args.r1_left_trimming}")
    print(f"[methy_caller] r1_right_trimming={args.r1_right_trimming}")
    print(f"[methy_caller] r2_left_trimming={args.r2_left_trimming}")
    print(f"[methy_caller] r2_right_trimming={args.r2_right_trimming}")
    if args.mode == "host":
        print(f"[methy_caller] barcodes={args.barcodes}")
        print(f"[methy_caller] spot_roster={args.spot_roster}")
        print(f"[methy_caller] out_dir={args.out_dir}")
    else:
        print(f"[methy_caller] output={args.output}")


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    print_resolved_parameters(args)
    if args.dry_run:
        print("[methy_caller] dry_run=1")
        return 0
    try:
        validate_runtime(args)
        if args.mode == "host":
            if not Path(args.barcodes).is_file():
                raise ValueError(f"barcode whitelist not found: {args.barcodes}")
            if not Path(args.spot_roster).is_file():
                raise ValueError(f"spot roster not found: {args.spot_roster}")
        bam_path = Path(args.bam)
        if not bam_path.is_file():
            raise ValueError(f"BAM not found: {bam_path}")
        if not bam_index_exists(bam_path):
            raise ValueError(f"BAM index not found: {bam_path}")
        if not Path(args.reference).is_file() or not Path(str(args.reference) + ".fai").is_file():
            raise ValueError(f"reference FASTA/.fai not found: {args.reference}")
        chromosomes = parse_chromosomes(args.chromosomes)
        validate_chromosomes(args.reference, chromosomes)
    except (ValueError, OSError) as exc:
        print(f"[methy_caller] error: {exc}", file=sys.stderr)
        return 1

    host = args.mode == "host"
    barcode_to_index: Optional[BarcodeMap] = None
    barcode_len = 0
    allowed_spots: Optional[Set[str]] = None
    skip_spots: Set[str] = set()
    roster: List[str] = []
    if host:
        try:
            _, barcode_to_index, barcode_len = load_barcode_whitelist(args.barcodes)
            roster = load_spot_roster(args.spot_roster)
        except (ValueError, OSError) as exc:
            print(f"[methy_caller] error: {exc}", file=sys.stderr)
            return 1
        allowed_spots = set(roster)
        skip_spots = existing_skip_spots(Path(args.out_dir), roster, args.context_mode)
        print(f"[methy_caller] roster_spots={len(roster)} skip_spots={len(skip_spots)}")

    suffixes = []
    if args.context_mode in {"cg", "both"}:
        suffixes.append("CG")
    if args.context_mode in {"ch", "both"}:
        suffixes.append("CH")

    if host:
        out_dir = Path(args.out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        combined_dir = out_dir
    else:
        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        combined_dir = out_path.parent
        if all(Path(str(args.output) + f".{suffix}.cov").is_file() for suffix in suffixes):
            print("[methy_caller] skip_existing_aggregate=1")
            return 0

    temp_suffix = f".parts.{os.getpid()}"
    combined_temps = {
        suffix: combined_dir / f".call.{suffix}{temp_suffix}"
        for suffix in suffixes
    }
    handles = {suffix: path.open("wb") for suffix, path in combined_temps.items()}
    try:
        with tempfile.TemporaryDirectory(prefix=".methy-call-parts.", dir=combined_dir) as part_dir:
            for chrom in chromosomes:
                process_chromosome(
                    args,
                    chrom,
                    Path(part_dir),
                    host,
                    barcode_to_index,
                    barcode_len,
                    allowed_spots,
                    skip_spots,
                    handles,
                )
        for handle in handles.values():
            handle.close()
        handles.clear()
        if host:
            shuffle_host_parts(
                combined_temps.get("CG", Path("/dev/null")),
                combined_temps.get("CH"),
                Path(args.out_dir),
                roster,
                skip_spots,
                args.context_mode,
            )
        else:
            for suffix, temp_path in combined_temps.items():
                os.replace(temp_path, Path(str(args.output) + f".{suffix}.cov"))
    except (ValueError, OSError, subprocess.CalledProcessError) as exc:
        print(f"[methy_caller] error: {exc}", file=sys.stderr)
        return 1
    finally:
        for handle in handles.values():
            handle.close()
        for temp_path in combined_temps.values():
            try:
                temp_path.unlink()
            except FileNotFoundError:
                pass
    print("[methy_caller] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
