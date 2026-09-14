#!/usr/bin/env python3
"""Shared methylation-call helpers (TAPS and EMSeq)."""

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
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, Dict, Iterable, List, Optional, Set, Tuple

import pysam

from scripts.split_bams import parse_cb_value, read_barcodes

BarcodeMap = Dict[str, int]
CgCounts = Dict[str, Dict[int, List[int]]]
ChCounts = Dict[str, Dict[int, List[object]]]
SMOKE_SPOT_LIMIT = 16

TOP_FLAGS: Set[int] = {99, 147}
BOT_FLAGS: Set[int] = {83, 163}
PAIRED_FLAGS: Set[int] = TOP_FLAGS | BOT_FLAGS
CH_FORWARD_CONTEXTS: Set[str] = {"CA", "CC", "CT"}
CH_REVERSE_MAP: Dict[str, str] = {"T": "CA", "G": "CC", "A": "CT"}
CH_REVERSE_NEIGHBOR: Dict[str, str] = {"CA": "T", "CC": "G", "CT": "A"}
HOST_CG_FIELDS = 7
HOST_CH_FIELDS = 9

WORKER_CHROMOSOME_SEQUENCE: Optional[str] = None


@dataclass(frozen=True)
class SpotCount:
    x_index: int
    y_index: int
    spot: str
    reads: int


def classify_cg_observation(
    dinucleotide: str, converted_means_methylated: bool
) -> Optional[bool]:
    if dinucleotide in {"TG", "CA"}:
        converted = True
    elif dinucleotide == "CG":
        converted = False
    else:
        return None
    return converted if converted_means_methylated else not converted


def classify_ch_observation(
    base: str, strand: str, converted_means_methylated: bool
) -> Optional[Tuple[int, int]]:
    if strand == "+":
        if base not in {"C", "T"}:
            return None
        retained = base == "C"
    else:
        if base not in {"G", "A"}:
            return None
        retained = base == "G"
    methylated = retained if not converted_means_methylated else not retained
    return (1, 0) if methylated else (0, 1)


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


def count_records_by_spot(
    records: Iterable[object],
    cb_tag: str,
    barcode_to_index: BarcodeMap,
    barcode_len: int,
) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for record in records:
        try:
            if hasattr(record, "has_tag") and not record.has_tag(cb_tag):
                continue
            cb = record.get_tag(cb_tag)
        except KeyError:
            continue
        if not isinstance(cb, str):
            cb = str(cb)
        spot = spot_id_from_cb(cb, barcode_to_index, barcode_len)
        if spot is None:
            continue
        counts[spot] = counts.get(spot, 0) + 1
    return counts


def write_counts_tsv(path: Path, counts: Dict[str, int]) -> None:
    rows = []
    for spot, reads in counts.items():
        x_text, y_text = spot.split("_", 1)
        rows.append((int(x_text), int(y_text), spot, int(reads)))
    rows.sort()
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", encoding="utf-8") as handle:
        handle.write("X_index\tY_index\tspot\treads\n")
        for x_index, y_index, spot, reads in rows:
            handle.write(f"{x_index}\t{y_index}\t{spot}\t{reads}\n")
    tmp.replace(path)


def load_counts_tsv(path: Path) -> List[SpotCount]:
    rows: List[SpotCount] = []
    with path.open(encoding="utf-8") as handle:
        header = handle.readline()
        if "spot" not in header:
            raise ValueError(f"counts TSV missing header: {path}")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            rows.append(
                SpotCount(
                    x_index=int(fields[0]),
                    y_index=int(fields[1]),
                    spot=fields[2],
                    reads=int(fields[3]),
                )
            )
    rows.sort(key=lambda row: (row.x_index, row.y_index))
    return rows


def load_roster_ids(path: Path) -> List[str]:
    text = path.read_text(encoding="utf-8")
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        return []
    if "\t" in lines[0] and "spot" in lines[0].split("\t"):
        header = lines[0].split("\t")
        spot_idx = header.index("spot")
        return [line.split("\t")[spot_idx] for line in lines[1:] if len(line.split("\t")) > spot_idx]
    return lines


def resolve_write_set(
    counts: Dict[str, int],
    smoke: bool,
    roster_ids: Optional[List[str]],
    limit: int = SMOKE_SPOT_LIMIT,
) -> List[str]:
    if roster_ids is not None:
        return list(roster_ids)
    spots = sorted(counts, key=lambda spot: (-counts[spot], spot))
    if smoke:
        return spots[:limit]
    return sorted(counts)


def counts_dict_from_rows(rows: List[SpotCount]) -> Dict[str, int]:
    return {row.spot: row.reads for row in rows}


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


def spot_is_allowed(spot_id: Optional[str], allowed: Optional[Set[str]]) -> bool:
    if spot_id is None:
        return False
    if allowed is None:
        return True
    return spot_id in allowed


def resolved_counts_path(args: argparse.Namespace) -> Optional[Path]:
    if args.counts_tsv:
        return Path(args.counts_tsv)
    if args.out_dir:
        return Path(args.out_dir).resolve().parent / "per_spot_read_counts.tsv"
    return None


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
    converted_means_methylated: bool,
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
        methylated = classify_cg_observation(
            (seq[qpos] + seq[qpos + 1]).upper(), converted_means_methylated
        )
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
    converted_means_methylated: bool,
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
        observation = classify_ch_observation(
            seq[qpos].upper(), strand, converted_means_methylated
        )
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
    converted_means_methylated: bool,
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
                    converted_means_methylated,
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
                    converted_means_methylated,
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
    converted_means_methylated: bool,
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
            converted_means_methylated,
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
        if not args.barcodes or not args.out_dir:
            raise ValueError("host mode requires --barcodes and --out-dir")
    elif not args.output:
        raise ValueError("aggregate mode requires --output")


def print_resolved_parameters(args: argparse.Namespace, log_prefix: str) -> None:
    print(f"[{log_prefix}] mode={args.mode}")
    print(f"[{log_prefix}] bam={args.bam}")
    print(f"[{log_prefix}] reference={args.reference}")
    print(f"[{log_prefix}] chromosomes={args.chromosomes}")
    print(f"[{log_prefix}] context_mode={args.context_mode}")
    print(f"[{log_prefix}] jobs={args.jobs}")
    print(f"[{log_prefix}] r1_left_trimming={args.r1_left_trimming}")
    print(f"[{log_prefix}] r1_right_trimming={args.r1_right_trimming}")
    print(f"[{log_prefix}] r2_left_trimming={args.r2_left_trimming}")
    print(f"[{log_prefix}] r2_right_trimming={args.r2_right_trimming}")
    if args.mode == "host":
        print(f"[{log_prefix}] barcodes={args.barcodes}")
        print(f"[{log_prefix}] spot_roster={args.spot_roster}")
        print(f"[{log_prefix}] out_dir={args.out_dir}")
        print(f"[{log_prefix}] smoke={int(bool(args.smoke))}")
        print(f"[{log_prefix}] counts_tsv={resolved_counts_path(args)}")
    else:
        print(f"[{log_prefix}] output={args.output}")


def run_caller(
    args: argparse.Namespace,
    converted_means_methylated: bool,
    log_prefix: str,
) -> int:
    print_resolved_parameters(args, log_prefix)
    if args.dry_run:
        if args.mode == "host":
            print(f"[{log_prefix}] roster_spots=0")
        print(f"[{log_prefix}] dry_run=1")
        return 0
    try:
        validate_runtime(args)
        if args.mode == "host":
            if not Path(args.barcodes).is_file():
                raise ValueError(f"barcode whitelist not found: {args.barcodes}")
            if args.spot_roster and not Path(args.spot_roster).is_file():
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
        print(f"[{log_prefix}] error: {exc}", file=sys.stderr)
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
            counts_path = resolved_counts_path(args)
            assert counts_path is not None
            if counts_path.is_file() and not args.force_recount:
                counts = counts_dict_from_rows(load_counts_tsv(counts_path))
            else:
                with pysam.AlignmentFile(args.bam, "rb") as bam_in:
                    counts = count_records_by_spot(
                        bam_in.fetch(until_eof=True),
                        args.cb_tag,
                        barcode_to_index,
                        barcode_len,
                    )
                write_counts_tsv(counts_path, counts)
            roster_ids = (
                load_roster_ids(Path(args.spot_roster)) if args.spot_roster else None
            )
            roster = resolve_write_set(counts, args.smoke, roster_ids)
        except (ValueError, OSError) as exc:
            print(f"[{log_prefix}] error: {exc}", file=sys.stderr)
            return 1
        allowed_spots = set(roster)
        skip_spots = existing_skip_spots(Path(args.out_dir), roster, args.context_mode)
        print(f"[{log_prefix}] roster_spots={len(roster)} skip_spots={len(skip_spots)}")

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
            print(f"[{log_prefix}] skip_existing_aggregate=1")
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
                    converted_means_methylated,
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
        print(f"[{log_prefix}] error: {exc}", file=sys.stderr)
        return 1
    finally:
        for handle in handles.values():
            handle.close()
        for temp_path in combined_temps.values():
            try:
                temp_path.unlink()
            except FileNotFoundError:
                pass
    print(f"[{log_prefix}] done")
    return 0
