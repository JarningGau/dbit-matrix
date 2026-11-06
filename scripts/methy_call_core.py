#!/usr/bin/env python3
"""Shared methylation-call helpers (TAPS and EMSeq)."""

from __future__ import annotations

import argparse
import json
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

from scripts.split_bams import normalize_cb_value, parse_cb_value, read_barcodes

BarcodeMap = Dict[str, int]
CgCounts = Dict[str, Dict[int, List[int]]]
ChCounts = Dict[str, Dict[int, List[object]]]
HostBamFunnel = Dict[str, object]
SMOKE_SPOT_LIMIT = 16

TOP_FLAGS: Set[int] = {99, 147}
BOT_FLAGS: Set[int] = {83, 163}
PAIRED_FLAGS: Set[int] = TOP_FLAGS | BOT_FLAGS
# Same proper-pair flags as scripts/summary.py host_valid_reads.
SUMMARY_VALID_FLAGS: Set[int] = PAIRED_FLAGS
CH_FORWARD_CONTEXTS: Set[str] = {"CA", "CC", "CT"}
CH_REVERSE_MAP: Dict[str, str] = {"T": "CA", "G": "CC", "A": "CT"}
CH_REVERSE_NEIGHBOR: Dict[str, str] = {"CA": "T", "CC": "G", "CT": "A"}
HOST_CG_FIELDS = 7
HOST_CH_FIELDS = 9
# Summary-compatible MAPQ gate: mapping_quality <= 10 excluded.
SUMMARY_MIN_MAPPING_QUALITY = 10
METH_STATS_FIELDS = [
    "barcode",
    "mean_5mCG",
    "mean_5mCH",
    "mean_5mCA",
    "mean_5mCC",
    "mean_5mCT",
    "cpg_site_count",
    "ch_site_count",
    "ca_site_count",
    "cc_site_count",
    "ct_site_count",
]

WORKER_CHROMOSOME_SEQUENCE: Optional[str] = None


@dataclass(frozen=True)
class SpotCount:
    barcode: str
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
    """Return normalized <bcB><bcA> when both halves are on the whitelist."""
    parsed = parse_cb_value(cb, barcode_len)
    if parsed is None:
        return None
    x_barcode, y_barcode = parsed
    if x_barcode not in barcode_to_index or y_barcode not in barcode_to_index:
        return None
    return normalize_cb_value(cb, barcode_len)


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


def empty_host_bam_funnel(mito_chromosomes: frozenset[str]) -> HostBamFunnel:
    return {
        "host_mapped_reads": 0,
        "host_valid_reads": 0,
        "host_valid_reads_mito": 0,
        "host_valid_reads_nuclear": 0,
        "mito_chromosomes": sorted(mito_chromosomes),
    }


def is_summary_confident_mapped(record: object) -> bool:
    """Match scripts/summary.py mapped filter (MAPQ>10, NH absent or ==1)."""
    if bool(getattr(record, "is_unmapped", False)):
        return False
    mapping_quality = int(getattr(record, "mapping_quality", 0) or 0)
    if mapping_quality <= SUMMARY_MIN_MAPPING_QUALITY:
        return False
    try:
        nh = record.get_tag("NH")  # type: ignore[attr-defined]
        if isinstance(nh, int) and nh > 1:
            return False
    except (KeyError, AttributeError):
        pass
    return True


def accumulate_host_bam_funnel(
    record: object,
    mito_chromosomes: frozenset[str],
    funnel: HostBamFunnel,
) -> None:
    if not is_summary_confident_mapped(record):
        return
    funnel["host_mapped_reads"] = int(funnel["host_mapped_reads"]) + 1
    flag = int(getattr(record, "flag", 0) or 0)
    if flag not in SUMMARY_VALID_FLAGS:
        return
    funnel["host_valid_reads"] = int(funnel["host_valid_reads"]) + 1
    reference_name = getattr(record, "reference_name", None)
    if isinstance(reference_name, str) and reference_name in mito_chromosomes:
        funnel["host_valid_reads_mito"] = int(funnel["host_valid_reads_mito"]) + 1
    else:
        funnel["host_valid_reads_nuclear"] = int(funnel["host_valid_reads_nuclear"]) + 1


def scan_host_bam_for_counts_and_funnel(
    records: Iterable[object],
    cb_tag: str,
    barcode_to_index: BarcodeMap,
    barcode_len: int,
    mito_chromosomes: frozenset[str],
) -> Tuple[Dict[str, int], HostBamFunnel]:
    """Single BAM pass: per-spot CB counts plus summary host funnel metrics."""
    counts: Dict[str, int] = {}
    funnel = empty_host_bam_funnel(mito_chromosomes)
    for record in records:
        accumulate_host_bam_funnel(record, mito_chromosomes, funnel)
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
    return counts, funnel


def write_host_bam_funnel_json(path: Path, funnel: HostBamFunnel) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(funnel, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    tmp.replace(path)


def load_host_bam_funnel_json(path: Path) -> HostBamFunnel | None:
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(payload, dict):
        return None
    required = (
        "host_mapped_reads",
        "host_valid_reads",
        "host_valid_reads_mito",
        "host_valid_reads_nuclear",
    )
    for key in required:
        if key not in payload:
            return None
        try:
            payload[key] = int(payload[key])
        except (TypeError, ValueError):
            return None
    mito = payload.get("mito_chromosomes")
    if mito is None:
        payload["mito_chromosomes"] = []
    elif isinstance(mito, list):
        payload["mito_chromosomes"] = [str(item) for item in mito]
    else:
        return None
    return payload


def write_counts_tsv(path: Path, counts: Dict[str, int]) -> None:
    rows = [(barcode, int(reads)) for barcode, reads in counts.items()]
    rows.sort(key=lambda item: item[0])
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", encoding="utf-8") as handle:
        handle.write("barcode\treads\n")
        for barcode, reads in rows:
            handle.write(f"{barcode}\t{reads}\n")
    tmp.replace(path)


def load_counts_tsv(path: Path) -> List[SpotCount]:
    rows: List[SpotCount] = []
    with path.open(encoding="utf-8") as handle:
        header = handle.readline()
        if "barcode" not in header:
            raise ValueError(f"counts TSV missing barcode header: {path}")
        header_fields = header.rstrip("\n").split("\t")
        barcode_idx = header_fields.index("barcode")
        reads_idx = header_fields.index("reads")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(barcode_idx, reads_idx):
                continue
            rows.append(
                SpotCount(
                    barcode=fields[barcode_idx],
                    reads=int(fields[reads_idx]),
                )
            )
    rows.sort(key=lambda row: row.barcode)
    return rows


def load_roster_ids(path: Path) -> List[str]:
    text = path.read_text(encoding="utf-8")
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        return []
    if "\t" in lines[0] and "barcode" in lines[0].split("\t"):
        header = lines[0].split("\t")
        barcode_idx = header.index("barcode")
        return [
            line.split("\t")[barcode_idx]
            for line in lines[1:]
            if len(line.split("\t")) > barcode_idx
        ]
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
    return {row.barcode: row.reads for row in rows}


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


def resolved_funnel_path(args: argparse.Namespace) -> Optional[Path]:
    counts_path = resolved_counts_path(args)
    if counts_path is None:
        return None
    return counts_path.parent / "host_bam_funnel.json"


def resolved_meth_stats_path(args: argparse.Namespace) -> Optional[Path]:
    counts_path = resolved_counts_path(args)
    if counts_path is None:
        return None
    return counts_path.parent / "per_spot_meth_stats.tsv"


def parse_mito_chromosomes(chromosome_csv: str) -> frozenset[str]:
    chromosomes = [item.strip() for item in chromosome_csv.split(",") if item.strip()]
    if not chromosomes:
        raise ValueError("mito chromosome list is empty")
    return frozenset(chromosomes)


def host_cov_path(out_dir: Path, spot_id: str, suffix: str) -> Path:
    """Host cov path: coverage/host/<first3>/<CB>.{CG,CH}.cov."""
    if len(spot_id) < 3:
        raise ValueError(f"CB spot id too short for prefix subdir: {spot_id}")
    return out_dir / spot_id[:3] / f"{spot_id}.{suffix}.cov"


def new_meth_stats_accumulator() -> Dict[str, Dict[str, float | int]]:
    """Per-spot running sums/counts for summary-compatible meth means."""
    return {}


def _spot_meth_bucket(acc: Dict[str, Dict[str, float | int]], spot_id: str) -> Dict[str, float | int]:
    bucket = acc.get(spot_id)
    if bucket is None:
        bucket = {
            "cg_sum": 0.0,
            "cg_count": 0,
            "ch_sum": 0.0,
            "ch_count": 0,
            "ca_sum": 0.0,
            "ca_count": 0,
            "cc_sum": 0.0,
            "cc_count": 0,
            "ct_sum": 0.0,
            "ct_count": 0,
        }
        acc[spot_id] = bucket
    return bucket


def accumulate_cov_line_for_meth_stats(
    acc: Dict[str, Dict[str, float | int]],
    suffix: str,
    fields: List[str],
    spot_id: Optional[str] = None,
) -> None:
    """Accumulate one part/cov line. Part lines include trailing spot id."""
    if suffix == "CG":
        min_fields = 6
        spot_idx = 6
        context_idx: Optional[int] = None
    elif suffix == "CH":
        min_fields = 6
        spot_idx = 8
        context_idx = 6
    else:
        return
    if spot_id is None:
        if len(fields) <= spot_idx:
            return
        spot_id = fields[spot_idx]
    if len(fields) < min_fields:
        return
    try:
        methylation = float(fields[3])
        methylated_count = int(fields[4])
        unmethylated_count = int(fields[5])
    except ValueError:
        return
    if methylated_count + unmethylated_count <= 0:
        return
    bucket = _spot_meth_bucket(acc, spot_id)
    if suffix == "CG":
        bucket["cg_sum"] = float(bucket["cg_sum"]) + methylation
        bucket["cg_count"] = int(bucket["cg_count"]) + 1
        return
    bucket["ch_sum"] = float(bucket["ch_sum"]) + methylation
    bucket["ch_count"] = int(bucket["ch_count"]) + 1
    if context_idx is None or len(fields) <= context_idx:
        return
    context = fields[context_idx].strip().upper()
    if context == "CA":
        bucket["ca_sum"] = float(bucket["ca_sum"]) + methylation
        bucket["ca_count"] = int(bucket["ca_count"]) + 1
    elif context == "CC":
        bucket["cc_sum"] = float(bucket["cc_sum"]) + methylation
        bucket["cc_count"] = int(bucket["cc_count"]) + 1
    elif context == "CT":
        bucket["ct_sum"] = float(bucket["ct_sum"]) + methylation
        bucket["ct_count"] = int(bucket["ct_count"]) + 1


def _format_mean_count(total: float, count: int) -> Tuple[str, str]:
    if count <= 0:
        return "NA", "NA"
    return f"{total / count:.2f}", str(count)


def finalize_meth_stats_rows(
    acc: Dict[str, Dict[str, float | int]],
    roster: List[str],
) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for spot_id in roster:
        bucket = acc.get(spot_id)
        if bucket is None:
            rows.append(
                {
                    "barcode": spot_id,
                    "mean_5mCG": "NA",
                    "mean_5mCH": "NA",
                    "mean_5mCA": "NA",
                    "mean_5mCC": "NA",
                    "mean_5mCT": "NA",
                    "cpg_site_count": "NA",
                    "ch_site_count": "NA",
                    "ca_site_count": "NA",
                    "cc_site_count": "NA",
                    "ct_site_count": "NA",
                }
            )
            continue
        mean_cg, cpg_count = _format_mean_count(
            float(bucket["cg_sum"]), int(bucket["cg_count"])
        )
        mean_ch, ch_count = _format_mean_count(
            float(bucket["ch_sum"]), int(bucket["ch_count"])
        )
        mean_ca, ca_count = _format_mean_count(
            float(bucket["ca_sum"]), int(bucket["ca_count"])
        )
        mean_cc, cc_count = _format_mean_count(
            float(bucket["cc_sum"]), int(bucket["cc_count"])
        )
        mean_ct, ct_count = _format_mean_count(
            float(bucket["ct_sum"]), int(bucket["ct_count"])
        )
        rows.append(
            {
                "barcode": spot_id,
                "mean_5mCG": mean_cg,
                "mean_5mCH": mean_ch,
                "mean_5mCA": mean_ca,
                "mean_5mCC": mean_cc,
                "mean_5mCT": mean_ct,
                "cpg_site_count": cpg_count,
                "ch_site_count": ch_count,
                "ca_site_count": ca_count,
                "cc_site_count": cc_count,
                "ct_site_count": ct_count,
            }
        )
    return rows


def write_per_spot_meth_stats_tsv(path: Path, rows: List[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(METH_STATS_FIELDS) + "\n")
        for row in rows:
            handle.write("\t".join(row.get(field, "NA") for field in METH_STATS_FIELDS) + "\n")
    tmp.replace(path)


def load_per_spot_meth_stats_tsv(path: Path) -> Dict[str, Dict[str, str]]:
    if not path.is_file():
        return {}
    rows: Dict[str, Dict[str, str]] = {}
    with path.open(encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if not header or header[0] != "barcode":
            return {}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < len(METH_STATS_FIELDS):
                continue
            row = {
                field: fields[idx] if idx < len(fields) else "NA"
                for idx, field in enumerate(METH_STATS_FIELDS)
            }
            barcode = row["barcode"].strip().upper()
            if barcode:
                rows[barcode] = row
    return rows


def parse_existing_cov_into_meth_stats(
    acc: Dict[str, Dict[str, float | int]],
    cov_path: Path,
    suffix: str,
    spot_id: str,
) -> None:
    if not cov_path.is_file():
        return
    with cov_path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split("\t")
            accumulate_cov_line_for_meth_stats(acc, suffix, fields, spot_id=spot_id)


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
    meth_acc: Optional[Dict[str, Dict[str, float | int]]] = None,
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
                if meth_acc is not None:
                    accumulate_cov_line_for_meth_stats(meth_acc, suffix, fields)
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


def merge_skip_spot_meth_stats(
    acc: Dict[str, Dict[str, float | int]],
    out_dir: Path,
    skip_spots: Set[str],
    context_mode: str,
    prior_rows: Dict[str, Dict[str, str]],
) -> None:
    """Fill accumulator for skipped spots from prior TSV or existing cov only."""
    for spot_id in sorted(skip_spots):
        if spot_id in acc:
            continue
        prior = prior_rows.get(spot_id.upper()) or prior_rows.get(spot_id)
        if prior is not None:
            # Rehydrate prior finalized means into accumulator via synthetic counts.
            # Prefer re-parse only when prior row is absent (below).
            _rehydrate_meth_stats_from_row(acc, prior)
            continue
        if context_mode in {"cg", "both"}:
            parse_existing_cov_into_meth_stats(
                acc, host_cov_path(out_dir, spot_id, "CG"), "CG", spot_id
            )
        if context_mode in {"ch", "both"}:
            parse_existing_cov_into_meth_stats(
                acc, host_cov_path(out_dir, spot_id, "CH"), "CH", spot_id
            )


def _rehydrate_meth_stats_from_row(
    acc: Dict[str, Dict[str, float | int]],
    row: Dict[str, str],
) -> None:
    spot_id = row["barcode"]
    bucket = _spot_meth_bucket(acc, spot_id)

    def apply(mean_key: str, count_key: str, sum_name: str, count_name: str) -> None:
        mean_text = row.get(mean_key, "NA")
        count_text = row.get(count_key, "NA")
        if mean_text == "NA" or count_text == "NA":
            return
        try:
            mean_val = float(mean_text)
            count_val = int(count_text)
        except ValueError:
            return
        if count_val <= 0:
            return
        bucket[sum_name] = mean_val * count_val
        bucket[count_name] = count_val

    apply("mean_5mCG", "cpg_site_count", "cg_sum", "cg_count")
    apply("mean_5mCH", "ch_site_count", "ch_sum", "ch_count")
    apply("mean_5mCA", "ca_site_count", "ca_sum", "ca_count")
    apply("mean_5mCC", "cc_site_count", "cc_sum", "cc_count")
    apply("mean_5mCT", "ct_site_count", "ct_sum", "ct_count")


def shuffle_host_parts(
    combined_cg: Path,
    combined_ch: Optional[Path],
    out_dir: Path,
    roster: List[str],
    skip_spots: Set[str],
    context_mode: str,
    meth_stats_path: Optional[Path] = None,
) -> None:
    allowed_spots = set(roster)
    meth_acc = new_meth_stats_accumulator()
    if context_mode in {"cg", "both"}:
        stream_sorted_parts(
            combined_cg,
            "CG",
            HOST_CG_FIELDS,
            out_dir,
            allowed_spots,
            skip_spots,
            meth_acc=meth_acc,
        )
    if context_mode in {"ch", "both"} and combined_ch is not None:
        stream_sorted_parts(
            combined_ch,
            "CH",
            HOST_CH_FIELDS,
            out_dir,
            allowed_spots,
            skip_spots,
            meth_acc=meth_acc,
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
    if meth_stats_path is not None:
        prior_rows = load_per_spot_meth_stats_tsv(meth_stats_path)
        merge_skip_spot_meth_stats(
            meth_acc, out_dir, skip_spots, context_mode, prior_rows
        )
        rows = finalize_meth_stats_rows(meth_acc, roster)
        write_per_spot_meth_stats_tsv(meth_stats_path, rows)


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
        print(f"[{log_prefix}] host_bam_funnel={resolved_funnel_path(args)}")
        print(f"[{log_prefix}] per_spot_meth_stats={resolved_meth_stats_path(args)}")
        print(
            f"[{log_prefix}] mito_chromosomes="
            f"{getattr(args, 'mito_chromosomes', 'chrM')}"
        )
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
            funnel_path = resolved_funnel_path(args)
            assert counts_path is not None
            assert funnel_path is not None
            mito_csv = getattr(args, "mito_chromosomes", None) or "chrM"
            mito_chromosomes = parse_mito_chromosomes(str(mito_csv))
            counts_ready = counts_path.is_file() and not args.force_recount
            funnel_ready = (
                load_host_bam_funnel_json(funnel_path) is not None
                and not args.force_recount
            )
            if counts_ready and funnel_ready:
                counts = counts_dict_from_rows(load_counts_tsv(counts_path))
            else:
                # One BAM pass refreshes counts and/or funnel together.
                with pysam.AlignmentFile(args.bam, "rb") as bam_in:
                    counts, funnel = scan_host_bam_for_counts_and_funnel(
                        bam_in.fetch(until_eof=True),
                        args.cb_tag,
                        barcode_to_index,
                        barcode_len,
                        mito_chromosomes,
                    )
                write_counts_tsv(counts_path, counts)
                write_host_bam_funnel_json(funnel_path, funnel)
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
                meth_stats_path=resolved_meth_stats_path(args),
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
