#!/usr/bin/env python3
"""Shared methylation-call helpers (TAPS and EMSeq)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from scripts.split_bams import parse_cb_value, read_barcodes

BarcodeMap = Dict[str, int]
SMOKE_SPOT_LIMIT = 16


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
