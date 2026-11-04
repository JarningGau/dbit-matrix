#!/usr/bin/env python3
"""Summarize per-spot outputs and render summary heatmaps (TAPS v2).

Extends the v1 summary with methylated-linker mC→T conversion QC columns
aggregated from demux ``*.stats.json``.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections.abc import Callable
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import theme_nature as tn

tn.apply_theme()


HEATMAP_SPECS = (
    {
        "field": "reads",
        "output_name": "reads_heatmap.png",
        "title": "Reads per Spot",
        "colorbar_label": "Reads",
        "cmap": "viridis",
        "vmin": 0.0,
        "vmax": None,
        "origin_corner": "upper_right",
    },
    {
        "field": "cpg_site_count",
        "output_name": "cpg_site_count_heatmap.png",
        "title": "CpG Sites per Spot",
        "colorbar_label": "CpG sites",
        "cmap": "magma",
        "vmin": 0.0,
        "vmax": None,
        "origin_corner": "upper_right",
    },
    {
        "field": "mean_5mCG",
        "output_name": "mean_5mCG_heatmap.png",
        "title": "Mean 5mCG per Spot",
        "colorbar_label": "Mean 5mCG (%)",
        "cmap": "coolwarm",
        "vmin": 50.0,
        "vmax": 80.0,
        "origin_corner": "upper_right",
        "hq_filter": True,
    },
    {
        "field": "mean_5mCH",
        "output_name": "mean_5mCH_heatmap.png",
        "title": "Mean 5mCH per Spot",
        "colorbar_label": "Mean 5mCH (%)",
        "cmap": "coolwarm",
        "vmin": 0.0,
        "vmax": 10.0,
        "origin_corner": "upper_right",
        "hq_filter": True,
    },
    {
        "field": "mean_5mCA",
        "output_name": "mean_5mCA_heatmap.png",
        "title": "Mean 5mCA per Spot",
        "colorbar_label": "Mean 5mCA (%)",
        "cmap": "coolwarm",
        "vmin": 0.0,
        "vmax": 10.0,
        "origin_corner": "upper_right",
        "hq_filter": True,
    },
    {
        "field": "mean_5mCC",
        "output_name": "mean_5mCC_heatmap.png",
        "title": "Mean 5mCC per Spot",
        "colorbar_label": "Mean 5mCC (%)",
        "cmap": "coolwarm",
        "vmin": 0.0,
        "vmax": 10.0,
        "origin_corner": "upper_right",
        "hq_filter": True,
    },
    {
        "field": "mean_5mCT",
        "output_name": "mean_5mCT_heatmap.png",
        "title": "Mean 5mCT per Spot",
        "colorbar_label": "Mean 5mCT (%)",
        "cmap": "coolwarm",
        "vmin": 0.0,
        "vmax": 10.0,
        "origin_corner": "upper_right",
        "hq_filter": True,
    },
)
CH_CONTEXTS = ("CA", "CC", "CT")
VALID_FLAGS = {99, 147, 83, 163}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate summary TSVs and per-spot heatmaps from call outputs "
            "under a sample work directory."
        )
    )
    parser.add_argument(
        "--work-path",
        required=True,
        help="Sample work directory containing coverage/ and split_bams/.",
    )
    parser.add_argument(
        "--spike-in-name",
        action="append",
        default=[],
        help=(
            "Expected spike-in name for stable sample summary columns "
            "(e.g. lambda). May be specified multiple times."
        ),
    )
    parser.add_argument(
        "--mito-chromosomes",
        default="chrM",
        help=(
            "Comma-separated host contigs treated as mitochondrial when "
            "splitting host_valid_reads into mito vs nuclear. Default: chrM."
        ),
    )
    parser.add_argument(
        "--barcode-positions",
        required=True,
        help=(
            "TSV with columns barcode,array_row,array_col "
            "(barcode = <bcB><bcA>); used for per-spot coords and heatmaps."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print resolved inputs/outputs and exit without writing files.",
    )
    return parser.parse_args()


def read_barcode_positions(path: str | Path) -> dict[str, tuple[str, str]]:
    """Map CB (<bcB><bcA>) -> (array_row, array_col) from a positions TSV."""
    positions: dict[str, tuple[str, str]] = {}
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"empty barcode positions file: {path}")
        required = {"barcode", "array_row", "array_col"}
        missing = required - {name.strip() for name in reader.fieldnames}
        if missing:
            raise ValueError(
                f"barcode positions missing columns {sorted(missing)}: {path}"
            )
        for row in reader:
            barcode = (row.get("barcode") or "").strip().upper()
            array_row = (row.get("array_row") or "").strip()
            array_col = (row.get("array_col") or "").strip()
            if not barcode:
                continue
            positions[barcode] = (array_row or "NA", array_col or "NA")
    if not positions:
        raise ValueError(f"empty barcode positions file: {path}")
    return positions


def summary_barcode_from_cb(cb: str) -> str:
    """Summary barcode column: <bcB><bcA>."""
    return cb


def parse_mito_chromosomes(chromosome_csv: str) -> frozenset[str]:
    chromosomes = [item.strip() for item in chromosome_csv.split(",") if item.strip()]
    if not chromosomes:
        raise ValueError("--mito-chromosomes resolved to an empty list")
    return frozenset(chromosomes)


def parse_cov_stats(cov_path: Path) -> tuple[float, int] | None:
    if not cov_path.exists():
        return None
    cpg_count = 0
    methylation_sum = 0.0
    with cov_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 6:
                continue
            try:
                methylation = float(fields[3])
                methylated_count = int(fields[4])
                unmethylated_count = int(fields[5])
            except ValueError:
                continue
            if methylated_count + unmethylated_count <= 0:
                continue
            methylation_sum += methylation
            cpg_count += 1
    if cpg_count == 0:
        return None
    return methylation_sum / cpg_count, cpg_count


def _finalize_mean_count(
    methylation_sum: float, site_count: int
) -> tuple[float, int] | None:
    if site_count <= 0:
        return None
    return methylation_sum / site_count, site_count


def parse_ch_cov_stats(
    cov_path: Path,
) -> tuple[tuple[float, int] | None, dict[str, tuple[float, int] | None]]:
    """Parse .CH.cov: overall CH mean plus CA/CC/CT means from the context column.

    Rows without a recognized context still count toward overall CH only.
    """
    empty_contexts: dict[str, tuple[float, int] | None] = {
        context: None for context in CH_CONTEXTS
    }
    if not cov_path.exists():
        return None, empty_contexts

    ch_sum = 0.0
    ch_count = 0
    context_sums = {context: 0.0 for context in CH_CONTEXTS}
    context_counts = {context: 0 for context in CH_CONTEXTS}
    with cov_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 6:
                continue
            try:
                methylation = float(fields[3])
                methylated_count = int(fields[4])
                unmethylated_count = int(fields[5])
            except ValueError:
                continue
            if methylated_count + unmethylated_count <= 0:
                continue
            ch_sum += methylation
            ch_count += 1
            if len(fields) < 7:
                continue
            context = fields[6].strip().upper()
            if context not in context_sums:
                continue
            context_sums[context] += methylation
            context_counts[context] += 1

    context_stats = {
        context: _finalize_mean_count(context_sums[context], context_counts[context])
        for context in CH_CONTEXTS
    }
    return _finalize_mean_count(ch_sum, ch_count), context_stats


def to_optional_int(value: object) -> int | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            return None
        return int(value)
    if isinstance(value, str):
        text = value.strip()
        if not text:
            return None
        try:
            return int(text)
        except ValueError:
            try:
                return int(float(text))
            except ValueError:
                return None
    return None


def format_optional_int(value: int | None) -> str:
    if value is None:
        return "NA"
    return f"{value:,}"


def format_percentage(numerator: int | None, denominator: int | None) -> str:
    if numerator is None or denominator is None or denominator <= 0:
        return "NA"
    return f"{(numerator / denominator) * 100:.2f}%"


def parse_fastp_raw_reads(fastp_json_path: Path) -> int | None:
    if not fastp_json_path.exists():
        return None
    try:
        payload = json.loads(fastp_json_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(payload, dict):
        return None

    summary = payload.get("summary")
    if isinstance(summary, dict):
        before_filtering = summary.get("before_filtering")
        if isinstance(before_filtering, dict):
            total_reads = to_optional_int(before_filtering.get("total_reads"))
            if total_reads is not None:
                return total_reads
        total_reads = to_optional_int(summary.get("total_reads"))
        if total_reads is not None:
            return total_reads

    return to_optional_int(payload.get("total_reads"))


def collect_barcoded_reads(demux_dir: Path) -> int | None:
    if not demux_dir.exists():
        return None
    stats_paths = sorted(demux_dir.glob("*.stats.json"))
    if not stats_paths:
        return None

    total = 0
    has_value = False
    for stats_path in stats_paths:
        try:
            payload = json.loads(stats_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if not isinstance(payload, dict):
            continue
        kept_reads = to_optional_int(payload.get("kept_reads"))
        if kept_reads is None:
            continue
        # demux stats report read pairs; summary read metrics are normalized to reads.
        total += kept_reads * 2
        has_value = True
    if not has_value:
        return None
    return total


def to_optional_float(value: object) -> float | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        if isinstance(value, float) and not math.isfinite(value):
            return None
        return float(value)
    if isinstance(value, str):
        text = value.strip()
        if not text or text.upper() == "NA":
            return None
        try:
            parsed = float(text)
        except ValueError:
            return None
        if not math.isfinite(parsed):
            return None
        return parsed
    return None


def format_optional_float(value: float | None, digits: int = 2) -> str:
    if value is None:
        return "NA"
    return f"{value:.{digits}f}"


def load_hmc_spike_summary(summary_json_path: Path) -> dict | None:
    if not summary_json_path.exists():
        return None
    try:
        payload = json.loads(summary_json_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(payload, dict):
        return None
    name = payload.get("name")
    if not isinstance(name, str) or not name.strip():
        return None
    return payload


def apply_hmc_spike_columns(
    row: dict[str, str],
    hmc_summary: dict | None,
    raw_reads: int | None,
) -> None:
    if hmc_summary is None:
        return
    name = str(hmc_summary["name"]).strip()
    matched = hmc_summary.get("reads_matched")
    matched_int: int | None
    try:
        matched_int = int(matched) if matched is not None else None
    except (TypeError, ValueError):
        matched_int = None
    mean_rate = hmc_summary.get("mean_conversion_rate")
    mean_float: float | None
    try:
        mean_float = float(mean_rate) if mean_rate is not None else None
    except (TypeError, ValueError):
        mean_float = None
    row[f"{name}_hmc_spike_reads"] = format_optional_int(matched_int)
    row[f"{name}_hmc_spike_reads_rate"] = format_percentage(matched_int, raw_reads)
    row[f"{name}_hmc_spike_mean_hmCtoT"] = format_optional_float(mean_float)


def collect_methylated_linker_qc(demux_dir: Path) -> dict[str, int | float | None]:
    """Aggregate methylated-linker conversion QC from demux stats.json chunks."""
    empty: dict[str, int | float | None] = {
        "spike_scored_reads": None,
        "scored_fully_converted": None,
        "spike_mean_conversion_rate": None,
        "spike_fraction_all_T": None,
        "spike_conversion_filtered_rejects": None,
    }
    if not demux_dir.exists():
        return empty
    stats_paths = sorted(demux_dir.glob("*.stats.json"))
    if not stats_paths:
        return empty

    scored_total = 0
    fully_converted_total = 0
    rate_weighted_sum = 0.0
    rate_weight = 0
    all_t_weighted_sum = 0.0
    all_t_weight = 0
    filtered_rejects = 0
    has_scored = False
    has_fully = False
    has_filter_reject = False

    for stats_path in stats_paths:
        try:
            payload = json.loads(stats_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if not isinstance(payload, dict):
            continue

        scored = to_optional_int(payload.get("spike_scored_reads"))
        if scored is not None:
            scored_total += scored
            has_scored = True

        fully = to_optional_int(payload.get("scored_fully_converted"))
        if fully is not None:
            fully_converted_total += fully
            has_fully = True

        mean_rate = to_optional_float(payload.get("spike_mCtoT_mean_rate_3_6_10"))
        if scored is not None and scored > 0 and mean_rate is not None:
            rate_weighted_sum += mean_rate * scored
            rate_weight += scored

        frac_all_t = to_optional_float(
            payload.get("spike_mCtoT_fraction_all_T_3_6_10")
        )
        if scored is not None and scored > 0 and frac_all_t is not None:
            all_t_weighted_sum += frac_all_t * scored
            all_t_weight += scored

        reject_counts = payload.get("reject_counts")
        if isinstance(reject_counts, dict):
            reject_n = to_optional_int(reject_counts.get("spikein_c_not_all_t"))
            if reject_n is not None:
                filtered_rejects += reject_n
                has_filter_reject = True

    return {
        "spike_scored_reads": scored_total if has_scored else None,
        "scored_fully_converted": fully_converted_total if has_fully else None,
        "spike_mean_conversion_rate": (
            rate_weighted_sum / rate_weight if rate_weight > 0 else None
        ),
        "spike_fraction_all_T": (
            all_t_weighted_sum / all_t_weight if all_t_weight > 0 else None
        ),
        "spike_conversion_filtered_rejects": (
            filtered_rejects if has_filter_reject else None
        ),
    }


def count_host_bam_metrics(
    host_bam_path: Path,
    mito_chromosomes: frozenset[str],
) -> tuple[int | None, int | None, int | None, int | None]:
    if not host_bam_path.exists():
        return None, None, None, None
    mapped_reads = 0
    valid_reads = 0
    valid_reads_mito = 0
    valid_reads_nuclear = 0
    try:
        with pysam.AlignmentFile(str(host_bam_path), "rb") as bam_handle:
            for read in bam_handle.fetch(until_eof=True):
                if read.is_unmapped:
                    continue
                # Only count reasonably confident unique alignments.
                if read.mapping_quality <= 10:
                    continue
                try:
                    nh = read.get_tag("NH")
                    if isinstance(nh, int) and nh > 1:
                        # Multiple alignments (NH>1) are excluded.
                        continue
                except KeyError:
                    # NH tag may be absent; treat as unique in that case.
                    pass
                mapped_reads += 1
                if read.flag not in VALID_FLAGS:
                    continue
                valid_reads += 1
                if read.reference_name in mito_chromosomes:
                    valid_reads_mito += 1
                else:
                    valid_reads_nuclear += 1
    except OSError:
        return None, None, None, None
    return mapped_reads, valid_reads, valid_reads_mito, valid_reads_nuclear


def count_mapped_reads(bam_path: Path) -> int | None:
    if not bam_path.exists():
        return None
    mapped_reads = 0
    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam_handle:
            for read in bam_handle.fetch(until_eof=True):
                if read.is_unmapped:
                    continue
                # Only count reasonably confident unique alignments (MAPQ>10, NH==1 if present).
                if read.mapping_quality <= 10:
                    continue
                try:
                    nh = read.get_tag("NH")
                    if isinstance(nh, int) and nh > 1:
                        continue
                except KeyError:
                    # No NH tag; treat as unique.
                    pass
                mapped_reads += 1
    except OSError:
        return None
    return mapped_reads


def read_per_spot_reads(read_counts_path: Path) -> dict[str, str]:
    """Map CB barcode -> reads string from per_spot_read_counts.tsv."""
    if not read_counts_path.exists():
        return {}
    rows: dict[str, str] = {}
    with read_counts_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            barcode = (row.get("barcode") or "").strip().upper()
            if not barcode:
                continue
            reads = (row.get("reads") or "").strip()
            rows[barcode] = reads or "NA"
    return rows


def format_float(value: float | None) -> str:
    if value is None:
        return "NA"
    return f"{value:.2f}"


def parse_optional_float(text: str | None) -> float | None:
    value = (text or "").strip()
    if not value or value == "NA":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def format_is_hq_spot(
    reads_raw: str | None, reads_threshold: float | None
) -> str:
    """HQ flag for per_spot_summary: 1/0, or NA if threshold or reads missing."""
    if reads_threshold is None:
        return "NA"
    reads = parse_optional_float(reads_raw)
    if reads is None:
        return "NA"
    return "1" if reads > reads_threshold else "0"


PER_SPOT_FIELDS = [
    "barcode",
    "array_row",
    "array_col",
    "reads",
    "is_HQ_spot",
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


def build_sample_summary_fieldnames(
    spike_names: list[str],
    sample_row: dict[str, str],
) -> list[str]:
    """Semantic column order for sample_summary.tsv (linker QC after funnel)."""
    fields = [
        "sample_id",
        "saturation_rate",
        "raw_reads",
        "barcoded_reads",
        "barcoded_reads_rate",
        "host_mapped_reads",
    ]
    for spike_name in spike_names:
        fields.append(f"{spike_name}_mapped_reads")
    fields.extend(
        [
            "host_valid_reads",
            "host_valid_reads_mito",
            "host_valid_reads_nuclear",
            "valid_reads_rate",
            "spike_scored_reads",
            "scored_fully_converted",
            "spike_mean_conversion_rate",
            "spike_fraction_all_T",
            "spike_conversion_filtered_rejects",
            "host_spot_mean_5mCG",
            "host_spot_mean_5mCH",
            "host_spot_mean_5mCA",
            "host_spot_mean_5mCC",
            "host_spot_mean_5mCT",
            "host_spot_median_cpg_sites",
            "host_mito_mean_5mCG",
            "host_mito_mean_5mCH",
            "host_mito_mean_5mCA",
            "host_mito_mean_5mCC",
            "host_mito_mean_5mCT",
        ]
    )
    for spike_name in spike_names:
        fields.extend(
            [
                f"{spike_name}_mean_5mCG",
                f"{spike_name}_mean_5mCH",
                f"{spike_name}_mean_5mCA",
                f"{spike_name}_mean_5mCC",
                f"{spike_name}_mean_5mCT",
            ]
        )
    hmc_prefixes = sorted(
        {
            key[: -len("_hmc_spike_reads")]
            for key in sample_row
            if key.endswith("_hmc_spike_reads")
        }
    )
    for prefix in hmc_prefixes:
        for suffix in (
            "_hmc_spike_reads",
            "_hmc_spike_reads_rate",
            "_hmc_spike_mean_hmCtoT",
        ):
            key = f"{prefix}{suffix}"
            if key in sample_row:
                fields.append(key)
    return fields


def parse_optional_coord(text: str | None) -> int | None:
    """Parse 0-based array_row/array_col (inclusive; 0 is valid)."""
    value = (text or "").strip()
    if not value or value == "NA":
        return None
    try:
        parsed = int(value)
    except ValueError:
        return None
    if parsed < 0:
        return None
    return parsed


def ch_cov_path_from_cg(cg_cov_path: Path) -> Path:
    """Sibling CH cov path for a CG cov path (suffix swap)."""
    name = cg_cov_path.name
    if name.endswith(".CG.cov"):
        return cg_cov_path.with_name(name[: -len(".CG.cov")] + ".CH.cov")
    return cg_cov_path.with_suffix("").with_suffix(".CH.cov")


def weighted_mean_from_rows(
    rows: list[dict[str, str]],
    mean_field: str,
    count_field: str,
) -> tuple[float | None, list[int]]:
    """Site-count-weighted mean and the list of positive site counts."""
    weighted_sum = 0.0
    total_sites = 0
    counts: list[int] = []
    for row in rows:
        mean_text = row[mean_field]
        count_text = row[count_field]
        if mean_text == "NA" or count_text == "NA":
            continue
        site_count = int(count_text)
        if site_count <= 0:
            continue
        weighted_sum += float(mean_text) * site_count
        total_sites += site_count
        counts.append(site_count)
    if total_sites <= 0:
        return None, counts
    return weighted_sum / total_sites, counts


def median_of_ints(values: list[int]) -> float | None:
    if not values:
        return None
    sorted_values = sorted(values)
    n = len(sorted_values)
    mid = n // 2
    if n % 2 == 1:
        return float(sorted_values[mid])
    return (sorted_values[mid - 1] + sorted_values[mid]) / 2.0


def _site_count_str(stats: tuple[float, int] | None) -> str:
    if stats is None:
        return "NA"
    return str(stats[1])


def _mean_from_stats(stats: tuple[float, int] | None) -> float | None:
    if stats is None:
        return None
    return stats[0]


def summarize_per_spot(
    host_cov_paths: list[Path],
    read_counts: dict[str, str],
    positions: dict[str, tuple[str, str]],
) -> list[dict[str, str]]:
    # cb -> (cg_stats, ch_stats, context_stats)
    by_cov: dict[
        str,
        tuple[
            tuple[float, int] | None,
            tuple[float, int] | None,
            dict[str, tuple[float, int] | None],
        ],
    ] = {}
    empty_contexts: dict[str, tuple[float, int] | None] = {
        context: None for context in CH_CONTEXTS
    }
    for cg_cov_path in host_cov_paths:
        cb = cg_cov_path.name[: -len(".CG.cov")].upper()
        cg_stats = parse_cov_stats(cg_cov_path)
        ch_stats, context_stats = parse_ch_cov_stats(ch_cov_path_from_cg(cg_cov_path))
        by_cov[cb] = (cg_stats, ch_stats, context_stats)

    all_cbs = sorted(set(by_cov.keys()) | set(read_counts.keys()))
    rows: list[dict[str, str]] = []
    for cb in all_cbs:
        cg_stats, ch_stats, context_stats = by_cov.get(
            cb, (None, None, empty_contexts)
        )
        array_row, array_col = positions.get(cb, ("NA", "NA"))
        reads = read_counts.get(cb, "NA") or "NA"
        rows.append(
            {
                "barcode": summary_barcode_from_cb(cb),
                "array_row": array_row if array_row else "NA",
                "array_col": array_col if array_col else "NA",
                "mean_5mCG": format_float(_mean_from_stats(cg_stats)),
                "mean_5mCH": format_float(_mean_from_stats(ch_stats)),
                "mean_5mCA": format_float(_mean_from_stats(context_stats["CA"])),
                "mean_5mCC": format_float(_mean_from_stats(context_stats["CC"])),
                "mean_5mCT": format_float(_mean_from_stats(context_stats["CT"])),
                "cpg_site_count": _site_count_str(cg_stats),
                "ch_site_count": _site_count_str(ch_stats),
                "ca_site_count": _site_count_str(context_stats["CA"]),
                "cc_site_count": _site_count_str(context_stats["CC"]),
                "ct_site_count": _site_count_str(context_stats["CT"]),
                "reads": reads,
            }
        )
    return rows


def build_sample_summary_row(
    sample_id: str,
    per_spot_rows: list[dict[str, str]],
    host_mito_cg_cov_path: Path,
    host_mito_ch_cov_path: Path,
    spike_cg_cov_paths: dict[str, Path],
    spike_ch_cov_paths: dict[str, Path],
    saturation_rate: str | None,
    raw_reads: int | None,
    barcoded_reads: int | None,
    host_mapped_reads: int | None,
    host_valid_reads: int | None,
    host_valid_reads_mito: int | None,
    host_valid_reads_nuclear: int | None,
    spike_mapped_reads: dict[str, int | None],
    methylated_linker_qc: dict[str, int | float | None] | None = None,
) -> dict[str, str]:
    host_spot_mean_5mcg, cpg_counts_for_median = weighted_mean_from_rows(
        per_spot_rows, "mean_5mCG", "cpg_site_count"
    )
    host_spot_mean_5mch, _ = weighted_mean_from_rows(
        per_spot_rows, "mean_5mCH", "ch_site_count"
    )
    host_spot_mean_5mca, _ = weighted_mean_from_rows(
        per_spot_rows, "mean_5mCA", "ca_site_count"
    )
    host_spot_mean_5mcc, _ = weighted_mean_from_rows(
        per_spot_rows, "mean_5mCC", "cc_site_count"
    )
    host_spot_mean_5mct, _ = weighted_mean_from_rows(
        per_spot_rows, "mean_5mCT", "ct_site_count"
    )
    host_spot_median_cpg = median_of_ints(cpg_counts_for_median)

    host_mito_cg_stats = parse_cov_stats(host_mito_cg_cov_path)
    host_mito_ch_stats, host_mito_context_stats = parse_ch_cov_stats(
        host_mito_ch_cov_path
    )
    host_mito_mean_5mcg = host_mito_cg_stats[0] if host_mito_cg_stats else None
    host_mito_mean_5mch = _mean_from_stats(host_mito_ch_stats)

    row: dict[str, str] = {"sample_id": sample_id}
    row.update({
        "saturation_rate": saturation_rate if saturation_rate else "NA",
        "host_spot_mean_5mCG": format_float(host_spot_mean_5mcg),
        "host_spot_mean_5mCH": format_float(host_spot_mean_5mch),
        "host_spot_mean_5mCA": format_float(host_spot_mean_5mca),
        "host_spot_mean_5mCC": format_float(host_spot_mean_5mcc),
        "host_spot_mean_5mCT": format_float(host_spot_mean_5mct),
        "host_spot_median_cpg_sites": format_float(host_spot_median_cpg),
        "host_mito_mean_5mCG": format_float(host_mito_mean_5mcg),
        "host_mito_mean_5mCH": format_float(host_mito_mean_5mch),
        "host_mito_mean_5mCA": format_float(
            _mean_from_stats(host_mito_context_stats["CA"])
        ),
        "host_mito_mean_5mCC": format_float(
            _mean_from_stats(host_mito_context_stats["CC"])
        ),
        "host_mito_mean_5mCT": format_float(
            _mean_from_stats(host_mito_context_stats["CT"])
        ),
    })
    for spike_name in sorted(spike_cg_cov_paths.keys()):
        spike_cg_stats = parse_cov_stats(spike_cg_cov_paths[spike_name])
        spike_ch_stats, spike_context_stats = parse_ch_cov_stats(
            spike_ch_cov_paths[spike_name]
        )
        row[f"{spike_name}_mean_5mCG"] = format_float(
            spike_cg_stats[0] if spike_cg_stats else None
        )
        row[f"{spike_name}_mean_5mCH"] = format_float(_mean_from_stats(spike_ch_stats))
        row[f"{spike_name}_mean_5mCA"] = format_float(
            _mean_from_stats(spike_context_stats["CA"])
        )
        row[f"{spike_name}_mean_5mCC"] = format_float(
            _mean_from_stats(spike_context_stats["CC"])
        )
        row[f"{spike_name}_mean_5mCT"] = format_float(
            _mean_from_stats(spike_context_stats["CT"])
        )
    row["raw_reads"] = format_optional_int(raw_reads)
    row["barcoded_reads"] = format_optional_int(barcoded_reads)
    row["barcoded_reads_rate"] = format_percentage(barcoded_reads, raw_reads)
    row["host_mapped_reads"] = format_optional_int(host_mapped_reads)
    for spike_name in sorted(spike_cg_cov_paths.keys()):
        row[f"{spike_name}_mapped_reads"] = format_optional_int(
            spike_mapped_reads.get(spike_name)
        )
    row["host_valid_reads"] = format_optional_int(host_valid_reads)
    row["host_valid_reads_mito"] = format_optional_int(host_valid_reads_mito)
    row["host_valid_reads_nuclear"] = format_optional_int(host_valid_reads_nuclear)
    row["valid_reads_rate"] = format_percentage(host_valid_reads, raw_reads)

    qc = methylated_linker_qc or {}
    row["spike_scored_reads"] = format_optional_int(
        to_optional_int(qc.get("spike_scored_reads"))
    )
    row["scored_fully_converted"] = format_optional_int(
        to_optional_int(qc.get("scored_fully_converted"))
    )
    row["spike_mean_conversion_rate"] = format_optional_float(
        to_optional_float(qc.get("spike_mean_conversion_rate"))
    )
    row["spike_fraction_all_T"] = format_optional_float(
        to_optional_float(qc.get("spike_fraction_all_T"))
    )
    row["spike_conversion_filtered_rejects"] = format_optional_int(
        to_optional_int(qc.get("spike_conversion_filtered_rejects"))
    )
    return row


def read_saturation_rate(summary_path: Path, sample_id: str) -> str | None:
    if not summary_path.exists():
        return None
    with summary_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if (row.get("sample_id") or "").strip() != sample_id:
                continue
            value = (row.get("saturation_rate") or "").strip()
            if not value:
                return None
            return value
    return None


def read_reads_threshold(summary_path: Path, sample_id: str) -> float | None:
    """HQ cutoff from saturation_summary.tsv (`reads > threshold`)."""
    if not summary_path.exists():
        return None
    with summary_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if (row.get("sample_id") or "").strip() != sample_id:
                continue
            return parse_optional_float(row.get("reads_threshold"))
    return None


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_tsv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def discover_spike_names(coverage_dir: Path) -> list[str]:
    names: list[str] = []
    for cov_path in sorted(coverage_dir.glob("*.CG.cov")):
        if cov_path.name == "host_mito.CG.cov":
            continue
        name = cov_path.name[: -len(".CG.cov")]
        if name.startswith("host"):
            continue
        if name and name not in names:
            names.append(name)
    return names


def compute_tick_values(max_index: int) -> list[int]:
    """Tick labels for 0-based inclusive indices 0..max_index."""
    if max_index < 0:
        return []
    if max_index <= 10:
        return list(range(0, max_index + 1))
    step = max(1, math.ceil((max_index + 1) / 10))
    ticks = list(range(0, max_index + 1, step))
    if ticks[-1] != max_index:
        ticks.append(max_index)
    return ticks


def write_empty_heatmap(output_path: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(2.2, 2.0))
    ax.axis("off")
    ax.text(0.5, 0.5, "No valid data", ha="center", va="center", fontsize=6)
    ax.set_title(title, fontsize=6)
    tn.style_axes(ax)
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def write_heatmap_from_rows(
    rows: list[dict[str, str]],
    field: str,
    output_path: Path,
    title: str,
    colorbar_label: str,
    cmap_name: str,
    vmin: float | None,
    vmax: float | None,
    origin_corner: str = "lower_left",
    include_row: Callable[[dict[str, str]], bool] | None = None,
) -> None:
    # x = array_col (bcB), y = array_row (bcA); both 0-based inclusive.
    # lower_left: (0,0) bottom-left; upper_right: (0,0) top-right.
    # Grid extent uses all rows with valid coords; values may be filtered.
    points: list[tuple[int, int, float]] = []
    max_col = -1
    max_row = -1
    has_coords = False
    for row in rows:
        array_row = parse_optional_coord(row.get("array_row"))
        array_col = parse_optional_coord(row.get("array_col"))
        if array_row is None or array_col is None:
            continue
        has_coords = True
        max_col = max(max_col, array_col)
        max_row = max(max_row, array_row)
        if include_row is not None and not include_row(row):
            continue
        value = parse_optional_float(row.get(field))
        if value is None:
            continue
        points.append((array_col, array_row, value))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    if not has_coords:
        write_empty_heatmap(output_path, title)
        return

    n_cols = max_col + 1
    n_rows = max_row + 1
    matrix = [[math.nan for _ in range(n_cols)] for _ in range(n_rows)]
    for array_col, array_row, value in points:
        matrix[array_row][array_col] = value

    cmap = plt.get_cmap(cmap_name).copy()
    cmap.set_bad("#f2f2f2")

    fig, ax = plt.subplots(figsize=(2.2, 2.0))
    image = ax.imshow(
        matrix,
        origin="lower",
        interpolation="nearest",
        cmap=cmap,
        aspect="equal",
        extent=(-0.5, max_col + 0.5, -0.5, max_row + 0.5),
        vmin=vmin,
        vmax=vmax,
    )
    ax.set_title(title, fontsize=6)
    ax.set_xlabel("array_col", fontsize=6)
    ax.set_ylabel("array_row", fontsize=6)

    x_ticks = compute_tick_values(max_col)
    y_ticks = compute_tick_values(max_row)
    ax.set_xticks(x_ticks)
    ax.set_yticks(y_ticks)
    if origin_corner == "upper_right":
        ax.set_xlim(max_col + 0.5, -0.5)
        ax.set_ylim(max_row + 0.5, -0.5)
    else:
        ax.set_xlim(-0.5, max_col + 0.5)
        ax.set_ylim(-0.5, max_row + 0.5)

    colorbar = fig.colorbar(image, ax=ax, shrink=0.9)
    colorbar.set_label(colorbar_label, fontsize=6)
    colorbar.ax.tick_params(labelsize=5, width=0.5, length=2)
    tn.style_axes(ax)
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def write_summary_heatmaps(
    per_spot_summary_path: Path,
    summary_dir: Path,
    reads_threshold: float | None = None,
) -> list[Path]:
    rows = read_tsv_rows(per_spot_summary_path)
    outputs: list[Path] = []

    def is_hq_spot(row: dict[str, str]) -> bool:
        return format_is_hq_spot(row.get("reads"), reads_threshold) == "1"

    for spec in HEATMAP_SPECS:
        output_path = summary_dir / spec["output_name"]
        include_row = None
        if spec.get("hq_filter") and reads_threshold is not None:
            include_row = is_hq_spot
        write_heatmap_from_rows(
            rows=rows,
            field=spec["field"],
            output_path=output_path,
            title=spec["title"],
            colorbar_label=spec["colorbar_label"],
            cmap_name=spec["cmap"],
            vmin=spec["vmin"],
            vmax=spec["vmax"],
            origin_corner=spec.get("origin_corner", "lower_left"),
            include_row=include_row,
        )
        outputs.append(output_path)
    return outputs


def main() -> int:
    args = parse_args()
    work_path = Path(args.work_path)
    coverage_dir = work_path / "coverage"
    host_cov_paths = sorted((coverage_dir / "host").rglob("*.CG.cov"))
    read_counts_path = work_path / "coverage" / "per_spot_read_counts.tsv"
    fastp_json_path = work_path / "shard_fastq" / "fastp.json"
    fastp_json_legacy_path = work_path / "fastp.json"
    demux_dir = work_path / "demux"
    saturation_summary_path = work_path / "qc" / "saturation" / "saturation_summary.tsv"
    pooled_dir = work_path / "pooled"
    host_bam_path = pooled_dir / "pooled.host.sorted.bam"
    host_mito_cg_cov_path = coverage_dir / "host_mito.CG.cov"
    host_mito_ch_cov_path = coverage_dir / "host_mito.CH.cov"
    summary_dir = work_path / "summary"
    per_spot_out = summary_dir / "per_spot_summary.tsv"
    sample_out = summary_dir / "sample_summary.tsv"
    heatmap_outputs = [summary_dir / spec["output_name"] for spec in HEATMAP_SPECS]
    mito_chromosomes = parse_mito_chromosomes(args.mito_chromosomes)

    requested_spikes = [item.strip() for item in args.spike_in_name if item.strip()]
    if requested_spikes:
        spike_names = []
        for item in requested_spikes:
            if item not in spike_names:
                spike_names.append(item)
    else:
        spike_names = discover_spike_names(coverage_dir)
    spike_cg_cov_paths = {
        spike_name: coverage_dir / f"{spike_name}.CG.cov" for spike_name in spike_names
    }
    spike_ch_cov_paths = {
        spike_name: coverage_dir / f"{spike_name}.CH.cov" for spike_name in spike_names
    }
    spike_bam_paths = {
        spike_name: pooled_dir / f"pooled.{spike_name}.sorted.bam"
        for spike_name in spike_names
    }

    print(f"[summary] work_path={work_path}")
    print(f"[summary] host_spot_cov_count={len(host_cov_paths)}")
    print(f"[summary] reads_table={read_counts_path}")
    print(f"[summary] barcode_positions={args.barcode_positions}")
    print(f"[summary] fastp_json={fastp_json_path}")
    print(f"[summary] fastp_json_legacy={fastp_json_legacy_path}")
    print(f"[summary] demux_dir={demux_dir}")
    print(f"[summary] saturation_summary={saturation_summary_path}")
    print(f"[summary] host_bam={host_bam_path}")
    print(f"[summary] host_mito_cg_cov={host_mito_cg_cov_path}")
    print(f"[summary] host_mito_ch_cov={host_mito_ch_cov_path}")
    print(f"[summary] mito_chromosomes={','.join(sorted(mito_chromosomes))}")
    print(f"[summary] spike_names={','.join(spike_names) if spike_names else 'none'}")
    for spike_name, spike_cov_path in spike_cg_cov_paths.items():
        print(f"[summary] spike_cg_cov[{spike_name}]={spike_cov_path}")
    for spike_name, spike_cov_path in spike_ch_cov_paths.items():
        print(f"[summary] spike_ch_cov[{spike_name}]={spike_cov_path}")
    for spike_name, spike_bam_path in spike_bam_paths.items():
        print(f"[summary] spike_bam[{spike_name}]={spike_bam_path}")
    print(f"[summary] output_per_spot={per_spot_out}")
    print(f"[summary] output_sample={sample_out}")
    for heatmap_path in heatmap_outputs:
        print(f"[summary] output_heatmap={heatmap_path}")

    sample_id = work_path.name
    reads_threshold = read_reads_threshold(saturation_summary_path, sample_id)
    if reads_threshold is None:
        print(
            "[summary] warning=missing_reads_threshold; "
            "mean_5mCG/mean_5mCH/mean_5mCA/mean_5mCC/mean_5mCT heatmaps use all spots"
        )
        print("[summary] hq_reads_threshold=NA")
        print("[summary] hq_filter=off")
    else:
        print(f"[summary] hq_reads_threshold={reads_threshold:g}")
        print("[summary] hq_filter=on")

    if args.dry_run:
        print("[summary] dry_run=1")
        return 0

    positions = read_barcode_positions(args.barcode_positions)
    read_counts = read_per_spot_reads(read_counts_path)
    per_spot_rows = summarize_per_spot(host_cov_paths, read_counts, positions)
    for row in per_spot_rows:
        row["is_HQ_spot"] = format_is_hq_spot(row.get("reads"), reads_threshold)
    write_tsv(per_spot_out, list(PER_SPOT_FIELDS), per_spot_rows)

    saturation_rate = read_saturation_rate(saturation_summary_path, sample_id)

    raw_reads = parse_fastp_raw_reads(fastp_json_path)
    if raw_reads is None:
        raw_reads = parse_fastp_raw_reads(fastp_json_legacy_path)
    barcoded_reads = collect_barcoded_reads(demux_dir)
    methylated_linker_qc = collect_methylated_linker_qc(demux_dir)
    (
        host_mapped_reads,
        host_valid_reads,
        host_valid_reads_mito,
        host_valid_reads_nuclear,
    ) = count_host_bam_metrics(host_bam_path, mito_chromosomes)
    spike_mapped_reads = {
        spike_name: count_mapped_reads(spike_bam_paths[spike_name])
        for spike_name in spike_names
    }

    sample_row = build_sample_summary_row(
        sample_id,
        per_spot_rows,
        host_mito_cg_cov_path,
        host_mito_ch_cov_path,
        spike_cg_cov_paths,
        spike_ch_cov_paths,
        saturation_rate,
        raw_reads,
        barcoded_reads,
        host_mapped_reads,
        host_valid_reads,
        host_valid_reads_mito,
        host_valid_reads_nuclear,
        spike_mapped_reads,
        methylated_linker_qc=methylated_linker_qc,
    )
    hmc_summary = load_hmc_spike_summary(
        work_path / "qc" / "hmc_spike" / "summary.json"
    )
    apply_hmc_spike_columns(sample_row, hmc_summary, raw_reads)
    sample_fields = build_sample_summary_fieldnames(
        sorted(spike_names), sample_row
    )
    write_tsv(sample_out, sample_fields, [sample_row])
    write_summary_heatmaps(per_spot_out, summary_dir, reads_threshold=reads_threshold)
    print(f"[summary] per_spot_rows={len(per_spot_rows)}")
    print("[summary] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
