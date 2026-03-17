#!/usr/bin/env python3
"""Summarize per-spot outputs and render summary heatmaps."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pysam


HEATMAP_SPECS = (
    {
        "field": "reads",
        "output_name": "reads_heatmap.png",
        "title": "Reads per Spot",
        "colorbar_label": "Reads",
        "cmap": "viridis",
        "vmin": 0.0,
        "vmax": None,
    },
    {
        "field": "cpg_site_count",
        "output_name": "cpg_site_count_heatmap.png",
        "title": "CpG Sites per Spot",
        "colorbar_label": "CpG sites",
        "cmap": "magma",
        "vmin": 0.0,
        "vmax": None,
    },
    {
        "field": "mean_methylation",
        "output_name": "mean_methylation_heatmap.png",
        "title": "Mean Methylation per Spot",
        "colorbar_label": "Mean methylation (%)",
        "cmap": "coolwarm",
        "vmin": 0.0,
        "vmax": 100.0,
    },
)
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
        "--dry-run",
        action="store_true",
        help="Print resolved inputs/outputs and exit without writing files.",
    )
    return parser.parse_args()


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
            if len(fields) < 4:
                continue
            try:
                methylation = float(fields[3])
            except ValueError:
                continue
            methylation_sum += methylation
            cpg_count += 1
    if cpg_count == 0:
        return None
    return methylation_sum / cpg_count, cpg_count


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


def count_host_bam_metrics(host_bam_path: Path) -> tuple[int | None, int | None]:
    if not host_bam_path.exists():
        return None, None
    mapped_reads = 0
    valid_reads = 0
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
                if read.flag in VALID_FLAGS:
                    valid_reads += 1
    except OSError:
        return None, None
    return mapped_reads, valid_reads


def count_mapped_reads(bam_path: Path) -> int | None:
    if not bam_path.exists():
        return None
    mapped_reads = 0
    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam_handle:
            for read in bam_handle.fetch(until_eof=True):
                if not read.is_unmapped:
                    mapped_reads += 1
    except OSError:
        return None
    return mapped_reads


def read_per_spot_reads(read_counts_path: Path) -> dict[str, tuple[str, str, str]]:
    if not read_counts_path.exists():
        return {}
    rows: dict[str, tuple[str, str, str]] = {}
    with read_counts_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            spot = (row.get("spot") or "").strip()
            if not spot:
                continue
            x_index = (row.get("X_index") or "").strip()
            y_index = (row.get("Y_index") or "").strip()
            reads = (row.get("reads") or "").strip()
            rows[spot] = (x_index, y_index, reads)
    return rows


def parse_xy_from_spot(spot: str) -> tuple[str, str]:
    if "_" not in spot:
        return "NA", "NA"
    left, right = spot.split("_", 1)
    left = left.strip()
    right = right.strip()
    if not left or not right:
        return "NA", "NA"
    try:
        return str(int(left)), str(int(right))
    except ValueError:
        return left, right


def format_float(value: float | None) -> str:
    if value is None:
        return "NA"
    return f"{value:.6f}"


def parse_optional_float(text: str | None) -> float | None:
    value = (text or "").strip()
    if not value or value == "NA":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_optional_index(text: str | None) -> int | None:
    value = (text or "").strip()
    if not value or value == "NA":
        return None
    try:
        parsed = int(value)
    except ValueError:
        return None
    if parsed <= 0:
        return None
    return parsed


def summarize_per_spot(
    host_cov_paths: list[Path],
    read_counts: dict[str, tuple[str, str, str]],
) -> list[dict[str, str]]:
    by_cov: dict[str, tuple[float | None, int | None]] = {}
    for cov_path in host_cov_paths:
        spot = cov_path.name[: -len(".CG.cov")]
        cov_stats = parse_cov_stats(cov_path)
        if cov_stats is None:
            by_cov[spot] = (None, None)
        else:
            mean_methylation, cpg_count = cov_stats
            by_cov[spot] = (mean_methylation, cpg_count)

    all_spots = sorted(set(by_cov.keys()) | set(read_counts.keys()))
    rows: list[dict[str, str]] = []
    for spot in all_spots:
        mean_methylation, cpg_count = by_cov.get(spot, (None, None))
        x_index, y_index = parse_xy_from_spot(spot)
        reads = "NA"
        if spot in read_counts:
            x_index, y_index, reads = read_counts[spot]
            if not x_index or not y_index:
                guessed_x, guessed_y = parse_xy_from_spot(spot)
                x_index = x_index or guessed_x
                y_index = y_index or guessed_y
            reads = reads or "NA"
        rows.append(
            {
                "X_index": x_index if x_index else "NA",
                "Y_index": y_index if y_index else "NA",
                "spot": spot,
                "mean_methylation": format_float(mean_methylation),
                "cpg_site_count": str(cpg_count) if cpg_count is not None else "NA",
                "reads": reads,
            }
        )
    return rows


def build_sample_summary_row(
    sample_id: str,
    per_spot_rows: list[dict[str, str]],
    host_mito_cov_path: Path,
    spike_cov_paths: dict[str, Path],
    raw_reads: int | None,
    barcoded_reads: int | None,
    host_mapped_reads: int | None,
    host_valid_reads: int | None,
    spike_mapped_reads: dict[str, int | None],
) -> dict[str, str]:
    weighted_sum = 0.0
    total_cpg = 0
    cpg_counts_for_median: list[int] = []
    for row in per_spot_rows:
        mean_text = row["mean_methylation"]
        cpg_text = row["cpg_site_count"]
        if mean_text == "NA" or cpg_text == "NA":
            continue
        cpg_count = int(cpg_text)
        if cpg_count <= 0:
            continue
        mean_value = float(mean_text)
        weighted_sum += mean_value * cpg_count
        total_cpg += cpg_count
        cpg_counts_for_median.append(cpg_count)

    host_spot_mean = (weighted_sum / total_cpg) if total_cpg > 0 else None
    host_spot_median_cpg: float | None
    if not cpg_counts_for_median:
        host_spot_median_cpg = None
    else:
        sorted_counts = sorted(cpg_counts_for_median)
        n = len(sorted_counts)
        mid = n // 2
        if n % 2 == 1:
            host_spot_median_cpg = float(sorted_counts[mid])
        else:
            host_spot_median_cpg = (sorted_counts[mid - 1] + sorted_counts[mid]) / 2.0

    host_mito_stats = parse_cov_stats(host_mito_cov_path)
    host_mito_mean = host_mito_stats[0] if host_mito_stats else None

    row: dict[str, str] = {"sample_id": sample_id}
    row.update({
        "host_spot_mean_methylation": format_float(host_spot_mean),
        "host_spot_median_cpg_sites": format_float(host_spot_median_cpg),
        "host_mito_mean_methylation": format_float(host_mito_mean),
    })
    for spike_name in sorted(spike_cov_paths.keys()):
        spike_stats = parse_cov_stats(spike_cov_paths[spike_name])
        spike_mean = spike_stats[0] if spike_stats else None
        row[f"{spike_name}_mean_methylation"] = format_float(spike_mean)
    row["raw_reads"] = format_optional_int(raw_reads)
    row["barcoded_reads"] = format_optional_int(barcoded_reads)
    row["barcoded_reads_rate"] = format_percentage(barcoded_reads, raw_reads)
    row["host_mapped_reads"] = format_optional_int(host_mapped_reads)
    for spike_name in sorted(spike_cov_paths.keys()):
        row[f"{spike_name}_mapped_reads"] = format_optional_int(
            spike_mapped_reads.get(spike_name)
        )
    row["host_valid_reads"] = format_optional_int(host_valid_reads)
    row["valid_reads_rate"] = format_percentage(host_valid_reads, raw_reads)
    return row


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
    if max_index <= 0:
        return []
    if max_index <= 10:
        return list(range(1, max_index + 1))
    step = max(1, math.ceil(max_index / 10))
    ticks = list(range(1, max_index + 1, step))
    if ticks[-1] != max_index:
        ticks.append(max_index)
    return ticks


def write_empty_heatmap(output_path: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.axis("off")
    ax.text(0.5, 0.5, "No valid data", ha="center", va="center", fontsize=12)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
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
) -> None:
    points: list[tuple[int, int, float]] = []
    max_x = 0
    max_y = 0
    for row in rows:
        x_index = parse_optional_index(row.get("X_index"))
        y_index = parse_optional_index(row.get("Y_index"))
        value = parse_optional_float(row.get(field))
        if x_index is None or y_index is None or value is None:
            continue
        points.append((x_index, y_index, value))
        max_x = max(max_x, x_index)
        max_y = max(max_y, y_index)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    if not points:
        write_empty_heatmap(output_path, title)
        return

    matrix = [[math.nan for _ in range(max_x)] for _ in range(max_y)]
    for x_index, y_index, value in points:
        matrix[y_index - 1][x_index - 1] = value

    cmap = plt.get_cmap(cmap_name).copy()
    cmap.set_bad("#f2f2f2")

    fig_width = min(max(6.0, max_x * 0.22), 12.0)
    fig_height = min(max(5.0, max_y * 0.22), 12.0)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    image = ax.imshow(
        matrix,
        origin="lower",
        interpolation="nearest",
        cmap=cmap,
        aspect="equal",
        extent=(0.5, max_x + 0.5, 0.5, max_y + 0.5),
        vmin=vmin,
        vmax=vmax,
    )
    ax.set_title(title)
    ax.set_xlabel("X index")
    ax.set_ylabel("Y index")

    x_ticks = compute_tick_values(max_x)
    y_ticks = compute_tick_values(max_y)
    ax.set_xticks(x_ticks)
    ax.set_yticks(y_ticks)
    ax.set_xlim(0.5, max_x + 0.5)
    ax.set_ylim(0.5, max_y + 0.5)

    colorbar = fig.colorbar(image, ax=ax, shrink=0.9)
    colorbar.set_label(colorbar_label)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def write_summary_heatmaps(per_spot_summary_path: Path, summary_dir: Path) -> list[Path]:
    rows = read_tsv_rows(per_spot_summary_path)
    outputs: list[Path] = []
    for spec in HEATMAP_SPECS:
        output_path = summary_dir / spec["output_name"]
        write_heatmap_from_rows(
            rows=rows,
            field=spec["field"],
            output_path=output_path,
            title=spec["title"],
            colorbar_label=spec["colorbar_label"],
            cmap_name=spec["cmap"],
            vmin=spec["vmin"],
            vmax=spec["vmax"],
        )
        outputs.append(output_path)
    return outputs


def main() -> int:
    args = parse_args()
    work_path = Path(args.work_path)
    coverage_dir = work_path / "coverage"
    host_cov_paths = sorted((coverage_dir / "host").rglob("*.CG.cov"))
    read_counts_path = work_path / "split_bams" / "per_spot_read_counts.tsv"
    fastp_json_path = work_path / "shard_fastq" / "fastp.json"
    fastp_json_legacy_path = work_path / "fastp.json"
    demux_dir = work_path / "demux"
    pooled_dir = work_path / "pooled"
    host_bam_path = pooled_dir / "pooled.byCB.bam"
    host_mito_cov_path = coverage_dir / "host_mito.CG.cov"
    summary_dir = work_path / "summary"
    per_spot_out = summary_dir / "per_spot_summary.tsv"
    sample_out = summary_dir / "sample_summary.tsv"
    heatmap_outputs = [summary_dir / spec["output_name"] for spec in HEATMAP_SPECS]

    requested_spikes = [item.strip() for item in args.spike_in_name if item.strip()]
    if requested_spikes:
        spike_names = []
        for item in requested_spikes:
            if item not in spike_names:
                spike_names.append(item)
    else:
        spike_names = discover_spike_names(coverage_dir)
    spike_cov_paths = {
        spike_name: coverage_dir / f"{spike_name}.CG.cov" for spike_name in spike_names
    }
    spike_bam_paths = {
        spike_name: pooled_dir / f"pooled.{spike_name}.sorted.bam"
        for spike_name in spike_names
    }

    print(f"[summary] work_path={work_path}")
    print(f"[summary] host_spot_cov_count={len(host_cov_paths)}")
    print(f"[summary] reads_table={read_counts_path}")
    print(f"[summary] fastp_json={fastp_json_path}")
    print(f"[summary] fastp_json_legacy={fastp_json_legacy_path}")
    print(f"[summary] demux_dir={demux_dir}")
    print(f"[summary] host_bam={host_bam_path}")
    print(f"[summary] host_mito_cov={host_mito_cov_path}")
    print(f"[summary] spike_names={','.join(spike_names) if spike_names else 'none'}")
    for spike_name, spike_cov_path in spike_cov_paths.items():
        print(f"[summary] spike_cov[{spike_name}]={spike_cov_path}")
    for spike_name, spike_bam_path in spike_bam_paths.items():
        print(f"[summary] spike_bam[{spike_name}]={spike_bam_path}")
    print(f"[summary] output_per_spot={per_spot_out}")
    print(f"[summary] output_sample={sample_out}")
    for heatmap_path in heatmap_outputs:
        print(f"[summary] output_heatmap={heatmap_path}")

    if args.dry_run:
        print("[summary] dry_run=1")
        return 0

    read_counts = read_per_spot_reads(read_counts_path)
    per_spot_rows = summarize_per_spot(host_cov_paths, read_counts)
    per_spot_fields = [
        "X_index",
        "Y_index",
        "spot",
        "mean_methylation",
        "cpg_site_count",
        "reads",
    ]
    write_tsv(per_spot_out, per_spot_fields, per_spot_rows)

    raw_reads = parse_fastp_raw_reads(fastp_json_path)
    if raw_reads is None:
        raw_reads = parse_fastp_raw_reads(fastp_json_legacy_path)
    barcoded_reads = collect_barcoded_reads(demux_dir)
    host_mapped_reads, host_valid_reads = count_host_bam_metrics(host_bam_path)
    spike_mapped_reads = {
        spike_name: count_mapped_reads(spike_bam_paths[spike_name])
        for spike_name in spike_names
    }

    sample_row = build_sample_summary_row(
        work_path.name,
        per_spot_rows,
        host_mito_cov_path,
        spike_cov_paths,
        raw_reads,
        barcoded_reads,
        host_mapped_reads,
        host_valid_reads,
        spike_mapped_reads,
    )
    sample_fields = list(sample_row.keys())
    write_tsv(sample_out, sample_fields, [sample_row])
    write_summary_heatmaps(per_spot_out, summary_dir)
    print(f"[summary] per_spot_rows={len(per_spot_rows)}")
    print("[summary] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
