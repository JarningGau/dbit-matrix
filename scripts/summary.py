#!/usr/bin/env python3
"""Summarize per-spot and sample-level methylation outputs after call step."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate summary TSVs from call outputs under a sample work directory."
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
    per_spot_rows: list[dict[str, str]],
    host_mito_cov_path: Path,
    spike_cov_paths: dict[str, Path],
) -> dict[str, str]:
    weighted_sum = 0.0
    total_cpg = 0
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

    host_spot_mean = (weighted_sum / total_cpg) if total_cpg > 0 else None
    host_mito_stats = parse_cov_stats(host_mito_cov_path)
    host_mito_mean = host_mito_stats[0] if host_mito_stats else None

    row = {
        "host_spot_mean_methylation": format_float(host_spot_mean),
        "host_spot_total_cpg_sites": str(total_cpg) if total_cpg > 0 else "NA",
        "host_mito_mean_methylation": format_float(host_mito_mean),
    }
    for spike_name in sorted(spike_cov_paths.keys()):
        spike_stats = parse_cov_stats(spike_cov_paths[spike_name])
        spike_mean = spike_stats[0] if spike_stats else None
        row[f"{spike_name}_mean_methylation"] = format_float(spike_mean)
    return row


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


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


def main() -> int:
    args = parse_args()
    work_path = Path(args.work_path)
    coverage_dir = work_path / "coverage"
    host_cov_paths = sorted((coverage_dir / "host").rglob("*.CG.cov"))
    read_counts_path = work_path / "split_bams" / "per_spot_read_counts.tsv"
    host_mito_cov_path = coverage_dir / "host_mito.CG.cov"
    summary_dir = work_path / "summary"
    per_spot_out = summary_dir / "per_spot_summary.tsv"
    sample_out = summary_dir / "sample_summary.tsv"

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

    print(f"[summary] work_path={work_path}")
    print(f"[summary] host_spot_cov_count={len(host_cov_paths)}")
    print(f"[summary] reads_table={read_counts_path}")
    print(f"[summary] host_mito_cov={host_mito_cov_path}")
    print(f"[summary] spike_names={','.join(spike_names) if spike_names else 'none'}")
    for spike_name, spike_cov_path in spike_cov_paths.items():
        print(f"[summary] spike_cov[{spike_name}]={spike_cov_path}")
    print(f"[summary] output_per_spot={per_spot_out}")
    print(f"[summary] output_sample={sample_out}")

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

    sample_row = build_sample_summary_row(per_spot_rows, host_mito_cov_path, spike_cov_paths)
    sample_fields = list(sample_row.keys())
    write_tsv(sample_out, sample_fields, [sample_row])
    print(f"[summary] per_spot_rows={len(per_spot_rows)}")
    print("[summary] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
