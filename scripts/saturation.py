#!/usr/bin/env python3
"""Run saturation QC from host CpG coverage outputs."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


DEFAULT_FRACTIONS = [
    0.01,
    0.02,
    0.05,
    0.10,
    0.20,
    0.30,
    0.40,
    0.50,
    0.60,
    0.70,
    0.80,
    0.90,
    1.00,
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Estimate host CpG saturation curve and write plot/summary under "
            "<work_path>/qc/saturation."
        )
    )
    parser.add_argument(
        "--work-path",
        required=True,
        help="Sample work directory containing split_bams/ and coverage/host/ outputs.",
    )
    parser.add_argument(
        "--reads-threshold",
        type=float,
        default=1_000_000.0,
        help="HQ spot threshold for per_spot_read_counts.tsv reads. Default: 1e6.",
    )
    parser.add_argument(
        "--pred-fraction",
        type=float,
        default=2.0,
        help="Coverage fraction used for forward prediction point. Default: 2.0.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory. Default: <work_path>/qc/saturation.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print resolved inputs/outputs and exit without writing files.",
    )
    return parser.parse_args()


def format_optional_int(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return "NA"
    return str(int(round(value)))


def format_optional_float(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return "NA"
    return f"{value:.6f}"


def parse_per_spot_reads(path: Path) -> dict[str, float]:
    if not path.exists():
        return {}
    rows: dict[str, float] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            spot = (row.get("spot") or "").strip()
            if not spot:
                continue
            reads_text = (row.get("reads") or "").strip()
            if not reads_text:
                continue
            try:
                reads = float(reads_text)
            except ValueError:
                continue
            if reads > 0:
                rows[spot] = reads
    return rows


def parse_cov_histogram(cov_path: Path) -> dict[int, int]:
    hist: dict[int, int] = {}
    with cov_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 6:
                continue
            try:
                methylated_count = int(fields[4])
                unmethylated_count = int(fields[5])
            except ValueError:
                continue
            depth = methylated_count + unmethylated_count
            if depth <= 0:
                continue
            hist[depth] = hist.get(depth, 0) + 1
    return hist


def expected_unique(hist: dict[int, int], fraction: float) -> float:
    total = 0.0
    for depth, count in hist.items():
        total += count * (1.0 - (1.0 - fraction) ** depth)
    return total


def sat_func(fraction: float, a: float, b: float) -> float:
    return a * (1.0 - math.exp(-b * fraction))


def fit_saturation_curve(fractions: list[float], y_values: list[float]) -> tuple[float, float]:
    if len(fractions) != len(y_values):
        raise ValueError("fractions and y_values length mismatch")
    if not fractions:
        raise ValueError("empty inputs for fitting")

    best_a = 0.0
    best_b = 1.0
    best_error = math.inf

    # Avoid extra dependencies: grid-search b and solve optimal a in closed-form.
    def search_b(log10_low: float, log10_high: float, count: int) -> None:
        nonlocal best_a, best_b, best_error
        if count <= 1:
            return
        step = (log10_high - log10_low) / float(count - 1)
        for index in range(count):
            b = 10 ** (log10_low + step * index)
            transformed = [1.0 - math.exp(-b * f) for f in fractions]
            denominator = sum(value * value for value in transformed)
            if denominator <= 0:
                continue
            a = sum(y * value for y, value in zip(y_values, transformed)) / denominator
            a = max(0.0, a)
            error = sum(
                (y - a * value) * (y - a * value)
                for y, value in zip(y_values, transformed)
            )
            if error < best_error:
                best_error = error
                best_a = a
                best_b = b

    search_b(log10_low=-4.0, log10_high=2.0, count=800)
    for _ in range(3):
        center = math.log10(best_b if best_b > 0 else 1.0)
        search_b(log10_low=center - 0.7, log10_high=center + 0.7, count=240)

    return best_a, best_b


def write_summary_tsv(path: Path, row: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sample_id",
        "observed_median_unique_cpgs",
        "theoretical_max_median_unique_cpgs",
        "predicted_median_unique_cpgs_at_2x",
        "saturation_rate",
        "hq_spot_count",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


def write_empty_plot(path: Path, sample_id: str, message: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    ax.axis("off")
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=12)
    ax.set_title(f"Saturation analysis ({sample_id})\nSaturation rate: NA")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def write_plot(
    path: Path,
    sample_id: str,
    fractions: list[float],
    median_uniques: list[float],
    a: float,
    b: float,
    pred_fraction: float,
    saturation_rate: float | None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    scale = 1e4

    max_x = max(max(fractions), pred_fraction, 1.0)
    fit_x = [max_x * index / 250.0 for index in range(251)]
    fit_y = [sat_func(value, a, b) / scale for value in fit_x]
    pred_y = sat_func(pred_fraction, a, b) / scale
    max_y = a / scale

    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    ax.plot(
        fractions,
        [value / scale for value in median_uniques],
        "o-",
        color="blue",
        linewidth=2,
        markersize=5,
        label="Observed median",
    )
    ax.plot(fit_x, fit_y, "r--", linewidth=2, label="Fitted curve")
    ax.scatter(
        [pred_fraction],
        [pred_y],
        color="green",
        s=55,
        zorder=4,
        label=f"Prediction at {pred_fraction:g}x",
    )
    ax.axhline(
        max_y,
        color="purple",
        linestyle="--",
        linewidth=1.8,
        label="Max CpG Reads",
    )
    if saturation_rate is None:
        sat_text = "NA"
    else:
        sat_text = f"{saturation_rate:.2f}%"
    ax.set_title(f"Saturation analysis ({sample_id})\nSaturation rate: {sat_text}")
    ax.set_xlabel("Coverage Fraction")
    ax.set_ylabel("Median Unique CpGs per Spot (x10^4)")
    ax.grid(True, alpha=0.35)
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    work_path = Path(args.work_path)
    sample_id = work_path.name
    read_counts_path = work_path / "split_bams" / "per_spot_read_counts.tsv"
    host_cov_paths = sorted((work_path / "coverage" / "host").rglob("*.CG.cov"))
    output_dir = Path(args.output_dir) if args.output_dir else work_path / "qc" / "saturation"
    plot_path = output_dir / "saturation_curve.png"
    summary_path = output_dir / "saturation_summary.tsv"

    print(f"[saturation] work_path={work_path}")
    print(f"[saturation] reads_table={read_counts_path}")
    print(f"[saturation] host_cov_count={len(host_cov_paths)}")
    print(f"[saturation] reads_threshold={args.reads_threshold}")
    print(f"[saturation] output_plot={plot_path}")
    print(f"[saturation] output_summary={summary_path}")

    if args.dry_run:
        print("[saturation] dry_run=1")
        return 0

    read_counts = parse_per_spot_reads(read_counts_path)
    spot_cov_map: dict[str, Path] = {
        cov_path.name[: -len(".CG.cov")]: cov_path for cov_path in host_cov_paths
    }
    hq_spots = [
        spot
        for spot, reads in read_counts.items()
        if reads > args.reads_threshold and spot in spot_cov_map
    ]
    hq_spots.sort()

    if not hq_spots:
        print("[saturation] warning=no_hq_spots_after_filter")
        row = {
            "sample_id": sample_id,
            "observed_median_unique_cpgs": "NA",
            "theoretical_max_median_unique_cpgs": "NA",
            "predicted_median_unique_cpgs_at_2x": "NA",
            "saturation_rate": "NA",
            "hq_spot_count": "0",
        }
        write_summary_tsv(summary_path, row)
        write_empty_plot(plot_path, sample_id, "No HQ spots after reads filter")
        print("[saturation] done")
        return 0

    per_fraction_values: dict[float, list[float]] = {fraction: [] for fraction in DEFAULT_FRACTIONS}
    for spot in hq_spots:
        hist = parse_cov_histogram(spot_cov_map[spot])
        if not hist:
            continue
        for fraction in DEFAULT_FRACTIONS:
            per_fraction_values[fraction].append(expected_unique(hist, fraction))

    usable_fractions: list[float] = []
    median_uniques: list[float] = []
    for fraction in DEFAULT_FRACTIONS:
        values = per_fraction_values[fraction]
        if not values:
            continue
        usable_fractions.append(fraction)
        median_uniques.append(float(statistics.median(values)))

    if not usable_fractions:
        print("[saturation] warning=no_coverage_values_for_hq_spots")
        row = {
            "sample_id": sample_id,
            "observed_median_unique_cpgs": "NA",
            "theoretical_max_median_unique_cpgs": "NA",
            "predicted_median_unique_cpgs_at_2x": "NA",
            "saturation_rate": "NA",
            "hq_spot_count": str(len(hq_spots)),
        }
        write_summary_tsv(summary_path, row)
        write_empty_plot(plot_path, sample_id, "No valid CG coverage rows")
        print("[saturation] done")
        return 0

    a, b = fit_saturation_curve(usable_fractions, median_uniques)
    observed = median_uniques[-1]
    theoretical = a if a > 0 else None
    predicted_2x = sat_func(args.pred_fraction, a, b) if a > 0 else None
    saturation_rate = None
    if theoretical is not None and theoretical > 0:
        saturation_rate = observed / theoretical * 100.0

    write_plot(
        path=plot_path,
        sample_id=sample_id,
        fractions=usable_fractions,
        median_uniques=median_uniques,
        a=a,
        b=b,
        pred_fraction=args.pred_fraction,
        saturation_rate=saturation_rate,
    )
    row = {
        "sample_id": sample_id,
        "observed_median_unique_cpgs": format_optional_int(observed),
        "theoretical_max_median_unique_cpgs": format_optional_int(theoretical),
        "predicted_median_unique_cpgs_at_2x": format_optional_int(predicted_2x),
        "saturation_rate": format_optional_float(saturation_rate),
        "hq_spot_count": str(len(hq_spots)),
    }
    write_summary_tsv(summary_path, row)
    print(f"[saturation] hq_spot_count={len(hq_spots)}")
    print("[saturation] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
