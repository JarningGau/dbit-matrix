#!/usr/bin/env python3
"""Run saturation QC from host CpG coverage outputs."""

from __future__ import annotations

import argparse
import csv
import json
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
DEFAULT_LINEAR_R2_THRESHOLD = 0.99


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
        "--linear-r2-threshold",
        type=float,
        default=DEFAULT_LINEAR_R2_THRESHOLD,
        help=(
            "If the linear (through-origin) fit reaches this R^2, use linear "
            "extrapolation; otherwise fall back to the saturation curve. "
            "Default: 0.99."
        ),
    )
    parser.add_argument(
        "--fastp-json",
        default=None,
        help=(
            "fastp JSON report for total sequencing depth on the plot x-axis. "
            "Default: auto-detect <work_path>/shard_fastq/fastp.json. If missing, "
            "the x-axis falls back to raw coverage fraction."
        ),
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


def discover_fastp_json(work_path: Path) -> Path | None:
    candidate = work_path / "shard_fastq" / "fastp.json"
    return candidate if candidate.is_file() else None


def resolve_fastp_json(work_path: Path, fastp_json: str | None) -> Path | None:
    if fastp_json:
        path = Path(fastp_json)
        if not path.is_file():
            raise FileNotFoundError(f"fastp JSON not found: {path}")
        return path
    return discover_fastp_json(work_path)


def load_fastp_sequencing_gbp(fastp_json: Path) -> float:
    data = json.loads(fastp_json.read_text(encoding="utf-8"))
    summary = data.get("summary", {})
    for section in ("after_filtering", "before_filtering"):
        total_bases = summary.get(section, {}).get("total_bases")
        if total_bases is not None:
            bases = int(total_bases)
            if bases > 0:
                return bases / 1e9
    raise ValueError(f"no total_bases found in fastp report: {fastp_json}")


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


def median_and_iqr(values: list[float]) -> tuple[float, float, float]:
    if not values:
        raise ValueError("empty values for median_and_iqr")
    median = float(statistics.median(values))
    if len(values) == 1:
        return median, 0.0, 0.0
    q1, _, q3 = statistics.quantiles(values, n=4, method="inclusive")
    return median, max(0.0, median - q1), max(0.0, q3 - median)


def sat_func(fraction: float, a: float, b: float) -> float:
    return a * (1.0 - math.exp(-b * fraction))


def fit_linear_through_origin(fractions: list[float], y_values: list[float]) -> float:
    denominator = sum(f * f for f in fractions)
    if denominator <= 0:
        return 0.0
    return sum(f * y for f, y in zip(fractions, y_values)) / denominator


def r_squared(y_values: list[float], predictions: list[float]) -> float:
    if not y_values:
        return 0.0
    mean_y = sum(y_values) / len(y_values)
    sst = sum((y - mean_y) ** 2 for y in y_values)
    sse = sum((y - p) ** 2 for y, p in zip(y_values, predictions))
    if sst <= 0:
        return 1.0 if sse <= 0 else 0.0
    return 1.0 - sse / sst


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
        "extrapolation_model",
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
    err_low: list[float],
    err_high: list[float],
    model: str,
    a: float,
    b: float,
    slope: float,
    pred_fraction: float,
    predicted: float | None,
    saturation_rate: float | None,
    sequencing_gbp: float | None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    scale = 1e4

    # When fastp depth is available the x-axis is sequencing depth in Gbp,
    # otherwise fall back to the raw coverage fraction.
    x_scale = sequencing_gbp if sequencing_gbp else 1.0
    x_label = "Sequencing depth (Gbp)" if sequencing_gbp else "Coverage Fraction"

    x_values = [fraction * x_scale for fraction in fractions]
    max_fraction = max(max(fractions), pred_fraction, 1.0)
    fit_fractions = [max_fraction * index / 250.0 for index in range(251)]
    fit_x = [fraction * x_scale for fraction in fit_fractions]
    if model == "linear":
        fit_y = [slope * value / scale for value in fit_fractions]
        fit_label = "Fitted line (linear)"
    else:
        fit_y = [sat_func(value, a, b) / scale for value in fit_fractions]
        fit_label = "Fitted curve (saturation)"
    pred_x = pred_fraction * x_scale
    pred_y = (predicted if predicted is not None else 0.0) / scale
    medians_scaled = [value / scale for value in median_uniques]
    yerr = [
        [value / scale for value in err_low],
        [value / scale for value in err_high],
    ]

    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    ax.errorbar(
        x_values,
        medians_scaled,
        yerr=yerr,
        fmt="o-",
        color="blue",
        linewidth=2,
        markersize=5,
        capsize=4,
        label="Observed median (IQR)",
    )
    ax.plot(fit_x, fit_y, "r--", linewidth=2, label=fit_label)
    ax.scatter(
        [pred_x],
        [pred_y],
        color="green",
        s=55,
        zorder=4,
        label=f"Prediction at {pred_fraction:g}x",
    )
    if model != "linear" and a > 0:
        ax.axhline(
            a / scale,
            color="purple",
            linestyle="--",
            linewidth=1.8,
            label="Max CpG Reads",
        )
    if saturation_rate is not None:
        sat_text = f"{saturation_rate:.2f}%"
    elif model == "linear":
        sat_text = "NA (linear, unsaturated)"
    else:
        sat_text = "NA"
    ax.set_title(f"Saturation analysis ({sample_id})\nSaturation rate: {sat_text}")
    ax.set_xlabel(x_label)
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

    fastp_json_path = resolve_fastp_json(work_path, args.fastp_json)
    sequencing_gbp: float | None = None
    if fastp_json_path is not None:
        sequencing_gbp = load_fastp_sequencing_gbp(fastp_json_path)
    else:
        print("[saturation] warning=fastp_json_not_found")

    print(f"[saturation] work_path={work_path}")
    print(f"[saturation] reads_table={read_counts_path}")
    print(f"[saturation] host_cov_count={len(host_cov_paths)}")
    print(f"[saturation] reads_threshold={args.reads_threshold}")
    print(f"[saturation] linear_r2_threshold={args.linear_r2_threshold}")
    print(f"[saturation] fastp_json={fastp_json_path}")
    print(
        "[saturation] sequencing_gbp="
        + (f"{sequencing_gbp:.6f}" if sequencing_gbp is not None else "NA")
    )
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
            "extrapolation_model": "NA",
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
    err_low: list[float] = []
    err_high: list[float] = []
    for fraction in DEFAULT_FRACTIONS:
        values = per_fraction_values[fraction]
        if not values:
            continue
        median, low, high = median_and_iqr(values)
        usable_fractions.append(fraction)
        median_uniques.append(median)
        err_low.append(low)
        err_high.append(high)

    if not usable_fractions:
        print("[saturation] warning=no_coverage_values_for_hq_spots")
        row = {
            "sample_id": sample_id,
            "observed_median_unique_cpgs": "NA",
            "theoretical_max_median_unique_cpgs": "NA",
            "predicted_median_unique_cpgs_at_2x": "NA",
            "saturation_rate": "NA",
            "extrapolation_model": "NA",
            "hq_spot_count": str(len(hq_spots)),
        }
        write_summary_tsv(summary_path, row)
        write_empty_plot(plot_path, sample_id, "No valid CG coverage rows")
        print("[saturation] done")
        return 0

    observed = median_uniques[-1]
    a, b = fit_saturation_curve(usable_fractions, median_uniques)
    slope = fit_linear_through_origin(usable_fractions, median_uniques)
    linear_r2 = r_squared(median_uniques, [slope * f for f in usable_fractions])
    exp_r2 = r_squared(median_uniques, [sat_func(f, a, b) for f in usable_fractions])

    # Prefer the linear (through-origin) extrapolation when the observed curve
    # is essentially linear (undersaturated): in that regime the exponential
    # asymptote/saturation rate is unidentifiable and unreliable. Only when the
    # linear fit is poor (clear curvature) do we trust the saturation curve.
    use_linear = linear_r2 >= args.linear_r2_threshold
    if use_linear:
        model = "linear"
        theoretical = None
        predicted_2x = slope * args.pred_fraction
        saturation_rate = None
    else:
        model = "saturation"
        theoretical = a if a > 0 else None
        predicted_2x = sat_func(args.pred_fraction, a, b) if a > 0 else None
        saturation_rate = None
        if theoretical is not None and theoretical > 0:
            saturation_rate = observed / theoretical * 100.0

    print(f"[saturation] linear_r2={linear_r2:.6f} exp_r2={exp_r2:.6f}")
    print(f"[saturation] extrapolation_model={model}")
    print(f"[saturation] saturation_rate={format_optional_float(saturation_rate)}")

    write_plot(
        path=plot_path,
        sample_id=sample_id,
        fractions=usable_fractions,
        median_uniques=median_uniques,
        err_low=err_low,
        err_high=err_high,
        model=model,
        a=a,
        b=b,
        slope=slope,
        pred_fraction=args.pred_fraction,
        predicted=predicted_2x,
        saturation_rate=saturation_rate,
        sequencing_gbp=sequencing_gbp,
    )
    row = {
        "sample_id": sample_id,
        "observed_median_unique_cpgs": format_optional_int(observed),
        "theoretical_max_median_unique_cpgs": format_optional_int(theoretical),
        "predicted_median_unique_cpgs_at_2x": format_optional_int(predicted_2x),
        "saturation_rate": format_optional_float(saturation_rate),
        "extrapolation_model": model,
        "hq_spot_count": str(len(hq_spots)),
    }
    write_summary_tsv(summary_path, row)
    print(f"[saturation] hq_spot_count={len(hq_spots)}")
    print("[saturation] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
