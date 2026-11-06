#!/usr/bin/env python3
"""Run saturation QC from host CpG coverage outputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path
from typing import NamedTuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import theme_nature as tn

tn.apply_theme()

_FIGSIZE = (2.5, 1.8)

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
KNEE_MIN_READS = 10.0
KNEE_SMOOTH_WINDOW = 51
KNEE_EXCLUDE_TOP_FRACTION = 0.05
KNEE_EXCLUDE_BOTTOM_FRACTION = 0.10
KNEE_HALF_STEEP_FRACTION = 0.5


class AutoKneeResult(NamedTuple):
    threshold: float
    knee_rank: int  # 1-based barcode rank on the truncated curve
    exclude_start: int  # 0-based index where slope search begins
    exclude_end: int  # 0-based exclusive end index of slope search


def _rolling_mean(values: list[float], window: int) -> list[float]:
    if window % 2 == 0:
        raise ValueError("smooth window must be odd")
    if window < 1:
        raise ValueError("smooth window must be >= 1")
    if not values:
        return []
    if len(values) == 1 or window == 1:
        return list(values)
    pad = window // 2
    padded = [values[0]] * pad + list(values) + [values[-1]] * pad
    inv = 1.0 / float(window)
    out: list[float] = []
    running = sum(padded[:window])
    out.append(running * inv)
    for index in range(window, len(padded)):
        running += padded[index] - padded[index - window]
        out.append(running * inv)
    return out


def auto_reads_threshold(
    reads: list[float],
    *,
    min_reads: float = KNEE_MIN_READS,
    smooth_window: int = KNEE_SMOOTH_WINDOW,
    exclude_top_fraction: float = KNEE_EXCLUDE_TOP_FRACTION,
    exclude_bottom_fraction: float = KNEE_EXCLUDE_BOTTOM_FRACTION,
    half_steep_fraction: float = KNEE_HALF_STEEP_FRACTION,
) -> AutoKneeResult:
    """Return HQ reads cutoff on (rank, log10(reads+1)).

    Spots are sorted descending; reads <= min_reads are dropped. The top
    ``exclude_top_fraction`` and bottom ``exclude_bottom_fraction`` of
    remaining ranks are excluded from the slope search. Among searchable
    ranks, candidates with slope at least ``half_steep_fraction`` as steep
    as the minimum slope are kept, and the leftmost such rank is the knee
    (avoids a steeper right-hand tail cliff).
    """
    sorted_reads = sorted((float(value) for value in reads), reverse=True)
    truncated = [value for value in sorted_reads if value > min_reads]
    if len(truncated) < 2:
        raise ValueError(
            "auto reads threshold needs at least 2 spots with reads > "
            f"{min_reads:g} (got {len(truncated)})"
        )
    n_truncated = len(truncated)
    exclude_start = int(math.floor(n_truncated * exclude_top_fraction))
    exclude_end = n_truncated - int(math.floor(n_truncated * exclude_bottom_fraction))
    if exclude_end - exclude_start < 2:
        raise ValueError(
            "auto reads threshold needs at least 2 searchable spots after "
            f"excluding top {exclude_top_fraction:.0%} and bottom "
            f"{exclude_bottom_fraction:.0%} "
            f"(n={n_truncated}, exclude_start={exclude_start}, "
            f"exclude_end={exclude_end})"
        )
    log_reads = [math.log10(value + 1.0) for value in truncated]
    window = smooth_window if smooth_window < n_truncated else (
        n_truncated if n_truncated % 2 == 1 else max(1, n_truncated - 1)
    )
    if window % 2 == 0:
        window = max(1, window - 1)
    smoothed = _rolling_mean(log_reads, window)

    # Ranks are unit-spaced, so dy/drank ~= dy via central differences.
    slopes: list[float] = []
    for index in range(len(smoothed)):
        if index == 0:
            slopes.append(smoothed[1] - smoothed[0])
        elif index == len(smoothed) - 1:
            slopes.append(smoothed[-1] - smoothed[-2])
        else:
            slopes.append(0.5 * (smoothed[index + 1] - smoothed[index - 1]))

    search = range(exclude_start, exclude_end)
    min_slope = min(slopes[index] for index in search)
    slope_cutoff = half_steep_fraction * min_slope
    candidates = [index for index in search if slopes[index] <= slope_cutoff]
    if not candidates:
        raise ValueError("auto reads threshold found no half-steep slope candidates")
    knee_index = min(candidates)
    return AutoKneeResult(
        threshold=float(truncated[knee_index]),
        knee_rank=knee_index + 1,
        exclude_start=exclude_start,
        exclude_end=exclude_end,
    )


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
        default="auto",
        help=(
            "HQ spot threshold for per_spot_read_counts.tsv reads. "
            "Use 'auto' (default) to pick the leftmost half-steep drop on "
            "barcode rank vs log10(reads+1) after dropping reads <= 10 and "
            "excluding the top 5%% and bottom 10%% of ranks, or a positive "
            "number for a manual cutoff."
        ),
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


def parse_reads_threshold_arg(raw: str | float | int) -> tuple[str, float | None]:
    """Return (mode, manual_threshold). mode is 'auto' or 'manual'."""
    text = str(raw).strip()
    if not text:
        raise ValueError("reads threshold must be 'auto' or a positive number")
    if text.lower() == "auto":
        return "auto", None
    try:
        value = float(text)
    except ValueError as exc:
        raise ValueError(
            "reads threshold must be 'auto' or a positive number"
        ) from exc
    if not math.isfinite(value) or value <= 0:
        raise ValueError("reads threshold must be > 0")
    return "manual", value


def resolve_reads_threshold(
    mode: str,
    manual_threshold: float | None,
    read_counts: dict[str, float],
) -> float:
    if mode == "auto":
        return auto_reads_threshold(list(read_counts.values())).threshold
    if manual_threshold is None:
        raise ValueError("manual reads threshold is missing")
    return float(manual_threshold)


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
            spot = (row.get("barcode") or "").strip().upper()
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
        "reads_threshold",
        "reads_threshold_mode",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


def write_reads_knee_plot(
    path: Path,
    sample_id: str,
    read_counts: dict[str, float],
    reads_threshold: float,
    mode: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    sorted_reads = sorted(read_counts.values(), reverse=True)
    ranks = list(range(1, len(sorted_reads) + 1))
    log_reads = [math.log10(value + 1.0) for value in sorted_reads]

    fig, ax = plt.subplots(figsize=_FIGSIZE)
    if ranks:
        ax.plot(ranks, log_reads, color="#333333", linewidth=0.8)
        ax.axhline(
            math.log10(reads_threshold + 1.0),
            color="#c0392b",
            linestyle="--",
            linewidth=0.7,
            label=f"threshold={reads_threshold:g} ({mode})",
        )
        # Mark the leftmost rank whose reads equal the threshold when present.
        mark_rank = None
        for rank, value in zip(ranks, sorted_reads):
            if value <= reads_threshold:
                mark_rank = rank
                break
        if mark_rank is not None:
            ax.scatter(
                [mark_rank],
                [math.log10(sorted_reads[mark_rank - 1] + 1.0)],
                color="#c0392b",
                s=12,
                zorder=4,
            )
        ax.legend(loc="upper right", fontsize=5)
    else:
        ax.text(0.5, 0.5, "No per-spot reads", ha="center", va="center", fontsize=6)
        ax.axis("off")
    ax.set_title(f"Reads knee ({sample_id})", fontsize=6)
    ax.set_xlabel("Barcode rank", fontsize=6)
    ax.set_ylabel("log10(reads + 1)", fontsize=6)
    tn.style_axes(ax)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def write_empty_plot(path: Path, sample_id: str, message: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.axis("off")
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=6)
    ax.set_title(f"Saturation analysis ({sample_id})\nSaturation rate: NA", fontsize=6)
    tn.style_axes(ax)
    fig.tight_layout()
    fig.savefig(path)
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
    y_ref_candidates = list(median_uniques)
    if predicted is not None and math.isfinite(predicted):
        y_ref_candidates.append(predicted)
    if model != "linear" and a > 0:
        y_ref_candidates.append(a)
    y_ref = max(y_ref_candidates) if y_ref_candidates else 0.0
    # Scale plotted y by 10^6 when values exceed 10^4, otherwise 10^3.
    if y_ref > 10_000:
        scale = 1e6
        scale_exponent = 6
    else:
        scale = 1e3
        scale_exponent = 3

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

    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.errorbar(
        x_values,
        medians_scaled,
        yerr=yerr,
        fmt="o-",
        color="blue",
        linewidth=0.8,
        markersize=2.5,
        capsize=1.5,
        elinewidth=0.5,
        label="Observed median (IQR)",
    )
    ax.plot(fit_x, fit_y, "r--", linewidth=0.8, label=fit_label)
    ax.scatter(
        [pred_x],
        [pred_y],
        color="green",
        s=12,
        zorder=4,
        label=f"Prediction at {pred_fraction:g}x",
    )
    if model != "linear" and a > 0:
        ax.axhline(
            a / scale,
            color="purple",
            linestyle="--",
            linewidth=0.7,
            label="Max CpG Reads",
        )
    if saturation_rate is not None:
        sat_text = f"{saturation_rate:.2f}%"
    elif model == "linear":
        sat_text = "NA (linear, unsaturated)"
    else:
        sat_text = "NA"
    ax.set_title(
        f"Saturation analysis ({sample_id})\nSaturation rate: {sat_text}",
        fontsize=6,
    )
    ax.set_xlabel(x_label, fontsize=6)
    ax.set_ylabel(
        rf"Median Unique CpGs per Spot ($\times 10^{{{scale_exponent}}}$)",
        fontsize=6,
    )
    ax.legend(loc="lower right", fontsize=5)
    tn.style_axes(ax)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    threshold_mode, manual_threshold = parse_reads_threshold_arg(args.reads_threshold)
    work_path = Path(args.work_path)
    sample_id = work_path.name
    read_counts_path = work_path / "coverage" / "per_spot_read_counts.tsv"
    host_cov_paths = sorted((work_path / "coverage" / "host").rglob("*.CG.cov"))
    output_dir = Path(args.output_dir) if args.output_dir else work_path / "qc" / "saturation"
    plot_path = output_dir / "saturation_curve.png"
    knee_plot_path = output_dir / "reads_knee.png"
    summary_path = output_dir / "saturation_summary.tsv"

    fastp_json_path = resolve_fastp_json(work_path, args.fastp_json)
    sequencing_gbp: float | None = None
    if fastp_json_path is not None:
        sequencing_gbp = load_fastp_sequencing_gbp(fastp_json_path)
    else:
        print("[saturation] warning=fastp_json_not_found")

    read_counts_preview = (
        parse_per_spot_reads(read_counts_path) if read_counts_path.exists() else {}
    )
    resolved_threshold: float | None = None
    knee_rank: int | None = None
    exclude_start: int | None = None
    exclude_end: int | None = None
    if threshold_mode == "manual":
        resolved_threshold = manual_threshold
    elif read_counts_preview:
        if threshold_mode == "auto":
            auto_knee = auto_reads_threshold(list(read_counts_preview.values()))
            resolved_threshold = auto_knee.threshold
            knee_rank = auto_knee.knee_rank
            exclude_start = auto_knee.exclude_start
            exclude_end = auto_knee.exclude_end
        else:
            resolved_threshold = resolve_reads_threshold(
                threshold_mode, manual_threshold, read_counts_preview
            )

    print(f"[saturation] work_path={work_path}")
    print(f"[saturation] reads_table={read_counts_path}")
    print(f"[saturation] host_cov_count={len(host_cov_paths)}")
    print(f"[saturation] reads_threshold_mode={threshold_mode}")
    print(
        "[saturation] reads_threshold="
        + (f"{resolved_threshold:g}" if resolved_threshold is not None else "unresolved")
    )
    if knee_rank is not None:
        print(f"[saturation] knee_rank={knee_rank}")
    if exclude_start is not None:
        print(f"[saturation] exclude_start={exclude_start}")
    if exclude_end is not None:
        print(f"[saturation] exclude_end={exclude_end}")
    print(f"[saturation] linear_r2_threshold={args.linear_r2_threshold}")
    print(f"[saturation] fastp_json={fastp_json_path}")
    print(
        "[saturation] sequencing_gbp="
        + (f"{sequencing_gbp:.6f}" if sequencing_gbp is not None else "NA")
    )
    print(f"[saturation] output_plot={plot_path}")
    print(f"[saturation] output_knee_plot={knee_plot_path}")
    print(f"[saturation] output_summary={summary_path}")

    if args.dry_run:
        print("[saturation] dry_run=1")
        return 0

    read_counts = parse_per_spot_reads(read_counts_path)
    if threshold_mode == "auto":
        reads_threshold = auto_reads_threshold(list(read_counts.values())).threshold
    else:
        reads_threshold = resolve_reads_threshold(
            threshold_mode, manual_threshold, read_counts
        )
    write_reads_knee_plot(
        knee_plot_path,
        sample_id,
        read_counts,
        reads_threshold,
        threshold_mode,
    )
    spot_cov_map: dict[str, Path] = {
        cov_path.name[: -len(".CG.cov")]: cov_path for cov_path in host_cov_paths
    }
    hq_spots = [
        spot
        for spot, reads in read_counts.items()
        if reads > reads_threshold and spot in spot_cov_map
    ]
    hq_spots.sort()

    def base_summary_fields() -> dict[str, str]:
        return {
            "reads_threshold": f"{reads_threshold:g}",
            "reads_threshold_mode": threshold_mode,
        }

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
            **base_summary_fields(),
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
            **base_summary_fields(),
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
        **base_summary_fields(),
    }
    write_summary_tsv(summary_path, row)
    print(f"[saturation] hq_spot_count={len(hq_spots)}")
    print("[saturation] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
