#!/usr/bin/env python3
"""Sliding 10-CpG span: scan +strand CG, span = p10 - p1 + 1, histograms by chr.

Requires a reference FASTA indexed with samtools faidx (same directory as .fa).

If ``<out-dir>/spans.tsv.gz`` exists (e.g. from a prior ``--write-spans-tsv`` run), the
reference is not rescanned: spans are loaded from that file.
``histogram_counts.tsv`` and ``cg_span_facet_by_chr.png`` are always refreshed together
(same IQR filtering and bin edges). ``summary_by_chr.tsv`` and ``spans.tsv.gz`` are not
overwritten in that mode. Use ``--recompute`` for a full scan and all outputs.

Example:

  pixi run python scripts/cg_span_histogram.py \\
    --reference /path/to/mm10.fa \\
    --out-dir .idea/cg_density_mm10
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pysam

_UCSC_AUTOSOME = re.compile(r"^chr(\d+)$")
_ENSEMBL_AUTOSOME = re.compile(r"^(\d+)$")


def parse_chromosome_csv(chromosome_csv: str) -> list[str]:
    chromosomes = [c.strip() for c in chromosome_csv.split(",") if c.strip()]
    if not chromosomes:
        raise ValueError("chromosome list is empty")
    return chromosomes


def validate_chromosomes_in_reference(reference_file: str, chromosomes: list[str]) -> None:
    with pysam.FastaFile(reference_file) as fasta:
        available = set(fasta.references)
    missing = [c for c in chromosomes if c not in available]
    if missing:
        raise ValueError(
            f"chromosome(s) not found in reference {reference_file!r}: {','.join(missing)}"
        )


def autosomes_and_x_from_reference(reference_path: str) -> list[str]:
    """Autosomes (numeric) + X; excludes Y, M/MT, unlocalized, alt contigs.

    UCSC (e.g. mm10, hg38): ``chr1``..``chrN`` then ``chrX`` if present.
    Ensembl primary: ``1``..``N`` then ``X`` if present.
    """
    with pysam.FastaFile(reference_path) as fasta:
        names = set(fasta.references)
    ucsc_autos = [n for n in names if _UCSC_AUTOSOME.match(n)]
    if ucsc_autos:
        pairs = sorted(
            (int(_UCSC_AUTOSOME.match(n).group(1)), n) for n in ucsc_autos
        )
        ordered = [n for _, n in pairs]
        if "chrX" in names:
            ordered.append("chrX")
        return ordered
    ens_autos = [n for n in names if _ENSEMBL_AUTOSOME.match(n)]
    pairs = sorted((int(_ENSEMBL_AUTOSOME.match(n).group(1)), n) for n in ens_autos)
    ordered = [n for _, n in pairs]
    if "X" in names:
        ordered.append("X")
    return ordered


def resolve_target_chromosomes(reference_path: str, chromosomes_csv: str | None) -> list[str]:
    if chromosomes_csv:
        chrs = parse_chromosome_csv(chromosomes_csv)
        validate_chromosomes_in_reference(reference_path, chrs)
        return chrs
    chrs = autosomes_and_x_from_reference(reference_path)
    if not chrs:
        raise ValueError(
            "Could not infer autosomes + X from reference FASTA index "
            "(expected UCSC chr1..chrN + chrX, or Ensembl 1..N + X). "
            "Pass --chromosomes explicitly."
        )
    return chrs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Scan reference for forward-strand CpGs on autosomes + X (default: inferred "
            "from the FASTA index, e.g. mm10 chr1–19+chrX, hg38 chr1–22+chrX); "
            "for each sliding window of 10 consecutive CpGs (step 1), "
            "compute span = p10 - p1 + 1 (0-based C positions). "
            "Write summary TSVs, binned counts, and a faceted histogram PNG."
        )
    )
    parser.add_argument(
        "-f",
        "--reference",
        required=True,
        help="Path to reference FASTA (must be indexed: samtools faidx <ref.fa>).",
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        type=Path,
        default=Path(".idea/cg_density_mm10"),
        help="Output directory (default: .idea/cg_density_mm10).",
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=80,
        help="Number of histogram bins (shared edges across chromosomes).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print per-chromosome CpG/window counts and exit without writing files.",
    )
    parser.add_argument(
        "--write-spans-tsv",
        action="store_true",
        help=(
            "Also write spans.tsv.gz (chr, span per window). "
            "Can be very large for a whole genome; use sparingly. "
            "Ignored if spans.tsv.gz already exists (use --recompute to replace)."
        ),
    )
    parser.add_argument(
        "--recompute",
        action="store_true",
        help=(
            "Ignore spans.tsv.gz in --out-dir and rescan the reference "
            "(default: if spans.tsv.gz exists, load it and skip scanning)."
        ),
    )
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable stderr progress output.",
    )
    parser.add_argument(
        "--chromosomes",
        metavar="LIST",
        help=(
            "Comma-separated chromosome names (order preserved). "
            "Default: autosomes + X from the reference (UCSC chrN + chrX, else Ensembl N + X)."
        ),
    )
    return parser.parse_args()


def _stderr_progress_line(
    label: str,
    current: int,
    total: int,
    detail: str,
    *,
    bar_width: int = 28,
) -> None:
    """Single-line progress on stderr (TTY: overwrite with \\r; else one line per update)."""
    if total <= 0:
        return
    frac = min(1.0, current / total)
    filled = min(bar_width, int(math.floor(frac * bar_width + 1e-9)))
    bar = "#" * filled + "-" * (bar_width - filled)
    tail = f"{current}/{total} ({100.0 * frac:.0f}%) {detail}"
    line = f"{label} [{bar}] {tail}"
    if hasattr(sys.stderr, "isatty") and sys.stderr.isatty():
        print(f"\r{line}", end="", file=sys.stderr, flush=True)
        if current >= total:
            print(file=sys.stderr)
    else:
        print(line, file=sys.stderr, flush=True)


def find_cpg_c_positions(reference_file: str, chromosome: str) -> list[int]:
    """0-based C coordinate for each forward CpG on one chromosome."""
    with pysam.FastaFile(reference_file) as fasta:
        sequence = fasta.fetch(chromosome).upper()
    positions: list[int] = []
    for start in range(len(sequence) - 1):
        if sequence[start] == "C" and sequence[start + 1] == "G":
            positions.append(start)
    return positions


def spans_for_chromosome(positions: np.ndarray) -> np.ndarray:
    """span[j] = pos[j+9] - pos[j] + 1."""
    if len(positions) < 10:
        return np.array([], dtype=np.int64)
    return positions[9:].astype(np.int64) - positions[:-9].astype(np.int64) + 1


def load_spans_from_tsv_gz(path: Path, chrom_order: list[str]) -> dict[str, np.ndarray]:
    """Load spans.tsv.gz into chr -> int64 array (keys follow chrom_order)."""
    buckets: dict[str, list[int]] = {c: [] for c in chrom_order}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as f:
        r = csv.reader(f, delimiter="\t")
        header = next(r, None)
        for row in r:
            if len(row) < 2:
                continue
            chrom, span_s = row[0], row[1]
            if chrom not in buckets:
                continue
            buckets[chrom].append(int(span_s))
    return {c: np.asarray(buckets[c], dtype=np.int64) for c in chrom_order}


def remove_outliers_iqr(x: np.ndarray) -> np.ndarray:
    """Tukey IQR rule; empty or tiny arrays unchanged."""
    if x.size == 0:
        return x
    if x.size < 4:
        return x
    q1, q3 = np.percentile(x, [25.0, 75.0])
    iqr = q3 - q1
    if iqr <= 0:
        return x
    lo = q1 - 1.5 * iqr
    hi = q3 + 1.5 * iqr
    return x[(x >= lo) & (x <= hi)]


def main() -> int:
    args = parse_args()
    ref = args.reference
    try:
        target_chromosomes = resolve_target_chromosomes(ref, args.chromosomes)
    except ValueError as e:
        print(e, file=sys.stderr)
        return 1
    out_dir: Path = args.out_dir
    spans_path = out_dir / "spans.tsv.gz"
    use_cache = spans_path.is_file() and not args.recompute

    spans_by_chr: dict[str, np.ndarray] = {}
    n_cpg_by_chr: dict[str, int | str] = {}
    available: set[str]

    if use_cache:
        if not args.no_progress:
            print(
                f"Using cached spans: {spans_path} "
                f"(skip scan; refresh histogram_counts.tsv + PNG; no summary/spans.gz)",
                file=sys.stderr,
            )
        spans_by_chr = load_spans_from_tsv_gz(spans_path, target_chromosomes)
        for chrom in target_chromosomes:
            n_cpg_by_chr[chrom] = ""
        available = set(target_chromosomes)
    else:
        with pysam.FastaFile(ref) as fasta:
            available = set(fasta.references)

        n_chr_target = len(target_chromosomes)
        for i, chrom in enumerate(target_chromosomes, start=1):
            if not args.no_progress:
                _stderr_progress_line(
                    "Scan reference (CpG)",
                    i,
                    n_chr_target,
                    chrom,
                )
            if chrom not in available:
                warnings.warn(f"Chromosome {chrom!r} not in reference; skipping.", stacklevel=1)
                n_cpg_by_chr[chrom] = 0
                continue
            pos_list = find_cpg_c_positions(ref, chrom)
            n_cpg_by_chr[chrom] = len(pos_list)
            pos = np.asarray(pos_list, dtype=np.int64)
            spans = spans_for_chromosome(pos)
            if len(spans) > 0:
                spans_by_chr[chrom] = spans

    if args.dry_run:
        print("dry-run: per-chromosome CpG counts and window counts", file=sys.stderr)
        for chrom in target_chromosomes:
            n_cpg = n_cpg_by_chr.get(chrom, 0)
            n_win = len(spans_by_chr.get(chrom, np.array([])))
            print(f"  {chrom}\tcpg={n_cpg}\twindows={n_win}", file=sys.stderr)
        return 0

    out_dir.mkdir(parents=True, exist_ok=True)

    if not spans_by_chr:
        print("No chromosome had >= 10 CpGs; nothing to write.", file=sys.stderr)
        return 1

    spans_plot_by_chr: dict[str, np.ndarray] = {}
    for chrom in target_chromosomes:
        raw = spans_by_chr.get(chrom)
        if raw is None or len(raw) == 0:
            spans_plot_by_chr[chrom] = np.array([], dtype=np.int64)
        else:
            spans_plot_by_chr[chrom] = remove_outliers_iqr(raw)

    # Histogram / PNG use IQR-filtered spans; bin edges from filtered data (fallback to raw if empty).
    valid_for_edges = [s for s in spans_plot_by_chr.values() if len(s) > 0]
    if not valid_for_edges:
        valid_for_edges = [s for s in spans_by_chr.values() if len(s) > 0]
    global_min = min(int(s.min()) for s in valid_for_edges)
    global_max = max(int(s.max()) for s in valid_for_edges)
    n_bins = max(1, args.bins)
    if global_min == global_max:
        edges = np.array([float(global_min), float(global_min) + 1.0], dtype=np.float64)
        n_bins = 1
    else:
        edges = np.linspace(float(global_min), float(global_max), n_bins + 1, dtype=np.float64)

    summary_path = out_dir / "summary_by_chr.tsv"
    hist_path = out_dir / "histogram_counts.tsv"

    if not use_cache:
        # summary_by_chr.tsv
        with summary_path.open("w", encoding="utf-8", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(
                [
                    "chr",
                    "n_cpg",
                    "n_windows",
                    "span_min",
                    "span_mean",
                    "span_median",
                    "span_max",
                ]
            )
            for chrom in target_chromosomes:
                if chrom not in available:
                    w.writerow([chrom, "", "", "", "", "", ""])
                    continue
                n_cpg = n_cpg_by_chr.get(chrom, 0)
                spans = spans_by_chr.get(chrom)
                if spans is None or len(spans) == 0:
                    w.writerow([chrom, n_cpg, 0, "", "", "", ""])
                    continue
                w.writerow(
                    [
                        chrom,
                        n_cpg,
                        len(spans),
                        int(spans.min()),
                        float(np.mean(spans)),
                        float(np.median(spans)),
                        int(spans.max()),
                    ]
                )

    # histogram_counts.tsv: same bins / IQR filter as the PNG (refreshed together)
    counts_by_chr: dict[str, np.ndarray] = {}
    for chrom in target_chromosomes:
        spans = spans_plot_by_chr.get(chrom)
        if spans is None or len(spans) == 0:
            counts_by_chr[chrom] = np.zeros(len(edges) - 1, dtype=np.int64)
        else:
            c, _ = np.histogram(spans.astype(np.float64), bins=edges)
            counts_by_chr[chrom] = c.astype(np.int64)

    with hist_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        header = ["bin_index", "bin_low", "bin_high"] + list(target_chromosomes)
        w.writerow(header)
        for b in range(len(edges) - 1):
            row = [
                b + 1,
                f"{edges[b]:.12g}",
                f"{edges[b + 1]:.12g}",
            ] + [int(counts_by_chr[c][b]) for c in target_chromosomes]
            w.writerow(row)

    # Faceted histogram PNG
    n_panel = len(target_chromosomes)
    fig, axes = plt.subplots(
        n_panel,
        1,
        figsize=(9, 0.55 * n_panel + 1),
        sharex=True,
        layout="constrained",
    )
    if n_panel == 1:
        axes = np.array([axes])
    for ax, chrom in zip(axes, target_chromosomes, strict=True):
        raw = spans_by_chr.get(chrom)
        spans = spans_plot_by_chr.get(chrom)
        if raw is None or len(raw) == 0:
            ax.text(0.5, 0.5, "no windows", ha="center", va="center", transform=ax.transAxes)
        elif spans is None or len(spans) == 0:
            ax.text(
                0.5,
                0.5,
                "all spans filtered (outliers)",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=7,
            )
        else:
            ax.hist(
                spans,
                bins=edges,
                density=True,
                color="steelblue",
                edgecolor="black",
                linewidth=0.2,
            )
        ax.set_ylabel(chrom, rotation=0, ha="right", va="center", fontsize=8)
        ax.tick_params(axis="y", labelsize=6)
    axes[-1].set_xlabel("Span (bp): p10 - p1 + 1 (0-based C)")
    fig.suptitle(
        "10-CpG span distribution (forward-strand CpGs, autosomes + X), facet by chromosome "
        "(Tukey IQR outliers removed for histogram)"
    )
    png_path = out_dir / "cg_span_facet_by_chr.png"
    fig.savefig(png_path, dpi=150)
    plt.close(fig)

    if args.write_spans_tsv and not use_cache:
        gz_path = out_dir / "spans.tsv.gz"
        total_spans = sum(
            len(spans_by_chr[c]) for c in target_chromosomes if c in spans_by_chr
        )
        step = max(1, total_spans // 100)
        written = 0
        with gzip.open(gz_path, "wt", encoding="utf-8", newline="") as gz:
            w = csv.writer(gz, delimiter="\t")
            w.writerow(["chr", "span"])
            for chrom in target_chromosomes:
                spans = spans_by_chr.get(chrom)
                if spans is None:
                    continue
                for s in spans:
                    w.writerow([chrom, int(s)])
                    written += 1
                    if (
                        not args.no_progress
                        and total_spans > 0
                        and (written == total_spans or written % step == 0)
                    ):
                        _stderr_progress_line(
                            "Write spans.tsv.gz",
                            written,
                            total_spans,
                            chrom,
                        )

    if not use_cache:
        print(f"Wrote: {summary_path}", file=sys.stderr)
    print(f"Wrote: {hist_path}", file=sys.stderr)
    print(f"Wrote: {png_path}", file=sys.stderr)
    if args.write_spans_tsv and not use_cache:
        print(f"Wrote: {out_dir / 'spans.tsv.gz'}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
