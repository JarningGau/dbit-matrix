#!/usr/bin/env python3
"""Run methscan steps via envs/methscan pixi workspace.

Reads per-cell coverage from work/coverage/host/; writes all methscan outputs under
work/methscan/ (see methscan_dirs). Context CG keeps the legacy layout; other contexts
nest under work/methscan/<CONTEXT>/.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

METHSCAN_DEFAULT_CHUNKSIZE = 10_000_000
DEFAULT_FILTER_MIN_SITES = 50_000
DEFAULT_PROFILE_STRAND_COLUMN = 6
DEFAULT_SCAN_THREADS = 10
DEFAULT_MATRIX_THREADS = 10
METHSCAN_CONTEXTS = ("CG", "CA", "CC", "CT")
DEFAULT_CONTEXT = "CG"


def repo_root_from_script() -> Path:
    return Path(__file__).resolve().parent.parent


def default_methscan_manifest() -> Path:
    return repo_root_from_script() / "envs" / "methscan"


def resolve_manifest(pixi_manifest: str | None) -> Path:
    if pixi_manifest:
        return Path(pixi_manifest).resolve()
    return default_methscan_manifest()


def validate_manifest(manifest: Path) -> str | None:
    pixi_toml = manifest / "pixi.toml"
    if not manifest.is_dir() or not pixi_toml.is_file():
        return f"methscan pixi workspace not found or missing pixi.toml: {manifest}"
    return None


def normalize_context(context: str | None) -> str:
    value = (context or DEFAULT_CONTEXT).upper()
    if value not in METHSCAN_CONTEXTS:
        raise ValueError(
            f"methscan context must be one of {', '.join(METHSCAN_CONTEXTS)}; got {context!r}"
        )
    return value


def coverage_host_dir(work: Path) -> Path:
    return work / "coverage" / "host"


def cov_glob_suffix(context: str) -> str:
    return f"*.{context}.cov"


def methscan_dirs(work: Path, context: str = DEFAULT_CONTEXT) -> dict[str, Path]:
    context = normalize_context(context)
    root = work / "methscan" if context == "CG" else work / "methscan" / context
    return {
        "root": root,
        "compact": root / "compact",
        "filter": root / "filter",
        "matrix": root / "matrix",
    }


def prepared_subdir_path(work: Path, name: str, context: str = DEFAULT_CONTEXT) -> Path:
    m = methscan_dirs(work, context)
    if name == "compact":
        return m["compact"]
    if name == "filter":
        return m["filter"]
    raise ValueError(f"unknown prepared subdir: {name}")


def run_methscan(
    manifest: Path,
    methscan_args: list[str],
    *,
    dry_run: bool,
    label: str = "",
) -> int:
    cmd = ["pixi", "run", "methscan", *methscan_args]
    prefix = f"{label}: " if label else ""
    if dry_run:
        print(f"{prefix}cwd={manifest}")
        print(f"{prefix}command=" + " ".join(cmd))
        return 0
    subprocess.run(cmd, cwd=str(manifest), check=True)
    return 0


def cmd_prepare(args: argparse.Namespace, work: Path, manifest: Path) -> int:
    context = normalize_context(getattr(args, "context", DEFAULT_CONTEXT))
    cov_dir = coverage_host_dir(work)
    out_dir = methscan_dirs(work, context)["compact"]
    pattern = cov_glob_suffix(context)
    cov_files = sorted(p for p in cov_dir.rglob(pattern) if p.is_file())
    if not cov_files and not args.dry_run:
        print(f"error: no {pattern} under {cov_dir}", file=sys.stderr)
        return 1
    if args.chunksize < 1:
        print("error: --chunksize must be >= 1", file=sys.stderr)
        return 1
    methscan_args = [
        "prepare",
        *[str(p) for p in cov_files],
        str(out_dir),
        "--input-format",
        "bismark",
        "--chunksize",
        str(args.chunksize),
    ]
    if args.dry_run:
        print(f"cwd={manifest}")
        if cov_files:
            display_args = methscan_args
        else:
            display_args = [
                "prepare",
                f"<INPUT: coverage/host/**/{pattern}>",
                str(out_dir),
                "--input-format",
                "bismark",
                "--chunksize",
                str(args.chunksize),
            ]
            print(
                f"note: no {pattern} under {cov_dir} yet; "
                f"prepare needs per-cell .{context}.cov files before this command is valid.",
                file=sys.stderr,
            )
        print("command=" + " ".join(["pixi", "run", "methscan", *display_args]))
        print(f"context={context}")
        print(f"input_count={len(cov_files)}")
        print(f"output_dir={out_dir}")
        print(f"chunksize={args.chunksize}")
        return 0
    if not cov_files:
        print(f"error: no {pattern} under {cov_dir}", file=sys.stderr)
        return 1
    out_dir.mkdir(parents=True, exist_ok=True)
    return run_methscan(manifest, methscan_args, dry_run=False)


def cmd_filter(args: argparse.Namespace, work: Path, manifest: Path) -> int:
    context = normalize_context(getattr(args, "context", DEFAULT_CONTEXT))
    paths = methscan_dirs(work, context)
    data_dir = paths["compact"]
    filtered_dir = paths["filter"]
    methscan_args = [
        "filter",
        str(data_dir),
        str(filtered_dir),
        "--min-sites",
        str(args.min_sites),
    ]
    return run_methscan(manifest, methscan_args, dry_run=args.dry_run)


def cmd_profile(args: argparse.Namespace, work: Path, manifest: Path) -> int:
    context = normalize_context(getattr(args, "context", DEFAULT_CONTEXT))
    bed = Path(args.tss_bed).resolve()
    if not bed.is_file() and not args.dry_run:
        print(f"error: TSS bed not found: {bed}", file=sys.stderr)
        return 1
    data_dir = prepared_subdir_path(work, args.prepared_dir, context)
    out_csv = (
        Path(args.profile_csv).resolve()
        if args.profile_csv
        else methscan_dirs(work, context)["root"] / "TSS_profile.csv"
    )
    methscan_args = [
        "profile",
        str(bed),
        str(data_dir),
        str(out_csv),
        "--strand-column",
        str(args.strand_column),
    ]
    return run_methscan(manifest, methscan_args, dry_run=args.dry_run)


def cmd_smooth(args: argparse.Namespace, work: Path, manifest: Path) -> int:
    context = normalize_context(getattr(args, "context", DEFAULT_CONTEXT))
    data_dir = methscan_dirs(work, context)["filter"]
    methscan_args = ["smooth", str(data_dir)]
    if args.bandwidth is not None:
        methscan_args.extend(["--bandwidth", str(args.bandwidth)])
    if args.use_weights:
        methscan_args.append("--use-weights")
    return run_methscan(manifest, methscan_args, dry_run=args.dry_run)


def cmd_scan(args: argparse.Namespace, work: Path, manifest: Path) -> int:
    context = normalize_context(getattr(args, "context", DEFAULT_CONTEXT))
    data_dir = methscan_dirs(work, context)["filter"]
    out_bed = (
        Path(args.vmrs_bed).resolve()
        if args.vmrs_bed
        else methscan_dirs(work, context)["root"] / "VMRs.bed"
    )
    methscan_args = [
        "scan",
        str(data_dir),
        str(out_bed),
        "--threads",
        str(args.threads),
    ]
    return run_methscan(manifest, methscan_args, dry_run=args.dry_run)


def cmd_matrix(args: argparse.Namespace, work: Path, manifest: Path) -> int:
    context = normalize_context(getattr(args, "context", DEFAULT_CONTEXT))
    regions = (
        Path(args.vmrs_bed).resolve()
        if args.vmrs_bed
        else methscan_dirs(work, context)["root"] / "VMRs.bed"
    )
    if not regions.is_file() and not args.dry_run:
        print(f"error: VMRs bed not found: {regions}", file=sys.stderr)
        return 1
    data_dir = methscan_dirs(work, context)["filter"]
    out_dir = (
        Path(args.matrix_prefix).resolve()
        if args.matrix_prefix
        else methscan_dirs(work, context)["matrix"]
    )
    methscan_args = [
        "matrix",
        str(regions),
        str(data_dir),
        str(out_dir),
        "--threads",
        str(args.threads),
    ]
    if args.sparse:
        methscan_args.append("--sparse")
    return run_methscan(manifest, methscan_args, dry_run=args.dry_run)


def add_common(p: argparse.ArgumentParser) -> None:
    p.add_argument(
        "--work-path",
        required=True,
        help="Sample work directory (contains coverage/ for input; methscan/ for outputs).",
    )
    p.add_argument(
        "--pixi-manifest",
        help="Directory containing pixi.toml for methscan (default: <repo>/envs/methscan).",
    )
    p.add_argument(
        "--context",
        choices=METHSCAN_CONTEXTS,
        default=DEFAULT_CONTEXT,
        help=f"Host cov context to prepare/analyze (default: {DEFAULT_CONTEXT}). "
        "CG uses methscan/; other contexts use methscan/<CONTEXT>/.",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Print resolved command; do not run pixi.",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run methscan via pixi (envs/methscan) with repo coverage/methscan paths.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_prep = sub.add_parser(
        "prepare",
        help="methscan prepare from coverage/host/*.<CONTEXT>.cov",
    )
    add_common(p_prep)
    p_prep.add_argument(
        "--chunksize",
        type=int,
        default=METHSCAN_DEFAULT_CHUNKSIZE,
        metavar="N",
        help=f"bp per chromosome chunk (default: {METHSCAN_DEFAULT_CHUNKSIZE}).",
    )

    p_filt = sub.add_parser("filter", help="methscan filter compact -> filter")
    add_common(p_filt)
    p_filt.add_argument(
        "--min-sites",
        type=int,
        default=DEFAULT_FILTER_MIN_SITES,
        metavar="N",
        help=f"Minimum methylation sites per cell (default: {DEFAULT_FILTER_MIN_SITES}).",
    )

    p_prof = sub.add_parser("profile", help="methscan profile TSS / regions")
    add_common(p_prof)
    p_prof.add_argument(
        "--tss-bed",
        required=True,
        help="Path to sorted .bed file of regions.",
    )
    p_prof.add_argument(
        "--strand-column",
        type=int,
        default=DEFAULT_PROFILE_STRAND_COLUMN,
        help=f"BED column (1-based) for strand (default: {DEFAULT_PROFILE_STRAND_COLUMN}).",
    )
    p_prof.add_argument(
        "--prepared-dir",
        choices=("compact", "filter"),
        default="compact",
        help="Which prepared directory under methscan/ to use (default: compact).",
    )
    p_prof.add_argument(
        "--profile-csv",
        help="Output .csv path (default: methscan[/CONTEXT]/TSS_profile.csv).",
    )

    p_smooth = sub.add_parser("smooth", help="methscan smooth on filter dir")
    add_common(p_smooth)
    p_smooth.add_argument(
        "--bandwidth",
        type=int,
        default=None,
        metavar="N",
        help="Smoothing bandwidth in bp (methscan default if omitted).",
    )
    p_smooth.add_argument(
        "--use-weights",
        action="store_true",
        help="Weigh each site by log1p(coverage).",
    )

    p_scan = sub.add_parser("scan", help="methscan scan VMRs from filter dir")
    add_common(p_scan)
    p_scan.add_argument(
        "--threads",
        type=int,
        default=DEFAULT_SCAN_THREADS,
        help=f"Thread count (default: {DEFAULT_SCAN_THREADS}).",
    )
    p_scan.add_argument(
        "--vmrs-bed",
        help="Output .bed path (default: methscan[/CONTEXT]/VMRs.bed).",
    )

    p_matrix = sub.add_parser("matrix", help="methscan matrix from regions bed + filter dir")
    add_common(p_matrix)
    p_matrix.add_argument(
        "--threads",
        type=int,
        default=DEFAULT_MATRIX_THREADS,
        help=f"Thread count (default: {DEFAULT_MATRIX_THREADS}).",
    )
    p_matrix.add_argument(
        "--sparse",
        action="store_true",
        help="Write sparse matrix output (methscan --sparse).",
    )
    p_matrix.add_argument(
        "--vmrs-bed",
        help="Input regions .bed (default: methscan[/CONTEXT]/VMRs.bed from scan).",
    )
    p_matrix.add_argument(
        "--matrix-prefix",
        help="Output directory (default: methscan[/CONTEXT]/matrix).",
    )

    p_all = sub.add_parser(
        "all",
        help="Run prepare, filter, profile, smooth, scan, matrix in order (sequential).",
    )
    add_common(p_all)
    p_all.add_argument(
        "--chunksize",
        type=int,
        default=METHSCAN_DEFAULT_CHUNKSIZE,
        metavar="N",
    )
    p_all.add_argument(
        "--min-sites",
        type=int,
        default=DEFAULT_FILTER_MIN_SITES,
    )
    p_all.add_argument(
        "--tss-bed",
        required=True,
        help="Path to sorted .bed for profile step.",
    )
    p_all.add_argument(
        "--strand-column",
        type=int,
        default=DEFAULT_PROFILE_STRAND_COLUMN,
    )
    p_all.add_argument(
        "--prepared-dir",
        choices=("compact", "filter"),
        default="compact",
        help="Profile step: which prepared directory (default: compact).",
    )
    p_all.add_argument(
        "--profile-csv",
        help="Profile output .csv (default: methscan[/CONTEXT]/TSS_profile.csv).",
    )
    p_all.add_argument(
        "--bandwidth",
        type=int,
        default=None,
    )
    p_all.add_argument(
        "--use-weights",
        action="store_true",
    )
    p_all.add_argument(
        "--scan-threads",
        type=int,
        default=DEFAULT_SCAN_THREADS,
    )
    p_all.add_argument(
        "--matrix-threads",
        type=int,
        default=DEFAULT_MATRIX_THREADS,
    )
    p_all.add_argument(
        "--sparse",
        action="store_true",
        help="Pass --sparse to methscan matrix.",
    )
    p_all.add_argument(
        "--vmrs-bed",
        help="VMRs .bed path for scan output / matrix input (default: methscan[/CONTEXT]/VMRs.bed).",
    )
    p_all.add_argument(
        "--matrix-prefix",
        help="Matrix output directory (default: methscan[/CONTEXT]/matrix).",
    )

    return parser


def _ns_with_context(args: argparse.Namespace, **kwargs) -> argparse.Namespace:
    return argparse.Namespace(
        work_path=str(kwargs.pop("work")),
        pixi_manifest=getattr(args, "pixi_manifest", None),
        dry_run=args.dry_run,
        context=normalize_context(getattr(args, "context", DEFAULT_CONTEXT)),
        **kwargs,
    )


def dispatch(args: argparse.Namespace) -> int:
    work = Path(args.work_path).resolve()
    try:
        normalize_context(getattr(args, "context", DEFAULT_CONTEXT))
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    manifest = resolve_manifest(getattr(args, "pixi_manifest", None))
    err = validate_manifest(manifest)
    if err:
        print(f"error: {err}", file=sys.stderr)
        return 1

    handlers = {
        "prepare": cmd_prepare,
        "filter": cmd_filter,
        "profile": cmd_profile,
        "smooth": cmd_smooth,
        "scan": cmd_scan,
        "matrix": cmd_matrix,
    }
    if args.command == "all":
        return cmd_all(args, work, manifest)

    return handlers[args.command](args, work, manifest)


def cmd_all(args: argparse.Namespace, work: Path, manifest: Path) -> int:
    """Run all six steps with shared flags; fail on first error."""
    steps = [
        (
            "prepare",
            _ns_with_context(args, work=work, chunksize=args.chunksize),
            cmd_prepare,
        ),
        (
            "filter",
            _ns_with_context(args, work=work, min_sites=args.min_sites),
            cmd_filter,
        ),
        (
            "profile",
            _ns_with_context(
                args,
                work=work,
                tss_bed=args.tss_bed,
                strand_column=args.strand_column,
                prepared_dir=args.prepared_dir,
                profile_csv=args.profile_csv,
            ),
            cmd_profile,
        ),
        (
            "smooth",
            _ns_with_context(
                args,
                work=work,
                bandwidth=args.bandwidth,
                use_weights=args.use_weights,
            ),
            cmd_smooth,
        ),
        (
            "scan",
            _ns_with_context(
                args,
                work=work,
                threads=args.scan_threads,
                vmrs_bed=args.vmrs_bed,
            ),
            cmd_scan,
        ),
        (
            "matrix",
            _ns_with_context(
                args,
                work=work,
                threads=args.matrix_threads,
                sparse=args.sparse,
                vmrs_bed=args.vmrs_bed,
                matrix_prefix=args.matrix_prefix,
            ),
            cmd_matrix,
        ),
    ]
    for _name, ns, fn in steps:
        rc = fn(ns, work, manifest)
        if rc != 0:
            return rc
    return 0


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    return dispatch(args)


if __name__ == "__main__":
    raise SystemExit(main())
