#!/usr/bin/env python3
"""Flatten host per-spot .CG.cov rows into two sorted headerless TSV tables.

This script is designed to keep Python memory usage roughly constant by:
- streaming `aggregated_cg_by_id.tsv` directly to disk (no full in-memory rows list)
- delegating coordinate-based ordering of `aggregated_cg_by_pos.tsv` to GNU `sort`
  (with a configurable memory cap via `-S`).
"""

from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path

SUFFIX = ".CG.cov"
OUT_BY_ID = "aggregated_cg_by_id.tsv"
OUT_BY_POS = "aggregated_cg_by_pos.tsv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read coverage/host/**/*.CG.cov (Bismark-like, from methy_caller), "
            "emit two headerless TSV files under coverage/."
        )
    )
    parser.add_argument(
        "--work-path",
        required=True,
        help="Sample work directory containing coverage/host/.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print resolved paths and row counts; do not write output files.",
    )
    parser.add_argument(
        "--sort-mem",
        default="8G",
        help=(
            "Memory limit passed to GNU `sort -S` when generating "
            "`aggregated_cg_by_pos.tsv`. Default: 8G."
        ),
    )
    return parser.parse_args()


def spot_id_from_cov_path(cov_path: Path) -> str:
    name = cov_path.name
    if not name.endswith(SUFFIX):
        raise ValueError(f"expected *{SUFFIX} file, got: {cov_path}")
    return name[: -len(SUFFIX)]


def iter_parsed_rows(cov_paths: list[Path]):
    """Yield parsed (spot_id, chrom, start, end, mc, c_unmeth) rows.

    Note: to keep memory usage constant, callers should consume this iterator
    in a streaming fashion (no `list()` / `sorted()` of all rows).
    """
    for cov_path in cov_paths:
        spot_id = spot_id_from_cov_path(cov_path)
        with cov_path.open("r", encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                fields = line.split("\t")
                if len(fields) < 6:
                    continue
                chrom = fields[0]
                try:
                    start = int(fields[1])
                    end = int(fields[2])
                    mc = int(fields[4])
                    c_unmeth = int(fields[5])
                except ValueError:
                    continue
                yield (spot_id, chrom, start, end, mc, c_unmeth)


def write_by_id(out_by_id: Path, cov_paths: list[Path]) -> int:
    out_by_id.parent.mkdir(parents=True, exist_ok=True)
    row_count = 0
    with out_by_id.open("w", encoding="utf-8") as handle:
        for spot_id, chrom, start, end, mc, c_unmeth in iter_parsed_rows(cov_paths):
            handle.write(
                f"{spot_id}\t{chrom}\t{start}\t{end}\t{mc}\t{c_unmeth}\n"
            )
            row_count += 1
    return row_count


def run_sort_by_pos(out_by_pos: Path, out_by_id: Path, sort_mem: str) -> None:
    """External sort for `aggregated_cg_by_pos.tsv`.

    Key: (chr, start, end, id) where:
    - id   = col1
    - chr  = col2
    - start= col3 (numeric)
    - end  = col4 (numeric)
    """
    sort_env = dict(os.environ)
    # Make sort order deterministic across locales; Python string ordering is codepoint-based.
    sort_env["LC_ALL"] = "C"

    sort_cmd = [
        "sort",
        "-t",
        "\t",
        "-k2,2",
        "-k3,3n",
        "-k4,4n",
        "-k1,1",
        "-S",
        sort_mem,
        "-o",
        str(out_by_pos),
        str(out_by_id),
    ]
    print(f"[aggregate] run: {' '.join(sort_cmd)}")
    subprocess.run(sort_cmd, check=True, env=sort_env)


def main() -> int:
    args = parse_args()
    work_path = Path(args.work_path)
    host_cov_root = work_path / "coverage" / "host"
    out_by_id = work_path / "coverage" / OUT_BY_ID
    out_by_pos = work_path / "coverage" / OUT_BY_POS

    print(f"[aggregate] work_path={work_path}")
    print(f"[aggregate] host_cov_glob={host_cov_root}/**/*{SUFFIX}")
    print(f"[aggregate] out_by_id={out_by_id}")
    print(f"[aggregate] out_by_pos={out_by_pos}")

    by_id_exists = out_by_id.exists()
    if by_id_exists:
        print(f"[aggregate] out_by_id_exists=1 path={out_by_id}")
    else:
        print(f"[aggregate] out_by_id_exists=0 path={out_by_id}")

    row_count = 0
    if args.dry_run:
        if by_id_exists:
            print("[aggregate] dry_run=1 skip_read_cov; will use existing out_by_id")
            row_count = 0
        else:
            cov_paths = sorted(
                host_cov_root.rglob(f"*{SUFFIX}"),
                key=lambda p: spot_id_from_cov_path(p),
            )
            print(f"[aggregate] cov_file_count={len(cov_paths)}")
            # Count rows without materializing them in memory.
            for _ in iter_parsed_rows(cov_paths):
                row_count += 1
    else:
        if by_id_exists:
            print("[aggregate] skip_write_by_id (existing file)")
            row_count = 0
        else:
            cov_paths = sorted(
                host_cov_root.rglob(f"*{SUFFIX}"),
                key=lambda p: spot_id_from_cov_path(p),
            )
            print(f"[aggregate] cov_file_count={len(cov_paths)}")
            row_count = write_by_id(out_by_id, cov_paths)
    print(f"[aggregate] row_count={row_count}")

    if args.dry_run:
        print("[aggregate] dry_run=1 skip_write")
        # Show the planned external sort invocation (still no file writes).
        print(
            "[aggregate] dry_run=1 skip sort; planned: "
            f"sort -t $'\\t' -k2,2 -k3,3n -k4,4n -k1,1 -S {args.sort_mem} "
            f"-o {out_by_pos} {out_by_id}"
        )
        return 0

    run_sort_by_pos(out_by_pos, out_by_id, sort_mem=args.sort_mem)
    print("[aggregate] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
