#!/usr/bin/env python3
"""Flatten host per-spot .CG.cov rows into two sorted headerless TSV tables."""

from __future__ import annotations

import argparse
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
    return parser.parse_args()


def spot_id_from_cov_path(cov_path: Path) -> str:
    name = cov_path.name
    if not name.endswith(SUFFIX):
        raise ValueError(f"expected *{SUFFIX} file, got: {cov_path}")
    return name[: -len(SUFFIX)]


def collect_rows(cov_paths: list[Path]) -> list[tuple[str, str, int, int, int, int]]:
    rows: list[tuple[str, str, int, int, int, int]] = []
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
                rows.append((spot_id, chrom, start, end, mc, c_unmeth))
    return rows


def write_tsv(path: Path, ordered: list[tuple[str, str, int, int, int, int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for spot_id, chrom, start, end, mc, c_unmeth in ordered:
            handle.write(
                f"{spot_id}\t{chrom}\t{start}\t{end}\t{mc}\t{c_unmeth}\n"
            )


def main() -> int:
    args = parse_args()
    work_path = Path(args.work_path)
    host_cov_root = work_path / "coverage" / "host"
    cov_paths = sorted(host_cov_root.rglob(f"*{SUFFIX}"))
    out_by_id = work_path / "coverage" / OUT_BY_ID
    out_by_pos = work_path / "coverage" / OUT_BY_POS

    print(f"[aggregate] work_path={work_path}")
    print(f"[aggregate] host_cov_glob={host_cov_root}/**/*{SUFFIX}")
    print(f"[aggregate] cov_file_count={len(cov_paths)}")
    print(f"[aggregate] out_by_id={out_by_id}")
    print(f"[aggregate] out_by_pos={out_by_pos}")

    rows = collect_rows(cov_paths)
    print(f"[aggregate] row_count={len(rows)}")

    by_id = sorted(rows, key=lambda r: (r[0], r[1], r[2], r[3]))
    by_pos = sorted(rows, key=lambda r: (r[1], r[2], r[3], r[0]))

    if args.dry_run:
        print("[aggregate] dry_run=1 skip_write")
        return 0

    write_tsv(out_by_id, by_id)
    write_tsv(out_by_pos, by_pos)
    print("[aggregate] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
