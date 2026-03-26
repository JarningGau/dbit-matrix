#!/usr/bin/env python3
"""Run methscan prepare on host per-spot *.CG.cov files via envs/methscan pixi workspace."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def repo_root_from_script() -> Path:
    return Path(__file__).resolve().parent.parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Collect coverage/host/**/*.CG.cov, run `pixi run methscan prepare` "
            "in the methscan pixi workspace (default: envs/methscan), writing to "
            "coverage/host_prepare/."
        )
    )
    parser.add_argument(
        "--work-path",
        required=True,
        help="Sample work directory (contains coverage/host/).",
    )
    parser.add_argument(
        "--pixi-manifest",
        help=(
            "Directory containing pixi.toml for methscan. "
            "Default: <repo>/envs/methscan."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print resolved paths and command; do not run pixi.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    work = Path(args.work_path).resolve()
    cov_dir = work / "coverage" / "host"
    out_dir = work / "coverage" / "host_prepare"
    if args.pixi_manifest:
        manifest = Path(args.pixi_manifest).resolve()
    else:
        manifest = repo_root_from_script() / "envs" / "methscan"

    pixi_toml = manifest / "pixi.toml"
    if not manifest.is_dir() or not pixi_toml.is_file():
        print(
            f"error: methscan pixi workspace not found or missing pixi.toml: {manifest}",
            file=sys.stderr,
        )
        return 1

    cov_files = sorted(p for p in cov_dir.rglob("*.CG.cov") if p.is_file())
    if not cov_files:
        print(
            f"error: no *.CG.cov under {cov_dir}",
            file=sys.stderr,
        )
        return 1

    cmd: list[str] = [
        "pixi",
        "run",
        "methscan",
        "prepare",
        *[str(p) for p in cov_files],
        str(out_dir),
        "--input-format",
        "bismark",
    ]
    if args.dry_run:
        print(f"cwd={manifest}")
        print("command=" + " ".join(cmd))
        print(f"input_count={len(cov_files)}")
        print(f"output_dir={out_dir}")
        return 0

    out_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(cmd, cwd=str(manifest), check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
