#!/usr/bin/env python3
"""Compare new EMSeq per-spot CG.cov files to a biscuit baseline tree."""

from __future__ import annotations

import argparse
import math
from pathlib import Path


def parse_cov(path: Path, pos_shift: int) -> dict[tuple[str, int], tuple[int, int, float]]:
    rows: dict[tuple[str, int], tuple[int, int, float]] = {}
    if not path.is_file():
        return rows
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if len(fields) < 6:
            continue
        chrom = fields[0]
        pos = int(fields[1]) + pos_shift
        pct = float(fields[3])
        meth = int(fields[4])
        unmeth = int(fields[5])
        rows[(chrom, pos)] = (meth, unmeth, pct)
    return rows


def pearson(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    if n < 2:
        return float("nan")
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    num = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    den_x = math.sqrt(sum((x - mean_x) ** 2 for x in xs))
    den_y = math.sqrt(sum((y - mean_y) ** 2 for y in ys))
    if den_x == 0 or den_y == 0:
        return float("nan")
    return num / (den_x * den_y)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-host", required=True, help="Baseline coverage/host dir")
    parser.add_argument("--new-host", required=True, help="New coverage/host dir")
    parser.add_argument(
        "--biscuit-pos-shift",
        type=int,
        default=-1,
        help="Added to biscuit pos to reach 0-based. Default: -1.",
    )
    args = parser.parse_args()
    baseline_root = Path(args.baseline_host)
    new_root = Path(args.new_host)
    spots = sorted({p.stem[: -len(".CG")] for p in baseline_root.rglob("*.CG.cov")})
    baseline_meth: list[float] = []
    new_meth: list[float] = []
    baseline_unmeth: list[float] = []
    new_unmeth: list[float] = []
    shared_abs_pct: list[float] = []
    disagree = 0
    shared = 0
    only_base = 0
    only_new = 0
    covered_shared = 0
    for spot in spots:
        base_matches = list(baseline_root.rglob(f"{spot}.CG.cov"))
        new_matches = list(new_root.rglob(f"{spot}.CG.cov"))
        base = parse_cov(base_matches[0], args.biscuit_pos_shift) if base_matches else {}
        new = parse_cov(new_matches[0], 0) if new_matches else {}
        only_base += len(set(base) - set(new))
        only_new += len(set(new) - set(base))
        for key in set(base) & set(new):
            shared += 1
            bm, bu, bp = base[key]
            nm, nu, np_ = new[key]
            baseline_meth.append(float(bm))
            new_meth.append(float(nm))
            baseline_unmeth.append(float(bu))
            new_unmeth.append(float(nu))
            shared_abs_pct.append(abs(bp - np_))
            if (bm + bu) > 0 and (nm + nu) > 0:
                covered_shared += 1
                if (bm > bu) != (nm > nu):
                    disagree += 1
    frac_disagree = (disagree / covered_shared) if covered_shared else float("nan")
    print(f"[compare] spots={len(spots)}")
    print(f"[compare] shared_sites={shared}")
    print(f"[compare] only_baseline={only_base}")
    print(f"[compare] only_new={only_new}")
    print(f"[compare] covered_shared={covered_shared}")
    print(f"[compare] pearson_meth={pearson(baseline_meth, new_meth)}")
    print(f"[compare] pearson_unmeth={pearson(baseline_unmeth, new_unmeth)}")
    print(
        f"[compare] mean_abs_pct_delta="
        f"{sum(shared_abs_pct) / len(shared_abs_pct) if shared_abs_pct else float('nan')}"
    )
    print(f"[compare] disagree_majority_fraction={frac_disagree}")
    if not covered_shared:
        print("[compare] warning: no covered shared sites")
        return 1
    if frac_disagree >= 0.5:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
