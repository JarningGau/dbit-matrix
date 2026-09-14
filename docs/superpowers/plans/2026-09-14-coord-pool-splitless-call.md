# Coordinate-sorted pool, splitless mainline, EMSeq caller Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Coordinate-sort and index the host pooled BAM, drop `split` from every `--stage all`, have `call` write `coverage/per_spot_read_counts.tsv`, and add an EMSeq pysam caller with inverted polarity.

**Architecture:** Extract shared pileup/CB/counts/smoke into `scripts/methy_call_core.py`. TAPS `scripts/methy_caller.py` and new `scripts-emseq/methy_caller.py` are thin entrypoints that set polarity. `scripts/pool.py` writes `pooled/pooled.host.sorted.bam` plus index. All three `make_cmd` drivers lose `split`. Downstream stages read the new BAM name and counts path.

**Tech Stack:** Python 3.11, pysam, samtools, pixi, unittest (stdlib; do not add pytest or change `pixi.lock`).

**Spec:** [`docs/superpowers/specs/2026-09-14-coord-pool-splitless-call-design.md`](../specs/2026-09-14-coord-pool-splitless-call-design.md)

**Branch / merge:** All implementation commits go on the existing
`feat/taps-pooled-methy-call` branch (continue the TAPS pooled-caller work).
Do **not** create `feat/coord-pool-splitless-call`. Do **not** merge to
`main` in this plan.

**Safety:** Do not modify `work/test-DNAme-EMSeq`. Biscuit compare uses a sibling work directory.

---

## File map

| File | Responsibility |
|------|----------------|
| `scripts/methy_call_core.py` | Polarity-aware classify, CB helpers, counts scan, smoke write-set, pileup/shuffle shared by both callers |
| `scripts/methy_caller.py` | TAPS CLI; `converted_means_methylated=True`; re-export wrappers so existing tests keep importing `scripts.methy_caller` |
| `scripts-emseq/methy_caller.py` | EMSeq CLI; polarity `False` |
| `scripts-emseq/compare_cov_to_biscuit.py` | Join new vs biscuit `.cov` with `pos - 1` |
| `scripts/pool.py` | Host coordinate-sort + index → `pooled.host.sorted.bam` |
| `scripts/call.py` | TAPS orchestrator: new BAM name, skip `pooled.host` as spike, optional smoke/roster/force-recount |
| `scripts-emseq/call.py` | Same orchestrator as TAPS; default caller `scripts-emseq/methy_caller.py`; no biscuit |
| `scripts/make_cmd.py`, `scripts-v2/make_cmd.py`, `scripts-emseq/make_cmd.py` | Drop `split`; `call_smoke`; new host checks |
| `scripts/saturation.py`, `scripts/summary.py`, `scripts-v2/summary.py` | Counts TSV path; host BAM name |
| `scripts/mbias.py`, `scripts-emseq/mbias.py` | Host BAM name |
| `tests/test_methy_call_core.py` | Polarity, counts TSV, smoke, write-set, roster parse |
| `tests/test_methy_caller.py` | Keep existing TAPS tests; update host-input test (roster no longer required) |
| `workflow/*.json` | Drop `split` slurm / `split_smoke`; add `call_smoke` where needed |
| Docs listed in the spec | BREAKING contracts |

Do not delete `scripts/split_bams.py`. Do not edit `doc/TODO.md`. Do not overwrite `work/test-DNAme-EMSeq`.

**Spike discovery trap:** `pooled.host.sorted.bam` matches `pooled.*.sorted.bam`. `discover_spike_bams` must skip stem `pooled.host` (and leftover `pooled.byCB`).

---

### Task 1: Branch and commit spec + plan

**Files:**
- Already present: `docs/superpowers/specs/2026-09-14-coord-pool-splitless-call-design.md`
- Create: `docs/superpowers/plans/2026-09-14-coord-pool-splitless-call.md` (this file)

- [ ] **Step 1: Stay on `feat/taps-pooled-methy-call`**

```bash
cd /home/jarning/project2026/dbit-matrix
git checkout feat/taps-pooled-methy-call
git branch --show-current
```

Expected: `feat/taps-pooled-methy-call`

If the current branch is anything else, stop. Do not create a new branch.
Do not commit these changes to `main`.

- [ ] **Step 2: Confirm spec status is approved**

First line block of the spec must contain `Status: approved`. If it still says draft, change it.

- [ ] **Step 3: Commit spec + plan only on `feat/taps-pooled-methy-call`**

Confirm `git branch --show-current` is `feat/taps-pooled-methy-call` before `git add`.

```bash
git add docs/superpowers/specs/2026-09-14-coord-pool-splitless-call-design.md \
        docs/superpowers/plans/2026-09-14-coord-pool-splitless-call.md
git commit -m "$(cat <<'EOF'
docs: add coord-pool splitless call spec and plan

EOF
)"
```

Do not add `doc/TODO.md` or work trees.

---

### Task 2: Failing unit tests for core polarity, counts, and smoke

**Files:**
- Create: `tests/test_methy_call_core.py`

- [ ] **Step 1: Write the test module**

```python
#!/usr/bin/env python3
from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts import methy_call_core as core


class ClassifyCgPolarityTests(unittest.TestCase):
    def test_taps_tg_methylated_cg_unmethylated(self) -> None:
        self.assertTrue(core.classify_cg_observation("TG", True))
        self.assertTrue(core.classify_cg_observation("CA", True))
        self.assertFalse(core.classify_cg_observation("CG", True))
        self.assertIsNone(core.classify_cg_observation("CC", True))

    def test_emseq_inverts(self) -> None:
        self.assertFalse(core.classify_cg_observation("TG", False))
        self.assertFalse(core.classify_cg_observation("CA", False))
        self.assertTrue(core.classify_cg_observation("CG", False))
        self.assertIsNone(core.classify_cg_observation("TT", False))


class ClassifyChPolarityTests(unittest.TestCase):
    def test_taps_plus_t_methylated(self) -> None:
        self.assertEqual(core.classify_ch_observation("T", "+", True), (1, 0))
        self.assertEqual(core.classify_ch_observation("C", "+", True), (0, 1))

    def test_emseq_plus_c_methylated(self) -> None:
        self.assertEqual(core.classify_ch_observation("C", "+", False), (1, 0))
        self.assertEqual(core.classify_ch_observation("T", "+", False), (0, 1))

    def test_emseq_minus_g_methylated(self) -> None:
        self.assertEqual(core.classify_ch_observation("G", "-", False), (1, 0))
        self.assertEqual(core.classify_ch_observation("A", "-", False), (0, 1))


class CountsTsvTests(unittest.TestCase):
    def test_write_and_load_sorted(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "per_spot_read_counts.tsv"
            core.write_counts_tsv(path, {"01_24": 70, "02_16": 550, "00_00": 3})
            text = path.read_text(encoding="utf-8")
            self.assertEqual(
                text.splitlines()[0],
                "X_index\tY_index\tspot\treads",
            )
            rows = core.load_counts_tsv(path)
            self.assertEqual(
                [(row.spot, row.reads) for row in rows],
                [("00_00", 3), ("01_24", 70), ("02_16", 550)],
            )


class SmokeWriteSetTests(unittest.TestCase):
    def test_top_16_by_reads_then_spot(self) -> None:
        counts = {f"{i:02d}_{i:02d}": i for i in range(20)}
        write_set = core.resolve_write_set(counts, smoke=True, roster_ids=None)
        self.assertEqual(len(write_set), 16)
        self.assertEqual(write_set[0], "19_19")
        self.assertEqual(write_set[-1], "04_04")
        self.assertNotIn("03_03", write_set)

    def test_roster_overrides_smoke(self) -> None:
        counts = {"00_00": 100, "01_01": 1}
        write_set = core.resolve_write_set(
            counts, smoke=True, roster_ids=["01_01"]
        )
        self.assertEqual(write_set, ["01_01"])

    def test_full_set_without_smoke(self) -> None:
        counts = {"02_02": 5, "00_00": 1}
        self.assertEqual(
            core.resolve_write_set(counts, smoke=False, roster_ids=None),
            ["00_00", "02_02"],
        )


class RosterParseTests(unittest.TestCase):
    def test_counts_tsv_spot_column(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "roster.tsv"
            path.write_text(
                "X_index\tY_index\tspot\treads\n1\t24\t01_24\t70\n",
                encoding="utf-8",
            )
            self.assertEqual(core.load_roster_ids(path), ["01_24"])

    def test_plain_spot_lines(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "spots.txt"
            path.write_text("01_24\n02_16\n", encoding="utf-8")
            self.assertEqual(core.load_roster_ids(path), ["01_24", "02_16"])


class FakeCbRecord:
    def __init__(self, cb: str | None) -> None:
        self._cb = cb

    def has_tag(self, name: str) -> bool:
        return name == "CB" and self._cb is not None

    def get_tag(self, name: str) -> str:
        if name != "CB" or self._cb is None:
            raise KeyError(name)
        return self._cb


class CountAssignedTests(unittest.TestCase):
    def test_assigned_and_skipped(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "barcodes.tsv"
            path.write_text("AAAA\nCCCC\n", encoding="utf-8")
            _, barcode_to_index, barcode_len = core.load_barcode_whitelist(str(path))
            records = [
                FakeCbRecord("AAAA+CCCC"),
                FakeCbRecord("AAAA+CCCC"),
                FakeCbRecord("TTTT+AAAA"),
                FakeCbRecord(None),
            ]
            counts = core.count_records_by_spot(
                records, "CB", barcode_to_index, barcode_len
            )
            self.assertEqual(counts, {"00_01": 2})


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pixi run python -m unittest tests.test_methy_call_core -v
```

Expected: `ModuleNotFoundError` or `ImportError` for `scripts.methy_call_core`.

- [ ] **Step 3: Commit the failing tests**

```bash
git add tests/test_methy_call_core.py
git commit -m "$(cat <<'EOF'
test: add failing core polarity counts and smoke tests

EOF
)"
```

---

### Task 3: Implement `scripts/methy_call_core.py` (classify, counts, smoke)

**Files:**
- Create: `scripts/methy_call_core.py`
- Test: `tests/test_methy_call_core.py`

- [ ] **Step 1: Create the module with the functions the tests import**

Do not move pileup yet. Only the symbols Task 2 imports.

```python
#!/usr/bin/env python3
"""Shared methylation-call helpers (TAPS and EMSeq)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from scripts.split_bams import parse_cb_value, read_barcodes

BarcodeMap = Dict[str, int]
SMOKE_SPOT_LIMIT = 16


@dataclass(frozen=True)
class SpotCount:
    x_index: int
    y_index: int
    spot: str
    reads: int


def classify_cg_observation(
    dinucleotide: str, converted_means_methylated: bool
) -> Optional[bool]:
    if dinucleotide in {"TG", "CA"}:
        converted = True
    elif dinucleotide == "CG":
        converted = False
    else:
        return None
    return converted if converted_means_methylated else not converted


def classify_ch_observation(
    base: str, strand: str, converted_means_methylated: bool
) -> Optional[Tuple[int, int]]:
    if strand == "+":
        if base not in {"C", "T"}:
            return None
        retained = base == "C"
    else:
        if base not in {"G", "A"}:
            return None
        retained = base == "G"
    methylated = retained if not converted_means_methylated else not retained
    return (1, 0) if methylated else (0, 1)


def load_barcode_whitelist(path: str) -> Tuple[List[str], BarcodeMap, int]:
    whitelist = Path(path)
    if not whitelist.is_file():
        raise ValueError(f"barcode whitelist not found: {whitelist}")
    barcodes, barcode_to_index = read_barcodes(path)
    if not barcodes:
        raise ValueError("barcode whitelist is empty")
    barcode_len = len(barcodes[0])
    if any(len(barcode) != barcode_len for barcode in barcodes):
        raise ValueError("barcode whitelist has unequal lengths")
    return barcodes, barcode_to_index, barcode_len


def spot_id_from_cb(
    cb: str, barcode_to_index: BarcodeMap, barcode_len: int
) -> Optional[str]:
    parsed = parse_cb_value(cb, barcode_len)
    if parsed is None:
        return None
    x_barcode, y_barcode = parsed
    x_index = barcode_to_index.get(x_barcode)
    y_index = barcode_to_index.get(y_barcode)
    if x_index is None or y_index is None:
        return None
    return f"{x_index:02d}_{y_index:02d}"


def count_records_by_spot(
    records: Iterable[object],
    cb_tag: str,
    barcode_to_index: BarcodeMap,
    barcode_len: int,
) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for record in records:
        try:
            if hasattr(record, "has_tag") and not record.has_tag(cb_tag):
                continue
            cb = record.get_tag(cb_tag)
        except KeyError:
            continue
        if not isinstance(cb, str):
            cb = str(cb)
        spot = spot_id_from_cb(cb, barcode_to_index, barcode_len)
        if spot is None:
            continue
        counts[spot] = counts.get(spot, 0) + 1
    return counts


def write_counts_tsv(path: Path, counts: Dict[str, int]) -> None:
    rows = []
    for spot, reads in counts.items():
        x_text, y_text = spot.split("_", 1)
        rows.append((int(x_text), int(y_text), spot, int(reads)))
    rows.sort()
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", encoding="utf-8") as handle:
        handle.write("X_index\tY_index\tspot\treads\n")
        for x_index, y_index, spot, reads in rows:
            handle.write(f"{x_index}\t{y_index}\t{spot}\t{reads}\n")
    tmp.replace(path)


def load_counts_tsv(path: Path) -> List[SpotCount]:
    rows: List[SpotCount] = []
    with path.open(encoding="utf-8") as handle:
        header = handle.readline()
        if "spot" not in header:
            raise ValueError(f"counts TSV missing header: {path}")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            rows.append(
                SpotCount(
                    x_index=int(fields[0]),
                    y_index=int(fields[1]),
                    spot=fields[2],
                    reads=int(fields[3]),
                )
            )
    rows.sort(key=lambda row: (row.x_index, row.y_index))
    return rows


def load_roster_ids(path: Path) -> List[str]:
    text = path.read_text(encoding="utf-8")
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        return []
    if "\t" in lines[0] and "spot" in lines[0].split("\t"):
        header = lines[0].split("\t")
        spot_idx = header.index("spot")
        return [line.split("\t")[spot_idx] for line in lines[1:] if len(line.split("\t")) > spot_idx]
    return lines


def resolve_write_set(
    counts: Dict[str, int],
    smoke: bool,
    roster_ids: Optional[List[str]],
    limit: int = SMOKE_SPOT_LIMIT,
) -> List[str]:
    if roster_ids is not None:
        return list(roster_ids)
    spots = sorted(counts, key=lambda spot: (-counts[spot], spot))
    if smoke:
        return spots[:limit]
    return sorted(counts)


def counts_dict_from_rows(rows: List[SpotCount]) -> Dict[str, int]:
    return {row.spot: row.reads for row in rows}
```

Fix `write_counts_tsv`: do **not** call `path.stat()` before the file exists. Use only `path.with_suffix(path.suffix + ".tmp")` or `path.with_name(path.name + ".tmp")`.

- [ ] **Step 2: Run the new tests**

```bash
pixi run python -m unittest tests.test_methy_call_core -v
```

Expected: all PASS.

- [ ] **Step 3: Commit**

```bash
git add scripts/methy_call_core.py
git commit -m "$(cat <<'EOF'
feat: add shared methy_call_core polarity and counts helpers

EOF
)"
```

---

### Task 4: Point TAPS `methy_caller.py` at core and add host scan CLI

**Files:**
- Modify: `scripts/methy_caller.py`
- Modify: `tests/test_methy_caller.py`

Keep existing TAPS wrappers so `tests/test_methy_caller.py` still imports `classify_cg_observation(dinuc)` with TAPS polarity.

- [ ] **Step 1: Change TAPS classify wrappers to delegate**

Replace the current `classify_cg_observation` / `classify_ch_observation` bodies:

```python
from scripts import methy_call_core as core

CONVERTED_MEANS_METHYLATED = True


def classify_cg_observation(dinucleotide: str) -> Optional[bool]:
    return core.classify_cg_observation(dinucleotide, CONVERTED_MEANS_METHYLATED)


def classify_ch_observation(base: str, strand: str) -> Optional[Tuple[int, int]]:
    return core.classify_ch_observation(base, strand, CONVERTED_MEANS_METHYLATED)
```

Keep `spot_id_from_cb` and `load_barcode_whitelist` as wrappers around `core.*` with the same signatures (existing tests call `mc.spot_id_from_cb` and `mc.load_barcode_whitelist`).

- [ ] **Step 2: Extend `parse_args`**

Replace `--spot-roster` help and add flags:

```python
    parser.add_argument(
        "--spot-roster",
        help="Optional .cov write-set (counts TSV or one XX_YY per line).",
    )
    parser.add_argument(
        "--counts-tsv",
        help="Host counts output. Default: <out-dir>/../per_spot_read_counts.tsv",
    )
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--force-recount", action="store_true")
```

Host mode no longer requires `--spot-roster` in `validate_runtime`. Require `--barcodes` and `--out-dir` only.

Default counts path when `--counts-tsv` is omitted: `Path(args.out_dir).resolve().parent / "per_spot_read_counts.tsv"` (so `--out-dir coverage/host` → `coverage/per_spot_read_counts.tsv`).

- [ ] **Step 3: Host main() builds roster from scan**

After whitelist load, replace `roster = load_spot_roster(args.spot_roster)` with:

```python
counts_path = Path(args.counts_tsv) if args.counts_tsv else Path(args.out_dir).resolve().parent / "per_spot_read_counts.tsv"
if counts_path.is_file() and not args.force_recount:
    counts = core.counts_dict_from_rows(core.load_counts_tsv(counts_path))
else:
    with pysam.AlignmentFile(args.bam, "rb") as bam_in:
        counts = core.count_records_by_spot(
            bam_in.fetch(until_eof=True),
            args.cb_tag,
            barcode_to_index,
            barcode_len,
        )
    core.write_counts_tsv(counts_path, counts)
roster_ids = core.load_roster_ids(Path(args.spot_roster)) if args.spot_roster else None
roster = core.resolve_write_set(counts, args.smoke, roster_ids)
```

Print `smoke`, `counts_tsv`, `roster_spots`. Dry-run must print these and write nothing (do not scan).

- [ ] **Step 4: Update `MissingHostInputsTests`**

`test_missing_barcodes_exits_1_without_traceback` must drop `--spot-roster` from argv. Host without `--barcodes` still exits 1.

- [ ] **Step 5: Run both test modules**

```bash
pixi run python -m unittest tests.test_methy_caller tests.test_methy_call_core -v
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add scripts/methy_caller.py tests/test_methy_caller.py
git commit -m "$(cat <<'EOF'
feat: scan pooled BAM for TAPS host counts and smoke write-set

EOF
)"
```

---

### Task 5: Add `scripts-emseq/methy_caller.py`

**Files:**
- Create: `scripts-emseq/methy_caller.py`

This file is an entrypoint, not a fork of pileup. Import core + TAPS module functions that are polarity-agnostic (`process_chromosome`, `shuffle_host_parts`, `validate_chromosomes`, …) **only if** those functions already take classify through the TAPS wrappers.

**Do not** import TAPS `classify_*` as the polarity source.

Cleaner approach that avoids hidden TAPS polarity:

1. Move remaining pileup/shuffle/`main` orchestration into `methy_call_core.py` as `run_caller(args, converted_means_methylated: bool) -> int`.
2. TAPS `main` becomes `return core.run_caller(parse_args(argv), True)`.
3. EMSeq `main` parses the same flags and calls `return core.run_caller(args, False)`.

If moving `run_caller` in this task would make the diff huge, do it here anyway. Do not leave pileup duplicated.

EMSeq argparse description: `"EMSeq methylation caller"`. Print `[emseq.methy_caller]` or reuse `[methy_caller]` — use `[emseq.methy_caller]` so logs are distinct.

- [ ] **Step 1: Move `run_caller` into core and keep TAPS CLI as a wrapper**

`run_caller` must use `classify_cg_observation(dinuc, polarity)` / `classify_ch_observation(base, strand, polarity)` inside the pileup loops. Update `count_cg_column` / `count_ch_column` signatures to take `converted_means_methylated: bool`.

- [ ] **Step 2: Write the EMSeq entry**

```python
#!/usr/bin/env python3
"""EMSeq methylation caller (pooled host or aggregate spike/mito)."""

from __future__ import annotations

import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from scripts import methy_call_core as core
from scripts.methy_caller import parse_args


def main(argv=None) -> int:
    args = parse_args(argv)
    return core.run_caller(args, converted_means_methylated=False, log_prefix="emseq.methy_caller")


if __name__ == "__main__":
    raise SystemExit(main())
```

TAPS `main` becomes `return core.run_caller(parse_args(argv), True, "methy_caller")`.

- [ ] **Step 3: Help / dry-run**

```bash
pixi run python scripts-emseq/methy_caller.py --help
pixi run python scripts/methy_caller.py --help
```

Expected: both show `--smoke`, `--force-recount`, `--counts-tsv`, optional `--spot-roster`.

```bash
pixi run python scripts-emseq/methy_caller.py --mode aggregate --bam x --reference y --chromosomes chr1 --output /tmp/x --dry-run
```

Expected: prints parameters and `dry_run=1`; exit 0.

- [ ] **Step 4: Re-run unit tests**

```bash
pixi run python -m unittest tests.test_methy_caller tests.test_methy_call_core -v
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add scripts/methy_call_core.py scripts/methy_caller.py scripts-emseq/methy_caller.py
git commit -m "$(cat <<'EOF'
feat: add EMSeq methy_caller using shared core polarity

EOF
)"
```

---

### Task 6: Coordinate-sort and index host pool

**Files:**
- Modify: `scripts/pool.py` host block (currently writes `pooled.cb.bam` / `pooled.byCB.bam` with `-t CB`, no index)

- [ ] **Step 1: Change the host outputs**

Replace the host block so it matches spike: unsorted temp, coordinate sort, index, delete temp.

```python
        host_tmp = pooled_dir / "pooled.host.unsorted.bam"
        host_sorted = pooled_dir / "pooled.host.sorted.bam"
        cat_cmd = [
            args.samtools_bin, "cat", "-o", str(host_tmp),
            *[str(path) for path in host_shards],
        ]
        sort_cmd = [
            args.samtools_bin, "sort",
            "-m", args.host_sort_mem,
            "-@", str(args.samtools_threads),
            "-o", str(host_sorted),
            str(host_tmp),
        ]
        index_cmd = [args.samtools_bin, "index", str(host_sorted)]
        print("[pool] target=host")
        print(f"[pool] host_shard_count={len(host_shards)}")
        print(f"[pool] command={quoted(cat_cmd)}")
        print(f"[pool] command={quoted(sort_cmd)}")
        print(f"[pool] command={quoted(index_cmd)}")
        if not args.dry_run:
            run_command(cat_cmd)
            run_command(sort_cmd)
            run_command(index_cmd)
            host_tmp.unlink(missing_ok=True)
```

No `-t CB`. Dry-run still prints all three commands.

- [ ] **Step 2: Dry-run**

```bash
pixi run python scripts/pool.py --help
pixi run python scripts/pool.py --work-path /tmp/missing-work --dry-run
```

`--help` must succeed. Dry-run without shards should still fail with `no host shards` (existing behavior) unless you pass a work dir that has shards. That is OK.

- [ ] **Step 3: Commit**

```bash
git add scripts/pool.py
git commit -m "$(cat <<'EOF'
feat: coordinate-sort and index pooled.host.sorted.bam

EOF
)"
```

---

### Task 7: TAPS `call.py` new BAM name and flags

**Files:**
- Modify: `scripts/call.py`

- [ ] **Step 1: Host BAM constant**

Use `pooled/pooled.host.sorted.bam` everywhere `pooled.byCB.bam` appears (`run_host_mode`, `ensure_host_mito_bam` fallback).

- [ ] **Step 2: Spike discovery**

In `discover_spike_bams`, skip stems `pooled.host` and `pooled.byCB`:

```python
        spike_name = stem.split(".", 1)[1]
        if spike_name in {"host", "byCB"}:
            continue
```

- [ ] **Step 3: Forward new caller flags**

Add `--smoke`, `--force-recount`, `--spot-roster` on `call.py`. `build_caller_command` for host mode:

- drop required `--spot-roster` pointing at `split_bams/`
- pass `--counts-tsv` `work/coverage/per_spot_read_counts.tsv`
- pass `--smoke` / `--force-recount` / `--spot-roster` when set

`run_host_mode` must not require `split_bams/per_spot_read_counts.tsv`.

- [ ] **Step 4: Help**

```bash
pixi run python scripts/call.py --help
```

Expected: mentions `pooled.host.sorted.bam`, `--smoke`, no `split_bams` requirement in the description.

- [ ] **Step 5: Commit**

```bash
git add scripts/call.py
git commit -m "$(cat <<'EOF'
feat: point TAPS call at pooled.host.sorted.bam and call-written counts

EOF
)"
```

---

### Task 8: Replace EMSeq biscuit `call.py`

**Files:**
- Modify: `scripts-emseq/call.py` (replace biscuit pileup with a TAPS-shaped orchestrator)

- [ ] **Step 1: Replace the file**

Copy `scripts/call.py` and change:

- module description: EMSeq, no biscuit
- `--caller-script` default: `scripts-emseq/methy_caller.py`
- keep the same host/mito/spike / trim / context-mode / smoke flags
- delete all biscuit pileup, `vcf2bed`, `mergecg`, tabix, `call_left_trimming`, `call_right_trimming`, `call_host_threads`, `discover_host_spot_bams`, mito-from-spot-cov stripping

Do not keep a biscuit fallback.

- [ ] **Step 2: Help**

```bash
pixi run python scripts-emseq/call.py --help
```

Expected: no `--biscuit-bin`, no `--call-left-trimming`. Has `--context-mode`, `--smoke`, `--barcodes`.

- [ ] **Step 3: Commit**

```bash
git add scripts-emseq/call.py
git commit -m "$(cat <<'EOF'
feat: replace EMSeq biscuit call with pooled pysam caller

EOF
)"
```

---

### Task 9: Saturation, summary, mbias host paths

**Files:**
- Modify: `scripts/saturation.py` (`read_counts_path = work_path / "coverage" / "per_spot_read_counts.tsv"`)
- Modify: `scripts/summary.py` (same TSV; `host_bam_path = pooled_dir / "pooled.host.sorted.bam"`)
- Modify: `scripts-v2/summary.py` (same two paths)
- Modify: `scripts/mbias.py` (host pooled BAM name)
- Modify: `scripts-emseq/mbias.py` (host pooled BAM name)

Grep each file for `pooled.byCB` and `split_bams/per_spot_read_counts` and replace. No other behavior changes.

- [ ] **Step 1: Apply the path replacements**

- [ ] **Step 2: Grep to confirm**

```bash
rg -n 'pooled\.byCB|split_bams/per_spot_read_counts' scripts scripts-v2 scripts-emseq
```

Expected: no matches except `scripts/split_bams.py` (legacy script) and comments if any. `discover_spike_bams` may still mention `byCB` as a skip name.

- [ ] **Step 3: Commit**

```bash
git add scripts/saturation.py scripts/summary.py scripts-v2/summary.py scripts/mbias.py scripts-emseq/mbias.py
git commit -m "$(cat <<'EOF'
feat: read call-written counts and pooled.host.sorted.bam downstream

EOF
)"
```

---

### Task 10: Drop `split` from TAPS `make_cmd.py`

**Files:**
- Modify: `scripts/make_cmd.py`

- [ ] **Step 1: Stage lists**

`STAGE_SEQUENCE` becomes:

```python
STAGE_SEQUENCE = [
    "fastp_split",
    "demux_extract_bc",
    "align",
    "pool",
    "mbias",
    "call",
    "saturation",
    "summary",
]
```

Delete `"split"` from `STAGE_REQUIRED_FIELDS`. `--stage split` then fails via `STAGE_CHOICES`.

- [ ] **Step 2: CLI / settings**

- Remove `--split-smoke` argparse. Add `--call-smoke` (`store_true`).
- Keep `--split-barcodes` as the CB whitelist passed to `call --barcodes`.
- Delete `build_split_command`, `build_split_sort_command`, `generate_split_slurm_submit_script`, and the `elif stage == "split":` command-generation branch (both local `.sh` and Slurm `.sbatch` pairs).
- Delete `split_bams_slurm_*` / `split_sort_slurm_*` settings.
- In `all` chaining, remove the special-case that submits `split_bams` then `sort`.

- [ ] **Step 3: Call command + workdir checks**

`build_call_command` / host invocation:

- `--barcodes` from `split_barcodes`
- `--smoke` when `call_smoke`
- do not pass `--spot-roster` by default

Workdir checks for `call`:

- require `pooled/pooled.host.sorted.bam` and `.bai` or `.csi`
- require barcode TSV
- do **not** require `split_bams/per_spot_read_counts.tsv`

Workdir checks for `pool` / `mbias` / `summary`: host BAM name `pooled.host.sorted.bam`.

- [ ] **Step 4: Dry-run**

```bash
pixi run python scripts/make_cmd.py --version
pixi run python scripts/make_cmd.py --workflow-config workflow/dbit_taps_test.json --stage all --runner local --dry-run
```

Expected: no `split_bams` / `05_split` in printed commands. Pool command contains `pooled.host.sorted.bam` and `samtools index`. Stage order `pool` then `mbias` then `call`.

- [ ] **Step 5: Commit**

```bash
git add scripts/make_cmd.py
git commit -m "$(cat <<'EOF'
feat: remove split from TAPS make_cmd mainline

EOF
)"
```

---

### Task 11: Drop `split` from TAPS v2 and EMSeq `make_cmd.py`

**Files:**
- Modify: `scripts-v2/make_cmd.py` (`STAGE_SEQUENCE` without `split`; `--call-smoke`; host checks for `pooled.host.sorted.bam` + index; delete split command/Slurm helpers)
- Modify: `scripts-emseq/make_cmd.py`

EMSeq extras:

- `build_call_command` must invoke `scripts-emseq/call.py` with TAPS-style `--r1-left-trimming` … and `--context-mode`, not `--call-left-trimming` / `--biscuit-bin` / `--host-threads`
- Add `call_context_mode` (default `cg`) and `call_trimming` nest like TAPS (`scripts/make_cmd.py` is the template)
- Delete `call_left_trimming`, `call_right_trimming`, `call_host_threads`, `call_spike_threads` settings
- Renumber generated scripts if prefixes assumed `05_split` (`mbias` was 06). After dropping split, keep existing numbers for mbias/call **or** compact them. Prefer **compact**: pool stays `04_`, mbias `05_`, call `06_`, saturation `07_`, summary `08_` so `--stage all` order matches filenames. Update any docstring that lists `05_split`.

- [ ] **Step 1: Apply the driver edits**

- [ ] **Step 2: Dry-run both**

```bash
pixi run python scripts-v2/make_cmd.py --workflow-config workflow/dbit_taps_v2_test.json --stage all --runner local --dry-run
pixi run python scripts-emseq/make_cmd.py --workflow-config workflow/dbit_emseq_test.json --stage all --runner local --dry-run
pixi run python scripts-emseq/make_cmd.py --workflow-config workflow/dbit_emseq_test.json --stage all --runner slurm --dry-run
```

Expected: no `split` stage; EMSeq call line has `scripts-emseq/call.py` and `scripts-emseq/methy_caller.py` or default caller; no biscuit pileup flags.

- [ ] **Step 3: Commit**

```bash
git add scripts-v2/make_cmd.py scripts-emseq/make_cmd.py
git commit -m "$(cat <<'EOF'
feat: remove split from v2 and EMSeq make_cmd mainlines

EOF
)"
```

---

### Task 12: Workflow JSON

**Files:** every `workflow/*.json` that contains `split_smoke`, `slurm.split`, `call_left_trimming`, or `call_right_trimming`.

Known TAPS e2e files from the previous caller plan plus:

- `workflow/dbit_taps_test.json`
- `workflow/dbit_taps_slurm.json`
- `workflow/dbit_taps_v2_test.json`
- `workflow/dbit_taps_v2_slurm.json`
- `workflow/dbit_emseq_test.json`
- `workflow/dbit_emseq_slurm.json`
- `workflow/dbit_emseq_slurm_HL.json`
- remaining `workflow/*TAPS*_e2e.json` / `*TAPSbeta*`

- [ ] **Step 1: Mechanical replacements**

For each file:

- Delete the `"split"` object under `"slurm"`
- Replace `"split_smoke": true` with `"call_smoke": true`
- Delete `"split_threads_read"`, `"split_threads_write"`, `"split_sort_jobs"` if present (unused)
- Keep `"split_barcodes"` (whitelist for call)
- EMSeq: delete `call_left_trimming`, `call_right_trimming`, `call_host_threads`, `call_spike_threads`
- EMSeq: add `"call_context_mode": "cg"` and, if those files had left/right 3, add:

```json
    "call_trimming": {
        "host": {
            "r1_left": 3,
            "r1_right": 3,
            "r2_left": 3,
            "r2_right": 3
        },
        "spike": {
            "r1_left": 3,
            "r1_right": 3,
            "r2_left": 3,
            "r2_right": 3
        }
    }
```

Copy the nest shape from `workflow/dbit_taps_test.json` if that file already has `call_trimming`.

- [ ] **Step 2: Re-dry-run EMSeq test config**

```bash
pixi run python scripts-emseq/make_cmd.py --workflow-config workflow/dbit_emseq_test.json --stage all --runner local --dry-run
```

Expected: success; `call_smoke` forwarded as `--smoke` on host call.

- [ ] **Step 3: Commit**

```bash
git add workflow/*.json
git commit -m "$(cat <<'EOF'
chore: drop split slurm blocks and add call_smoke

EOF
)"
```

---

### Task 13: Docs

**Files:**
- `AGENTS.md` — mainline without `split`; pool writes `pooled.host.sorted.bam` plus index; delete required Slurm `split` layout
- `docs/developers/contracts.md` — pool / split / call / saturation / summary / mbias paths
- `docs/developers/config-reference.md` — `call_smoke`; drop biscuit call keys; `split_*` not required for mainline
- `docs/developers/architecture.md` and `docs/developers/runner-reference.md` — stage list
- `docs/users/taps.md` and `docs/users/emseq.md`
- `docs/maintenance/taps-regression.md` and `docs/maintenance/emseq-regression.md` — stage order `pool -> mbias -> call`; biscuit compare pointer
- `docs/maintenance/changelog.md` under `## Current` with `BREAKING` bullets from the spec
- `README.md` only if it still lists `split` in the mainline sentence

Do not edit `doc/`. Do not duplicate defaults; link to contracts / config-reference.

`contracts.md` `split` section: keep a short “legacy script” note that `scripts/split_bams.py` exists but is not a workflow stage. Do not document per-spot BAM as a mainline output.

- [ ] **Step 1: Edit the pages**

- [ ] **Step 2: Commit**

```bash
git add AGENTS.md README.md docs
git commit -m "$(cat <<'EOF'
docs: record splitless mainline and pooled.host.sorted.bam

EOF
)"
```

Do not `git add` `docs/superpowers` if those files were already committed in Task 1 unless this task edits them.

---

### Task 14: Biscuit compare script

**Files:**
- Create: `scripts-emseq/compare_cov_to_biscuit.py`

- [ ] **Step 1: Implement**

```python
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
    shared_meth: list[float] = []
    shared_unmeth: list[float] = []
    shared_abs_pct: list[float] = []
    disagree = 0
    shared = 0
    only_base = 0
    only_new = 0
    covered_shared = 0
    for spot in spots:
        base = parse_cov(next(baseline_root.rglob(f"{spot}.CG.cov")), args.biscuit_pos_shift)
        new_matches = list(new_root.rglob(f"{spot}.CG.cov"))
        new = parse_cov(new_matches[0], 0) if new_matches else {}
        keys = set(base) | set(new)
        only_base += len(set(base) - set(new))
        only_new += len(set(new) - set(base))
        for key in set(base) & set(new):
            shared += 1
            bm, bu, bp = base[key]
            nm, nu, np_ = new[key]
            shared_meth.append(float(bm))
            shared_unmeth.append(float(nm))
            shared_abs_pct.append(abs(bp - np_))
            if (bm + bu) > 0 and (nm + nu) > 0:
                covered_shared += 1
                base_m = bm > bu
                new_m = nm > nu
                if base_m != new_m:
                    disagree += 1
    frac_disagree = (disagree / covered_shared) if covered_shared else float("nan")
    print(f"[compare] spots={len(spots)}")
    print(f"[compare] shared_sites={shared}")
    print(f"[compare] only_baseline={only_base}")
    print(f"[compare] only_new={only_new}")
    print(f"[compare] covered_shared={covered_shared}")
    print(f"[compare] disagree_majority_fraction={frac_disagree}")
    print(f"[compare] pearson_meth={pearson(shared_meth, [float(parse_cov(next(new_root.rglob(f'{s}.CG.cov')), 0).get(k, (0,0,0))[0]) for s in spots for k in []])}")
    return 0
```

**Do not ship that last `pearson_meth` line.** It is wrong. Print Pearson as:

```python
    meth_corr = pearson(shared_meth, shared_new_meth)
```

Collect `shared_new_meth` / `shared_new_unmeth` in the same loop as `bm, bu` / `nm, nu`. The implementation in the repo must be the loop version below, not the broken one-liner.

Correct loop body:

```python
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
    print(f"[compare] pearson_meth={pearson(baseline_meth, new_meth)}")
    print(f"[compare] pearson_unmeth={pearson(baseline_unmeth, new_unmeth)}")
    print(f"[compare] mean_abs_pct_delta={sum(shared_abs_pct)/len(shared_abs_pct) if shared_abs_pct else float('nan')}")
    print(f"[compare] disagree_majority_fraction={frac_disagree}")
    if covered_shared and frac_disagree >= 0.5:
        return 1
    return 0
```

Exit 1 if `disagree_majority_fraction >= 0.5`. Exit 0 if there are no covered shared sites **only** after printing a warning; treat that as fail (`return 1`) so an empty compare cannot pass.

- [ ] **Step 2: Help**

```bash
pixi run python scripts-emseq/compare_cov_to_biscuit.py --help
```

- [ ] **Step 3: Commit**

```bash
git add scripts-emseq/compare_cov_to_biscuit.py
git commit -m "$(cat <<'EOF'
feat: add biscuit vs pooled EMSeq coverage compare script

EOF
)"
```

---

### Task 15: Dry-run gate (required)

- [ ] **Step 1: Version and help**

```bash
pixi run python scripts/make_cmd.py --version
pixi run python scripts-v2/make_cmd.py --version
pixi run python scripts-emseq/make_cmd.py --version
pixi run python scripts/pool.py --help
pixi run python scripts/call.py --help
pixi run python scripts-emseq/call.py --help
pixi run python scripts/methy_caller.py --help
pixi run python scripts-emseq/methy_caller.py --help
```

Expected: versions print; no `split` in EMSeq/TAPS call help as a required input.

- [ ] **Step 2: Workflow dry-runs**

```bash
pixi run python scripts/make_cmd.py --workflow-config workflow/dbit_taps_test.json --stage all --runner local --dry-run
pixi run python scripts-v2/make_cmd.py --workflow-config workflow/dbit_taps_v2_test.json --stage all --runner local --dry-run
pixi run python scripts-emseq/make_cmd.py --workflow-config workflow/dbit_emseq_test.json --stage all --runner local --dry-run
pixi run python scripts-emseq/make_cmd.py --workflow-config workflow/dbit_emseq_test.json --stage all --runner slurm --dry-run
```

Expected: each succeeds; printed scripts have no `split_bams.py` and no `05_split`; pool indexes `pooled.host.sorted.bam`.

- [ ] **Step 3: Unit tests**

```bash
pixi run python -m unittest tests.test_methy_call_core tests.test_methy_caller -v
```

Expected: PASS.

If anything fails, fix it and make a **new** commit (`fix: ...`). Do not amend.

---

### Task 16: Biscuit regression on frozen work tree

**Do not write into `work/test-DNAme-EMSeq`.**

- [ ] **Step 1: Sibling directory and host BAM**

```bash
mkdir -p work/test-DNAme-EMSeq-pooled-call/pooled \
         work/test-DNAme-EMSeq-pooled-call/coverage
pixi run samtools sort -o work/test-DNAme-EMSeq-pooled-call/pooled/pooled.host.sorted.bam \
  work/test-DNAme-EMSeq/pooled/pooled.byCB.bam
pixi run samtools index work/test-DNAme-EMSeq-pooled-call/pooled/pooled.host.sorted.bam
```

Write a roster of the 16 existing biscuit spots:

```bash
find work/test-DNAme-EMSeq/coverage/host -name '*.CG.cov' -printf '%f\n' \
  | sed 's/\.CG\.cov$//' | sort > work/test-DNAme-EMSeq-pooled-call/biscuit_spots.txt
```

- [ ] **Step 2: Run new EMSeq host call**

Use `workflow/dbit_emseq_test.json` for reference, barcodes (`configs/barcodes_50a_CtoT.tsv`), chromosomes, and trims (map old left/right 3 into four trim knobs = 3). Example:

```bash
pixi run python scripts-emseq/call.py \
  --work-path work/test-DNAme-EMSeq-pooled-call \
  --mode host \
  --reference-file /mnt/e/LiLab_HL/resource/biscuit/mm10/genome.fa \
  --chromosomes chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chrX \
  --barcodes configs/barcodes_50a_CtoT.tsv \
  --spot-roster work/test-DNAme-EMSeq-pooled-call/biscuit_spots.txt \
  --r1-left-trimming 3 --r1-right-trimming 3 \
  --r2-left-trimming 3 --r2-right-trimming 3 \
  --jobs 8
```

If the reference path is missing on the machine, stop and tell the user. Do not invent another FASTA.

- [ ] **Step 3: Compare**

```bash
pixi run python scripts-emseq/compare_cov_to_biscuit.py \
  --baseline-host work/test-DNAme-EMSeq/coverage/host \
  --new-host work/test-DNAme-EMSeq-pooled-call/coverage/host
```

Expected: exit 0; `disagree_majority_fraction < 0.5`.

- [ ] **Step 4: Record numbers**

Append the printed metrics under `docs/maintenance/changelog.md` `## Current` as one bullet, then commit:

```bash
git add docs/maintenance/changelog.md
git commit -m "$(cat <<'EOF'
docs: record EMSeq biscuit compare metrics

EOF
)"
```

If the compare fails the polarity bar, stop. Do not “fix” biscuit files. Debug the EMSeq classifier.

---

## Self-review

**Spec coverage**

| Spec item | Task |
|-----------|------|
| `pooled.host.sorted.bam` + index | 6 |
| Drop `split` from three drivers | 10, 11 |
| Call-written counts TSV + smoke/roster | 2, 3, 4 |
| Core + two entrypoints, no `--assay` on TAPS | 3, 4, 5 |
| EMSeq polarity / no biscuit call | 5, 8 |
| EMSeq `call_trimming` | 11, 12 |
| Downstream path updates | 9 |
| Docs + BREAKING | 13 |
| Unit tests | 2, 3, 4, 15 |
| Biscuit compare, frozen baseline | 14, 16 |
| Skip `pooled.host` as spike | 7, 8 |

**Placeholders:** none. Task 14 spells the compare loop; do not implement the broken Pearson one-liner.

**Types:** `converted_means_methylated: bool`, `SpotCount`, `resolve_write_set(counts, smoke, roster_ids)`, counts path `coverage/per_spot_read_counts.tsv`, host BAM `pooled.host.sorted.bam`.
