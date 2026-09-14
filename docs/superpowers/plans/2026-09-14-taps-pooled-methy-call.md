# TAPS pooled methylation caller Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace TAPS/TAPS v2 spot-BAM calling with a template-aligned pooled caller that still writes the current per-spot and aggregate `.cov` files.

**Architecture:** One `scripts/methy_caller.py` with `--mode host|aggregate`. Host pileups `pooled/pooled.byCB.bam` in genomic intervals, maps `CB` with `split_bams.py` rules, and shuffles part files into `coverage/host/<XX>/<XX>_<YY>.*.cov`. Spike/mito use the same classifiers without barcodes. `scripts/call.py` orchestrates; `methy_caller_CH.py` is deleted.

**Tech Stack:** Python 3.11, pysam, pixi, unittest (stdlib; do not add pytest / change `pixi.lock`).

**Spec:** [`docs/superpowers/specs/2026-09-14-taps-call-alignment-design.md`](../specs/2026-09-14-taps-call-alignment-design.md)

**Branch / merge:** All implementation commits go on `feat/taps-pooled-methy-call`. Do **not** merge to `main` in this plan. Merge only after Task 9 (`work/L126-brain-S2-TAPS_1M`) passes and the user asks.

**Sample:** `work/L126-brain-S2-TAPS_1M` already has pooled BAMs, smoke split (16 spots in `split_bams/per_spot_read_counts.tsv`), mbias subsample BAM, and old `coverage/` from the previous caller. Back up that `coverage/` before overwriting. `pooled/pooled.byCB.bam` currently has **no index**; Task 9 indexes it.

---

## File map

| File | Responsibility |
|------|----------------|
| `scripts/methy_caller.py` | TAPS caller: classification, interval pileup, host shuffle, aggregate write |
| `scripts/call.py` | Stage wrapper: host pooled BAM + barcodes + roster; mito/spike aggregate |
| `scripts/methy_caller_CH.py` | Delete |
| `scripts/make_cmd.py` | TAPS driver: `--barcodes`, new defaults, drop CH script |
| `scripts-v2/make_cmd.py` | Same driver changes for TAPS v2 |
| `scripts/workflow_input_checks.py` | Host-call workdir checks for pooled BAM/index/roster |
| `tests/test_methy_caller.py` | Unit tests for classifiers, neighbors, CB, dry-run |
| `docs/developers/contracts.md` | `call` inputs |
| `docs/developers/config-reference.md` | defaults, `call_jobs` meaning |
| `docs/maintenance/changelog.md` | BREAKING |
| `docs/maintenance/taps-regression.md` | drop CH-script checks |
| `docs/users/taps.md` | only if user-facing call inputs change (pooled BAM + barcodes) |
| `workflow/dbit_taps_test.json`, `workflow/dbit_taps_slurm.json`, `workflow/dbit_taps_v2_test.json`, `workflow/dbit_taps_v2_slurm.json`, `workflow/L126-brain-S2-TAPS_1M_e2e.json` | drop `call_ch_caller_script`; `call_max_depth` 50000 |
| remaining `workflow/*TAPS*_e2e.json` | drop `call_ch_caller_script`; set `call_max_depth` 50000 |

Do not edit `scripts-emseq/*`. Do not edit `doc/TODO.md`.

---

### Task 1: Create feature branch and commit the spec

**Files:**
- Already present: `docs/superpowers/specs/2026-09-14-taps-call-alignment-design.md`
- Create: `docs/superpowers/plans/2026-09-14-taps-pooled-methy-call.md` (this file)

- [ ] **Step 1: Create the branch from current `main`**

```bash
cd /home/jarning/project2026/dbit-matrix
git checkout main
git checkout -b feat/taps-pooled-methy-call
```

Expected: `On branch feat/taps-pooled-methy-call`

- [ ] **Step 2: Mark the spec approved**

In `docs/superpowers/specs/2026-09-14-taps-call-alignment-design.md`, change:

```markdown
Status: draft, pending user review
```

to:

```markdown
Status: approved
```

- [ ] **Step 3: Commit spec + plan only**

```bash
git add docs/superpowers/specs/2026-09-14-taps-call-alignment-design.md \
        docs/superpowers/plans/2026-09-14-taps-pooled-methy-call.md
git commit -m "$(cat <<'EOF'
docs: add TAPS pooled caller spec and implementation plan

EOF
)"
```

Do not add `doc/TODO.md` or unrelated dirty files.

---

### Task 2: Failing unit tests for classification and CB mapping

**Files:**
- Create: `tests/__init__.py`
- Create: `tests/test_methy_caller.py`

- [ ] **Step 1: Create the test module**

`tests/__init__.py` is empty.

`tests/test_methy_caller.py`:

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

from scripts import methy_caller as mc
from scripts.split_bams import parse_cb_value, read_barcodes


class ClassifyCgTests(unittest.TestCase):
    def test_tg_and_ca_are_methylated(self) -> None:
        self.assertTrue(mc.classify_cg_observation("TG"))
        self.assertTrue(mc.classify_cg_observation("CA"))

    def test_cg_is_unmethylated(self) -> None:
        self.assertFalse(mc.classify_cg_observation("CG"))

    def test_other_dinucleotides_are_ignored(self) -> None:
        self.assertIsNone(mc.classify_cg_observation("CC"))
        self.assertIsNone(mc.classify_cg_observation("TT"))


class ClassifyChTests(unittest.TestCase):
    def test_plus_c_unmethylated_t_methylated(self) -> None:
        self.assertEqual(mc.classify_ch_observation("C", "+"), (0, 1))
        self.assertEqual(mc.classify_ch_observation("T", "+"), (1, 0))

    def test_minus_g_unmethylated_a_methylated(self) -> None:
        self.assertEqual(mc.classify_ch_observation("G", "-"), (0, 1))
        self.assertEqual(mc.classify_ch_observation("A", "-"), (1, 0))

    def test_wrong_base_ignored(self) -> None:
        self.assertIsNone(mc.classify_ch_observation("G", "+"))
        self.assertIsNone(mc.classify_ch_observation("C", "-"))


class NeighborTests(unittest.TestCase):
    def test_cc_allows_converted_neighbor(self) -> None:
        self.assertTrue(mc.ch_neighbor_matches("CC", "+", "C"))
        self.assertTrue(mc.ch_neighbor_matches("CC", "+", "T"))
        self.assertTrue(mc.ch_neighbor_matches("CC", "-", "G"))
        self.assertTrue(mc.ch_neighbor_matches("CC", "-", "A"))
        self.assertFalse(mc.ch_neighbor_matches("CC", "+", "A"))

    def test_ca_ct_require_exact_neighbor(self) -> None:
        self.assertTrue(mc.ch_neighbor_matches("CA", "+", "A"))
        self.assertFalse(mc.ch_neighbor_matches("CA", "+", "T"))
        self.assertTrue(mc.ch_neighbor_matches("CT", "+", "T"))
        self.assertTrue(mc.ch_neighbor_matches("CA", "-", "T"))
        self.assertTrue(mc.ch_neighbor_matches("CT", "-", "A"))


class YdFilterTests(unittest.TestCase):
    def test_missing_yd_is_kept(self) -> None:
        record = SimpleNamespace(
            get_tag=lambda name: (_ for _ in ()).throw(KeyError(name))
        )
        self.assertFalse(mc.has_unclear_conversion_strand(record))

    def test_yd_u_is_skipped(self) -> None:
        record = SimpleNamespace(get_tag=lambda name: "u")
        self.assertTrue(mc.has_unclear_conversion_strand(record))

    def test_yd_f_is_kept(self) -> None:
        record = SimpleNamespace(get_tag=lambda name: "f")
        self.assertFalse(mc.has_unclear_conversion_strand(record))


class SpotIdTests(unittest.TestCase):
    def test_plus_and_concat_match_split_bams(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "barcodes.tsv"
            path.write_text("AAAA\nCCCC\nGGGG\n", encoding="utf-8")
            _, barcode_to_index = read_barcodes(str(path))
            barcode_len = 4
            self.assertEqual(parse_cb_value("AAAA+CCCC", barcode_len), ("AAAA", "CCCC"))
            self.assertEqual(parse_cb_value("AAAACCCC", barcode_len), ("AAAA", "CCCC"))
            self.assertEqual(
                mc.spot_id_from_cb("AAAA+CCCC", barcode_to_index, barcode_len),
                "00_01",
            )
            self.assertEqual(
                mc.spot_id_from_cb("AAAACCCC", barcode_to_index, barcode_len),
                "00_01",
            )
            self.assertIsNone(
                mc.spot_id_from_cb("TTTT+AAAA", barcode_to_index, barcode_len)
            )

    def test_roster_filters_unlisted_spots(self) -> None:
        allowed = {"00_01"}
        self.assertTrue(mc.spot_is_allowed("00_01", allowed))
        self.assertFalse(mc.spot_is_allowed("02_03", allowed))
        self.assertFalse(mc.spot_is_allowed(None, allowed))


class ReferenceClassifyTests(unittest.TestCase):
    def test_cpg_and_ch_sites(self) -> None:
        # 012345678
        # ACGCATGCT
        seq = "ACGCATGCT"
        self.assertEqual(mc.classify_reference_position(seq, 1), (True, None))
        self.assertEqual(mc.classify_reference_position(seq, 3), (False, ("CA", "+")))
        self.assertEqual(mc.classify_reference_position(seq, 6), (False, ("CA", "-")))


class DryRunTests(unittest.TestCase):
    def test_dry_run_writes_nothing(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "out.CG.cov"
            argv = [
                "--mode",
                "aggregate",
                "--bam",
                str(Path(tmp) / "missing.bam"),
                "--reference",
                str(Path(tmp) / "missing.fa"),
                "--chromosomes",
                "chr1",
                "--output",
                str(out),
                "--dry-run",
            ]
            self.assertEqual(mc.main(argv), 0)
            self.assertFalse(out.exists())


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd /home/jarning/project2026/dbit-matrix
pixi run python -m unittest tests.test_methy_caller -v
```

Expected: `FAIL` or `ERROR` with `cannot import name 'classify_cg_observation'` (or equivalent). Do not implement yet.

- [ ] **Step 3: Commit the failing tests**

```bash
git add tests/__init__.py tests/test_methy_caller.py
git commit -m "$(cat <<'EOF'
test: add failing TAPS caller classification tests

EOF
)"
```

---

### Task 3: Replace `scripts/methy_caller.py`

**Files:**
- Replace: `scripts/methy_caller.py`

This file is the only caller. Keep it in one module. Import barcode parsing from `scripts/split_bams.py` so spot IDs cannot drift.

- [ ] **Step 1: Write the new caller**

Replace `scripts/methy_caller.py` entirely with:

```python
#!/usr/bin/env python3
"""TAPS methylation caller (pooled host or aggregate spike/mito)."""

from __future__ import annotations

import argparse
import multiprocessing
import os
import shutil
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import BinaryIO, Dict, Iterable, List, Optional, Set, Tuple

import pysam

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from scripts.split_bams import parse_cb_value, read_barcodes

TOP_FLAGS: Set[int] = {99, 147}
BOT_FLAGS: Set[int] = {83, 163}
PAIRED_FLAGS: Set[int] = TOP_FLAGS | BOT_FLAGS
CH_FORWARD_CONTEXTS: Set[str] = {"CA", "CC", "CT"}
CH_REVERSE_MAP: Dict[str, str] = {"T": "CA", "G": "CC", "A": "CT"}
CH_REVERSE_NEIGHBOR: Dict[str, str] = {"CA": "T", "CC": "G", "CT": "A"}

WORKER_CHROMOSOME_SEQUENCE: Optional[str] = None
BarcodeMap = Dict[str, int]
CgCounts = Dict[str, Dict[int, List[int]]]
ChCounts = Dict[str, Dict[int, List[object]]]


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="TAPS methylation caller")
    parser.add_argument("--mode", choices=("host", "aggregate"), required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("-c", "--chromosomes", required=True)
    parser.add_argument(
        "--context-mode", choices=("cg", "ch", "both"), default="cg"
    )
    parser.add_argument("--min-base-quality", type=int, default=30)
    parser.add_argument("--min-mapping-quality", type=int, default=10)
    parser.add_argument("--max-depth", type=int, default=50_000)
    parser.add_argument("--batch-size", type=int, default=10_000_000)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--sample-size", type=int, default=None)
    parser.add_argument("--r1-left-trimming", type=int, default=0)
    parser.add_argument("--r1-right-trimming", type=int, default=0)
    parser.add_argument("--r2-left-trimming", type=int, default=0)
    parser.add_argument("--r2-right-trimming", type=int, default=0)
    parser.add_argument("--cb-tag", default="CB")
    parser.add_argument("--barcodes", help="Whitelist TSV for host mode.")
    parser.add_argument(
        "--spot-roster",
        help="split_bams/per_spot_read_counts.tsv for host mode.",
    )
    parser.add_argument("--out-dir", help="Host output directory (coverage/host).")
    parser.add_argument("--output", help="Aggregate CG/CH stem path without suffix.")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args(argv)


def parse_chromosomes(csv: str) -> List[str]:
    chromosomes = [item.strip() for item in csv.split(",") if item.strip()]
    if not chromosomes:
        raise ValueError("chromosome list is empty")
    return chromosomes


def validate_chromosomes(reference: str, chromosomes: List[str]) -> None:
    with pysam.FastaFile(reference) as fasta:
        available = set(fasta.references)
    missing = [chrom for chrom in chromosomes if chrom not in available]
    if missing:
        raise ValueError(
            f"chromosomes not found in reference '{reference}': {','.join(missing)}"
        )


def bam_index_exists(path: Path) -> bool:
    return any(
        candidate.is_file()
        for candidate in (
            Path(str(path) + ".bai"),
            path.with_suffix(".bai"),
            Path(str(path) + ".csi"),
        )
    )


def classify_cg_observation(dinucleotide: str) -> Optional[bool]:
    if dinucleotide in {"TG", "CA"}:
        return True
    if dinucleotide == "CG":
        return False
    return None


def classify_ch_observation(base: str, strand: str) -> Optional[Tuple[int, int]]:
    if strand == "+":
        if base not in {"C", "T"}:
            return None
        methylated = base == "T"
    else:
        if base not in {"G", "A"}:
            return None
        methylated = base == "A"
    return (1, 0) if methylated else (0, 1)


def ch_neighbor_matches(context: str, strand: str, base: str) -> bool:
    observed = base.upper()
    if context == "CC":
        return observed in (("C", "T") if strand == "+" else ("G", "A"))
    expected = context[1] if strand == "+" else CH_REVERSE_NEIGHBOR[context]
    return observed == expected


def classify_reference_position(
    sequence: str, pos: int
) -> Tuple[bool, Optional[Tuple[str, str]]]:
    base = sequence[pos]
    next_base = sequence[pos + 1] if pos + 1 < len(sequence) else ""
    if base == "C":
        dinucleotide = base + next_base
        if dinucleotide == "CG":
            return True, None
        if dinucleotide in CH_FORWARD_CONTEXTS:
            return False, (dinucleotide, "+")
        return False, None
    if base == "G" and pos > 0:
        context = CH_REVERSE_MAP.get(sequence[pos - 1])
        if context is not None:
            return False, (context, "-")
    return False, None


def has_unclear_conversion_strand(record: pysam.AlignedSegment) -> bool:
    try:
        return record.get_tag("YD") == "u"
    except KeyError:
        return False


def is_trimmed(
    record: pysam.AlignedSegment,
    query_pos: int,
    read_len: int,
    r1_left: int,
    r1_right: int,
    r2_left: int,
    r2_right: int,
) -> bool:
    if record.is_read1:
        left_trim, right_trim = r1_left, r1_right
    elif record.is_read2:
        left_trim, right_trim = r2_left, r2_right
    else:
        return False
    if left_trim == 0 and right_trim == 0:
        return False
    if record.is_reverse:
        left_cycle = read_len - query_pos
        right_cycle = query_pos + 1
    else:
        left_cycle = query_pos + 1
        right_cycle = read_len - query_pos
    return left_cycle <= left_trim or right_cycle <= right_trim


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


def spot_is_allowed(spot_id: Optional[str], allowed: Optional[Set[str]]) -> bool:
    if spot_id is None:
        return False
    if allowed is None:
        return True
    return spot_id in allowed


def load_spot_roster(path: str) -> List[str]:
    spots: List[str] = []
    with Path(path).open(encoding="utf-8") as handle:
        header = handle.readline()
        if not header:
            raise ValueError(f"spot roster is empty: {path}")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            spots.append(fields[2])
    if not spots:
        raise ValueError(f"spot roster has no spots: {path}")
    return spots


def host_cov_path(out_dir: Path, spot_id: str, suffix: str) -> Path:
    row = spot_id.split("_", 1)[0]
    return out_dir / row / f"{spot_id}.{suffix}.cov"


def format_cg_line(chrom: str, pos: int, meth: int, unmeth: int) -> str:
    total = meth + unmeth
    pct = round((meth / total) * 100, 2) if total > 0 else 0.0
    return f"{chrom}\t{pos}\t{pos}\t{pct:.2f}\t{meth}\t{unmeth}\n"


def format_ch_line(
    chrom: str, pos: int, meth: int, unmeth: int, context: str, strand: str
) -> str:
    total = meth + unmeth
    pct = round((meth / total) * 100, 2) if total > 0 else 0.0
    return (
        f"{chrom}\t{pos}\t{pos}\t{pct:.2f}\t{meth}\t{unmeth}\t{context}\t{strand}\n"
    )


def effective_batch_size(
    chromosome_length: int, configured_batch_size: int, jobs: int
) -> int:
    return min(configured_batch_size, max(1, chromosome_length // max(jobs, 1)))


def count_cg_column(
    pileup_col,
    dest: Dict[int, List[int]] | CgCounts,
    host: bool,
    cb_tag: str,
    barcode_to_index: Optional[BarcodeMap],
    barcode_len: int,
    allowed_spots: Optional[Set[str]],
    skip_spots: Set[str],
    r1_left: int,
    r1_right: int,
    r2_left: int,
    r2_right: int,
) -> None:
    pos = pileup_col.pos
    for pileup_read in pileup_col.pileups:
        record = pileup_read.alignment
        if record.flag not in PAIRED_FLAGS:
            continue
        if has_unclear_conversion_strand(record):
            continue
        if pileup_read.is_del or pileup_read.is_refskip:
            continue
        qpos = pileup_read.query_position
        if qpos is None:
            continue
        seq = record.query_sequence
        if not seq or qpos + 1 >= len(seq):
            continue
        if is_trimmed(record, qpos, len(seq), r1_left, r1_right, r2_left, r2_right):
            continue
        if is_trimmed(
            record, qpos + 1, len(seq), r1_left, r1_right, r2_left, r2_right
        ):
            continue
        methylated = classify_cg_observation((seq[qpos] + seq[qpos + 1]).upper())
        if methylated is None:
            continue
        key = 0 if methylated else 1
        if not host:
            dest.setdefault(pos, [0, 0])[key] += 1
            continue
        try:
            cb = record.get_tag(cb_tag)
        except KeyError:
            continue
        assert barcode_to_index is not None
        spot = spot_id_from_cb(str(cb), barcode_to_index, barcode_len)
        if not spot_is_allowed(spot, allowed_spots) or spot in skip_spots:
            continue
        assert isinstance(spot, str)
        slot = dest.setdefault(spot, {}).setdefault(pos, [0, 0])
        slot[key] += 1


def count_ch_column(
    pileup_col,
    context: str,
    strand: str,
    dest,
    host: bool,
    cb_tag: str,
    barcode_to_index: Optional[BarcodeMap],
    barcode_len: int,
    allowed_spots: Optional[Set[str]],
    skip_spots: Set[str],
    r1_left: int,
    r1_right: int,
    r2_left: int,
    r2_right: int,
) -> None:
    pos = pileup_col.pos
    strand_flags = TOP_FLAGS if strand == "+" else BOT_FLAGS
    for pileup_read in pileup_col.pileups:
        record = pileup_read.alignment
        if record.flag not in strand_flags:
            continue
        if has_unclear_conversion_strand(record):
            continue
        if pileup_read.is_del or pileup_read.is_refskip:
            continue
        qpos = pileup_read.query_position
        if qpos is None:
            continue
        seq = record.query_sequence
        if not seq or qpos >= len(seq):
            continue
        if is_trimmed(record, qpos, len(seq), r1_left, r1_right, r2_left, r2_right):
            continue
        if strand == "+":
            if qpos + 1 >= len(seq):
                continue
            neighbor_qpos = qpos + 1
            neighbor = seq[neighbor_qpos]
        else:
            if qpos == 0:
                continue
            neighbor_qpos = qpos - 1
            neighbor = seq[neighbor_qpos]
        if is_trimmed(
            record, neighbor_qpos, len(seq), r1_left, r1_right, r2_left, r2_right
        ):
            continue
        if not ch_neighbor_matches(context, strand, neighbor):
            continue
        observation = classify_ch_observation(seq[qpos].upper(), strand)
        if observation is None:
            continue
        meth, unmeth = observation
        if not host:
            values = dest.setdefault((pos, context, strand), [0, 0])
            values[0] += meth
            values[1] += unmeth
            continue
        try:
            cb = record.get_tag(cb_tag)
        except KeyError:
            continue
        assert barcode_to_index is not None
        spot = spot_id_from_cb(str(cb), barcode_to_index, barcode_len)
        if not spot_is_allowed(spot, allowed_spots) or spot in skip_spots:
            continue
        assert isinstance(spot, str)
        slot = dest.setdefault(spot, {}).setdefault(pos, [0, 0, context, strand])
        slot[0] += meth
        slot[1] += unmeth


def process_interval(
    bam_path: str,
    chrom: str,
    batch_index: int,
    start: int,
    stop: int,
    part_dir: str,
    context_mode: str,
    host: bool,
    cb_tag: str,
    barcode_to_index: Optional[BarcodeMap],
    barcode_len: int,
    allowed_spots: Optional[Set[str]],
    skip_spots: Set[str],
    min_base_quality: int,
    min_mapping_quality: int,
    max_depth: int,
    r1_left: int,
    r1_right: int,
    r2_left: int,
    r2_right: int,
    sample_size: Optional[int],
) -> Tuple[int, Dict[str, str]]:
    sequence = WORKER_CHROMOSOME_SEQUENCE
    if sequence is None:
        raise RuntimeError("worker chromosome sequence was not initialized")
    cg_host: CgCounts = defaultdict(dict)
    ch_host: ChCounts = defaultdict(dict)
    cg_agg: Dict[int, List[int]] = defaultdict(lambda: [0, 0])
    ch_agg: Dict[Tuple[int, str, str], List[int]] = defaultdict(lambda: [0, 0])
    seen_sites = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for pileup_col in bam.pileup(
            chrom,
            start,
            stop,
            truncate=True,
            ignore_overlaps=True,
            min_base_quality=min_base_quality,
            stepper="samtools",
            redo_baq=False,
            max_depth=max_depth,
            ignore_orphans=True,
            compute_baq=True,
            min_mapping_quality=min_mapping_quality,
        ):
            is_cpg, ch_target = classify_reference_position(sequence, pileup_col.pos)
            if sample_size is not None and (is_cpg or ch_target is not None):
                seen_sites += 1
                if seen_sites > sample_size:
                    continue
            if is_cpg and context_mode in {"cg", "both"}:
                count_cg_column(
                    pileup_col,
                    cg_host if host else cg_agg,
                    host,
                    cb_tag,
                    barcode_to_index,
                    barcode_len,
                    allowed_spots,
                    skip_spots,
                    r1_left,
                    r1_right,
                    r2_left,
                    r2_right,
                )
            if ch_target is not None and context_mode in {"ch", "both"}:
                ctx, strand = ch_target
                count_ch_column(
                    pileup_col,
                    ctx,
                    strand,
                    ch_host if host else ch_agg,
                    host,
                    cb_tag,
                    barcode_to_index,
                    barcode_len,
                    allowed_spots,
                    skip_spots,
                    r1_left,
                    r1_right,
                    r2_left,
                    r2_right,
                )
    part_root = Path(part_dir)
    paths: Dict[str, str] = {}
    if context_mode in {"cg", "both"}:
        cg_path = part_root / f"batch_{batch_index:06d}.CG.part"
        with cg_path.open("w", encoding="utf-8") as handle:
            if host:
                for spot_id, pos_map in cg_host.items():
                    for pos, (meth, unmeth) in pos_map.items():
                        if meth + unmeth <= 0:
                            continue
                        handle.write(
                            format_cg_line(chrom, pos, meth, unmeth).rstrip("\n")
                            + f"\t{spot_id}\n"
                        )
            else:
                for pos, (meth, unmeth) in sorted(cg_agg.items()):
                    if meth + unmeth <= 0:
                        continue
                    handle.write(format_cg_line(chrom, pos, meth, unmeth))
        paths["CG"] = str(cg_path)
    if context_mode in {"ch", "both"}:
        ch_path = part_root / f"batch_{batch_index:06d}.CH.part"
        with ch_path.open("w", encoding="utf-8") as handle:
            if host:
                for spot_id, pos_map in ch_host.items():
                    for pos, values in pos_map.items():
                        meth, unmeth, ctx, strand = values
                        if int(meth) + int(unmeth) <= 0:
                            continue
                        handle.write(
                            format_ch_line(
                                chrom, pos, int(meth), int(unmeth), str(ctx), str(strand)
                            ).rstrip("\n")
                            + f"\t{spot_id}\n"
                        )
            else:
                for (pos, ctx, strand), (meth, unmeth) in sorted(ch_agg.items()):
                    if meth + unmeth <= 0:
                        continue
                    handle.write(format_ch_line(chrom, pos, meth, unmeth, ctx, strand))
        paths["CH"] = str(ch_path)
    return batch_index, paths


def append_parts(
    part_paths: Dict[str, str],
    handles: Dict[str, BinaryIO],
) -> None:
    for context, path_str in part_paths.items():
        path = Path(path_str)
        if path.stat().st_size == 0:
            continue
        with path.open("rb") as part_handle:
            shutil.copyfileobj(part_handle, handles[context], length=16 * 1024 * 1024)


def shuffle_host_parts(
    combined_cg: Path,
    combined_ch: Optional[Path],
    out_dir: Path,
    roster: List[str],
    skip_spots: Set[str],
    context_mode: str,
) -> None:
    handles: Dict[Tuple[str, str], BinaryIO] = {}
    temps: Dict[Tuple[str, str], Path] = {}

    def get_handle(spot_id: str, suffix: str) -> BinaryIO:
        key = (spot_id, suffix)
        if key not in handles:
            final_path = host_cov_path(out_dir, spot_id, suffix)
            final_path.parent.mkdir(parents=True, exist_ok=True)
            temp_path = final_path.with_name(final_path.name + f".tmp.{os.getpid()}")
            temps[key] = temp_path
            handles[key] = temp_path.open("wb")
        return handles[key]

    def consume(part_path: Path, suffix: str, extra_cols: int) -> None:
        if not part_path.is_file() or part_path.stat().st_size == 0:
            return
        with part_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                spot_id = fields[-1]
                body = "\t".join(fields[:-1]) + "\n"
                if extra_cols and len(fields) < 7 + extra_cols:
                    continue
                get_handle(spot_id, suffix).write(body.encode("utf-8"))

    try:
        if context_mode in {"cg", "both"}:
            consume(combined_cg, "CG", 0)
        if context_mode in {"ch", "both"} and combined_ch is not None:
            consume(combined_ch, "CH", 1)
    finally:
        for handle in handles.values():
            handle.close()
    for (spot_id, suffix), temp_path in temps.items():
        os.replace(temp_path, host_cov_path(out_dir, spot_id, suffix))
    suffixes = []
    if context_mode in {"cg", "both"}:
        suffixes.append("CG")
    if context_mode in {"ch", "both"}:
        suffixes.append("CH")
    for spot_id in roster:
        if spot_id in skip_spots:
            continue
        for suffix in suffixes:
            path = host_cov_path(out_dir, spot_id, suffix)
            if not path.exists():
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("", encoding="utf-8")


def process_chromosome(
    args: argparse.Namespace,
    chrom: str,
    part_root: Path,
    host: bool,
    barcode_to_index: Optional[BarcodeMap],
    barcode_len: int,
    allowed_spots: Optional[Set[str]],
    skip_spots: Set[str],
    combined_handles: Dict[str, BinaryIO],
) -> None:
    global WORKER_CHROMOSOME_SEQUENCE
    with pysam.FastaFile(args.reference) as fasta:
        sequence = fasta.fetch(chrom).upper()
    WORKER_CHROMOSOME_SEQUENCE = sequence
    chrom_len = len(sequence)
    batch = effective_batch_size(chrom_len, args.batch_size, args.jobs)
    intervals = [
        (index, start, min(start + batch, chrom_len))
        for index, start in enumerate(range(0, chrom_len, batch))
    ]
    chrom_dir = part_root / chrom
    chrom_dir.mkdir(parents=True, exist_ok=True)
    worker_kwargs = dict(
        bam_path=args.bam,
        chrom=chrom,
        part_dir=str(chrom_dir),
        context_mode=args.context_mode,
        host=host,
        cb_tag=args.cb_tag,
        barcode_to_index=barcode_to_index,
        barcode_len=barcode_len,
        allowed_spots=allowed_spots,
        skip_spots=skip_spots,
        min_base_quality=args.min_base_quality,
        min_mapping_quality=args.min_mapping_quality,
        max_depth=args.max_depth,
        r1_left=args.r1_left_trimming,
        r1_right=args.r1_right_trimming,
        r2_left=args.r2_left_trimming,
        r2_right=args.r2_right_trimming,
        sample_size=args.sample_size,
    )
    results = []
    if args.jobs == 1 or len(intervals) == 1:
        for batch_index, start, stop in intervals:
            results.append(
                process_interval(batch_index=batch_index, start=start, stop=stop, **{
                    k: v for k, v in worker_kwargs.items() if k != "part_dir"
                } | {"part_dir": str(chrom_dir)})
            )
    else:
        fork_context = multiprocessing.get_context("fork")
        with ProcessPoolExecutor(
            max_workers=min(args.jobs, len(intervals)),
            mp_context=fork_context,
        ) as executor:
            futures = [
                executor.submit(
                    process_interval,
                    args.bam,
                    chrom,
                    batch_index,
                    start,
                    stop,
                    str(chrom_dir),
                    args.context_mode,
                    host,
                    args.cb_tag,
                    barcode_to_index,
                    barcode_len,
                    allowed_spots,
                    skip_spots,
                    args.min_base_quality,
                    args.min_mapping_quality,
                    args.max_depth,
                    args.r1_left_trimming,
                    args.r1_right_trimming,
                    args.r2_left_trimming,
                    args.r2_right_trimming,
                    args.sample_size,
                )
                for batch_index, start, stop in intervals
            ]
            for future in as_completed(futures):
                results.append(future.result())
    for _, part_paths in sorted(results):
        append_parts(part_paths, combined_handles)
    shutil.rmtree(chrom_dir)
    WORKER_CHROMOSOME_SEQUENCE = None


def existing_skip_spots(
    out_dir: Path, roster: List[str], context_mode: str
) -> Set[str]:
    skip: Set[str] = set()
    suffixes = []
    if context_mode in {"cg", "both"}:
        suffixes.append("CG")
    if context_mode in {"ch", "both"}:
        suffixes.append("CH")
    for spot_id in roster:
        if all(host_cov_path(out_dir, spot_id, suffix).is_file() for suffix in suffixes):
            skip.add(spot_id)
    return skip


def validate_runtime(args: argparse.Namespace) -> None:
    if args.jobs <= 0:
        raise ValueError("--jobs must be > 0")
    if args.batch_size <= 0 or args.max_depth <= 0:
        raise ValueError("--batch-size and --max-depth must be > 0")
    if args.min_base_quality < 0 or args.min_mapping_quality < 0:
        raise ValueError("quality filters must be >= 0")
    for name in (
        "r1_left_trimming",
        "r1_right_trimming",
        "r2_left_trimming",
        "r2_right_trimming",
    ):
        if getattr(args, name) < 0:
            raise ValueError(f"--{name.replace('_', '-')} must be >= 0")
    if args.sample_size is not None and args.sample_size <= 0:
        raise ValueError("--sample-size must be > 0 when provided")
    if args.mode == "host":
        if not args.barcodes or not args.spot_roster or not args.out_dir:
            raise ValueError("host mode requires --barcodes, --spot-roster, and --out-dir")
    elif not args.output:
        raise ValueError("aggregate mode requires --output")


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    print(f"[methy_caller] mode={args.mode}")
    print(f"[methy_caller] bam={args.bam}")
    print(f"[methy_caller] reference={args.reference}")
    print(f"[methy_caller] chromosomes={args.chromosomes}")
    print(f"[methy_caller] context_mode={args.context_mode}")
    print(f"[methy_caller] jobs={args.jobs}")
    if args.mode == "host":
        print(f"[methy_caller] barcodes={args.barcodes}")
        print(f"[methy_caller] spot_roster={args.spot_roster}")
        print(f"[methy_caller] out_dir={args.out_dir}")
    else:
        print(f"[methy_caller] output={args.output}")
    if args.dry_run:
        print("[methy_caller] dry_run=1")
        return 0
    try:
        validate_runtime(args)
        chromosomes = parse_chromosomes(args.chromosomes)
        validate_chromosomes(args.reference, chromosomes)
        bam_path = Path(args.bam)
        if not bam_path.is_file():
            raise ValueError(f"BAM not found: {bam_path}")
        if not bam_index_exists(bam_path):
            raise ValueError(f"BAM index not found: {bam_path}")
        if not Path(args.reference).is_file() or not Path(str(args.reference) + ".fai").is_file():
            raise ValueError(f"reference FASTA/.fai not found: {args.reference}")
    except ValueError as exc:
        print(f"[methy_caller] error: {exc}", file=sys.stderr)
        return 1

    host = args.mode == "host"
    barcode_to_index: Optional[BarcodeMap] = None
    barcode_len = 0
    allowed_spots: Optional[Set[str]] = None
    skip_spots: Set[str] = set()
    roster: List[str] = []
    if host:
        barcodes, barcode_to_index = read_barcodes(args.barcodes)
        if not barcodes:
            print("[methy_caller] error: barcode whitelist is empty", file=sys.stderr)
            return 1
        barcode_len = len(barcodes[0])
        roster = load_spot_roster(args.spot_roster)
        allowed_spots = set(roster)
        skip_spots = existing_skip_spots(Path(args.out_dir), roster, args.context_mode)
        print(f"[methy_caller] roster_spots={len(roster)} skip_spots={len(skip_spots)}")

    suffixes = []
    if args.context_mode in {"cg", "both"}:
        suffixes.append("CG")
    if args.context_mode in {"ch", "both"}:
        suffixes.append("CH")

    if host:
        out_dir = Path(args.out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        combined_dir = out_dir
    else:
        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        combined_dir = out_path.parent
        if all(Path(str(args.output) + f".{suffix}.cov").is_file() for suffix in suffixes):
            print("[methy_caller] skip_existing_aggregate=1")
            return 0

    temp_suffix = f".parts.{os.getpid()}"
    combined_temps = {
        suffix: combined_dir / f".call.{suffix}{temp_suffix}"
        for suffix in suffixes
    }
    handles = {suffix: path.open("wb") for suffix, path in combined_temps.items()}
    try:
        with tempfile.TemporaryDirectory(prefix=".methy-call-parts.", dir=combined_dir) as part_dir:
            for chrom in chromosomes:
                process_chromosome(
                    args,
                    chrom,
                    Path(part_dir),
                    host,
                    barcode_to_index,
                    barcode_len,
                    allowed_spots,
                    skip_spots,
                    handles,
                )
        for handle in handles.values():
            handle.close()
        handles.clear()
        if host:
            shuffle_host_parts(
                combined_temps.get("CG", Path("/dev/null")),
                combined_temps.get("CH"),
                Path(args.out_dir),
                roster,
                skip_spots,
                args.context_mode,
            )
        else:
            for suffix, temp_path in combined_temps.items():
                os.replace(temp_path, Path(str(args.output) + f".{suffix}.cov"))
    finally:
        for handle in handles.values():
            handle.close()
        for temp_path in combined_temps.values():
            try:
                temp_path.unlink()
            except FileNotFoundError:
                pass
    print("[methy_caller] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

Fix the `jobs==1` `process_interval` call if the dict-merge form is awkward: call `process_interval` with explicit positional args matching the parallel path.

**Sample-size note:** the spec is “first N reference sites per chromosome”. The snippet above counts sites inside each interval, which is wrong across intervals. Implement sample-size by pre-scanning the chromosome sequence for CG/CH positions, keeping `positions[:sample_size]`, and skipping pileup columns not in that set. Do not use a per-interval counter.

- [ ] **Step 2: Run unit tests**

```bash
pixi run python -m unittest tests.test_methy_caller -v
```

Expected: all tests `OK`.

- [ ] **Step 3: CLI help**

```bash
pixi run python scripts/methy_caller.py --help
```

Expected: shows `--mode host|aggregate`, `--barcodes`, `--spot-roster`.

- [ ] **Step 4: Commit**

```bash
git add scripts/methy_caller.py
git commit -m "$(cat <<'EOF'
feat: call TAPS methylation from pooled BAM with template rules

BREAKING: host calling reads pooled.byCB.bam; MAPQ default 10;
max_depth default 50000; CH neighbor and strand flags match the
template caller.

EOF
)"
```

---

### Task 4: Rewire `scripts/call.py`

**Files:**
- Modify: `scripts/call.py`

- [ ] **Step 1: Change the CLI**

Remove `--ch-caller-script`. Add:

```python
    parser.add_argument(
        "--barcodes",
        help="Barcode whitelist TSV (required for host calling).",
    )
```

Keep `--caller-script` default `scripts/methy_caller.py`. Change `--min-mapping-quality` default to `10` and `--max-depth` default to `50000`.

- [ ] **Step 2: Replace host/spike command builders**

Delete `discover_host_spot_bams`, `run_host_spot`, and `build_context_caller_command`. Add:

```python
def build_caller_command(
    args: argparse.Namespace,
    *,
    mode: str,
    bam: Path,
    reference: str,
    chromosomes: str,
    output: Path | None = None,
    out_dir: Path | None = None,
    barcodes: Path | None = None,
    spot_roster: Path | None = None,
) -> list[str]:
    command = [
        sys.executable,
        args.caller_script,
        "--mode",
        mode,
        "--bam",
        str(bam),
        "--reference",
        reference,
        "--chromosomes",
        chromosomes,
        "--context-mode",
        args.context_mode,
        "--min-base-quality",
        str(args.min_base_quality),
        "--min-mapping-quality",
        str(args.min_mapping_quality),
        "--max-depth",
        str(args.max_depth),
        "--batch-size",
        str(args.batch_size),
        "--jobs",
        str(args.jobs),
        "--r1-left-trimming",
        str(args.r1_left_trimming),
        "--r1-right-trimming",
        str(args.r1_right_trimming),
        "--r2-left-trimming",
        str(args.r2_left_trimming),
        "--r2-right-trimming",
        str(args.r2_right_trimming),
    ]
    if args.sample_size is not None:
        command.extend(["--sample-size", str(args.sample_size)])
    if mode == "host":
        assert out_dir is not None and barcodes is not None and spot_roster is not None
        command.extend(
            [
                "--out-dir",
                str(out_dir),
                "--barcodes",
                str(barcodes),
                "--spot-roster",
                str(spot_roster),
            ]
        )
    else:
        assert output is not None
        command.extend(["--output", str(output)])
    return command
```

Host run:

```python
def run_host_mode(args: argparse.Namespace, work_path: Path) -> None:
    if not args.barcodes:
        raise ValueError("host calling requires --barcodes")
    pooled = work_path / "pooled" / "pooled.byCB.bam"
    roster = work_path / "split_bams" / "per_spot_read_counts.tsv"
    out_dir = work_path / "coverage" / "host"
    if not pooled.is_file():
        raise ValueError(f"missing pooled host BAM: {pooled}")
    if not roster.is_file():
        raise ValueError(f"missing spot roster: {roster}")
    command = build_caller_command(
        args,
        mode="host",
        bam=pooled,
        reference=args.reference_file,
        chromosomes=args.chromosomes,
        out_dir=out_dir,
        barcodes=Path(args.barcodes),
        spot_roster=roster,
    )
    run_command(command, args.dry_run)
    mito_outputs = run_host_mito_aggregate(args, work_path)
    for mito_out in mito_outputs:
        print(f"[call] output_host_mito={mito_out}")
```

Mito/spike: one aggregate invocation per target with `--output` set to `coverage/host_mito` or `coverage/<spike>` (caller appends `.CG.cov` / `.CH.cov`). Skip if all required suffixes exist.

`run_or_skip` for aggregate: skip when every required `.{CG,CH}.cov` exists.

- [ ] **Step 3: Host dry-run should not require split BAMs**

`--dry-run` still prints the host command using `pooled/pooled.byCB.bam`.

- [ ] **Step 4: Help check**

```bash
pixi run python scripts/call.py --help
```

Expected: `--barcodes` present, `--ch-caller-script` absent.

- [ ] **Step 5: Commit**

```bash
git add scripts/call.py
git commit -m "$(cat <<'EOF'
feat: drive TAPS call from pooled host BAM

BREAKING: call.py host mode no longer reads split_bams/*.sorted.bam.

EOF
)"
```

---

### Task 5: TAPS and TAPS v2 `make_cmd` + input checks

**Files:**
- Modify: `scripts/make_cmd.py`
- Modify: `scripts-v2/make_cmd.py`
- Modify: `scripts/workflow_input_checks.py`

- [ ] **Step 1: Add pooled-BAM host-call checks**

In `scripts/workflow_input_checks.py` add:

```python
def require_bam_index(label: str, bam: Path) -> None:
    require_file(label, bam)
    indexed = (
        Path(str(bam) + ".bai").is_file()
        or bam.with_suffix(".bai").is_file()
        or Path(str(bam) + ".csi").is_file()
    )
    if not indexed:
        raise ValueError(f"missing BAM index ({label}): {bam}")
```

- [ ] **Step 2: `scripts/make_cmd.py` (mirror in `scripts-v2/make_cmd.py`)**

1. Remove `--call-ch-caller-script` argparse, `settings["call_ch_caller_script"]`, and `build_call_command` `--ch-caller-script`.
2. Add to `build_call_command` when `mode in ("all", "host")`:

```python
        "--barcodes",
        args.split_barcodes,
```

3. Change help/defaults:

- `call_min_mapping_quality` default `10` (`else 10`)
- `call_max_depth` default `50000` (`else 50_000`)

4. In `elif stage == "call":` when `call_mode in ("all", "host")` and not `skip_workdir_inputs`:

```python
            wic.require_bam_index(
                "pooled.byCB.bam",
                sample_work / "pooled" / "pooled.byCB.bam",
            )
            wic.require_file(
                "per_spot_read_counts.tsv",
                sample_work / "split_bams" / "per_spot_read_counts.tsv",
            )
            wic.require_file(
                "split_barcodes",
                wic.resolve_config_path(settings["split_barcodes"]),
            )
```

Keep existing `call_reference_file` / spike index checks. Stop requiring `call_ch_caller_script`.

- [ ] **Step 3: Dry-run drivers**

```bash
pixi run python scripts/make_cmd.py --version
pixi run python scripts-v2/make_cmd.py --version
pixi run python scripts/make_cmd.py --workflow-config workflow/dbit_taps_test.json --stage call --runner local --dry-run --skip-workdir-input-checks
pixi run python scripts-v2/make_cmd.py --workflow-config workflow/dbit_taps_v2_test.json --stage call --runner local --dry-run --skip-workdir-input-checks
```

Expected: generated `call.py` lines include `--barcodes` and `scripts/methy_caller.py --mode host`; no `methy_caller_CH.py`.

- [ ] **Step 4: Commit**

```bash
git add scripts/make_cmd.py scripts-v2/make_cmd.py scripts/workflow_input_checks.py
git commit -m "$(cat <<'EOF'
feat: wire pooled TAPS caller through make_cmd

BREAKING: drop call_ch_caller_script; default MAPQ 10 and max_depth 50000.

EOF
)"
```

---

### Task 6: Delete CH script and update workflow JSON

**Files:**
- Delete: `scripts/methy_caller_CH.py`
- Modify: every `workflow/*.json` that contains `call_ch_caller_script`

- [ ] **Step 1: Remove the CH script**

```bash
git rm scripts/methy_caller_CH.py
```

- [ ] **Step 2: Update JSON**

For each matching workflow file:

- delete the `call_ch_caller_script` line
- set `"call_max_depth": 50000` if it is `250`

Required files include:

- `workflow/dbit_taps_test.json`
- `workflow/dbit_taps_slurm.json`
- `workflow/dbit_taps_v2_test.json`
- `workflow/dbit_taps_v2_slurm.json`
- `workflow/L126-brain-S2-TAPS_1M_e2e.json`

and the other `workflow/*TAPS*_e2e.json` hits from `rg call_ch_caller_script workflow`.

- [ ] **Step 3: Confirm no remaining references**

```bash
rg -n "methy_caller_CH|call_ch_caller_script" --glob '!docs/superpowers/**' --glob '!docs/maintenance/changelog.md'
```

Expected: no code/workflow hits. Changelog is updated in Task 7.

- [ ] **Step 4: Commit**

```bash
git add -u scripts/methy_caller_CH.py workflow
git commit -m "$(cat <<'EOF'
chore: remove TAPS CH caller script and config key

EOF
)"
```

---

### Task 7: Documentation

**Files:**
- Modify: `docs/developers/contracts.md`
- Modify: `docs/developers/config-reference.md`
- Modify: `docs/users/taps.md`
- Modify: `docs/maintenance/changelog.md`
- Modify: `docs/maintenance/taps-regression.md`

- [ ] **Step 1: Contracts — replace `call` inputs**

```markdown
Inputs:

- host: `pooled/pooled.byCB.bam` (indexed) plus `split_barcodes` and
  `split_bams/per_spot_read_counts.tsv`
- host mito fallback source: `qc/mbias/host.subsampled.sorted.bam`
- spike-in: `pooled/pooled.<spike_name>.sorted.bam`
```

Add: host calling does not read `split_bams/**/*.sorted.bam`. Spot filenames still `{X_index}/{X_index}_{Y_index}`. Smoke calling uses the split roster.

- [ ] **Step 2: Config reference**

- Remove `call_ch_caller_script` from the field list and notes.
- Document `call_jobs` as interval workers inside a chromosome (default 8).
- Document defaults `call_min_mapping_quality=10`, `call_max_depth=50000`.
- Note host call needs `split_barcodes`.

- [ ] **Step 3: Changelog `## Current` (top)**

```markdown
- **BREAKING TAPS/TAPS v2 call:** host methylation calling reads
  `pooled/pooled.byCB.bam` with interval-parallel pileup instead of
  per-spot split BAMs. CG/CH share one pass. CH neighbor/strand rules
  match the template TAPS caller. Defaults: MAPQ 10, `max_depth` 50000.
  Removed `scripts/methy_caller_CH.py` and `call_ch_caller_script`.
  Per-spot `.cov` paths are unchanged. See [contracts](../developers/contracts.md#call).
```

- [ ] **Step 4: TAPS user guide**

Under Common Mistakes add: host `call` needs an indexed `pooled/pooled.byCB.bam` and `split_bams/per_spot_read_counts.tsv`.

- [ ] **Step 5: Regression doc**

Keep the existing `make_cmd --dry-run` commands. Add:

```bash
pixi run python scripts/methy_caller.py --help
pixi run python scripts/call.py --help
```

State that generated call commands must not reference `methy_caller_CH.py`.

- [ ] **Step 6: Commit**

```bash
git add docs/developers/contracts.md docs/developers/config-reference.md \
        docs/users/taps.md docs/maintenance/changelog.md \
        docs/maintenance/taps-regression.md
git commit -m "$(cat <<'EOF'
docs: record TAPS pooled caller contract and breaking defaults

EOF
)"
```

---

### Task 8: Driver dry-run regression

- [ ] **Step 1: TAPS local + slurm dry-run**

```bash
pixi run python scripts/make_cmd.py --workflow-config workflow/dbit_taps_test.json --stage all --runner local --dry-run
pixi run python scripts/make_cmd.py --workflow-config workflow/dbit_taps_test.json --stage all --runner slurm --dry-run
```

Expected: stage order unchanged; call command uses `methy_caller.py --mode host` and `--barcodes`.

- [ ] **Step 2: TAPS v2 dry-run**

```bash
pixi run python scripts-v2/make_cmd.py --workflow-config workflow/dbit_taps_v2_test.json --stage all --runner local --dry-run
```

Expected: same call wiring; v2 still uses `scripts/call.py`.

- [ ] **Step 3: Unit tests still pass**

```bash
pixi run python -m unittest tests.test_methy_caller -v
```

Expected: `OK`.

---

### Task 9: Sample run on `work/L126-brain-S2-TAPS_1M`

This sample is smoke-split (16 spots), 1M-read TAPS, `call_context_mode` currently `cg`. Use it as the real check. **Do not merge** after this task.

- [ ] **Step 1: Backup old coverage**

```bash
cd /home/jarning/project2026/dbit-matrix
test -d work/L126-brain-S2-TAPS_1M/coverage
mv work/L126-brain-S2-TAPS_1M/coverage \
   work/L126-brain-S2-TAPS_1M/coverage.pre-pooled-caller
mkdir -p work/L126-brain-S2-TAPS_1M/coverage
```

- [ ] **Step 2: Index the pooled host BAM**

```bash
pixi run samtools index work/L126-brain-S2-TAPS_1M/pooled/pooled.byCB.bam
ls work/L126-brain-S2-TAPS_1M/pooled/pooled.byCB.bam.bai
```

Expected: `.bai` exists.

- [ ] **Step 3: Confirm reference FASTA and `.fai`**

```bash
test -f /mnt/e/LiLab_HL/resource/bwa/mm10/genome.fa
test -f /mnt/e/LiLab_HL/resource/bwa/mm10/genome.fa.fai
test -f /mnt/e/LiLab_HL/resource/bwa/lambda/genome.fa.fai
test -f /mnt/e/LiLab_HL/resource/bwa/puc19/genome.fa.fai
```

If mm10 FASTA is missing, stop and report; do not fake a genome.

- [ ] **Step 4: Smoke host call on `chr1` with `--sample-size 5000` and `both` contexts**

```bash
pixi run python scripts/call.py \
  --work-path work/L126-brain-S2-TAPS_1M \
  --mode host \
  --reference-file /mnt/e/LiLab_HL/resource/bwa/mm10/genome.fa \
  --chromosomes chr1 \
  --barcodes configs/barcodes_50a.tsv \
  --context-mode both \
  --jobs 2 \
  --sample-size 5000 \
  --r2-left-trimming 10
```

Expected:

- 16 `coverage/host/<XX>/<XX>_<YY>.CG.cov` matching roster spots in `split_bams/per_spot_read_counts.tsv`
- 16 matching `.CH.cov` files
- no extra host spots
- CG lines have 6 columns; CH lines have 8 (`context` in `CA|CC|CT`, `strand` in `+| -`)
- `host_mito.CG.cov` and `host_mito.CH.cov` exist (mito still runs inside `run_host_mode`; mito uses `--mito-chromosomes` default `chrM` and should **not** inherit `--sample-size` restriction across host nuclear chromosomes — pass mito chromosomes separately in `call.py` as today)

If mito would also get `--sample-size 5000` on `chrM` only, that is acceptable.

- [ ] **Step 5: Spike aggregate call (full small genomes)**

```bash
pixi run python scripts/call.py \
  --work-path work/L126-brain-S2-TAPS_1M \
  --mode spike \
  --reference-file /mnt/e/LiLab_HL/resource/bwa/mm10/genome.fa \
  --chromosomes chr1 \
  --context-mode both \
  --jobs 1 \
  --r1-left-trimming 30 \
  --r2-left-trimming 10 \
  --spike-reference lambda=/mnt/e/LiLab_HL/resource/bwa/lambda/genome.fa \
  --spike-reference puc19=/mnt/e/LiLab_HL/resource/bwa/puc19/genome.fa
```

`--chromosomes chr1` must not be applied to spike; existing `call.py` spike path loads all contigs from each spike FASTA. After Task 4 that behavior must remain. Expected files: `coverage/lambda.CG.cov`, `coverage/lambda.CH.cov`, `coverage/puc19.CG.cov`, `coverage/puc19.CH.cov`.

- [ ] **Step 6: Sanity-check counts**

```bash
wc -l work/L126-brain-S2-TAPS_1M/coverage/host/*/*.CG.cov \
      work/L126-brain-S2-TAPS_1M/coverage/lambda.CG.cov \
      work/L126-brain-S2-TAPS_1M/coverage/puc19.CG.cov
```

Compare roughly to `coverage.pre-pooled-caller` (old caller, CG-only, MAPQ 10, max_depth 250). Numbers will not match exactly. Fail the task if:

- roster spots are missing
- extra spots appear
- CH files have the wrong column count
- lambda/puc19 CG files are empty when the old ones were not

- [ ] **Step 7: Do not merge, do not push unless asked**

Leave `feat/taps-pooled-methy-call` unmerged. Report:

- unit test result
- dry-run result
- L126 file list and `wc -l`
- any remaining issues

Merge to `main` only when the user explicitly asks after this report.

---

## Spec coverage

| Spec item | Task |
|-----------|------|
| TAPS/v2 only; EMSeq untouched | File map / Task 5 |
| Pooled host + interval jobs | Tasks 3–4 |
| Per-spot `.cov` paths unchanged | Tasks 3–4, 9 |
| Template TAPS CG/CH rules, YD, MAPQ 10, max_depth 50000 | Tasks 2–3 |
| One pileup for `both` | Task 3 |
| CB mapping = `split_bams.py` | Tasks 2–3 |
| Smoke roster filter + empty files | Tasks 3, 9 |
| Spike/mito aggregate same classifiers | Tasks 3–4, 9 |
| Delete CH script / config key | Task 6 |
| Docs + BREAKING changelog | Task 7 |
| `--help` / `--dry-run` | Tasks 3, 4, 8 |
| Unit tests | Task 2 |
| L126 sample, new branch, merge later | Tasks 1, 9 |

## Placeholder / type check

- No TBD. `classify_cg_observation`, `spot_id_from_cb`, `ch_neighbor_matches`, `has_unclear_conversion_strand`, `spot_is_allowed` names are shared by tests and Task 3.
- `--output` for aggregate is a stem: `coverage/lambda` → `lambda.CG.cov`. `call.py` must pass that stem, not a `.cov` path.
- `sample_size` is per-chromosome site cap after a full-chromosome site scan, not a per-interval counter.
