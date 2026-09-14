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
