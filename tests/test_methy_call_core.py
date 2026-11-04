#!/usr/bin/env python3
from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

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
            core.write_counts_tsv(
                path,
                {
                    "CCCCAAAA": 70,
                    "GGGGAAAA": 550,
                    "AAAAAAAA": 3,
                },
            )
            text = path.read_text(encoding="utf-8")
            self.assertEqual(text.splitlines()[0], "barcode\treads")
            rows = core.load_counts_tsv(path)
            self.assertEqual(
                [(row.barcode, row.reads) for row in rows],
                [("AAAAAAAA", 3), ("CCCCAAAA", 70), ("GGGGAAAA", 550)],
            )


class SmokeWriteSetTests(unittest.TestCase):
    def test_top_16_by_reads_then_spot(self) -> None:
        counts = {f"B{i:02d}A{i:02d}": i for i in range(20)}
        write_set = core.resolve_write_set(counts, smoke=True, roster_ids=None)
        self.assertEqual(len(write_set), 16)
        self.assertEqual(write_set[0], "B19A19")
        self.assertEqual(write_set[-1], "B04A04")
        self.assertNotIn("B03A03", write_set)

    def test_roster_overrides_smoke(self) -> None:
        counts = {"AAAAAAAA": 100, "CCCCAAAA": 1}
        write_set = core.resolve_write_set(
            counts, smoke=True, roster_ids=["CCCCAAAA"]
        )
        self.assertEqual(write_set, ["CCCCAAAA"])

    def test_full_set_without_smoke(self) -> None:
        counts = {"GGGGAAAA": 5, "AAAAAAAA": 1}
        self.assertEqual(
            core.resolve_write_set(counts, smoke=False, roster_ids=None),
            ["AAAAAAAA", "GGGGAAAA"],
        )


class RosterParseTests(unittest.TestCase):
    def test_counts_tsv_barcode_column(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "roster.tsv"
            path.write_text(
                "barcode\treads\nCCCCAAAA\t70\n",
                encoding="utf-8",
            )
            self.assertEqual(core.load_roster_ids(path), ["CCCCAAAA"])

    def test_plain_barcode_lines(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "spots.txt"
            path.write_text("CCCCAAAA\nGGGGAAAA\n", encoding="utf-8")
            self.assertEqual(core.load_roster_ids(path), ["CCCCAAAA", "GGGGAAAA"])


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
                FakeCbRecord("CCCCAAAA"),
                FakeCbRecord("CCCCAAAA"),
                FakeCbRecord("TTTTAAAA"),
                FakeCbRecord(None),
            ]
            counts = core.count_records_by_spot(
                records, "CB", barcode_to_index, barcode_len
            )
            self.assertEqual(counts, {"CCCCAAAA": 2})


class HostCovPathTests(unittest.TestCase):
    def test_prefix_subdir_from_cb(self) -> None:
        path = core.host_cov_path(Path("/tmp/host"), "CCCCAAAA", "CG")
        self.assertEqual(path, Path("/tmp/host/CCC/CCCCAAAA.CG.cov"))


if __name__ == "__main__":
    unittest.main()
