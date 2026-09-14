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
