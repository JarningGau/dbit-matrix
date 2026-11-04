#!/usr/bin/env python3
from __future__ import annotations

import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
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
    def test_concat_bcb_bca_spot_id(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "barcodes.tsv"
            path.write_text("AAAA\nCCCC\nGGGG\n", encoding="utf-8")
            _, barcode_to_index = read_barcodes(str(path))
            barcode_len = 4
            # CB = <bcB><bcA>; parse returns (bcA/X, bcB/Y)
            self.assertEqual(parse_cb_value("CCCCAAAA", barcode_len), ("AAAA", "CCCC"))
            self.assertIsNone(parse_cb_value("AAAA+CCCC", barcode_len))
            self.assertEqual(
                mc.spot_id_from_cb("CCCCAAAA", barcode_to_index, barcode_len),
                "CCCCAAAA",
            )
            self.assertEqual(
                mc.spot_id_from_cb("AAAACCCC", barcode_to_index, barcode_len),
                "AAAACCCC",
            )
            self.assertIsNone(
                mc.spot_id_from_cb("TTTTAAAA", barcode_to_index, barcode_len)
            )

    def test_roster_filters_unlisted_spots(self) -> None:
        allowed = {"CCCCAAAA"}
        self.assertTrue(mc.spot_is_allowed("CCCCAAAA", allowed))
        self.assertFalse(mc.spot_is_allowed("GGGGAAAA", allowed))
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

    def test_host_dry_run_prints_scan_keys_without_writing(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out_dir = root / "coverage" / "host"
            expected_counts = out_dir.resolve().parent / "per_spot_read_counts.tsv"
            argv = [
                "--mode",
                "host",
                "--bam",
                str(root / "missing.bam"),
                "--reference",
                str(root / "missing.fa"),
                "--chromosomes",
                "chr1",
                "--barcodes",
                str(root / "barcodes.tsv"),
                "--out-dir",
                str(out_dir),
                "--dry-run",
            ]
            stdout = StringIO()
            with redirect_stdout(stdout):
                self.assertEqual(mc.main(argv), 0)
            text = stdout.getvalue()
            self.assertIn("smoke=", text)
            self.assertIn("counts_tsv=", text)
            self.assertIn("roster_spots=", text)
            self.assertIn(str(expected_counts), text)
            self.assertFalse(out_dir.exists())
            self.assertFalse(expected_counts.exists())


class SamplePositionTests(unittest.TestCase):
    SEQ = "ACGCATGCT"

    def test_cg_first_site(self) -> None:
        self.assertEqual(mc.select_sample_positions(self.SEQ, "cg", 1), {1})

    def test_ch_first_two_sites(self) -> None:
        self.assertEqual(mc.select_sample_positions(self.SEQ, "ch", 2), {3, 6})

    def test_both_first_three_sites(self) -> None:
        self.assertEqual(mc.select_sample_positions(self.SEQ, "both", 3), {1, 3, 6})

    def test_none_disables_filter(self) -> None:
        self.assertIsNone(mc.select_sample_positions(self.SEQ, "cg", None))


class BamIndexTests(unittest.TestCase):
    def test_bam_bai_exists(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            bam = Path(tmp) / "pooled.byCB.bam"
            bam.write_bytes(b"")
            self.assertFalse(mc.bam_index_exists(bam))
            Path(str(bam) + ".bai").write_bytes(b"")
            self.assertTrue(mc.bam_index_exists(bam))

    def test_missing_index(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            bam = Path(tmp) / "pooled.byCB.bam"
            bam.write_bytes(b"")
            self.assertFalse(mc.bam_index_exists(bam))


class BarcodeWhitelistTests(unittest.TestCase):
    def test_unequal_lengths_raise(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "barcodes.tsv"
            path.write_text("AAAA\nCCC\n", encoding="utf-8")
            with self.assertRaises(ValueError) as ctx:
                mc.load_barcode_whitelist(str(path))
            self.assertIn("unequal", str(ctx.exception))

    def test_missing_file_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "missing_barcodes.tsv"
            with self.assertRaises(ValueError) as ctx:
                mc.load_barcode_whitelist(str(path))
            self.assertIn("not found", str(ctx.exception))


class MissingHostInputsTests(unittest.TestCase):
    def test_missing_barcodes_exits_1_without_traceback(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            argv = [
                "--mode",
                "host",
                "--bam",
                str(root / "missing.bam"),
                "--reference",
                str(root / "missing.fa"),
                "--chromosomes",
                "chr1",
                "--barcodes",
                str(root / "missing_barcodes.tsv"),
                "--out-dir",
                str(root / "out"),
            ]
            stderr = StringIO()
            with redirect_stderr(stderr):
                code = mc.main(argv)
            text = stderr.getvalue()
            self.assertEqual(code, 1)
            self.assertIn("[methy_caller] error:", text)
            self.assertIn("barcode", text.lower())
            self.assertNotIn("Traceback", text)


class ShuffleHostPartsTests(unittest.TestCase):
    def test_streams_roster_spots_and_ignores_unknown(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out_dir = root / "out"
            combined_cg = root / "combined.CG"
            combined_cg.write_text(
                "chr1\t11\t11\t50.00\t1\t1\tAAAAAAAA\n"
                "chr1\t10\t10\t100.00\t1\t0\tCCCCAAAA\n"
                "chr1\t20\t20\t0.00\t0\t1\tZZZZZZZZ\n"
                "badline\n"
                "chr1\t12\t12\t0.00\t0\t2\tAAAAAAAA\n",
                encoding="utf-8",
            )
            roster = ["AAAAAAAA", "CCCCAAAA", "GGGGAAAA"]
            mc.shuffle_host_parts(
                combined_cg, None, out_dir, roster, set(), "cg"
            )
            path_aa = out_dir / "AAA" / "AAAAAAAA.CG.cov"
            path_ca = out_dir / "CCC" / "CCCCAAAA.CG.cov"
            path_ga = out_dir / "GGG" / "GGGGAAAA.CG.cov"
            self.assertTrue(path_aa.is_file())
            self.assertTrue(path_ca.is_file())
            self.assertTrue(path_ga.is_file())
            self.assertEqual(path_ga.read_text(encoding="utf-8"), "")
            text_aa = path_aa.read_text(encoding="utf-8")
            self.assertEqual(
                text_aa,
                "chr1\t11\t11\t50.00\t1\t1\nchr1\t12\t12\t0.00\t0\t2\n",
            )
            self.assertNotIn("AAAAAAAA", text_aa)
            self.assertEqual(
                path_ca.read_text(encoding="utf-8"),
                "chr1\t10\t10\t100.00\t1\t0\n",
            )
            self.assertFalse((out_dir / "ZZZ").exists())

    def test_skip_spots_are_not_overwritten(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out_dir = root / "out"
            existing = out_dir / "AAA" / "AAAAAAAA.CG.cov"
            existing.parent.mkdir(parents=True, exist_ok=True)
            existing.write_text("KEEP\n", encoding="utf-8")
            combined_cg = root / "combined.CG"
            combined_cg.write_text(
                "chr1\t10\t10\t100.00\t1\t0\tAAAAAAAA\n"
                "chr1\t11\t11\t0.00\t0\t1\tCCCCAAAA\n",
                encoding="utf-8",
            )
            mc.shuffle_host_parts(
                combined_cg,
                None,
                out_dir,
                ["AAAAAAAA", "CCCCAAAA"],
                {"AAAAAAAA"},
                "cg",
            )
            self.assertEqual(existing.read_text(encoding="utf-8"), "KEEP\n")
            self.assertEqual(
                (out_dir / "CCC" / "CCCCAAAA.CG.cov").read_text(encoding="utf-8"),
                "chr1\t11\t11\t0.00\t0\t1\n",
            )

    def test_ch_requires_nine_fields(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out_dir = root / "out"
            combined_ch = root / "combined.CH"
            combined_ch.write_text(
                "chr1\t10\t10\t100.00\t1\t0\tCA\t+\tAAAAAAAA\n"
                "chr1\t11\t11\t0.00\t0\t1\tAAAAAAAA\n",
                encoding="utf-8",
            )
            mc.shuffle_host_parts(
                Path("/dev/null"),
                combined_ch,
                out_dir,
                ["AAAAAAAA"],
                set(),
                "ch",
            )
            text = (out_dir / "AAA" / "AAAAAAAA.CH.cov").read_text(encoding="utf-8")
            self.assertEqual(text, "chr1\t10\t10\t100.00\t1\t0\tCA\t+\n")
            self.assertNotIn("00_00", text)


if __name__ == "__main__":
    unittest.main()
