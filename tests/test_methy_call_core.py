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


class FakeFunnelRecord:
    def __init__(
        self,
        *,
        cb: str | None = "CCCCAAAA",
        is_unmapped: bool = False,
        mapping_quality: int = 60,
        nh: int | None = 1,
        flag: int = 99,
        reference_name: str = "chr1",
    ) -> None:
        self._cb = cb
        self.is_unmapped = is_unmapped
        self.mapping_quality = mapping_quality
        self._nh = nh
        self.flag = flag
        self.reference_name = reference_name

    def has_tag(self, name: str) -> bool:
        if name == "CB":
            return self._cb is not None
        if name == "NH":
            return self._nh is not None
        return False

    def get_tag(self, name: str):
        if name == "CB":
            if self._cb is None:
                raise KeyError(name)
            return self._cb
        if name == "NH":
            if self._nh is None:
                raise KeyError(name)
            return self._nh
        raise KeyError(name)


class HostBamFunnelTests(unittest.TestCase):
    def test_scan_counts_and_funnel_same_pass(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "barcodes.tsv"
            path.write_text("AAAA\nCCCC\n", encoding="utf-8")
            _, barcode_to_index, barcode_len = core.load_barcode_whitelist(str(path))
            mito = frozenset({"chrM"})
            records = [
                FakeFunnelRecord(flag=99, reference_name="chr1"),
                FakeFunnelRecord(flag=99, reference_name="chrM"),
                FakeFunnelRecord(flag=0, reference_name="chr1"),  # mapped, not valid
                FakeFunnelRecord(mapping_quality=5),  # excluded from mapped
                FakeFunnelRecord(nh=2),  # excluded from mapped
                FakeFunnelRecord(is_unmapped=True),  # still counted by CB
                FakeFunnelRecord(cb="TTTTAAAA"),  # funnel only (not whitelist)
                FakeFunnelRecord(cb=None, flag=147, reference_name="chr2"),
            ]
            counts, funnel = core.scan_host_bam_for_counts_and_funnel(
                records,
                "CB",
                barcode_to_index,
                barcode_len,
                mito,
            )
            # Six records carry whitelist CB CCCCCAAA (including low-MAPQ / NH / unmapped).
            self.assertEqual(counts, {"CCCCAAAA": 6})
            # Mapped: two valid host + one invalid-flag + off-whitelist + no-CB.
            self.assertEqual(funnel["host_mapped_reads"], 5)
            self.assertEqual(funnel["host_valid_reads"], 4)
            self.assertEqual(funnel["host_valid_reads_mito"], 1)
            self.assertEqual(funnel["host_valid_reads_nuclear"], 3)
            self.assertEqual(funnel["mito_chromosomes"], ["chrM"])

    def test_write_and_load_funnel_json(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "host_bam_funnel.json"
            payload = {
                "host_mapped_reads": 10,
                "host_valid_reads": 8,
                "host_valid_reads_mito": 1,
                "host_valid_reads_nuclear": 7,
                "mito_chromosomes": ["chrM"],
            }
            core.write_host_bam_funnel_json(path, payload)
            loaded = core.load_host_bam_funnel_json(path)
            self.assertEqual(loaded, payload)


class PerSpotMethStatsTests(unittest.TestCase):
    def test_accumulate_cg_and_ch_contexts(self) -> None:
        acc = core.new_meth_stats_accumulator()
        core.accumulate_cov_line_for_meth_stats(
            acc,
            "CG",
            ["chr1", "10", "10", "80.00", "8", "2", "CCCCAAAA"],
        )
        core.accumulate_cov_line_for_meth_stats(
            acc,
            "CG",
            ["chr1", "20", "20", "60.00", "3", "2", "CCCCAAAA"],
        )
        core.accumulate_cov_line_for_meth_stats(
            acc,
            "CH",
            ["chr1", "30", "30", "10.00", "1", "9", "CA", "+", "CCCCAAAA"],
        )
        core.accumulate_cov_line_for_meth_stats(
            acc,
            "CH",
            ["chr1", "40", "40", "20.00", "2", "8", "CC", "+", "CCCCAAAA"],
        )
        core.accumulate_cov_line_for_meth_stats(
            acc,
            "CH",
            ["chr1", "50", "50", "0.00", "0", "5", "XX", "+", "CCCCAAAA"],
        )
        rows = core.finalize_meth_stats_rows(acc, ["CCCCAAAA", "GGGGAAAA"])
        by_bc = {row["barcode"]: row for row in rows}
        self.assertEqual(by_bc["CCCCAAAA"]["mean_5mCG"], "70.00")
        self.assertEqual(by_bc["CCCCAAAA"]["cpg_site_count"], "2")
        self.assertEqual(by_bc["CCCCAAAA"]["mean_5mCH"], "10.00")
        self.assertEqual(by_bc["CCCCAAAA"]["ch_site_count"], "3")
        self.assertEqual(by_bc["CCCCAAAA"]["mean_5mCA"], "10.00")
        self.assertEqual(by_bc["CCCCAAAA"]["ca_site_count"], "1")
        self.assertEqual(by_bc["CCCCAAAA"]["mean_5mCC"], "20.00")
        self.assertEqual(by_bc["CCCCAAAA"]["cc_site_count"], "1")
        self.assertEqual(by_bc["CCCCAAAA"]["mean_5mCT"], "NA")
        self.assertEqual(by_bc["CCCCAAAA"]["ct_site_count"], "NA")
        self.assertEqual(by_bc["GGGGAAAA"]["mean_5mCG"], "NA")
        self.assertEqual(by_bc["GGGGAAAA"]["cpg_site_count"], "NA")

    def test_write_and_load_meth_stats_tsv(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "per_spot_meth_stats.tsv"
            rows = [
                {
                    "barcode": "AAAAAAAA",
                    "mean_5mCG": "50.00",
                    "mean_5mCH": "NA",
                    "mean_5mCA": "NA",
                    "mean_5mCC": "NA",
                    "mean_5mCT": "NA",
                    "cpg_site_count": "1",
                    "ch_site_count": "NA",
                    "ca_site_count": "NA",
                    "cc_site_count": "NA",
                    "ct_site_count": "NA",
                }
            ]
            core.write_per_spot_meth_stats_tsv(path, rows)
            loaded = core.load_per_spot_meth_stats_tsv(path)
            self.assertEqual(loaded["AAAAAAAA"]["mean_5mCG"], "50.00")


if __name__ == "__main__":
    unittest.main()
