#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
for path in (ROOT, SCRIPTS):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import summary as summary_mod  # noqa: E402


class SummarySidecarTests(unittest.TestCase):
    def test_per_spot_from_meth_stats_joins_reads_and_positions(self) -> None:
        meth_rows = {
            "CCCCAAAA": {
                "barcode": "CCCCAAAA",
                "mean_5mCG": "70.00",
                "mean_5mCH": "10.00",
                "mean_5mCA": "10.00",
                "mean_5mCC": "NA",
                "mean_5mCT": "NA",
                "cpg_site_count": "2",
                "ch_site_count": "1",
                "ca_site_count": "1",
                "cc_site_count": "NA",
                "ct_site_count": "NA",
            }
        }
        read_counts = {"CCCCAAAA": "100", "GGGGAAAA": "5"}
        positions = {"CCCCAAAA": ("0", "1"), "GGGGAAAA": ("2", "3")}
        rows = summary_mod.per_spot_rows_from_meth_stats(
            meth_rows, read_counts, positions
        )
        by_bc = {row["barcode"]: row for row in rows}
        self.assertEqual(by_bc["CCCCAAAA"]["mean_5mCG"], "70.00")
        self.assertEqual(by_bc["CCCCAAAA"]["reads"], "100")
        self.assertEqual(by_bc["CCCCAAAA"]["array_row"], "0")
        self.assertEqual(by_bc["CCCCAAAA"]["array_col"], "1")
        self.assertEqual(by_bc["GGGGAAAA"]["mean_5mCG"], "NA")
        self.assertEqual(by_bc["GGGGAAAA"]["reads"], "5")

    def test_load_sidecars_prefer_or_none(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            coverage = Path(tmp) / "coverage"
            coverage.mkdir()
            funnel_path = coverage / "host_bam_funnel.json"
            meth_path = coverage / "per_spot_meth_stats.tsv"
            self.assertIsNone(
                summary_mod.load_summary_sidecars(funnel_path, meth_path)
            )
            funnel_path.write_text(
                json.dumps(
                    {
                        "host_mapped_reads": 10,
                        "host_valid_reads": 8,
                        "host_valid_reads_mito": 1,
                        "host_valid_reads_nuclear": 7,
                        "mito_chromosomes": ["chrM"],
                    }
                ),
                encoding="utf-8",
            )
            meth_path.write_text(
                "barcode\tmean_5mCG\tmean_5mCH\tmean_5mCA\tmean_5mCC\tmean_5mCT\t"
                "cpg_site_count\tch_site_count\tca_site_count\tcc_site_count\tct_site_count\n"
                "AAAAAAAA\t50.00\tNA\tNA\tNA\tNA\t1\tNA\tNA\tNA\tNA\n",
                encoding="utf-8",
            )
            sidecars = summary_mod.load_summary_sidecars(funnel_path, meth_path)
            assert sidecars is not None
            funnel, meth_rows = sidecars
            self.assertEqual(funnel["host_mapped_reads"], 10)
            self.assertEqual(meth_rows["AAAAAAAA"]["mean_5mCG"], "50.00")


if __name__ == "__main__":
    unittest.main()
