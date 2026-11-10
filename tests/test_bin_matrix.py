#!/usr/bin/env python3
from __future__ import annotations

import gzip
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from scripts import bin_matrix as bm
import make_cmd as taps_make_cmd


def _write_fai(path: Path, contigs: list[tuple[str, int]]) -> None:
    lines = []
    offset = 0
    for name, length in contigs:
        # name length offset linebases linewidth
        lines.append(f"{name}\t{length}\t{offset}\t60\t61\n")
        offset += length + (length // 60) + 1
    path.write_text("".join(lines), encoding="utf-8")


def _write_cg(path: Path, rows: list[tuple[str, int, int, int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for chrom, pos, meth, unmeth in rows:
            total = meth + unmeth
            pct = round((meth / total) * 100, 2) if total else 0.0
            handle.write(f"{chrom}\t{pos}\t{pos}\t{pct:.2f}\t{meth}\t{unmeth}\n")


def _write_ch(
    path: Path, rows: list[tuple[str, int, int, int, str, str]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for chrom, pos, meth, unmeth, context, strand in rows:
            total = meth + unmeth
            pct = round((meth / total) * 100, 2) if total else 0.0
            handle.write(
                f"{chrom}\t{pos}\t{pos}\t{pct:.2f}\t{meth}\t{unmeth}\t{context}\t{strand}\n"
            )


def _read_mtx(path: Path) -> tuple[int, int, int, list[tuple[int, int, int, int, float]]]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        lines = [ln.strip() for ln in handle if ln.strip() and not ln.startswith("%")]
    nrows, ncols, nnz = (int(x) for x in lines[0].split())
    entries = []
    for line in lines[1:]:
        parts = line.split()
        entries.append(
            (int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3]), float(parts[4]))
        )
    return nrows, ncols, nnz, entries


class BuildBinsTests(unittest.TestCase):
    def test_bins_include_short_last_bin(self) -> None:
        bins, chrom_offsets = bm.build_bins(
            [("chr1", 250)], bin_size=100, chromosomes=None
        )
        self.assertEqual(len(bins), 3)
        self.assertEqual(bins[0], ("chr1:0-100", "chr1", 0, 100))
        self.assertEqual(bins[2], ("chr1:200-250", "chr1", 200, 250))
        self.assertEqual(chrom_offsets["chr1"], 0)

    def test_zero_based_site_maps_to_bin(self) -> None:
        bins, chrom_offsets = bm.build_bins(
            [("chr1", 250)], bin_size=100, chromosomes=None
        )
        self.assertEqual(bm.bin_index_for_site(chrom_offsets, "chr1", 0, 100), 0)
        self.assertEqual(bm.bin_index_for_site(chrom_offsets, "chr1", 99, 100), 0)
        self.assertEqual(bm.bin_index_for_site(chrom_offsets, "chr1", 100, 100), 1)
        self.assertEqual(bm.bin_index_for_site(chrom_offsets, "chr1", 249, 100), 2)
        self.assertEqual(len(bins), 3)


class BinMatrixIntegrationTests(unittest.TestCase):
    def test_cg_and_ca_matrices(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            ref = work / "ref.fa"
            ref.write_text(">chr1\n" + ("N" * 250) + "\n", encoding="utf-8")
            _write_fai(Path(str(ref) + ".fai"), [("chr1", 250)])

            host = work / "coverage" / "host" / "AAA"
            _write_cg(
                host / "AAAABBBB.CG.cov",
                [
                    ("chr1", 10, 3, 1),  # bin 0
                    ("chr1", 150, 1, 1),  # bin 1
                ],
            )
            _write_cg(
                host / "CCCCDDDD.CG.cov",
                [
                    ("chr1", 20, 0, 2),  # bin 0 only
                ],
            )
            _write_ch(
                host / "AAAABBBB.CA.cov",
                [
                    ("chr1", 5, 2, 2, "CA", "+"),
                    ("chr1", 130, 1, 1, "CA", "-"),
                ],
            )
            _write_ch(
                host / "AAAABBBB.CC.cov",
                [
                    ("chr1", 120, 4, 0, "CC", "+"),
                ],
            )
            _write_ch(
                host / "CCCCDDDD.CT.cov",
                [
                    ("chr1", 8, 1, 0, "CT", "+"),
                ],
            )

            rc = bm.run_bin_matrix(
                work_path=work,
                reference=ref,
                bin_size=100,
                contexts=["CG", "CA", "CH"],
                chromosomes=["chr1"],
                barcodes=None,
                out_root="matrix",
                dry_run=False,
            )
            self.assertEqual(rc, 0)

            cg_dir = work / "matrix" / "100bp" / "CG"
            self.assertTrue((cg_dir / "matrix.mtx.gz").is_file())
            with gzip.open(cg_dir / "barcodes.tsv.gz", "rt", encoding="utf-8") as handle:
                barcodes = handle.read().splitlines()
            self.assertEqual(barcodes, ["AAAABBBB", "CCCCDDDD"])
            with gzip.open(cg_dir / "features.tsv.gz", "rt", encoding="utf-8") as handle:
                features = handle.read().splitlines()
            self.assertEqual(features[0], "bin_id\tchrom\tstart\tend")
            self.assertEqual(len(features) - 1, 3)

            nrows, ncols, nnz, entries = _read_mtx(cg_dir / "matrix.mtx.gz")
            self.assertEqual((nrows, ncols), (3, 2))
            # AAAABBBB col1: bins 0 and 1; CCCCDDDD col2: bin 0 only → 3 nonzeros
            self.assertEqual(nnz, 3)
            self.assertEqual(len(entries), 3)
            by_rc = {(r, c): (m, cov, f) for r, c, m, cov, f in entries}
            self.assertEqual(by_rc[(1, 1)], (3, 4, 0.75))
            self.assertEqual(by_rc[(2, 1)], (1, 2, 0.5))
            self.assertEqual(by_rc[(1, 2)], (0, 2, 0.0))
            self.assertNotIn((2, 2), by_rc)  # uncovered omitted
            self.assertNotIn((3, 1), by_rc)

            ca_dir = work / "matrix" / "100bp" / "CA"
            _, _, ca_nnz, ca_entries = _read_mtx(ca_dir / "matrix.mtx.gz")
            # AAAABBBB: CA at bin0 (2+2) and bin1 (1+1); CCCCDDDD: no CA
            self.assertEqual(ca_nnz, 2)
            ca_by = {(r, c): (m, cov, f) for r, c, m, cov, f in ca_entries}
            self.assertEqual(ca_by[(1, 1)], (2, 4, 0.5))
            self.assertEqual(ca_by[(2, 1)], (1, 2, 0.5))

            ch_dir = work / "matrix" / "100bp" / "CH"
            _, _, ch_nnz, ch_entries = _read_mtx(ch_dir / "matrix.mtx.gz")
            # AAAABBBB bin0 CA; bin1 CC+CA; CCCCDDDD bin0 CT
            self.assertEqual(ch_nnz, 3)
            ch_by = {(r, c): (m, cov, f) for r, c, m, cov, f in ch_entries}
            self.assertEqual(ch_by[(1, 1)][:2], (2, 4))
            self.assertAlmostEqual(ch_by[(1, 1)][2], 0.5, places=5)
            self.assertEqual(ch_by[(2, 1)][:2], (5, 6))
            self.assertAlmostEqual(ch_by[(2, 1)][2], 5 / 6, places=5)
            self.assertEqual(ch_by[(1, 2)][:2], (1, 1))
            self.assertAlmostEqual(ch_by[(1, 2)][2], 1.0, places=5)

    def test_dry_run_writes_nothing(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            ref = work / "ref.fa"
            ref.write_text(">chr1\n" + ("N" * 100) + "\n", encoding="utf-8")
            _write_fai(Path(str(ref) + ".fai"), [("chr1", 100)])
            host = work / "coverage" / "host" / "AAA"
            _write_cg(host / "AAAABBBB.CG.cov", [("chr1", 10, 1, 1)])
            rc = bm.run_bin_matrix(
                work_path=work,
                reference=ref,
                bin_size=100,
                contexts=["CG"],
                chromosomes=None,
                barcodes=None,
                out_root="matrix",
                dry_run=True,
            )
            self.assertEqual(rc, 0)
            self.assertFalse((work / "matrix").exists())


class FormatBinSizeLabelTests(unittest.TestCase):
    def test_compact_units(self) -> None:
        self.assertEqual(bm.format_bin_size_label(500), "500bp")
        self.assertEqual(bm.format_bin_size_label(100), "100bp")
        self.assertEqual(bm.format_bin_size_label(999), "999bp")
        self.assertEqual(bm.format_bin_size_label(1_000), "1k")
        self.assertEqual(bm.format_bin_size_label(10_000), "10k")
        self.assertEqual(bm.format_bin_size_label(100_000), "100k")
        self.assertEqual(bm.format_bin_size_label(1_500), "1.5k")
        self.assertEqual(bm.format_bin_size_label(1_000_000), "1M")
        self.assertEqual(bm.format_bin_size_label(2_500_000), "2.5M")

    def test_rejects_non_positive(self) -> None:
        with self.assertRaises(ValueError):
            bm.format_bin_size_label(0)


class AdaptWhitelistTests(unittest.TestCase):
    def test_expand_short_barcodes_to_concat_cb(self) -> None:
        wl = {"AAAA", "TTTT"}
        barcodes = {"AAAATTTT", "TTTTAAAA", "GGGGCCCC"}
        adapted = bm.adapt_whitelist_to_barcodes(wl, barcodes)
        self.assertEqual(adapted, {"AAAAAAAA", "AAAATTTT", "TTTTAAAA", "TTTTTTTT"})
        self.assertEqual(barcodes & adapted, {"AAAATTTT", "TTTTAAAA"})

    def test_keep_direct_match_whitelist(self) -> None:
        wl = {"AAAATTTT", "OTHER"}
        barcodes = {"AAAATTTT", "TTTTAAAA"}
        adapted = bm.adapt_whitelist_to_barcodes(wl, barcodes)
        self.assertIs(adapted, wl)


class BinMatrixDriverConfigTests(unittest.TestCase):
    def test_validate_defaults(self) -> None:
        cfg = taps_make_cmd.validate_bin_matrix_config({})
        self.assertFalse(cfg["enabled"])
        self.assertEqual(cfg["bin_size"], 100_000)
        self.assertEqual(cfg["contexts"], "CG,CA")

    def test_validate_rejects_bad_bin_size(self) -> None:
        with self.assertRaises(ValueError):
            taps_make_cmd.validate_bin_matrix_config({"bin_size": 0})

    def test_validate_rejects_bad_contexts(self) -> None:
        with self.assertRaises(ValueError):
            taps_make_cmd.validate_bin_matrix_config({"contexts": ["CG", "XX"]})

    def test_stages_for_all_include_when_enabled(self) -> None:
        stages = taps_make_cmd.stages_for_all(
            {"bin_matrix": {"enabled": True}, "hmc_spike": {"enabled": False}}
        )
        self.assertIn("bin_matrix", stages)
        call_i = stages.index("call")
        bin_i = stages.index("bin_matrix")
        sat_i = stages.index("saturation")
        self.assertEqual(call_i + 1, bin_i)
        self.assertEqual(bin_i + 1, sat_i)
        self.assertNotIn("hmc_spike", stages)

    def test_stages_for_all_skip_when_disabled(self) -> None:
        stages = taps_make_cmd.stages_for_all(
            {"bin_matrix": {"enabled": False}, "hmc_spike": {"enabled": False}}
        )
        self.assertNotIn("bin_matrix", stages)
        self.assertEqual(
            stages.index("saturation"),
            stages.index("call") + 1,
        )

    def test_slurm_deps_call_only(self) -> None:
        self.assertEqual(
            taps_make_cmd.STAGE_SLURM_DEPENDENCIES["bin_matrix"], ["call"]
        )
        self.assertEqual(
            taps_make_cmd.STAGE_SLURM_DEPENDENCIES["saturation"], ["call"]
        )
        self.assertNotIn(
            "bin_matrix", taps_make_cmd.STAGE_SLURM_DEPENDENCIES["summary"]
        )


if __name__ == "__main__":
    unittest.main()
