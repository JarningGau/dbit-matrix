#!/usr/bin/env python3
from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import hmc_spike as hs


CLAI = "AAGGTGCAATCAGGTGTCAGTCTTAAAGGATCGATAAGGAAGCGGAA"


class HmcMotifMatchTests(unittest.TestCase):
    def test_exact_backbone_with_ct_hit(self) -> None:
        text = "TTT" + CLAI.replace("C", "T") + "GGG"
        regex = hs.compile_hmc_motif_regex(CLAI)
        got = hs.find_hmc_motif(text, CLAI, 0, len(text), regex)
        self.assertEqual(got, (3, "hit"))

    def test_retained_c_sites_hit(self) -> None:
        text = CLAI  # all C sites still C
        got = hs.find_hmc_motif(text, CLAI, 0, len(text))
        self.assertEqual(got, (0, "hit"))

    def test_reject_backbone_mismatch(self) -> None:
        oligo = list(CLAI.replace("C", "T"))
        oligo[0] = "G"  # backbone A→G
        text = "".join(oligo)
        got = hs.find_hmc_motif(text, CLAI, 0, len(text))
        self.assertEqual(got, (None, "none"))

    def test_reject_non_ct_at_sequence_c(self) -> None:
        oligo = list(CLAI.replace("C", "T"))
        oligo[6] = "A"  # sequence C site must be C/T
        text = "".join(oligo)
        got = hs.find_hmc_motif(text, CLAI, 0, len(text))
        self.assertEqual(got, (None, "none"))

    def test_ambiguous_two_exact(self) -> None:
        unit = CLAI.replace("C", "T")
        text = unit + "ACGT" + unit
        got = hs.find_hmc_motif(text, CLAI, 0, len(text))
        self.assertEqual(got, (None, "ambiguous"))

    def test_pattern_uses_ct_class(self) -> None:
        pattern = hs.compile_hmc_motif_regex(CLAI).pattern
        self.assertIn("[CT]", pattern)
        self.assertNotIn(".", pattern.replace("[CT]", ""))


if __name__ == "__main__":
    unittest.main()
