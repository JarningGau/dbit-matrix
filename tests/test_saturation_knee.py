#!/usr/bin/env python3
from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import saturation


class AutoReadsThresholdTests(unittest.TestCase):
    def test_rejects_too_few_points_after_tail_filter(self) -> None:
        with self.assertRaises(ValueError):
            saturation.auto_reads_threshold([5.0, 3.0, 1.0])
        with self.assertRaises(ValueError):
            saturation.auto_reads_threshold([20.0])

    def test_steepest_slope_on_synthetic_cliff(self) -> None:
        # Plateau of high counts, then a sharp drop, then a low tail (<=10 dropped).
        plateau = [1000.0 - i * 2.0 for i in range(80)]
        cliff = [800.0 * math.exp(-0.08 * i) + 15.0 for i in range(60)]
        tail = [8.0, 6.0, 4.0, 2.0]
        reads = plateau + cliff + tail
        threshold = saturation.auto_reads_threshold(reads)
        # Steepest drop should land on the cliff, not the plateau head or <=10 tail.
        self.assertGreater(threshold, 10.0)
        self.assertLess(threshold, plateau[-1])
        self.assertGreaterEqual(threshold, min(cliff))
        self.assertLessEqual(threshold, max(cliff))


if __name__ == "__main__":
    unittest.main()
