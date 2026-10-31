"""Closed-form checks on reporting metrics, separate from descriptor chemistry."""

import importlib.util
import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
HAS_REFERENCE = importlib.util.find_spec("morfeus") is not None
if HAS_REFERENCE:
    from study_conformer_ensembles import metrics


@unittest.skipUnless(
    HAS_REFERENCE, "requires the scientific Morfeus reference environment"
)
class ConformerStudyMetrics(unittest.TestCase):
    def test_closed_form_residuals_not_squared_correlation(self):
        rows = [{"reference": x, "native": 2 * x + 1} for x in [1, 2, 3]]
        rows.append({"reference": None, "native": 99})
        result = metrics(rows)
        self.assertEqual(result["n"], 3)
        self.assertEqual(result["missing_pairs"], 1)
        self.assertEqual(result["mae"], 3)
        self.assertAlmostEqual(result["rmse"], math.sqrt(29 / 3))
        self.assertEqual(result["median_absolute_error"], 3)
        self.assertEqual(result["maximum_absolute_error"], 4)
        self.assertEqual(result["r_squared"], -13.5)
        self.assertEqual(result["slope"], 2)
        self.assertEqual(result["intercept"], 1)

    def test_constant_and_missing_reference_are_not_perfect_agreement(self):
        result = metrics([{"reference": 1, "native": y} for y in [1, 2]])
        self.assertIsNone(result["r_squared"])
        self.assertIsNone(result["slope"])
        self.assertIsNone(result["intercept"])
        self.assertEqual(result["mae"], 0.5)
        self.assertEqual(
            metrics([{"reference": None, "native": None}]), {"n": 0, "missing_pairs": 1}
        )


if __name__ == "__main__":
    unittest.main()
