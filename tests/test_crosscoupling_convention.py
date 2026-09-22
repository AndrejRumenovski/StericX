"""Primary-table boundary witnesses must use the official inclusive labels."""

import csv
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "studies"))
try:
    import study_007_crosscoupling as nickel
    import study_009_pd_crosscoupling as palladium
except ImportError as error:
    raise unittest.SkipTest("science extra required for classifier checks") from error


class ThresholdConventionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        result = ROOT / "docs/scientific_accuracy_audit/models/results"
        with (result / "crosscoupling_complete_inputs.csv").open() as stream:
            cls.rows = list(csv.DictReader(stream))
        cls.expected = json.loads(
            (result / "crosscoupling_reference_models.json").read_text()
        )

    def test_all_seven_frozen_boundary_labels_are_active(self):
        boundaries = set()
        for module in (nickel, palladium):
            for reaction, paper in module.PAPER.items():
                source = [row for row in self.rows if row["reaction"] == reaction]
                rows = [
                    {"id": int(r["id"]), "yield": float(r["yield"])} for r in source
                ]
                descriptors = {r["id"]: float(r["id"]) for r in rows}
                _, labels = module.reaction_arrays(
                    {reaction: rows}, descriptors, reaction
                )
                for row, label in zip(rows, labels, strict=True):
                    self.assertEqual(label, int(row["yield"] >= paper["y_cut"]))
                    if row["yield"] == paper["y_cut"]:
                        boundaries.add((reaction, row["id"]))
        self.assertEqual(
            boundaries,
            {
                ("I", 162),
                ("II", 88),
                ("II", 179),
                ("IV", 310),
                ("V", 88),
                ("VII", 84),
                ("VIII", 11),
            },
        )

    def test_corrected_classifiers_match_immutable_inclusive_reference_fits(self):
        # The same frozen descriptors isolate the label correction from geometry
        # corrections. These are resubstitution fits, not external validation.
        for module in (nickel, palladium):
            fit = (
                module.single_node_threshold
                if module is nickel
                else module.fit_single_node
            )
            for reaction in module.PAPER:
                source = [row for row in self.rows if row["reaction"] == reaction]
                rows = [
                    {"id": int(r["id"]), "yield": float(r["yield"])} for r in source
                ]
                descriptors = {
                    int(r["id"]): float(r["native_percent_buried_volume_min"])
                    for r in source
                }
                x, y = module.reaction_arrays({reaction: rows}, descriptors, reaction)
                actual = fit(x, y)
                expected = next(
                    r
                    for r in self.expected
                    if r["reaction"] == reaction
                    and r["descriptor"] == "current_stericx_fullprecision"
                    and r["label_convention"] == "official_ge"
                )
                self.assertEqual(actual["n_active"], expected["n_active"])
                self.assertEqual(actual["threshold"], expected["sklearn_threshold"])
                self.assertEqual(actual["accuracy"], expected["sklearn_accuracy"])
                self.assertEqual(actual["mcc"], expected["sklearn_mcc"])


if __name__ == "__main__":
    unittest.main()
