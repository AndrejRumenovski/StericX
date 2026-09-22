"""Regression tests for the corrected-model reference adapter's admission rules."""

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

REPO = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "remediation_model_recheck",
    REPO / "docs/scientific_remediation/models/recheck.py",
)
recheck = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(recheck)


class ModelReferenceAdmissionTests(unittest.TestCase):
    def test_nonfinite_and_missing_values_never_pass_numeric_checks(self):
        check = recheck.Checks()
        for value in (None, float("nan"), float("inf"), True):
            check.close("invalid", value, 1.0)
        self.assertFalse(any(item["passed"] for item in check.items))

    def test_caught_observer_panic_is_not_an_expected_invalid_input(self):
        check = recheck.Checks()
        recheck.extra_references(
            [{"id": "malformed_feature_index", "op": "model"}],
            [{"error": "observer caught panic; see stderr"}],
            check,
        )
        self.assertEqual(len(check.items), 2)
        self.assertTrue(all(not item["passed"] for item in check.items))

    def test_extra_observer_missing_row_is_rejected(self):
        with self.assertRaises(ValueError):
            recheck.extra_references([{"id": "missing"}], [], recheck.Checks())

    def test_unchanged_valid_model_number_requires_exact_equality(self):
        check = recheck.Checks()
        recheck.exact_historical_numbers(
            check, "fit", {"weights": [1.000000000000001]}, {"weights": [1.0]}
        )
        self.assertFalse(check.items[-1]["passed"])

    def test_missing_historical_numeric_field_is_rejected(self):
        check = recheck.Checks()
        recheck.exact_historical_numbers(check, "fit", {}, {"weights": [1.0]})
        self.assertEqual(len(check.items), 1)
        self.assertFalse(check.items[0]["passed"])

    def test_numpy_scalar_reference_comparisons_serialize(self):
        import numpy as np

        check = recheck.Checks()
        check.close("descriptor", 1.0, np.float64(1.0), atol=0, rtol=0)
        self.assertTrue(check.items[0]["passed"])
        json.dumps(check.items, allow_nan=False)

    def test_linked_raw_file_mutation_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "raw.json"
            path.write_text('{"prediction":1}\n')
            record = recheck.oracle.file_record(path)
            recheck.verify_records({"files": [record]})
            path.write_text('{"prediction":2}\n')
            with self.assertRaises(ValueError):
                recheck.verify_records({"files": [record]})

    def test_cli_receipt_cannot_substitute_another_build(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            build = root / "build.json"
            build.write_text(json.dumps({"identity": "new"}))
            receipt = {
                "kind": "corrected_cli_capture",
                "complete_capture": True,
                "build_manifest": {"path": str(build), "sha256": "wrong"},
            }
            with mock.patch.object(
                recheck.replay, "load_manifest", return_value=receipt
            ):
                with self.assertRaisesRegex(ValueError, "build identities differ"):
                    recheck.verify_cli(root, build, {})


if __name__ == "__main__":
    unittest.main()
