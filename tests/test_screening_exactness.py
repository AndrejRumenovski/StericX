"""Validate the strict oracle independently of either StericX binary."""

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location(
    "check_screening_exactness",
    Path(__file__).resolve().parents[1] / "scripts/check_screening_exactness.py",
)
oracle = importlib.util.module_from_spec(spec)
spec.loader.exec_module(oracle)


class ScreeningOracleTests(unittest.TestCase):
    def test_every_observable_must_match_exactly(self):
        expected = {
            "returncode": 0,
            "stdout_sha256": "science",
            "stderr_sha256": "warn",
        }
        oracle.assert_exact(expected, dict(expected, argv=["different_binary"]), "same")
        for field, changed in (
            ("returncode", 2),
            ("stdout_sha256", "rounding changed"),
            ("stderr_sha256", "different error"),
        ):
            with (
                self.subTest(field=field),
                self.assertRaisesRegex(ValueError, "changed"),
            ):
                oracle.assert_exact(
                    expected, dict(expected, **{field: changed}), "case"
                )

    def test_mutating_a_frozen_fixture_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "input"
            path.write_text("valid")
            manifest = {"inputs": [{"path": str(path), "sha256": oracle.sha(path)}]}
            oracle.verify_inputs(manifest)
            path.write_text("other")
            with self.assertRaisesRegex(ValueError, "input changed"):
                oracle.verify_inputs(manifest)

    def test_matching_failures_cannot_masquerade_as_valid_coverage(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            with self.assertRaisesRegex(ValueError, "unexpected baseline status"):
                oracle.validate_coverage(
                    {"name": "must succeed", "expected_status": 0},
                    {"returncode": 2},
                    directory,
                )

    def test_missing_intended_exclusion_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            (directory / "stdout.txt").write_text(
                json.dumps({"screened": 1, "skipped": 0})
            )
            (directory / "stderr.txt").write_text("")
            with self.assertRaisesRegex(ValueError, "intended exclusions"):
                oracle.validate_coverage(
                    {"name": "mixed", "expected_status": 0, "minimum_skipped": 1},
                    {"returncode": 0},
                    directory,
                )


if __name__ == "__main__":
    unittest.main()
