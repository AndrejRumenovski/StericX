"""Admission tests for the additive corrected CLI comparison gate."""

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
SPEC = importlib.util.spec_from_file_location(
    "corrected_cli_compare", REPO / "scripts/compare_corrected_cli.py"
)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


class CorrectedCliCompareTests(unittest.TestCase):
    def test_changed_case_metadata_is_rejected_even_when_argv_matches(self):
        root = Path("/new-corrected-cli")
        old = {
            "harness_sha256": "normalizer",
            "cases": [
                {
                    "id": f"case_{i}",
                    "argv": ["command", "--output", f"/old/{i}/output"],
                    "outputs": [f"/old/{i}/output"],
                    "expected_status": 0,
                    "expected_artifacts": {},
                }
                for i in range(93)
            ],
        }
        plan = {
            "kind": "corrected_cli_requests",
            "helper": {"sha256": "capture"},
            "historical_manifest": {"path": "/old/manifest"},
            "historical_helper_sha256": "normalizer",
            "cases": module.canonical_cases(root, old),
        }
        plan["cases"][0]["expected_status"] = 2
        with (
            mock.patch.object(module.capture, "verify_plan", return_value=plan),
            mock.patch.object(module.oracle, "sha", return_value="capture"),
            mock.patch.object(
                module.oracle, "verify", return_value=Path("/old/manifest")
            ),
            mock.patch.object(module.oracle, "read", return_value=old),
            mock.patch.object(module.cli, "verify"),
        ):
            with self.assertRaisesRegex(ValueError, "full canonical"):
                module.verify_plan(root)

    def test_new_success_requires_outputs_even_if_historical_case_failed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "stdout").write_text("records_processed=1\n")
            (root / "artifacts").mkdir()
            case = {
                "id": "new_success",
                "expected_status": 2,
                "expected_artifacts": {},
                "outputs": ["/scratch/parsed.sigpack"],
            }
            with self.assertRaisesRegex(ValueError, "nonempty output"):
                module.required_outputs(case, root, 0)

    def test_only_named_metrics_are_ignored_and_error_order_is_exact(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "artifacts").mkdir()
            (root / "result.json").write_text(json.dumps({"returncode": 0}))
            (root / "stdout").write_text("value=1\ntotal_ms=2\n")
            (root / "stderr").write_text("first\nsecond\n")
            original = module.cli.fingerprint(root)
            (root / "stdout").write_text("value=1\ntotal_ms=999\n")
            self.assertEqual(original, module.cli.fingerprint(root))
            (root / "stdout").write_text("value=2\ntotal_ms=999\n")
            self.assertNotEqual(original, module.cli.fingerprint(root))
            (root / "stdout").write_text("value=1\ntotal_ms=2\n")
            (root / "stderr").write_text("second\nfirst\n")
            self.assertNotEqual(original, module.cli.fingerprint(root))


if __name__ == "__main__":
    unittest.main()
