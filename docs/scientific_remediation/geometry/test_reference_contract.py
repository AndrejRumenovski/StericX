"""Guard the seal and complete-stream claims of the independent replay helper."""

import argparse
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import recompute_volume_means as replay
import sealed_reference as sealed


class ReferenceContract(unittest.TestCase):
    def test_changed_reference_is_rejected_before_execution(self):
        with tempfile.TemporaryDirectory() as directory:
            audit = Path(directory)
            reference = audit / "scripts/geometry_reference.py"
            reference.parent.mkdir()
            reference.write_text("raise AssertionError('tampered code executed')\n")
            with (
                patch.object(sealed, "AUDIT", audit),
                patch.object(sealed, "REFERENCE", reference),
                self.assertRaisesRegex(ValueError, "differs from immutable"),
            ):
                sealed.load_reference()

    def test_changed_audit_seal_and_original_output_path_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            manifest = Path(directory) / "manifest.json"
            manifest.write_text("{}\n")
            with (
                patch.object(sealed, "MANIFEST", manifest),
                self.assertRaisesRegex(ValueError, "pre-remediation seal"),
            ):
                sealed.load_reference()
        with self.assertRaisesRegex(ValueError, "immutable"):
            sealed.refuse_audit_output(sealed.AUDIT / "must_not_be_created")

    def fixture(self, root, request_ids, observation_ids, expected, limit=None):
        request = root / "requests.jsonl"
        observed = root / "observations.jsonl"
        request.write_text("".join(json.dumps({"id": x}) + "\n" for x in request_ids))
        observed.write_text(
            "".join(json.dumps({"id": x}) + "\n" for x in observation_ids)
        )
        return argparse.Namespace(
            requests=request,
            observations=observed,
            expected_rows=expected,
            limit=limit,
        )

    def test_wrong_expected_count_or_truncated_scope_cannot_be_complete(self):
        with tempfile.TemporaryDirectory() as directory:
            args = self.fixture(Path(directory), ["a"], ["a"], 2)
            coverage = {}
            with self.assertRaisesRegex(ValueError, "Expected 2"):
                list(replay.requests(args, coverage))
            self.assertNotIn("complete_source_stream_validated", coverage)

    def test_limited_compute_still_checks_entire_source_stream(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            args = self.fixture(root, ["a", "b"], ["a", "wrong"], 2, limit=1)
            with self.assertRaisesRegex(ValueError, "ID order"):
                list(replay.requests(args, {}))
            args = self.fixture(root, ["a", "b"], ["a", "b"], 2, limit=1)
            coverage = {}
            self.assertEqual(len(list(replay.requests(args, coverage))), 1)
            self.assertEqual(coverage["source_rows"], 2)
            self.assertTrue(coverage["complete_source_stream_validated"])


if __name__ == "__main__":
    unittest.main()
