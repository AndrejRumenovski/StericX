"""Reject scientific changes even when ordinary numeric equality hides them."""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import scientifically_exact_optimization as exact


class ExactOptimizationOracleTests(unittest.TestCase):
    def compare(self, left, right, *, stderr=b"", request=b"identical"):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            observations = {}
            for name, rows in (("baseline", left), ("candidate", right)):
                directory = root / name
                directory.mkdir()
                (directory / "requests").write_bytes(
                    request if name == "candidate" else b"identical"
                )
                (directory / "stdout").write_text(
                    "".join(json.dumps(row) + "\n" for row in rows)
                )
                (directory / "stderr").write_bytes(
                    stderr if name == "candidate" else b""
                )
                path = directory / "manifest.json"
                path.write_text("{}")
                observations[path] = {
                    "lanes": {
                        "full_bins": {
                            **{
                                key: exact.oracle.file_record(directory / key)
                                for key in ("requests", "stdout", "stderr")
                            },
                            "returncode": 0,
                        }
                    }
                }
            args = argparse.Namespace(
                snapshot=root / "snapshot",
                baseline=root / "baseline",
                candidate=root / "candidate",
                output=root / "comparison",
            )
            with (
                patch.object(
                    exact,
                    "verify_frozen",
                    return_value={
                        "admitted_observations": exact.oracle.file_record(
                            args.baseline / "manifest.json"
                        )
                    },
                ),
                patch.object(
                    exact.replay,
                    "verify_observations",
                    side_effect=lambda path: (observations[path], {}),
                ),
                contextlib.redirect_stdout(io.StringIO()),
            ):
                exact.compare_observations(args)
            result = json.loads((args.output / "comparison.json").read_text())
            self.assertTrue(result["passed"])
            self.assertEqual(result["tolerance"], 0)

    def test_exact_complete_observation_passes(self):
        row = {"id": "x", "public": 1.0, "bins": [0.0, 2.0], "bits": [0, 1]}
        self.compare([row], [row])

    def test_private_bin_change_fails_with_identical_public_aggregate(self):
        with self.assertRaisesRegex(ValueError, "equivalence FAILED"):
            self.compare(
                [{"id": "x", "public": 1.0, "bins": [0.0, 2.0]}],
                [{"id": "x", "public": 1.0, "bins": [2.0, 0.0]}],
            )

    def test_signed_zero_change_is_not_hidden_by_python_equality(self):
        with self.assertRaisesRegex(ValueError, "equivalence FAILED"):
            self.compare([{"id": "x", "v": -0.0}], [{"id": "x", "v": 0.0}])

    def test_missing_error_case_fails(self):
        with self.assertRaisesRegex(ValueError, "equivalence FAILED"):
            self.compare([{"id": "x", "error": "ambiguous"}], [])

    def test_changed_stderr_fails(self):
        with self.assertRaisesRegex(ValueError, "equivalence FAILED"):
            self.compare([{"id": "x"}], [{"id": "x"}], stderr=b"new warning")

    def test_changed_request_fails_even_when_outputs_match(self):
        with self.assertRaisesRegex(ValueError, "requests differ"):
            self.compare([{"id": "x"}], [{"id": "x"}], request=b"changed")

    def test_frozen_copy_survives_live_edit_but_rejects_snapshot_edit(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / "live.rs"
            source.write_text("baseline")
            copied = exact.bound_copy(source, root / "frozen.rs")
            source.write_text("candidate")
            exact.oracle.verify(copied["snapshot"])
            (root / "frozen.rs").write_text("tampered")
            with self.assertRaises(ValueError):
                exact.oracle.verify(copied["snapshot"])

    def test_unadmitted_baseline_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "manifest.json").write_text("{}")
            args = argparse.Namespace(snapshot=root, baseline=root)
            with (
                patch.object(
                    exact,
                    "verify_frozen",
                    return_value={
                        "admitted_observations": {"sha256": "other baseline"}
                    },
                ),
                self.assertRaisesRegex(ValueError, "not the frozen admitted"),
            ):
                exact.compare_observations(args)

    def test_historical_or_other_build_workloads_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            build = root / "build.json"
            build.write_text("{}")
            workloads = root / "workloads.json"
            for value in (
                {"kind": "old_workloads"},
                {"kind": "corrected_science_profile_inputs", "build_manifest": {}},
            ):
                workloads.write_text(json.dumps(value))
                with self.assertRaisesRegex(ValueError, "admitted build"):
                    exact.corrected_workloads(workloads, build)

    def test_raw_observation_tampering_invalidates_frozen_baseline(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            raw = root / "stdout.jsonl"
            raw.write_text('{"id":"baseline"}\n')
            frozen = {
                "kind": "reviewed_corrected_optimization_baseline",
                "copies": [],
                "evidence": [],
                "raw_evidence": [exact.oracle.file_record(raw)],
            }
            raw.write_text('{"id":"tampered"}\n')
            with (
                patch.object(exact.replay, "load_manifest", return_value=frozen),
                self.assertRaises(ValueError),
            ):
                exact.verify_frozen(root / "manifest.json")


if __name__ == "__main__":
    unittest.main()
