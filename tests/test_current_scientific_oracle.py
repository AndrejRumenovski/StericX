"""Corrupt or incomplete evidence must never pass the exact scientific gate."""

from __future__ import annotations

import copy
import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import check_current_scientific_equivalence as oracle


def first(path):
    with path.open() as stream:
        return json.loads(next(stream))


class CurrentScientificOracleTests(unittest.TestCase):
    def geometry(self):
        audit = oracle.DEFAULT_AUDIT
        return (
            first(audit / "geometry/requests.jsonl"),
            first(audit / "geometry/sut_outputs.jsonl"),
        )

    def test_private_dump_error_cannot_replace_successful_public_bins(self):
        request, result = self.geometry()
        self.assertTrue(oracle.validate_dump(request, result, {}))
        result["dump"] = {"error": "broken private observer"}
        for allowed in ({}, {request["id"]: result["dump"]}):
            with self.assertRaisesRegex(
                ValueError, "successful public BV without bins"
            ):
                oracle.validate_dump(request, result, allowed)

    def test_public_error_does_not_waive_private_success_or_allow_new_error(self):
        request, result = self.geometry()
        result["buried_volume"] = {"error": "public asymmetry guard"}
        self.assertTrue(oracle.validate_dump(request, result, {}))
        result["dump"] = {"error": "public asymmetry guard"}
        with self.assertRaisesRegex(ValueError, "Unreviewed private dump error"):
            oracle.validate_dump(request, result, {})
        self.assertFalse(
            oracle.validate_dump(request, result, {request["id"]: result["dump"]})
        )

    def test_missing_frame_bin_basis_atom_or_requested_points_fails(self):
        request, good = self.geometry()
        for field in ("quadrants", "octants", "basis", "aligned_atoms", "near_vbur"):
            with self.subTest(field=field):
                result = copy.deepcopy(good)
                del result["dump"]["orientations"][0][field]
                with self.assertRaises(ValueError):
                    oracle.validate_dump(request, result, {})
        result = copy.deepcopy(good)
        result["dump"]["orientations"].pop()
        with self.assertRaisesRegex(ValueError, "all three frames"):
            oracle.validate_dump(request, result, {})
        request["include_points"] = True
        with self.assertRaisesRegex(ValueError, "grid points"):
            oracle.validate_dump(request, good, {})

    def test_missing_public_section_and_one_bit_tamper_fail(self):
        request, result = self.geometry()
        del result["sterimol_bond"]
        with self.assertRaisesRegex(ValueError, "required observation fields"):
            oracle.validate_result(request, result)
        request, result = self.geometry()
        bits = result["buried_volume"]["_bits"]["buried_volume"]
        result["buried_volume"]["_bits"]["buried_volume"] = f"{int(bits, 16) ^ 1:08x}"
        with self.assertRaisesRegex(ValueError, "Number/bit observation mismatch"):
            oracle.validate_result(request, result)

    def test_all_bins_require_finite_numbers_and_positive_grid(self):
        request, good = self.geometry()
        self.assertTrue(oracle.validate_dump(request, good, {}, finite=True))
        for field in ("quadrants", "octants", "basis", "aligned_atoms", "near_vbur"):
            with self.subTest(field=field):
                result = copy.deepcopy(good)
                frame = result["dump"]["orientations"][0]
                if field in ("quadrants", "octants"):
                    frame[field][0] = None
                elif field == "basis":
                    frame[field][0][0] = None
                elif field == "aligned_atoms":
                    frame[field][0]["radius_squared"] = None
                else:
                    frame[field] = None
                with self.assertRaisesRegex(ValueError, "finite private numeric"):
                    oracle.validate_dump(request, result, {}, finite=True)
        result = copy.deepcopy(good)
        result["dump"]["grid_count"] = 0
        with self.assertRaisesRegex(ValueError, "grid count"):
            oracle.validate_dump(request, result, {}, finite=True)
        # Other lanes preserve their exact reviewed nonfinite serialization;
        # the all-bin finite policy is explicit and cannot silently broaden.
        oracle.numeric_array([None], 1, "documented numerical lane")

    def test_private_frame_indices_must_be_in_bounds_integers(self):
        request, good = self.geometry()
        for indices in (
            ["a", "b", "c"],
            [True, 2, 3],
            [-1, 2, 3],
            [len(request["atoms"]), 2, 3],
        ):
            with self.subTest(indices=indices):
                result = copy.deepcopy(good)
                result["dump"]["neighbors"] = indices
                for frame, plane in zip(
                    result["dump"]["orientations"], indices, strict=True
                ):
                    frame["plane"] = plane
                with self.assertRaisesRegex(ValueError, "frame neighbors"):
                    oracle.validate_dump(request, result, {})
        result = copy.deepcopy(good)
        frame = result["dump"]["orientations"][0]
        frame["plane"] = float(frame["plane"])
        with self.assertRaisesRegex(ValueError, "reordered private frames"):
            oracle.validate_dump(request, result, {})

    def fixture(self, root):
        audit = oracle.DEFAULT_AUDIT
        request = json.loads((audit / "kinetics/inputs/kinetics.json").read_text())[0]
        result = first(audit / "kinetics/frozen_outputs/kinetics.jsonl")
        for name, value in (("request.jsonl", request), ("expected.jsonl", result)):
            (root / name).write_bytes(oracle.canonical(value))
        (root / "empty").write_bytes(b"")
        (root / "helper").write_bytes(b"fixed helper")
        (root / "binary").write_bytes(b"fixed executable")
        (root / "oracle.json").write_bytes(b"{}")
        (root / "build.json").write_bytes(b"{}")
        (root / "evidence.json").write_bytes(b"{}")
        lane = {
            "requests": oracle.file_record(root / "request.jsonl"),
            "expected_stdout": oracle.file_record(root / "expected.jsonl"),
            "expected_stderr": oracle.file_record(root / "empty"),
            "expected_rows": 1,
            "expected_returncode": 0,
            "allowed_dump_errors": {},
            "expected_coverage": oracle.coverage(
                root / "request.jsonl", root / "expected.jsonl"
            ),
        }
        prepared = {
            "helper": oracle.file_record(root / "helper"),
            "evidence_lock": oracle.file_record(root / "evidence.json"),
            "excluded_scope": [],
            "audit": str(audit),
            "lanes": {"kinetics": lane},
        }
        compiled = {"binary": oracle.file_record(root / "binary")}
        directory = root / "observed/kinetics"
        directory.mkdir(parents=True)
        (directory / "stdout.jsonl").write_bytes(oracle.canonical(result))
        (directory / "stderr.txt").write_bytes(b"")
        row = {
            "requests": lane["requests"],
            "stdout": oracle.file_record(directory / "stdout.jsonl"),
            "stderr": oracle.file_record(directory / "stderr.txt"),
            "returncode": 0,
            "timed_out": False,
            "coverage": lane["expected_coverage"],
            "frozen_fidelity": True,
        }
        oracle.write_new(directory / "result.json", row)
        manifest = {
            "schema_version": oracle.SCHEMA_VERSION,
            "kind": "live_scientific_observations",
            "complete": True,
            "error": None,
            "oracle_manifest": oracle.file_record(root / "oracle.json"),
            "build_manifest": oracle.file_record(root / "build.json"),
            **compiled,
            "helper": prepared["helper"],
            "evidence_lock": prepared["evidence_lock"],
            "classification_counts": {},
            "excluded_scope": [],
            "requested_lanes": ["kinetics"],
            "full_observer_scope": False,
            "frozen_fidelity": True,
            "lanes": {"kinetics": row},
        }
        oracle.write_new(root / "observed/manifest.json", manifest)
        return manifest, prepared, compiled

    def test_empty_scope_forged_coverage_fidelity_and_stitched_identity_fail(self):
        mutations = (
            lambda m: m.update(lanes={}, requested_lanes=[], full_observer_scope=True),
            lambda m: m.update(requested_lanes=["geometry"]),
            lambda m: m.update(full_observer_scope=True),
            lambda m: m["lanes"]["kinetics"]["coverage"].update(rows=195),
            lambda m: m["lanes"]["kinetics"].update(frozen_fidelity=False),
            lambda m: m.update(frozen_fidelity=False),
            lambda m: m.update(binary=m["helper"]),
            lambda m: m.update(helper=m["binary"]),
        )
        for mutation in mutations:
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as name:
                root = Path(name)
                manifest, prepared, compiled = self.fixture(root)
                with (
                    patch.object(
                        oracle,
                        "load_oracle",
                        return_value=(prepared, {"classification_counts": {}}),
                    ),
                    patch.object(oracle, "load_build", return_value=compiled),
                ):
                    oracle.validate_observations(root / "observed")
                    mutation(manifest)
                    (root / "observed/manifest.json").write_text(json.dumps(manifest))
                    with self.assertRaises(ValueError):
                        oracle.validate_observations(root / "observed")

    def test_build_rejects_mixed_oracle_helper_binary_and_changed_source(self):
        for mutation, message in (
            ("oracle_manifest", "identical prepared oracle"),
            ("helper", "helper identities differ"),
            ("binary", "recorded build output"),
            ("source", "Locked file changed"),
        ):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as name:
                root = Path(name)
                (root / "target/release").mkdir(parents=True)
                for name in (
                    "manifest.json",
                    "helper",
                    "source.rs",
                    "target/release/stericx-audit-observer",
                ):
                    (root / name).write_bytes(b"{}")
                prepared = {"helper": oracle.file_record(root / "helper")}
                compiled = {
                    "kind": "live_scientific_observer_build",
                    "schema_version": oracle.SCHEMA_VERSION,
                    "returncode": 0,
                    "oracle_manifest": oracle.file_record(root / "manifest.json"),
                    "helper": prepared["helper"],
                    "binary": oracle.file_record(
                        root / "target/release/stericx-audit-observer"
                    ),
                    "sources": [oracle.file_record(root / "source.rs")],
                    "adapter_files": [],
                }
                if mutation == "source":
                    (root / "source.rs").write_text("changed")
                else:
                    compiled[mutation] = oracle.file_record(root / "source.rs")
                with (
                    patch.object(oracle, "read", return_value=compiled),
                    self.assertRaisesRegex(ValueError, message),
                ):
                    oracle.load_build(root, root, prepared)

    def test_equal_unanchored_pair_is_not_admitted(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            manifest, prepared, _ = self.fixture(root)
            manifest["frozen_fidelity"] = False
            with patch.object(
                oracle, "validate_observations", return_value=(manifest, prepared)
            ):
                output = root / "comparison.json"
                with (
                    redirect_stdout(io.StringIO()),
                    self.assertRaisesRegex(ValueError, "anchored equivalence"),
                ):
                    oracle.comparison(
                        SimpleNamespace(
                            baseline=root / "observed",
                            candidate=root / "observed",
                            output=output,
                        )
                    )
                result = json.loads(output.read_text())
                self.assertTrue(result["equivalent"])
                self.assertFalse(result["admitted"])

    def test_changed_stream_missing_rows_or_reordered_ids_fail(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            _, prepared, compiled = self.fixture(root)
            stdout = root / "observed/kinetics/stdout.jsonl"
            good = stdout.read_bytes()
            stdout.write_bytes(good + b"\n")
            with (
                patch.object(
                    oracle,
                    "load_oracle",
                    return_value=(prepared, {"classification_counts": {}}),
                ),
                patch.object(oracle, "load_build", return_value=compiled),
                self.assertRaisesRegex(ValueError, "stream identity changed"),
            ):
                oracle.validate_observations(root / "observed")
            stdout.write_bytes(b"")
            with self.assertRaisesRegex(ValueError, "count mismatch"):
                oracle.coverage(root / "request.jsonl", stdout)
            row = json.loads(good)
            row["id"] = "another conformer"
            stdout.write_bytes(oracle.canonical(row))
            with self.assertRaisesRegex(ValueError, "out-of-order ID"):
                oracle.coverage(root / "request.jsonl", stdout)


if __name__ == "__main__":
    unittest.main()
