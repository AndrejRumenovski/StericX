"""The performance gate must detect scientific differences, even one output bit."""

from __future__ import annotations

import json
import shutil
import struct
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import check_scientific_cli_equivalence as oracle


class ScientificCliOracleTests(unittest.TestCase):
    def prepare_capture(self, root):
        run = root / "run"
        case = {
            "id": "export",
            "expected_status": 0,
            "outputs": ["record.sigpack"],
            "expected_artifacts": {},
        }
        directory = run / case["id"]
        (directory / "artifacts").mkdir(parents=True)
        (directory / "artifacts/record.sigpack").write_bytes(b"record")
        (directory / "result.json").write_text('{"returncode": 0}')
        (directory / "stdout").write_bytes(b"value=1\ntotal_ms=2\n")
        (directory / "stderr").write_bytes(b"")
        (run / "bin").mkdir()
        (run / "bin/stericx").write_bytes(b"executable")
        (run / "python").mkdir()
        (run / "python/python.json").write_bytes(b"{}")
        (run / "python_source").mkdir()
        shutil.copyfile(oracle.__file__, root / "harness.py")
        manifest = {
            "cases": [case],
            "inputs": {},
            "harness_sha256": oracle.digest(root / "harness.py"),
        }
        oracle.save(root / "manifest.json", manifest)
        identity = {
            "manifest_sha256": oracle.digest(root / "manifest.json"),
            "harness_sha256": manifest["harness_sha256"],
            "binary_sha256": oracle.digest(run / "bin/stericx"),
            "python_sut": {},
        }
        oracle.save(run / "identity.json", identity)
        oracle.save(
            run / "complete.json",
            {
                "results": {
                    case["id"]: oracle.fingerprint(directory),
                    "python": oracle.digest(run / "python/python.json"),
                },
                "raw_hashes": oracle.raw_hashes(run, manifest),
            },
        )
        return run, manifest

    def test_capture_rejects_changed_manifest_executable_and_raw_metrics(self):
        for target, value, message in (
            ("manifest.json", b"{}", "manifest changed"),
            ("run/bin/stericx", b"another build", "executable changed"),
            ("run/export/stdout", b"value=1\ntotal_ms=3\n", "observation changed"),
            ("harness.py", b"different helper", "harness identity differs"),
        ):
            with self.subTest(target=target), tempfile.TemporaryDirectory() as name:
                root = Path(name)
                run, manifest = self.prepare_capture(root)
                oracle.verify_run(root, manifest, run)
                (root / target).write_bytes(value)
                with self.assertRaisesRegex(ValueError, message):
                    oracle.verify_run(root, manifest, run)

    def test_two_missing_success_outputs_cannot_be_an_equivalent_pair(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            run, manifest = self.prepare_capture(root)
            path = run / "export/artifacts/record.sigpack"
            path.unlink()
            with self.assertRaisesRegex(ValueError, "missing or empty required output"):
                oracle.verify_run(root, manifest, run)
            path.write_bytes(b"")
            with self.assertRaisesRegex(ValueError, "missing or empty required output"):
                oracle.verify_run(root, manifest, run)

    def test_only_explicit_process_metrics_are_ignored(self):
        raw = b"total_ms=1.2\nburied_volume=4.2\nother_ms=7\nerror: total_ms=0\n"
        self.assertEqual(
            oracle.scientific_stdout(raw),
            b"buried_volume=4.2\nother_ms=7\nerror: total_ms=0\n",
        )
        self.assertNotEqual(
            oracle.scientific_stdout(b"value=1.0000000000000000\n"),
            oracle.scientific_stdout(b"value=1.0000000000000002\n"),
        )

    def test_complete_fingerprint_detects_bits_errors_and_missing_artifacts(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            artifacts = root / "artifacts"
            artifacts.mkdir()
            (root / "result.json").write_text('{"returncode": 0}')
            (root / "stdout").write_bytes(b"value=1\ntotal_ms=2\n")
            (root / "stderr").write_bytes(b"")
            packed = artifacts / "record.sigpack"
            packed.write_bytes(struct.pack("=I", 0x3F800000))
            original = oracle.fingerprint(root)
            (root / "stdout").write_bytes(b"value=1\ntotal_ms=3\n")
            self.assertEqual(original, oracle.fingerprint(root))
            packed.write_bytes(struct.pack("=I", 0x3F800001))
            self.assertNotEqual(original, oracle.fingerprint(root))
            packed.write_bytes(struct.pack("=I", 0x3F800000))
            (root / "stderr").write_bytes(b"warning\n")
            self.assertNotEqual(original, oracle.fingerprint(root))
            (root / "stderr").write_bytes(b"")
            (root / "result.json").write_text('{"returncode": 2}')
            self.assertNotEqual(original, oracle.fingerprint(root))
            (root / "result.json").write_text('{"returncode": 0}')
            packed.unlink()
            self.assertNotEqual(original, oracle.fingerprint(root))

    def test_only_portable_creation_timestamp_is_normalized(self):
        with tempfile.TemporaryDirectory() as name:
            path = Path(name) / "portable.json"
            document = {
                "created": {"created_utc": "2026-09-21T00:00:00Z"},
                "weight": 1.0,
            }
            path.write_text(json.dumps(document))
            first = oracle.artifact_bytes(path)
            document["created"]["created_utc"] = "2026-09-22T00:00:00Z"
            path.write_text(json.dumps(document))
            self.assertEqual(first, oracle.artifact_bytes(path))
            document["weight"] = 1.0000000000000002
            path.write_text(json.dumps(document))
            self.assertNotEqual(first, oracle.artifact_bytes(path))


if __name__ == "__main__":
    unittest.main()
