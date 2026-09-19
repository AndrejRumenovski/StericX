"""Reject incomplete or approximately equal scientific-volume oracle outputs."""

import csv
import importlib.util
import struct
import tempfile
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location(
    "check_volume_exactness",
    Path(__file__).resolve().parents[1] / "scripts/check_volume_exactness.py",
)
oracle = importlib.util.module_from_spec(spec)
spec.loader.exec_module(oracle)


class VolumeOracleTests(unittest.TestCase):
    def fixture(self, directory):
        packed = struct.pack(
            "=8sIIQII32s32f",
            b"SIGPKV2\0",
            2,
            0x01020304,
            1,
            128,
            32,
            bytes(32),
            *[1.0] * 32,
        )
        (directory / "volumes.sigpack").write_bytes(packed)
        with (directory / "conformers.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(oracle.FIELDS)
            writer.writerow(["1.0"] * len(oracle.FIELDS))
        (directory / "stdout.txt").write_bytes(b"records_processed=1\ntotal_ms=4.0\n")
        (directory / "stderr.txt").write_bytes(b"")
        return packed

    def test_truncated_header_or_record_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            packed = self.fixture(directory)
            oracle.validate_artifacts(directory, 1, 1)
            for incomplete in (packed[:20], packed[:-1], packed + b"\0"):
                (directory / "volumes.sigpack").write_bytes(incomplete)
                with self.assertRaisesRegex(ValueError, "invalid packed"):
                    oracle.validate_artifacts(directory, 1, 1)

    def test_wrong_record_schema_or_nonfinite_values_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            packed = self.fixture(directory)
            invalid_header = bytearray(packed)
            struct.pack_into("=I", invalid_header, 28, 31)
            (directory / "volumes.sigpack").write_bytes(invalid_header)
            with self.assertRaisesRegex(ValueError, "invalid packed"):
                oracle.validate_artifacts(directory, 1, 1)
            for invalid in (float("nan"), float("inf")):
                invalid_numeric = bytearray(packed)
                struct.pack_into("=f", invalid_numeric, 64 + 31 * 4, invalid)
                (directory / "volumes.sigpack").write_bytes(invalid_numeric)
                with self.assertRaisesRegex(ValueError, "non-finite packed"):
                    oracle.validate_artifacts(directory, 1, 1)

    def test_every_conformer_and_all_nine_fields_are_required(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            self.fixture(directory)
            with self.assertRaisesRegex(ValueError, "incomplete conformer"):
                oracle.validate_artifacts(directory, 1, 2)
            path = directory / "conformers.csv"
            path.write_text(",".join(oracle.FIELDS[:-1]) + "\n" + ",".join(["1"] * 8))
            with self.assertRaisesRegex(ValueError, "missing conformer"):
                oracle.validate_artifacts(directory, 1, 1)
            path.write_text(
                ",".join(oracle.FIELDS) + "\n" + ",".join(["1"] * 8 + ["nan"])
            )
            with self.assertRaisesRegex(ValueError, "non-finite conformer"):
                oracle.validate_artifacts(directory, 1, 1)

    def test_single_bit_and_unrounded_csv_changes_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            packed = self.fixture(directory)
            expected = {"returncode": 0, "fingerprints": oracle.fingerprints(directory)}
            changes = (
                ("volumes.sigpack", packed[:-4] + struct.pack("=I", 0x3F800001)),
                (
                    "conformers.csv",
                    (directory / "conformers.csv")
                    .read_bytes()
                    .replace(b"1.0", b"1.0000001", 1),
                ),
                ("stderr.txt", b"unexpected warning\n"),
                ("stdout.txt", b"records_processed=2\ntotal_ms=4.0\n"),
            )
            for name, content in changes:
                self.fixture(directory)
                (directory / name).write_bytes(content)
                actual = {
                    "returncode": 0,
                    "fingerprints": oracle.fingerprints(directory),
                }
                with (
                    self.subTest(name=name),
                    self.assertRaisesRegex(ValueError, "changed"),
                ):
                    oracle.assert_exact(expected, actual, "fixture")
            with self.assertRaisesRegex(ValueError, "changed"):
                oracle.assert_exact(expected, dict(expected, returncode=2), "fixture")

    def test_stdout_normalization_excludes_only_named_metadata(self):
        raw = b"total_ms=1\nvbur=1.2345678\nnear_vbur=2\nother_ms=3\ntotal_ms\n"
        self.assertEqual(
            oracle.normalized_stdout(raw),
            b"vbur=1.2345678\nnear_vbur=2\nother_ms=3\ntotal_ms\n",
        )

    def test_mutating_frozen_inputs_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "fixture.xyz"
            path.write_text("first")
            manifest = {"inputs": [{"path": str(path), "sha256": oracle.sha(path)}]}
            oracle.verify_inputs(manifest)
            path.write_text("second")
            with self.assertRaisesRegex(ValueError, "input changed"):
                oracle.verify_inputs(manifest)


if __name__ == "__main__":
    unittest.main()
