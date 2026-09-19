"""Protect measurement validity and exact scientific comparison contracts."""

import importlib.util
import json
import os
import tempfile
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location(
    "profile_stericx",
    Path(__file__).resolve().parents[1] / "scripts/profile_stericx.py",
)
profile = importlib.util.module_from_spec(spec)
spec.loader.exec_module(profile)


class ProfileHarnessTests(unittest.TestCase):
    def test_sdf_preserves_all_coordinate_decimal_tokens(self):
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "ligand.xyz"
            source.write_text(
                "2\nexample\nP -1.123456789 2.987654321 0.000000001\nH 0 0 1\n"
            )
            sdf = profile.xyz_to_sdf([source, source])
            self.assertEqual(sdf.count("$$$$\n"), 2)
            self.assertEqual(sdf.count("-1.123456789 2.987654321 0.000000001 P"), 2)

    def test_repeated_sigpack_preserves_record_bytes_and_partial_cycle(self):
        with tempfile.TemporaryDirectory() as temporary:
            source, target = Path(temporary) / "in", Path(temporary) / "out"
            source.write_bytes(bytes(range(128)))
            self.assertEqual(profile.repeat_sigpack(source, target, 5), 2)
            self.assertEqual(
                target.read_bytes(), source.read_bytes() * 2 + bytes(range(64))
            )

    def test_only_named_metrics_and_output_paths_are_excluded(self):
        source = (
            b"total_ms=2.001\nmse_kcal2_per_mol2=0.12345678\nunknown_ms=123\n"
            b"rayon_threads=1\nsigpack_output=/a\n"
        )
        normalized = profile.normalized_stdout(source, "predict")
        self.assertEqual(
            normalized,
            b"mse_kcal2_per_mol2=0.12345678\nunknown_ms=123\nrayon_threads=1\n",
        )
        self.assertEqual(profile.normalized_stdout(source, "descriptors"), source)

    def test_descriptor_success_cannot_hide_skipped_geometry(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            (folder / "stdout.txt").write_text(json.dumps([{"conformers": 1}]))
            (folder / "stderr.txt").write_text("skipped invalid.xyz: no donor\n")
            with self.assertRaisesRegex(RuntimeError, "skipped inputs"):
                profile.check_result(
                    {"kind": "descriptors", "expected_records": 1}, folder
                )

    def test_frozen_input_modification_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "input"
            path.write_bytes(b"fixed")
            manifest = {
                "files": [
                    {"path": str(path), "bytes": 5, "sha256": profile.digest(path)}
                ]
            }
            profile.verify_manifest(manifest)
            path.write_bytes(b"other")
            with self.assertRaisesRegex(RuntimeError, "frozen input changed"):
                profile.verify_manifest(manifest)

    @unittest.skipUnless(
        hasattr(os, "wait4") and hasattr(os, "sched_getaffinity"), "Linux measurement"
    )
    def test_wait4_reaps_exact_child_and_restores_controller_affinity(self):
        with tempfile.TemporaryDirectory() as temporary:
            before = os.sched_getaffinity(0)
            for index in range(2):
                metrics = profile.measure(
                    ["/bin/true"],
                    dict(os.environ),
                    {min(before)},
                    Path(temporary) / str(index),
                )
                self.assertEqual(metrics["returncode"], 0)
                self.assertGreater(metrics["wall_ns"], 0)
                self.assertGreaterEqual(metrics["peak_rss_bytes"], 0)
                with self.assertRaises(ChildProcessError):
                    os.waitpid(metrics["pid"], os.WNOHANG)
                self.assertEqual(os.sched_getaffinity(0), before)

    @unittest.skipUnless(
        hasattr(os, "wait4") and hasattr(os, "sched_getaffinity"), "Linux measurement"
    )
    def test_signaled_child_is_rejected_and_actual_signal_retained(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / "signal"
            with self.assertRaisesRegex(RuntimeError, "child failed"):
                profile.measure(
                    ["/bin/sh", "-c", "kill -TERM $$"],
                    dict(os.environ),
                    {min(os.sched_getaffinity(0))},
                    directory,
                )
            metrics = json.loads((directory / "metrics.json").read_text())
            self.assertEqual(metrics["returncode"], -15)


if __name__ == "__main__":
    unittest.main()
