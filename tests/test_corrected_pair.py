"""Paired benchmark admission/order tests; no native benchmark is executed."""

import argparse
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import benchmark_corrected_pair as pair


class CorrectedPairTests(unittest.TestCase):
    def test_all_pairs_alternate_and_both_binaries_have_one_warmup(self):
        plan = list(pair.schedule())
        self.assertEqual(
            plan[:2], [("warmup", 0, "baseline"), ("warmup", 0, "candidate")]
        )
        self.assertEqual(len(plan), 18)
        for index in range(8):
            expected = (
                ["baseline", "candidate"]
                if index % 2 == 0
                else ["candidate", "baseline"]
            )
            self.assertEqual(
                plan[2 + 2 * index : 4 + 2 * index],
                [("measured", index, side) for side in expected],
            )
        self.assertEqual(len(plan) * len(pair.THREADS) * len(pair.WORKLOADS), 720)

    def test_settings_preserve_frozen_affinity_and_disable_diagnostic_output(self):
        affinities = {
            str(t): value for t, value in pair.exact.inventory.AFFINITIES.items()
        }
        with (
            mock.patch.dict(
                os.environ, {"STERICX_PROFILE_PATH": "/bad", "RAYON_NUM_THREADS": "99"}
            ),
            mock.patch.object(pair.os, "sched_getaffinity", return_value=set(range(6))),
        ):
            for thread in pair.THREADS:
                env, affinity = pair.settings(thread, affinities)
                self.assertEqual(env["RAYON_NUM_THREADS"], str(thread))
                self.assertEqual(
                    affinity, set(map(int, affinities[str(thread)].split(",")))
                )
                self.assertNotIn("STERICX_PROFILE_PATH", env)
                self.assertEqual(env["LC_ALL"], "C")
                self.assertEqual(env["OMP_NUM_THREADS"], "1")
        with mock.patch.object(pair.os, "sched_getaffinity", return_value={5}):
            with self.assertRaisesRegex(ValueError, "CPUs unavailable"):
                pair.settings(1, affinities)

    def test_rejected_gate_cannot_launch_even_a_warmup_and_failure_is_retained(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "new"
            args = argparse.Namespace(
                snapshot=Path(temporary) / "snapshot",
                candidate_build=Path(temporary) / "candidate/manifest.json",
                gate=Path(temporary) / "gate",
                output=output,
            )
            with (
                mock.patch.object(
                    pair, "verify_gate", side_effect=ValueError("Wrong candidate build")
                ),
                mock.patch.object(pair.profile, "measure") as measure,
                mock.patch.object(pair.profile, "build_launcher") as launcher,
            ):
                with self.assertRaisesRegex(ValueError, "Wrong candidate build"):
                    pair.run(args)
                measure.assert_not_called()
                launcher.assert_not_called()
            summary = json.loads((output / "summary.json").read_text())
            self.assertEqual(summary["status"], "failed")
            self.assertEqual(summary["attempts"], [])
            self.assertFalse(
                json.loads((output / "manifest.json").read_text())["complete"]
            )

    def test_passing_boolean_does_not_replace_domain_evidence(self):
        value = dict(
            complete=True,
            passed=True,
            baseline_build="baseline",
            candidate_build="candidate",
            domains={},
        )
        with mock.patch.object(pair, "receipt", return_value=value):
            with self.assertRaisesRegex(ValueError, "All three"):
                pair.verify_reference_review({}, "baseline", "candidate")
        with self.assertRaisesRegex(ValueError, "Wrong candidate"):
            pair.bound_builds(value, "baseline", "other")

    def test_reference_result_failure_is_not_hidden_by_review_pass(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "actual.json"
            candidate = {"build": "candidate"}
            path.write_text(json.dumps({"failures": 1, "candidate_build": candidate}))
            evidence = pair.oracle.file_record(path)
            domain = dict(
                passed=True,
                evidence=[evidence],
                checks=[dict(evidence=evidence, field=["failures"], equals=0)],
            )
            value = dict(
                complete=True,
                passed=True,
                baseline_build="a",
                candidate_build=candidate,
                domains={
                    name: domain for name in ("geometry", "thermodynamics", "models")
                },
            )
            with mock.patch.object(pair, "receipt", return_value=value):
                with self.assertRaisesRegex(ValueError, "reference check failed"):
                    pair.verify_reference_review({}, "a", candidate)
            # A true-looking old result without candidate provenance is not
            # admitted even when every outer wrapper claims success.
            path.write_text('{"failures": 0}')
            evidence = pair.oracle.file_record(path)
            domain.update(
                evidence=[evidence],
                checks=[dict(evidence=evidence, field=["failures"], equals=0)],
            )
            with mock.patch.object(pair, "receipt", return_value=value):
                with self.assertRaisesRegex(ValueError, "candidate build provenance"):
                    pair.verify_reference_review({}, "a", candidate)

    def test_legacy_reference_raw_hash_maps_are_verified(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            raw = root / "raw.json"
            raw.write_text("original\n")
            result = root / "result.json"
            result.write_text(
                json.dumps(
                    {
                        "inputs": {str(raw): pair.oracle.sha(raw)},
                        "artifacts": {"raw.json": pair.oracle.sha(raw)},
                    }
                )
            )
            record = pair.oracle.file_record(result)
            pair.verify_result(record)
            raw.write_text("changed\n")
            with self.assertRaisesRegex(ValueError, "raw evidence changed"):
                pair.verify_result(record)

    def test_missing_prediction_exporter_provenance_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "exporter build receipts"):
            pair.verify_prediction_exporters({}, {}, {})

    def test_top_gate_cannot_substitute_another_candidate_build(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name in (
                "baseline.json",
                "candidate.json",
                "wrong.json",
                "manifest.json",
            ):
                (root / name).write_text(json.dumps({"name": name}))
            baseline = pair.oracle.file_record(root / "baseline.json")
            frozen = {
                "copies": [
                    {
                        "original": baseline,
                        "snapshot": {"path": str(root / "receipts/build.json")},
                    }
                ]
            }
            gate = root / "gate.json"
            pair.replay.manifest(
                gate,
                dict(
                    kind="corrected_candidate_performance_gate",
                    complete=True,
                    baseline_snapshot=pair.oracle.file_record(root / "manifest.json"),
                    baseline_build=baseline,
                    candidate_build=pair.oracle.file_record(root / "wrong.json"),
                ),
            )
            with (
                mock.patch.object(pair.exact, "verify_frozen", return_value=frozen),
                mock.patch.object(pair.replay, "verify_build") as build,
            ):
                with self.assertRaisesRegex(ValueError, "Wrong candidate build"):
                    pair.verify_gate(root, root / "candidate.json", gate)
                build.assert_not_called()

    def test_prediction_gate_rejects_truncated_or_mismatching_full_stream(self):
        value = dict(
            complete=True,
            exact_equivalence=True,
            baseline_build="a",
            candidate_build="b",
            encoding="ieee754-f32-le",
            records=1_000_000,
            data="data",
            weights="weights",
            threads={
                str(t): {
                    "baseline": dict(bytes=4_000_000, sha256="same"),
                    "candidate": dict(bytes=4_000_000, sha256="same"),
                }
                for t in pair.THREADS
            },
        )
        workloads = {
            "workloads": {
                "predict_1000000": {
                    "argv": ["predict", "--data", "data", "--weights", "weights"]
                }
            }
        }
        with (
            mock.patch.object(pair, "receipt", return_value=value),
            mock.patch.object(pair, "verify_prediction_exporters"),
            mock.patch.object(
                pair.oracle, "file_record", side_effect=lambda path: str(path)
            ),
        ):
            pair.verify_prediction({}, "a", "b", workloads)
            value["threads"]["4"]["candidate"]["sha256"] = "different"
            with self.assertRaisesRegex(ValueError, "bitstreams differ"):
                pair.verify_prediction({}, "a", "b", workloads)
            value["threads"]["4"]["candidate"].update(sha256="same", bytes=4)
            with self.assertRaisesRegex(ValueError, "truncated"):
                pair.verify_prediction({}, "a", "b", workloads)

    def test_stdout_or_stderr_mismatch_stops_and_retains_both_attempts(self):
        for changed in ("stdout", "stderr"):
            with (
                self.subTest(changed=changed),
                tempfile.TemporaryDirectory() as temporary,
            ):
                output = Path(temporary)
                baseline = output / "expected"
                baseline.mkdir()
                (baseline / "stdout.txt").write_text('[{"B1":1}]')
                (baseline / "stderr.txt").write_text("")
                workload = dict(
                    kind="descriptors", expected_records=1, argv=["descriptors"]
                )
                expected = pair.full_fingerprint(workload, baseline)
                context = dict(
                    frozen={"thread_affinities": {"1": "2"}},
                    expected={"1": {"case": expected}},
                    binaries={
                        side: {"path": side} for side in ("baseline", "candidate")
                    },
                )
                journal = {"attempts": []}

                def fake_measure(
                    argv, env, affinity, directory, launcher, changed=changed
                ):
                    directory.mkdir(parents=True)
                    value = (
                        '[{"B1":2}]'
                        if argv[0] == "candidate" and changed == "stdout"
                        else '[{"B1":1}]'
                    )
                    (directory / "stdout.txt").write_text(value)
                    (directory / "stderr.txt").write_text(
                        "new warning"
                        if argv[0] == "candidate" and changed == "stderr"
                        else ""
                    )
                    return {"wall_ns": 100, "peak_rss_bytes": 4096, "returncode": 0}

                with (
                    mock.patch.object(
                        pair.profile, "measure", side_effect=fake_measure
                    ),
                    mock.patch.object(pair.os, "sched_getaffinity", return_value={2}),
                ):
                    with self.assertRaisesRegex(ValueError, "Exact output mismatch"):
                        pair.measure_configuration(
                            "case",
                            workload,
                            1,
                            context,
                            output,
                            Path("launcher"),
                            journal,
                        )
                self.assertEqual(
                    [row["status"] for row in journal["attempts"]],
                    ["complete", "failed"],
                )
                self.assertTrue(
                    Path(journal["attempts"][1]["directory"])
                    .joinpath("metrics.json")
                    .is_file()
                )

    def test_paired_statistics_keep_an_extreme_sample_and_reject_missing_pairs(self):
        samples = []
        for phase, index, side in pair.schedule():
            samples.append(
                dict(
                    phase=phase,
                    pair=index,
                    side=side,
                    wall_ns=(1000 if index == 7 else 100) if side == "baseline" else 50,
                    cpu_user_s=0.1,
                    cpu_system_s=0.01,
                    cpu_percent_one_core=90,
                    peak_rss_bytes=4096,
                    minor_faults=1,
                    major_faults=0,
                    voluntary_context_switches=0,
                    involuntary_context_switches=0,
                )
            )
        result = pair.paired_summary(samples)
        self.assertEqual(result["paired_wall_speedup"]["values"], [2.0] * 7 + [20.0])
        self.assertEqual(result["paired_wall_speedup"]["maximum"], 20.0)
        self.assertEqual(result["exclusions"], [])
        with self.assertRaisesRegex(ValueError, "Incomplete measured pairs"):
            pair.paired_summary(samples[:-1])


if __name__ == "__main__":
    unittest.main()
