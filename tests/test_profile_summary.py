"""Check profiling attribution and comparability with tiny synthetic runs."""

from __future__ import annotations

import copy
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "summarize_profile.py"
SPEC = importlib.util.spec_from_file_location("stericx_profile_summary", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
profile_summary = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(profile_summary)


def function(
    name: str,
    stage: str,
    exclusive: int,
    *,
    inclusive: int | None = None,
    thread: int = 0,
) -> dict:
    return {
        "function": name,
        "stage": stage,
        "thread_index": thread,
        "calls": 1,
        "exclusive_ns": exclusive,
        "inclusive_ns": exclusive if inclusive is None else inclusive,
    }


def report(functions: list[dict]) -> dict:
    return {
        "dropped_scopes": 0,
        "observed_threads": len({item["thread_index"] for item in functions}),
        "functions": functions,
        "allocations": {
            "allocation_calls": 10,
            "requested_bytes": 128,
            "peak_live_requested_bytes": 64,
        },
    }


class AmdahlTests(unittest.TestCase):
    def test_zero_fraction_and_complete_removal_boundaries(self) -> None:
        self.assertEqual(profile_summary.amdahl(0), 1)
        self.assertIsNone(profile_summary.amdahl(1))
        self.assertEqual(profile_summary.amdahl(1, 2), 2)
        self.assertEqual(profile_summary.amdahl(0.75), 4)
        self.assertEqual(profile_summary.amdahl(0.75, 2), 1.6)
        self.assertEqual(profile_summary.amdahl(0.75, 1), 1)

    def test_invalid_fractions_or_slowdown_factors_are_rejected(self) -> None:
        for fraction in (-0.01, 1.01, float("nan")):
            with self.subTest(fraction=fraction), self.assertRaises(ValueError):
                profile_summary.amdahl(fraction)
        for improvement in (-1, 0, 0.99):
            with self.subTest(improvement=improvement), self.assertRaises(ValueError):
                profile_summary.amdahl(0.5, improvement)


class SummaryTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.native = self.series("native-binary")
        self.diagnostic = self.series("diagnostic-binary")
        self.sample = report(
            [
                function("cli::run", "orchestration", 100, inclusive=800),
                function("pipeline", "conformer_processing", 100, inclusive=700),
                function("occupied", "buried_volume", 500),
                function("parse", "file_parsing", 100),
            ]
        )

    @staticmethod
    def series(binary: str) -> dict:
        return {
            "status": "complete",
            "manifest_sha256": "fixed-inputs",
            "affinity": [2],
            "environment": {"RAYON_NUM_THREADS": "1"},
            "binary_sha256": binary,
            "results": {
                "workload": {
                    "fingerprints": {"stdout_sha256": "same-scientific-output"},
                    "summary": {"wall_ns": {"median": 1000}},
                    "measured_runs": 1,
                    "exactness_scope": "all emitted values byte-for-byte",
                }
            },
        }

    def save_series(self) -> None:
        for name, value in (("native", self.native), ("diagnostic", self.diagnostic)):
            path = self.root / "runs" / name / "summary.json"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(value))

    def save_sample(self, value: dict, *, index: int = 0, wall: int = 1000) -> None:
        directory = (
            self.root / "runs" / "diagnostic" / "workload" / f"measured_{index:03d}"
        )
        directory.mkdir(parents=True)
        (directory / "stages.json").write_text(json.dumps(value))
        (directory / "metrics.json").write_text(json.dumps({"wall_ns": wall}))

    def summarize(self) -> dict:
        self.save_series()
        return profile_summary.summarize(self.root, "native", "diagnostic")["results"][
            "workload"
        ]

    def test_exclusive_and_inclusive_bounds_use_external_process_wall(self) -> None:
        self.save_sample(self.sample)
        result = self.summarize()
        buried = result["main_thread_stage_wall"]["buried_volume"]
        self.assertEqual(buried["median_diagnostic_wall_fraction"], 0.5)
        self.assertEqual(buried["maximum_end_to_end_speedup"], 2)
        self.assertAlmostEqual(
            buried["end_to_end_speedup_if_stage_twice_as_fast"], 4 / 3
        )
        pipeline = next(
            item
            for item in result["main_thread_functions"]
            if item["function"] == "pipeline"
        )
        self.assertAlmostEqual(
            pipeline["maximum_end_to_end_speedup_removing_self_work"], 10 / 9
        )
        self.assertAlmostEqual(
            pipeline["maximum_end_to_end_speedup_removing_function_and_callees"],
            10 / 3,
        )
        self.assertEqual(result["outside_cli_scope_median_ns"], 200)
        self.assertAlmostEqual(
            sum(
                item["median_diagnostic_wall_fraction"]
                for item in result["main_thread_stage_wall"].values()
            ),
            1,
        )

    def test_worker_elapsed_is_kept_out_of_main_wall_shares(self) -> None:
        # A long worker span overlaps the main thread's wait; adding it would
        # incorrectly produce more than one process wall's worth of stages.
        self.sample["functions"].append(
            function("worker", "buried_volume", 900, thread=1)
        )
        self.sample["observed_threads"] = 2
        self.save_sample(self.sample)
        result = self.summarize()
        self.assertEqual(
            result["worker_summed_elapsed_ns_by_stage"], {"buried_volume": 900}
        )
        buried = result["main_thread_stage_wall"]["buried_volume"]
        self.assertEqual(buried["median_exclusive_ns"], 500)
        self.assertEqual(buried["maximum_end_to_end_speedup"], 2)
        self.assertNotIn(
            "worker", {item["function"] for item in result["main_thread_functions"]}
        )

    def test_cli_root_identifies_main_thread_without_assuming_index_zero(self) -> None:
        for item in self.sample["functions"]:
            item["thread_index"] = 1
        self.sample["functions"].append(
            function("worker", "buried_volume", 900, thread=0)
        )
        self.sample["observed_threads"] = 2
        self.save_sample(self.sample)
        result = self.summarize()
        self.assertEqual(result["outside_cli_scope_median_ns"], 200)
        self.assertEqual(
            result["worker_summed_elapsed_ns_by_stage"]["buried_volume"], 900
        )

    def test_missing_function_stage_and_worker_samples_are_zero_filled(self) -> None:
        self.diagnostic["results"]["workload"]["measured_runs"] = 3
        for index in range(3):
            functions = [function("cli::run", "orchestration", 800, inclusive=800)]
            if index == 1:
                functions[0]["exclusive_ns"] -= 200
                functions.extend(
                    [
                        function("rare", "rare_stage", 200),
                        function("rare_worker", "worker_stage", 400, thread=1),
                    ]
                )
            self.save_sample(report(functions), index=index)
        result = self.summarize()
        rare = next(
            item
            for item in result["main_thread_functions"]
            if item["function"] == "rare"
        )
        self.assertEqual(rare["median_calls"], 0)
        self.assertEqual(rare["median_exclusive_ns"], 0)
        self.assertEqual(rare["median_inclusive_ns"], 0)
        self.assertEqual(
            rare["maximum_end_to_end_speedup_removing_function_and_callees"], 1
        )
        self.assertEqual(
            result["main_thread_stage_wall"]["rare_stage"]["median_exclusive_ns"], 0
        )
        self.assertEqual(result["worker_summed_elapsed_ns_by_stage"]["worker_stage"], 0)

    def test_fractions_are_calculated_per_run_before_taking_median(self) -> None:
        self.diagnostic["results"]["workload"]["measured_runs"] = 3
        for index, (wall, buried) in enumerate(((100, 10), (1000, 900), (10000, 1000))):
            root = wall * 9 // 10
            self.save_sample(
                report(
                    [
                        function(
                            "cli::run", "orchestration", root - buried, inclusive=root
                        ),
                        function("occupied", "buried_volume", buried),
                    ]
                ),
                index=index,
                wall=wall,
            )
        result = self.summarize()
        stage = result["main_thread_stage_wall"]["buried_volume"]
        self.assertEqual(stage["median_exclusive_ns"], 900)
        self.assertEqual(stage["median_diagnostic_wall_fraction"], 0.1)
        self.assertAlmostEqual(stage["maximum_end_to_end_speedup"], 10 / 9)

    def test_input_execution_setting_and_output_mismatches_are_rejected(self) -> None:
        self.save_sample(self.sample)
        original = copy.deepcopy(self.diagnostic)
        cases = (
            ("manifest_sha256", "changed-inputs", "inputs differ"),
            ("affinity", [3], "settings differ"),
            ("environment", {"RAYON_NUM_THREADS": "2"}, "settings differ"),
        )
        for key, value, error in cases:
            with self.subTest(key=key):
                self.diagnostic = copy.deepcopy(original)
                self.diagnostic[key] = value
                with self.assertRaisesRegex(ValueError, error):
                    self.summarize()
        self.diagnostic = copy.deepcopy(original)
        self.diagnostic["results"]["workload"]["fingerprints"]["stdout_sha256"] = (
            "different-output"
        )
        with self.assertRaisesRegex(ValueError, "scientific output changed"):
            self.summarize()

    def test_incomplete_series_or_missing_stage_reports_are_rejected(self) -> None:
        self.save_sample(self.sample)
        for series in (self.native, self.diagnostic):
            series["status"] = "running"
            with self.assertRaisesRegex(ValueError, "completed measurement series"):
                self.summarize()
            series["status"] = "complete"
        self.diagnostic["results"]["workload"]["measured_runs"] = 2
        with self.assertRaisesRegex(ValueError, "missing stage reports"):
            self.summarize()

    def test_dropped_scopes_are_rejected(self) -> None:
        self.sample["dropped_scopes"] = 1
        self.save_sample(self.sample)
        with self.assertRaisesRegex(ValueError, "scope capacity exceeded"):
            self.summarize()

    def test_invalid_main_wall_partition_is_rejected(self) -> None:
        self.sample["functions"][0]["exclusive_ns"] += 1
        self.save_sample(self.sample)
        with self.assertRaisesRegex(ValueError, "do not partition CLI wall"):
            self.summarize()

    def test_internal_root_cannot_exceed_external_wall(self) -> None:
        self.save_sample(self.sample, wall=799)
        with self.assertRaisesRegex(ValueError, "scope exceeds measured process wall"):
            self.summarize()

    def test_cli_root_must_be_unique(self) -> None:
        self.sample["functions"].append(function("cli::run", "orchestration", 1))
        self.save_sample(self.sample)
        with self.assertRaisesRegex(ValueError, "exactly one CLI wall root"):
            self.summarize()

    def test_partial_numerical_comparison_scope_is_preserved(self) -> None:
        # Prediction stdout only contains a preview and aggregate; this report
        # must retain that limitation until a separate all-predictions audit.
        limitation = "CLI preview and aggregate only; complete f32 audit required"
        self.native["results"]["workload"]["exactness_scope"] = limitation
        self.save_sample(self.sample)
        self.assertEqual(self.summarize()["native"]["exactness_scope"], limitation)


if __name__ == "__main__":
    unittest.main()
