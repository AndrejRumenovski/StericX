#!/usr/bin/env python3
"""Recheck every new performance sample without changing measurement rules."""

from __future__ import annotations

import argparse
import json
import statistics
from pathlib import Path

import benchmark_paton_extensions as benchmark
import report_paton_performance as retained


def read(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def distribution(values):
    median = statistics.median(values)
    return {
        "values_ns": values,
        "median_ns": median,
        "mad_ns": statistics.median(abs(v - median) for v in values),
        "minimum_ns": min(values),
        "maximum_ns": max(values),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.input = args.input.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    journal = read(args.input / "journal.json")
    retained.require(journal["complete"], "Incomplete measurement campaign")
    retained.require(
        journal["quiet_window"]["user_confirmed"]
        and journal["quiet_window"]["agent_heavy_work_complete"],
        "Unconfirmed window",
    )
    for item in journal["bindings"]:
        benchmark.verify(item)
    for command in journal["commands"]:
        retained.require(command["returncode"] == 0, "Failed measured command")
        benchmark.verify(command["log"])
        for item in command.get("produced_records", []):
            benchmark.verify(item)
    native_summary = read(args.input / "native/summary.json")
    workload_path = next(
        Path(item["path"])
        for item in native_summary["bindings"]
        if item["path"].endswith("workloads/manifest.json")
    )
    workloads = read(workload_path)
    benchmark.profile.verify_manifest(workloads)
    # Absolute path intentionally selects this new campaign while reusing the
    # unchanged 504-attempt fingerprint/recomputation/investigation logic.
    native = retained.review_campaign(
        str(args.input / "native"), workloads["workloads"]
    )
    native["confirmed_quiet_window"] = True
    save(args.output / "native_review.json", native)
    (args.output / "TIMING_TABLE.md").write_text(retained.timing_table(native))
    allocations = {}
    for row in native_summary["allocations"]:
        benchmark.verify(row["report"])
        actual = read(row["report"]["path"])["allocations"]
        retained.require(actual == row["allocations"], "Allocation record changed")
        allocations.setdefault((row["threads"], row["workload"]), {})[row["side"]] = (
            actual
        )
    retained.require(
        len(allocations) == 14 and all(len(v) == 2 for v in allocations.values()),
        "Incomplete allocation comparisons",
    )
    allocation_rows = [
        {
            "threads": key[0],
            "workload": key[1],
            "baseline": value["baseline"],
            "candidate": value["candidate"],
            "differences": {
                field: {
                    "baseline": value["baseline"].get(field),
                    "candidate": value["candidate"].get(field),
                }
                for field in set(value["baseline"]) | set(value["candidate"])
                if value["baseline"].get(field) != value["candidate"].get(field)
            },
        }
        for key, value in sorted(allocations.items())
    ]
    python = {
        side: {method: [] for method in ["native", "bridge", "public"]}
        for side in ["baseline", "candidate"]
    }
    expected_result = None
    retained.require(
        journal["python_campaign_order"]
        == ["baseline", "candidate", "candidate", "baseline"],
        "Python A/B/B/A order changed",
    )
    for index, side in enumerate(journal["python_campaign_order"]):
        directory = args.input / f"python_{index}_{side}"
        result = read(directory / "summary.json")
        retained.require(
            result["complete"] and result["scientific_results_exact"],
            "Incomplete Python timing",
        )
        retained.require(
            result["iterations_per_block"] == 100 and len(result["attempts"]) == 24,
            "Python timing schedule changed",
        )
        for item in result["bindings"]:
            benchmark.verify(item)
        command = next(
            c for c in journal["commands"] if c["name"] == f"python_{index}_{side}"
        )
        retained.require(
            len(command["produced_records"]) == 11,
            "Python raw-output hash coverage changed",
        )
        expected_schedule = [
            (pair, method)
            for pair in range(8)
            for method in (
                ["native", "bridge", "public"]
                if pair % 2 == 0
                else ["public", "bridge", "native"]
            )
        ]
        retained.require(
            [(a["pair"], a["method"]) for a in result["attempts"]] == expected_schedule,
            "Python within-campaign schedule changed",
        )
        retained.require(
            {p.name for p in directory.glob("native_*.json")}
            == {"native_preflight.json", *[f"native_{i}.json" for i in range(8)]},
            "Native Python-timing evidence coverage changed",
        )
        current = read(directory / "bridge_preflight.json")
        if expected_result is None:
            expected_result = current
        retained.require(
            current == expected_result, "Python baseline/candidate result changed"
        )
        campaign_samples = {method: [] for method in ["native", "bridge", "public"]}
        for attempt in result["attempts"]:
            retained.require(
                attempt["status"] == "complete", "Failed Python timing block"
            )
            python[side][attempt["method"]].append(attempt["elapsed_ns"] / 100)
            campaign_samples[attempt["method"]].append(attempt["elapsed_ns"] / 100)
            if attempt["method"] == "native":
                raw = read(directory / f"native_{attempt['pair']}.json")
                retained.require(
                    raw["iterations"] == 100
                    and raw["elapsed_ns"] == attempt["elapsed_ns"],
                    "Native raw elapsed time differs from Python timing record",
                )
        for method, values in campaign_samples.items():
            calculated = distribution(values)
            calculated["values_ns_per_call"] = calculated.pop("values_ns")
            retained.require(
                calculated == result["results"][method],
                "Python timing summary failed recomputation",
            )
        for path in directory.glob("native_*.json"):
            value = read(path)
            retained.require(
                value["result"] == expected_result, "Native ensemble output changed"
            )
    python_rows = []
    for method in ["native", "bridge", "public"]:
        retained.require(
            all(len(python[side][method]) == 16 for side in python),
            "Missing Python samples",
        )
        before, after = (
            distribution(python[side][method]) for side in ["baseline", "candidate"]
        )
        python_rows.append(
            {
                "method": method,
                "baseline": before,
                "candidate": after,
                "candidate_over_baseline_medians": after["median_ns"]
                / before["median_ns"],
            }
        )
    ensemble_allocations = journal["ensemble_allocation_runs"]
    retained.require(len(ensemble_allocations) == 2, "Missing ensemble diagnostics")
    retained.require(
        {r["side"] for r in ensemble_allocations} == {"baseline", "candidate"},
        "Wrong ensemble diagnostic sides",
    )
    for row in ensemble_allocations:
        benchmark.verify(row["metrics"])
        benchmark.verify(row["stdout"])
        retained.require(
            read(row["stdout"]["path"]) == row["values"], "Ensemble allocation changed"
        )
        row["process_metrics"] = read(row["metrics"]["path"])
        retained.require(
            row["values"]["conformers"] == 6
            and row["values"]["measured_calls"] == 1
            and row["values"]["failed_allocation_calls"] == 0,
            "Invalid ensemble allocation diagnostic",
        )
    result = {
        "kind": "dbstep_resolution_matched_performance_review",
        "measurement_journal": benchmark.record(args.input / "journal.json"),
        "reviewer": benchmark.record(__file__),
        "native_review": benchmark.record(args.output / "native_review.json"),
        "all_raw_fingerprints_and_summaries_rechecked": True,
        "native_attempts": native["attempts"],
        "native_regression_triggers": native["regression_triggers"],
        "maximum_median_slowdown_percent": native["maximum_median_slowdown_percent"],
        "native_allocation_comparisons": allocation_rows,
        "ensemble_allocation_runs": ensemble_allocations,
        "python_campaign_order": journal["python_campaign_order"],
        "python_methods": python_rows,
        "python_regression_triggers": [
            r["method"]
            for r in python_rows
            if r["candidate_over_baseline_medians"] > 1.1
        ],
        "native_binary_bytes_identical": native_summary["binaries"]["baseline"][
            "sha256"
        ]
        == native_summary["binaries"]["candidate"]["sha256"],
        "release_decision": (
            "Separate scientific/release review; triggers require investigation."
        ),
    }
    save(args.output / "review.json", result)
    print(
        json.dumps(
            {
                "native_attempts": native["attempts"],
                "regression_triggers": native["regression_triggers"],
                "python_regression_triggers": result["python_regression_triggers"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
