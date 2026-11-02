#!/usr/bin/env python3
"""Paired descriptor timings against the Paton task's frozen native baseline."""

from __future__ import annotations

import argparse
import json
import platform
import sys
from pathlib import Path

import benchmark_corrected_pair as paired
import profile_stericx as profile

WORKLOADS = (
    "single_small",
    "single_large",
    "ensemble_sdf",
    "conformers_56",
    "descriptors_10000",
    "parse_ensembles_1000",
    "db_build_ensembles",
)
AFFINITIES = {"1": "2", "2": "2,3", "4": "0,1,2,3", "6": "0,1,2,3,4,5"}


def record(path):
    path = Path(path).resolve()
    return {
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": profile.digest(path),
    }


def verify(item):
    if record(item["path"]) != item:
        raise ValueError(f"Evidence changed: {item['path']}")


def load_build(path):
    path = path.resolve() / "manifest.json"
    result = json.loads(path.read_text())
    for item in [*result["sources"], *result["binaries"].values()]:
        verify(item)
    return record(path), result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--workloads", required=True, type=Path)
    parser.add_argument("--review", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--allocation-baseline", type=Path)
    parser.add_argument("--allocation-candidate", type=Path)
    args = parser.parse_args()
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    baseline_record, baseline = load_build(args.baseline)
    candidate_record, candidate = load_build(args.candidate)
    review = json.loads(args.review.read_text())
    if not review.get("descriptor_performance_admitted"):
        raise ValueError("Independent scientific review has not admitted timing")
    if review["baseline_build"] != baseline_record:
        raise ValueError("Review binds a different baseline")
    if review["candidate_build"] != candidate_record:
        raise ValueError("Review binds a different candidate")
    for item in review["evidence"]:
        verify(item)
    workloads = json.loads(args.workloads.read_text())
    profile.verify_manifest(workloads)
    bindings = [
        baseline_record,
        candidate_record,
        record(args.workloads),
        record(args.review),
        record(__file__),
        record(profile.__file__),
        record(paired.__file__),
        record(Path(__file__).with_name("profile_child.c")),
        record(
            Path(__file__).resolve().parents[1]
            / "docs/paton_recommendations/performance/PROTOCOL.md"
        ),
    ]
    binaries = {
        "baseline": baseline["binaries"]["native"],
        "candidate": candidate["binaries"]["native"],
    }
    launcher = profile.build_launcher(args.output / "launcher")
    journal = {
        "kind": "paton_descriptor_paired_performance",
        "created_utc": profile.now(),
        "command": sys.argv,
        "bindings": bindings,
        "binaries": binaries,
        "machine": platform.uname()._asdict(),
        "affinities": AFFINITIES,
        "controlled_environment": profile.CONTROLLED_ENV,
        "cache_policy": "warm; no flushing",
        "launcher": record(launcher),
        "attempts": [],
        "results": {},
        "complete": False,
    }
    summary_path = args.output / "summary.json"
    profile.write_json(summary_path, journal)
    for thread in paired.THREADS:
        env, affinity = paired.settings(thread, AFFINITIES)
        journal["results"][str(thread)] = {}
        for name in WORKLOADS:
            workload = workloads["workloads"][name]
            expected = None
            samples = []
            for phase, pair, side in paired.schedule():
                directory = (
                    args.output / f"t{thread}" / name / f"{phase}_{pair:02d}" / side
                )
                argv = [
                    binaries[side]["path"],
                    *[x.replace("{run_dir}", str(directory)) for x in workload["argv"]],
                ]
                attempt = dict(
                    workload=name,
                    threads=thread,
                    phase=phase,
                    pair=pair,
                    side=side,
                    directory=str(directory),
                )
                journal["attempts"].append(attempt)
                profile.write_json(summary_path, journal)
                try:
                    metrics = profile.measure(argv, env, affinity, directory, launcher)
                    fingerprint = paired.full_fingerprint(workload, directory)
                    if expected is None:
                        if side != "baseline" or phase != "warmup":
                            raise ValueError("Baseline warmup must establish identity")
                        expected = fingerprint
                    if fingerprint != expected:
                        raise ValueError(f"Scientific output changed: {directory}")
                    metrics.update(
                        phase=phase, pair=pair, side=side, fingerprints=fingerprint
                    )
                    profile.write_json(directory / "metrics.json", metrics)
                    attempt.update(
                        status="complete", metrics=record(directory / "metrics.json")
                    )
                    samples.append(metrics)
                except BaseException as error:
                    attempt.update(status="failed", error=str(error))
                    raise
                finally:
                    profile.write_json(summary_path, journal)
            result = paired.paired_summary(samples)
            for side in binaries:
                result[side]["records_per_second"] = (
                    workload["expected_records"]
                    * 1e9
                    / result[side]["wall_ns"]["median"]
                )
            result["fingerprints"] = expected
            journal["results"][str(thread)][name] = result
            profile.write_json(summary_path, journal)
            print(f"Completed t{thread} {name}", flush=True)

    # Instrumented diagnostics are separate and never contribute timing samples.
    diagnostics = {
        "baseline": args.allocation_baseline,
        "candidate": args.allocation_candidate,
    }
    if any(diagnostics.values()) and not all(diagnostics.values()):
        raise ValueError("Both allocation executables are required")
    journal["allocations"] = []
    if all(diagnostics.values()):
        journal["allocation_binaries"] = {
            side: record(path) for side, path in diagnostics.items()
        }
        for thread in (1, 6):
            env, affinity = paired.settings(thread, AFFINITIES)
            for name in WORKLOADS:
                workload = workloads["workloads"][name]
                for side, binary in diagnostics.items():
                    directory = args.output / "allocations" / f"t{thread}" / name / side
                    report = directory / "profile.json"
                    env = {**env, "STERICX_PROFILE_PATH": str(report)}
                    argv = [
                        str(binary.resolve()),
                        *[
                            x.replace("{run_dir}", str(directory))
                            for x in workload["argv"]
                        ],
                    ]
                    profile.measure(argv, env, affinity, directory, launcher)
                    fingerprint = paired.full_fingerprint(workload, directory)
                    if (
                        fingerprint
                        != journal["results"][str(thread)][name]["fingerprints"]
                    ):
                        raise ValueError("Instrumented scientific output changed")
                    report_value = json.loads(report.read_text())
                    journal["allocations"].append(
                        dict(
                            threads=thread,
                            workload=name,
                            side=side,
                            report=record(report),
                            allocations=report_value["allocations"],
                        )
                    )
                    profile.write_json(summary_path, journal)
    for item in bindings + list(binaries.values()):
        verify(item)
    profile.verify_manifest(workloads)
    load_build(args.baseline)
    load_build(args.candidate)
    journal.update(complete=True, completed_utc=profile.now())
    profile.write_json(summary_path, journal)


if __name__ == "__main__":
    main()
