#!/usr/bin/env python3
"""Verify and exactly compare corrected CLI captures without modifying old receipts.

This gate covers all 93 CLI cases. Python scientific behavior requires the
separate corrected thermodynamics gate. Exact behavior does not itself establish
scientific validity; the accepted baseline must also carry independent receipts.
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import capture_remediation_cli as capture
import replay_scientific_remediation as replay

oracle = replay.oracle
cli = capture.cli


def canonical_cases(root, historical):
    result = []
    for old in historical["cases"]:
        replacements = {
            path: str(root / "scratch" / old["id"] / Path(path).name)
            for path in old["outputs"]
        }
        result.append(
            {
                **old,
                "argv": [replacements.get(arg, arg) for arg in old["argv"]],
                "outputs": [replacements[path] for path in old["outputs"]],
                "historical_argv": old["argv"],
            }
        )
    return result


def verify_plan(root):
    plan = capture.verify_plan(root)
    if plan["kind"] != "corrected_cli_requests":
        raise ValueError("Wrong corrected CLI plan kind")
    if plan["helper"]["sha256"] != oracle.sha(Path(capture.__file__)):
        raise ValueError("Executing capture verifier differs from frozen plan helper")
    historical = oracle.read(oracle.verify(plan["historical_manifest"]))
    cli.verify(historical)
    if plan["historical_helper_sha256"] != historical["harness_sha256"]:
        raise ValueError("Historical normalization helper identity changed")
    if plan["cases"] != canonical_cases(root, historical):
        raise ValueError(
            "Corrected case differs from full canonical historical request"
        )
    ids = [case["id"] for case in plan["cases"]]
    if len(ids) != 93 or len(set(ids)) != 93:
        raise ValueError("All 93 distinct corrected CLI cases are required")
    return plan


def required_outputs(case, directory, actual_status):
    if actual_status != 0:
        return
    if (directory / "stdout").stat().st_size == 0:
        raise ValueError("Successful CLI case has empty stdout: " + case["id"])
    required = set(case["expected_artifacts"]) | {
        Path(path).name for path in case["outputs"]
    }
    for name in required:
        artifact = directory / "artifacts" / name
        if not artifact.is_file() or not artifact.stat().st_size:
            raise ValueError(
                "Successful case lacks nonempty output: " + case["id"] + "/" + name
            )
    if case["id"].startswith("deck_"):
        cli.validate_deck(directory)


def verify_run(root, label):
    if Path(label).name != label:
        raise ValueError("Run label must be a directory name")
    plan = verify_plan(root)
    run = root / "runs" / label
    receipt_identity = oracle.file_record(run / "manifest.json")
    receipt = replay.load_manifest(run / "manifest.json")
    if receipt.get("kind") != "corrected_cli_capture" or not receipt.get(
        "complete_capture"
    ):
        raise ValueError("Complete corrected CLI capture required")
    if receipt["plan"] != oracle.file_record(root / "manifest.json"):
        raise ValueError("Run plan identity differs")
    if receipt["helper"] != plan["helper"]:
        raise ValueError("Capture helper differs from request plan")
    oracle.verify(receipt["helper"])
    build_path = oracle.verify(receipt["build_manifest"])
    built, _ = replay.verify_build(build_path)
    if receipt["binary"] != built["binaries"]["native"]:
        raise ValueError("Captured native binary is not the built frozen executable")
    oracle.verify(receipt["binary"])
    if receipt["threads"] not in (1, 2, 4, 6):
        raise ValueError("Unexpected thread configuration")
    expected_ids = [case["id"] for case in plan["cases"]]
    if set(receipt["cases"]) != set(expected_ids):
        raise ValueError("Incomplete or extra CLI case capture")
    actual_files = {
        str(path.resolve())
        for path in run.rglob("*")
        if path.is_file()
        and path not in (run / "manifest.json", run / "manifest.json.sha256")
    }
    declared_files = [item["path"] for item in receipt["files"]]
    if (
        len(set(declared_files)) != len(declared_files)
        or set(declared_files) != actual_files
    ):
        raise ValueError("Incomplete, duplicate, or extra raw file inventory")
    for item in receipt["files"]:
        oracle.verify(item)
    fingerprints = {}
    for case in plan["cases"]:
        directory = run / case["id"]
        observed = cli.fingerprint(directory)
        recorded = receipt["cases"][case["id"]]
        expected_record = dict(
            fingerprint=observed,
            historical_status=case["expected_status"],
            status_changed=observed["returncode"] != case["expected_status"],
        )
        if recorded != expected_record:
            raise ValueError(
                "Raw output differs from captured case receipt: " + case["id"]
            )
        required_outputs(case, directory, observed["returncode"])
        fingerprints[case["id"]] = observed
    oracle.verify(receipt_identity)
    return receipt, fingerprints, receipt_identity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--mode", choices=("threads", "candidate"), required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    output = replay.new_output(args.output)
    helpers = [
        Path(__file__),
        Path(capture.__file__),
        Path(cli.__file__),
        Path(replay.__file__),
        Path(oracle.__file__),
        Path(replay.refs.__file__),
    ]
    executing = [oracle.file_record(path) for path in helpers]
    for path in helpers:
        shutil.copyfile(path, output / path.name)
    try:
        baseline, before, baseline_identity = verify_run(root, args.baseline)
        candidate, after, candidate_identity = verify_run(root, args.candidate)
        if baseline["plan"] != candidate["plan"]:
            raise ValueError("Runs used different exact inputs/command plans")
        if args.mode == "threads":
            if (
                baseline["binary"] != candidate["binary"]
                or baseline["build_manifest"] != candidate["build_manifest"]
            ):
                raise ValueError(
                    "Thread comparison requires the identical source build/binary"
                )
        elif baseline["threads"] != candidate["threads"]:
            raise ValueError("Candidate comparison requires matching thread counts")
        differences = [
            {"case": name, "baseline": before[name], "candidate": after[name]}
            for name in before
            if before[name] != after[name]
        ]
        for item in executing:
            oracle.verify(item)
        # Check all raw bytes once more after comparison, not merely the stored
        # fingerprint maps. Changes during verification cannot pass silently.
        for receipt in (baseline, candidate):
            for item in receipt["files"]:
                oracle.verify(item)
            oracle.verify(receipt["plan"])
            replay.verify_build(oracle.verify(receipt["build_manifest"]))
            oracle.verify(receipt["binary"])
        oracle.verify(baseline_identity)
        oracle.verify(candidate_identity)
        replay.manifest(
            output / "manifest.json",
            dict(
                kind="corrected_cli_exact_comparison",
                complete=True,
                exact_equivalence=not differences,
                scientific_global_pass=False,
                baseline=baseline_identity,
                candidate=candidate_identity,
                mode=args.mode,
                baseline_threads=baseline["threads"],
                candidate_threads=candidate["threads"],
                cases=len(before),
                differences=differences,
                plan=baseline["plan"],
                executing_helpers=executing,
                comparison=(
                    "Exact return codes, ordered stdout/stderr bytes, "
                    "and every artifact; "
                    "only declared process metrics and portable.created_utc omitted."
                ),
                python_scope=(
                    "Separately checked by corrected thermodynamics receipt; "
                    "not included in these 93 CLI cases."
                ),
                files=[
                    oracle.file_record(path)
                    for path in sorted(output.iterdir())
                    if path.is_file()
                ],
            ),
        )
        print(
            json.dumps(
                dict(
                    cases=len(before), differences=len(differences), output=str(output)
                )
            )
        )
        if differences:
            raise SystemExit(2)
    except Exception as exc:
        oracle.write_new(
            output / "failure.json", dict(type=type(exc).__name__, detail=str(exc))
        )
        raise


if __name__ == "__main__":
    main()
