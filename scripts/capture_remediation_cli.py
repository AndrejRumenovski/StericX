#!/usr/bin/env python3
"""Freeze corrected CLI behavior without blessing historical failures.

Replays every historical CLI request, remapping only output destinations. Records
all status/output changes for review. Python scientific behavior is covered by
the separate thermodynamics replay, whose driver handles newly rejected inputs.
Subsequent optimization comparisons must use the same prepared output paths.
"""

from __future__ import annotations

import argparse
import fcntl
import json
import os
import shutil
import subprocess
from pathlib import Path

import check_scientific_cli_equivalence as cli
import replay_scientific_remediation as replay


def prepare(args):
    original = args.historical.resolve() / "manifest.json"
    source = json.loads(original.read_text())
    cli.verify(source)
    root = replay.new_output(args.root)
    (root / "scratch").mkdir()
    cases = []
    for old in source["cases"]:
        replacements = {
            output: str(root / "scratch" / old["id"] / Path(output).name)
            for output in old["outputs"]
        }
        cases.append(
            {
                **old,
                "argv": [replacements.get(arg, arg) for arg in old["argv"]],
                "outputs": [replacements[path] for path in old["outputs"]],
                "historical_argv": old["argv"],
            }
        )
    replay.manifest(
        root / "manifest.json",
        {
            "kind": "corrected_cli_requests",
            "historical_manifest": replay.oracle.file_record(original),
            "inputs": source["inputs"],
            "historical_helper_sha256": source["harness_sha256"],
            "helper": replay.oracle.file_record(Path(__file__)),
            "cases": cases,
        },
    )


def verify_plan(root):
    plan = replay.load_manifest(root / "manifest.json")
    replay.oracle.verify(plan["historical_manifest"])
    replay.oracle.verify(plan["helper"])
    old = json.loads(Path(plan["historical_manifest"]["path"]).read_text())
    cli.verify(old)
    if plan["inputs"] != old["inputs"] or len(plan["cases"]) != len(old["cases"]):
        raise ValueError("CLI input inventory changed")
    for original, case in zip(old["cases"], plan["cases"], strict=True):
        if case["id"] != original["id"]:
            raise ValueError("CLI case identity/order changed")
        replacements = {
            output: str(root / "scratch" / original["id"] / Path(output).name)
            for output in original["outputs"]
        }
        if case["argv"] != [replacements.get(a, a) for a in original["argv"]]:
            raise ValueError("CLI arguments differ beyond output relocation")
    return plan


def observe(args):
    root = args.root.resolve()
    plan = verify_plan(root)
    build_path = args.build.resolve() / "manifest.json"
    build_record = replay.oracle.file_record(build_path)
    built, _ = replay.verify_build(build_path)
    helper_record = replay.oracle.file_record(Path(__file__))
    binary = replay.oracle.verify(built["binaries"]["native"])
    destination = replay.new_output(root / "runs" / args.label)
    env = dict(os.environ, LC_ALL="C", TZ="UTC", RAYON_NUM_THREADS=str(args.threads))
    env.pop("STERICX_PROFILE_PATH", None)
    cases = {}
    with (root / "capture.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        for case in plan["cases"]:
            output = destination / case["id"]
            output.mkdir()
            scratch = root / "scratch" / case["id"]
            scratch.mkdir(exist_ok=False)
            try:
                with (
                    (output / "stdout").open("xb") as stdout,
                    (output / "stderr").open("xb") as stderr,
                ):
                    run = subprocess.run(
                        [str(binary), *case["argv"]],
                        cwd=cli.REPO,
                        env=env,
                        stdout=stdout,
                        stderr=stderr,
                        timeout=600,
                    )
                shutil.copytree(scratch, output / "artifacts")
                cli.save(output / "result.json", {"returncode": run.returncode})
                if case["id"].startswith("deck_") and run.returncode == 0:
                    cli.validate_deck(output)
                cases[case["id"]] = {
                    "fingerprint": cli.fingerprint(output),
                    "historical_status": case["expected_status"],
                    "status_changed": run.returncode != case["expected_status"],
                }
            finally:
                shutil.rmtree(scratch)
    verify_plan(root)
    replay.verify_build(build_path)
    replay.oracle.verify(build_record)
    replay.oracle.verify(helper_record)
    replay.manifest(
        destination / "manifest.json",
        {
            "kind": "corrected_cli_capture",
            "complete_capture": True,
            "scientific_pass": False,
            "build_manifest": build_record,
            "plan": replay.oracle.file_record(root / "manifest.json"),
            "helper": helper_record,
            "binary": built["binaries"]["native"],
            "threads": args.threads,
            "cases": cases,
            "files": [
                replay.oracle.file_record(path)
                for path in sorted(destination.rglob("*"))
                if path.is_file()
            ],
        },
    )
    print(json.dumps({"rows": len(cases), "output": str(destination)}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    commands = parser.add_subparsers(dest="action", required=True)
    command = commands.add_parser("prepare")
    command.add_argument("--historical", required=True, type=Path)
    command.set_defaults(run=prepare)
    command = commands.add_parser("observe")
    command.add_argument("--build", required=True, type=Path)
    command.add_argument("--label", required=True)
    command.add_argument("--threads", type=int, choices=[1, 2, 4, 6], default=1)
    command.set_defaults(run=observe)
    args = parser.parse_args()
    args.run(args)


if __name__ == "__main__":
    main()
