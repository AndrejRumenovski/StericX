#!/usr/bin/env python3
"""Retain bounded verification of the completed, separately prepared inputs."""

import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "scripts"))
import prepare_corrected_profile as prep  # noqa: E402


def main():
    destination = Path(__file__).resolve().parent
    manifest_path = REPO / (
        ".stericx/profiling/scientifically_validated_optimization/"
        "corrected_workloads_v2/workloads/manifest.json"
    )
    manifest = json.loads(manifest_path.read_text())
    prep.profile.verify_manifest(manifest)
    prep.replay.verify_build(Path(manifest["build_manifest"]["path"]))
    commands = []
    for name, argv in (
        (
            "tests",
            [
                sys.executable,
                "-m",
                "unittest",
                "discover",
                "-s",
                "tests",
                "-p",
                "test_corrected_profile_preparation.py",
            ],
        ),
        (
            "lint",
            [
                str(REPO / ".venv/bin/ruff"),
                "check",
                "scripts/prepare_corrected_profile.py",
                "tests/test_corrected_profile_preparation.py",
            ],
        ),
        (
            "format",
            [
                str(REPO / ".venv/bin/ruff"),
                "format",
                "--check",
                "scripts/prepare_corrected_profile.py",
                "tests/test_corrected_profile_preparation.py",
            ],
        ),
    ):
        result = subprocess.run(argv, cwd=REPO, capture_output=True, check=False)
        log = destination / (name + ".log")
        with log.open("xb") as stream:
            stream.write(result.stdout + result.stderr)
        commands.append(
            {"argv": argv, "returncode": result.returncode, "log": prep.record(log)}
        )
    selection = json.loads(Path(manifest["corpus_selection"]["path"]).read_text())
    receipt = {
        "manifest": prep.record(manifest_path),
        "build": manifest["build_manifest"],
        "input_bundle": manifest["input_manifest"],
        "helper": prep.record(Path(prep.__file__)),
        "verifier": prep.record(Path(__file__)),
        "tests": prep.record(REPO / "tests/test_corrected_profile_preparation.py"),
        "commands": commands,
        "all_checks_passed": all(row["returncode"] == 0 for row in commands),
        "manifest_files_reverified": len(manifest["files"]),
        "source_copy_pairs": len(manifest["copied_inputs"]),
        "native_preparation_returncodes": [
            row["returncode"] for row in manifest["preparation_commands"]
        ],
        "selected_molecules": selection["selected_molecules"],
        "selected_conformers": selection["selected_conformers"],
        "excluded_molecules": len(selection["excluded"]),
        "workload_expected_records": {
            key: value["expected_records"]
            for key, value in manifest["workloads"].items()
        },
        "ensemble_molecule": manifest["ensemble_molecule_id"],
        "ensemble_conformers": manifest["workloads"]["ensemble_sdf"][
            "expected_conformers"
        ],
        "benchmarks_run": False,
        "scientific_acceptance": manifest["scientific_acceptance"],
    }
    with (destination / "receipt.json").open("x") as stream:
        json.dump(receipt, stream, indent=2)
        stream.write("\n")
    if not receipt["all_checks_passed"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
