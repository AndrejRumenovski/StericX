#!/usr/bin/env python3
"""Freeze advanced replay evidence and explicitly separate the historical gates."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import re
import shutil
import tarfile
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source, output = args.input.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    references = source / "references"
    lane_paths = {
        "sterimol": "sterimol/summary.json",
        "volume": "volume/summary.json",
        "density": "density/summary.json",
        "profiles": "profiles/summary.json",
        "populations": "populations/report.json",
        "ensemble_reference": "ensemble_reference/summary.json",
    }
    lanes = {key: read(references / path) for key, path in lane_paths.items()}
    legacy_paths = {
        "sterimol": "docs/validation/dbstep/sterimol_v1/summary.json",
        "volume": "docs/validation/dbstep/volume_v2/summary.json",
        "density": "docs/validation/dbstep/density_v2/summary.json",
        "profiles": "docs/validation/dbstep/profiles_v1/summary.json",
        "ensemble_reference": (
            "docs/paton_recommendations/conformers/ensemble_study_v1/summary.json"
        ),
    }
    comparisons = {}
    for key, old_path in legacy_paths.items():
        old, new = ROOT / old_path, references / lane_paths[key]
        comparisons[key] = {
            "historical_path": old_path,
            "historical_sha256": sha(old),
            "new_sha256": sha(new),
            "json_equal": read(old) == lanes[key],
        }
    installs = {}
    for name in ["wheel-3.11", "wheel-3.12", "wheel-3.13", "wheel-3.14", "sdist-3.12"]:
        directory = source / "installs" / name
        receipt = read(directory / "receipt.json")
        tests = (directory / "tests.log").read_text()
        installs[name] = {
            "passed": receipt.get("passed", False),
            "artifact_unchanged": receipt.get("artifact_unchanged", False),
            "test_outcome": re.findall(r"\d+ passed[^\n]*", tests),
            "artifact": receipt["artifact"],
            "all_commands_exit_zero": all(
                c["returncode"] == 0 for c in receipt["commands"]
            ),
        }
    phase_receipts = {
        phase: read(source / phase / "receipt.json")
        for phase in ["references", "distribution", "installs"]
    }
    failed_commands = [
        {"phase": phase, **command}
        for phase, receipt in phase_receipts.items()
        for command in receipt["commands"]
        if command["returncode"] != 0
    ]
    summary = {
        "created_utc": datetime.now(UTC).isoformat(),
        "reviewer_sha256": sha(__file__),
        "source_commits": {p: r["source_commit"] for p, r in phase_receipts.items()},
        "source_changes": {p: r["source_changes"] for p, r in phase_receipts.items()},
        "execution_failures_preserved": failed_commands,
        "historical_summary_comparison": comparisons,
        "scientific_gates": {
            "sterimol_conditional": lanes["sterimol"][
                "conditional_numerical_gate_pass"
            ],
            "hard_sphere_original_strict": lanes["volume"]["finite_lattice_gate_pass"],
            "density_declared_grid_scope": lanes["density"]["admitted"],
            "profile_preservation": lanes["profiles"]["profile_preservation_pass"],
            "profile_hard_sphere_original_strict": lanes["profiles"][
                "hard_sphere_strict_gate_all_positions"
            ],
            "populations": not lanes["populations"]["failures"],
            "ensemble_reference": not any(
                lanes["ensemble_reference"][key]
                for key in [
                    "reference_failures",
                    "conformer_investigations",
                    "ensemble_investigations",
                ]
            ),
        },
        "local_installations": installs,
        "cross_platform_ci": (
            "Separate fresh remote run required; local Linux installs "
            "are not cross-platform evidence."
        ),
        "limitations": [
            "Original hard-sphere strict comparison remains false; "
            "no new policy is applied here.",
            "Density agreement covers declared grid integration, "
            "not physical choice of density threshold.",
            "Ensemble study recomputes references against frozen outputs; "
            "installed wheels separately verify all 56 native goldens.",
            "Build-package pip inventory failed because the build environment "
            "has no pip; importlib.metadata supplement retained.",
            "Wheel build used Cargo incremental artifacts whose source hashes "
            "are unchanged; install environments are new.",
            "No performance conclusion follows from these concurrent "
            "validation durations.",
        ],
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    for key, path in lane_paths.items():
        shutil.copyfile(references / path, output / (key + "_summary.json"))
    files = []
    for phase in ["references", "distribution", "installs"]:
        for path in sorted((source / phase).rglob("*")):
            rel = path.relative_to(source)
            if not path.is_file() or set(rel.parts).intersection(
                {"environments", "wheels", "sdist", "__pycache__"}
            ):
                continue
            files.append(path)
    manifest = {
        str(p.relative_to(source)): {"sha256": sha(p), "bytes": p.stat().st_size}
        for p in files
    }
    (output / "raw_evidence_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    archive = output / "raw_evidence.tar.gz"
    with (
        archive.open("wb") as stream,
        gzip.GzipFile(fileobj=stream, mode="wb", filename="", mtime=0) as gz,
        tarfile.open(fileobj=gz, mode="w|") as tar,
    ):
        for path in files:
            data = path.read_bytes()
            info = tarfile.TarInfo(str(path.relative_to(source)))
            info.size = len(data)
            info.mode = 0o644
            tar.addfile(info, io.BytesIO(data))
    (output / "archive_receipt.json").write_text(
        json.dumps(
            {
                "path": archive.name,
                "sha256": sha(archive),
                "bytes": archive.stat().st_size,
                "entries": len(files),
                "raw_bytes": sum(r["bytes"] for r in manifest.values()),
            },
            indent=2,
        )
        + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
