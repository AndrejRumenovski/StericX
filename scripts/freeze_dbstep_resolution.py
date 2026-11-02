#!/usr/bin/env python3
"""Freeze the pre-resolution implementation and all original DBSTEP failures."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASELINE = "f04b156f7b8b31700e1c6ab978e337b00cc734e0"
WITNESS = "SIG-NIHDA-401/0/sphere_radius=2.5"


def record(path):
    path = Path(path).resolve()
    with path.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    return {
        "path": str(path.relative_to(ROOT)),
        "bytes": path.stat().st_size,
        "sha256": digest,
    }


def command(argv):
    result = subprocess.run(argv, cwd=ROOT, capture_output=True, text=True, check=True)
    return {"argv": argv, "stdout": result.stdout, "stderr": result.stderr}


def read(path):
    return json.loads(Path(path).read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if (output / "manifest.json").exists():
        raise FileExistsError("Never overwrite a baseline freeze")
    if command(["git", "rev-parse", "HEAD"])["stdout"].strip() != BASELINE:
        raise ValueError(
            "This baseline must be frozen at the recorded pre-resolution HEAD"
        )
    command(["git", "diff", "--exit-code", BASELINE])
    build_path = args.build.resolve() / "manifest.json"
    build = read(build_path)
    for item in build["sources"]:
        relative = Path(item["path"]).relative_to(build["source"])
        if record(ROOT / relative)["sha256"] != item["sha256"]:
            raise ValueError(
                f"Production source differs from baseline build: {relative}"
            )
    for item in build["binaries"].values():
        if record(item["path"])["sha256"] != item["sha256"]:
            raise ValueError("Baseline executable changed")

    tracked = command(["git", "ls-files", "-z"])["stdout"].split("\0")
    sources, evidence, reports = [], [], []
    for name in filter(None, tracked):
        path = ROOT / name
        if name.startswith(
            ("src/", "python/", "scripts/", "tests/", "examples/", ".github/")
        ) or name in {"Cargo.toml", "Cargo.lock", "pyproject.toml", "uv.lock"}:
            sources.append(record(path))
        if name.startswith("docs/validation/dbstep/"):
            evidence.append(record(path))
        if path.suffix == ".md" and name.startswith(
            (
                "docs/scientific_accuracy_audit/",
                "docs/scientific_remediation/",
                "docs/kraken_residual_forensics/",
                "docs/paton_recommendations/",
            )
        ):
            reports.append(record(path))

    original = ROOT / "docs/validation/dbstep/volume_v2"
    retained = output / "original_failure"
    retained.mkdir(exist_ok=False)
    for name in ("requests.jsonl", "native_stdout.jsonl"):
        lines = [
            line
            for line in (original / name).read_bytes().splitlines(keepends=True)
            if json.loads(line).get("id") == WITNESS
        ]
        if len(lines) != 1:
            raise ValueError("Original witness coverage changed")
        (retained / name).write_bytes(lines[0])
    request = read(retained / "requests.jsonl")
    source = (
        ROOT
        / "docs/paton_recommendations/conformers/sensitivity_v1"
        / "SIG-NIHDA-401/report.json"
    )
    conformer = read(source)["descriptors"]["conformers"][0]
    (retained / "source_conformer.json").write_text(
        json.dumps(conformer, indent=2) + "\n"
    )

    original_summary = read(original / "summary.json")
    failures = []
    for category in ("failed_matched_frames", "failed_matched_summaries"):
        for row in original_summary[category]:
            failures.append(
                {
                    "source": record(original / "summary.json"),
                    "category": category,
                    **row,
                    "absolute_difference": abs(row["residual"]),
                    "relative_difference_to_dbstep": abs(
                        row["residual"] / row["reference"]
                    ),
                    "status_at_freeze": "failed_original_strict_gate",
                }
            )
    profiles = read(ROOT / "docs/validation/dbstep/profiles_v1/reference_raw.json")
    volume_v1 = read(ROOT / "docs/validation/dbstep/volume_v1/summary.json")
    ledger = {
        "strict_numeric_failures": failures,
        "strict_denominator_failures": original_summary["errors"],
        "original_volume_v1_failures": {
            key: volume_v1[key]
            for key in ("errors", "failed_matched_frames", "failed_matched_summaries")
        },
        "profile_repeated_witness": [p for p in profiles if p.get("strict_errors")],
        "retained_setup_failures": [
            {
                "source": "docs/validation/dbstep/OBSERVER_BUILD.md",
                "reason": (
                    "Observer build initially lacked a direct serde dependency; "
                    "corrected build preserved locked crates."
                ),
            },
            {
                "source": "docs/validation/dbstep/DENSITY_RESULTS.md",
                "reason": (
                    "Initial reference controller failed to serialize a NumPy boolean "
                    "after calculation; original attempt retained."
                ),
            },
        ],
        "classification_at_freeze": (
            "Original f32/f64 sphere-membership difference explained previously; "
            "release policy still unresolved; no new correction or pass declared."
        ),
        "no_other_original_numerical_failures": {
            "sterimol": read(ROOT / "docs/validation/dbstep/sterimol_v1/summary.json")[
                "conditional_numerical_gate_pass"
            ],
            "density": read(ROOT / "docs/validation/dbstep/density_v2/summary.json")[
                "admitted"
            ],
            "profiles_preservation": read(
                ROOT / "docs/validation/dbstep/profiles_v1/summary.json"
            )["profile_preservation_pass"],
        },
    }
    (output / "original_failures.json").write_text(json.dumps(ledger, indent=2) + "\n")
    (output / "build_manifest.json").write_bytes(build_path.read_bytes())
    manifest = {
        "kind": "dbstep_resolution_baseline",
        "created_utc": datetime.now(UTC).isoformat(),
        "command": sys.argv,
        "git_commit": BASELINE,
        "git_tree": command(["git", "rev-parse", "HEAD^{tree}"])["stdout"].strip(),
        "rustc": command(["rustc", "-Vv"]),
        "cargo": command(["cargo", "-V"]),
        "research_python": sys.version,
        "build": record(build_path),
        "binaries": build["binaries"],
        "reference_environment": record(output / "dbstep/environment.json"),
        "controller": record(__file__),
        "sources": sources,
        "original_dbstep_evidence": evidence,
        "scientific_report_inventory": reports,
        "witness_id": WITNESS,
        "witness_configuration": request["config"],
        "witness_center": request["center"],
        "witness_donor": request["donor"],
        "witness_neighbors": request["neighbors"],
        "retained_outputs": [record(p) for p in sorted(retained.iterdir())],
        "failure_ledger": record(output / "original_failures.json"),
        "scientific_pass": False,
        "scope": (
            "Immutable pre-resolution source, outputs and failures. "
            "No production correction or new acceptance policy evaluated."
        ),
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(record(output / "manifest.json")))


if __name__ == "__main__":
    main()
