#!/usr/bin/env python3
"""Retain final v4 native observations for the existing large-value witnesses."""

import csv
import hashlib
import json
import math
import struct
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parent
BUILD = REPO / ".stericx/scientific_remediation/validated_build_v4"
OUT = HERE / "final_v4"


def record(path):
    return {
        "path": str(path.resolve()),
        "bytes": path.stat().st_size,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


def f32(value):
    return struct.unpack("<f", struct.pack("<f", value))[0]


def main():
    OUT.mkdir()
    build = json.loads((BUILD / "manifest.json").read_text())
    native = Path(build["binaries"]["native"]["path"])
    if record(native) != build["binaries"]["native"]:
        raise ValueError("frozen native identity mismatch")
    inputs = [
        BUILD / "manifest.json",
        native,
        Path(__file__),
        HERE / "ungrouped.csv",
        HERE / "large_varying.csv",
        *sorted((HERE / "inputs").rglob("*.sdf")),
    ]
    bindings = [record(path) for path in inputs]
    commands = []

    def run(name, arguments):
        argv = [str(native), *map(str, arguments)]
        result = subprocess.run(argv, capture_output=True, check=False)
        stdout, stderr = OUT / (name + ".stdout"), OUT / (name + ".stderr")
        stdout.write_bytes(result.stdout)
        stderr.write_bytes(result.stderr)
        commands.append(
            {
                "name": name,
                "argv": argv,
                "returncode": result.returncode,
                "stdout": record(stdout),
                "stderr": record(stderr),
            }
        )
        return result

    db = run(
        "database",
        [
            "db",
            "build",
            "--source",
            HERE / "inputs",
            "--output",
            OUT / "database.csv",
            "--group-by-parent",
            "--label-from",
            "parent",
            "--sterimol-axis",
            "coordination",
            "--density",
            "1",
        ],
    )
    search_args = [
        "--similar-to",
        HERE / "inputs/ligand/conformer_0.sdf",
        "--features",
        "sterimol_l",
        "--sterimol-axis",
        "coordination",
        "--density",
        "1",
        "--format",
        "json",
    ]
    identical = run(
        "search_identical",
        ["search", "--database", HERE / "ungrouped.csv", *search_args],
    )
    varying = run(
        "search_varying",
        ["search", "--database", HERE / "large_varying.csv", *search_args],
    )
    with (HERE / "ungrouped.csv").open() as stream:
        reference = list(csv.DictReader(stream))
    with (OUT / "database.csv").open() as stream:
        grouped = list(csv.DictReader(stream))
    checks = {
        "database_success": db.returncode == 0,
        "four_source_observations": len(reference) == 4,
        "one_group": len(grouped) == 1,
        "grouped_exact_mean": all(
            row["sterimol_l"] == grouped[0]["sterimol_l"] for row in reference
        ),
        "grouped_finite": math.isfinite(float(grouped[0]["sterimol_l"])),
        "identical_axis_rejected": identical.returncode != 0,
        "identical_reason": b"constant" in identical.stderr,
        "identical_no_successful_json": not identical.stdout,
        "varying_success": varying.returncode == 0,
    }
    expected = {}
    query = f32(float(reference[0]["sterimol_l"]))
    for hit in json.loads(varying.stdout)["hits"]:
        candidate = {"low": 2.0**126, "high": 3 * 2.0**126}[hit["ligand"]]
        # Two equally weighted endpoints: mean=2**127, population SD=2**126.
        distance = f32(abs(query - candidate) / 2.0**126)
        expected[hit["ligand"]] = distance
        checks["distance_exact_" + hit["ligand"]] = hit["distance"] == distance
    checks["exact_hit_order"] = [
        hit["ligand"] for hit in json.loads(varying.stdout)["hits"]
    ] == ["high", "low"]
    checks["bindings_unchanged"] = bindings == [record(path) for path in inputs]
    receipt = {
        "scope": "final v4 same-input DB/search finite-value correction witnesses",
        "bindings": bindings,
        "commands": commands,
        "checks": checks,
        "expected_distances_f32": expected,
        "all_checks_passed": all(checks.values()),
        "no_tolerance": True,
        "limitations": "Focused arithmetic checks, not global scientific validation.",
        "files": [record(path) for path in sorted(OUT.iterdir()) if path.is_file()],
    }
    (OUT / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    if not receipt["all_checks_passed"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
