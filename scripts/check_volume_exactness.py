#!/usr/bin/env python3
"""Freeze complete v2 packed records and all nine unrounded conformer volumes.

This is a correctness oracle, never a performance measurement. Artifacts are
compared byte for byte. Raw stdout remains available; only the explicitly named
timing, RSS and destination-path lines are excluded from its comparison.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import struct
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULT_ROOT = REPO / ".stericx/profiling/optimization/volume_oracles"
FIELDS = (
    "vbur",
    "percent_vbur",
    "qvbur_min",
    "qvbur_max",
    "max_delta_qvbur",
    "ovbur_min",
    "ovbur_max",
    "near_vbur",
    "far_vbur",
)
VOLATILE = {
    "binary_export_ms",
    "total_ms",
    "throughput_conformers_per_second",
    "rss_start_bytes",
    "rss_end_bytes",
    "rss_delta_bytes",
    "sigpack_v2_output",
    "per_conformer_output",
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def normalized_stdout(raw: bytes) -> bytes:
    excluded = {key.encode() for key in VOLATILE}
    return b"".join(
        line
        for line in raw.splitlines(keepends=True)
        if b"=" not in line or line.split(b"=", 1)[0] not in excluded
    )


def validate_artifacts(directory: Path, records: int, conformers: int) -> None:
    packed = (directory / "volumes.sigpack").read_bytes()
    if len(packed) != 64 + records * 128:
        raise ValueError("incomplete or invalid packed v2 output")
    header = struct.unpack("=8sIIQII32s", packed[:64])
    expected = (b"SIGPKV2\0", 2, 0x01020304, records, 128, 32, bytes(32))
    if header != expected or len(packed) != 64 + records * 128:
        raise ValueError("incomplete or invalid packed v2 output")
    if not all(
        math.isfinite(v) for v in struct.unpack(f"={records * 32}f", packed[64:])
    ):
        raise ValueError("non-finite packed numeric output")
    with (directory / "conformers.csv").open(newline="") as stream:
        reader = csv.DictReader(stream)
        if not set(FIELDS).issubset(reader.fieldnames or ()):
            raise ValueError("missing conformer volume fields")
        rows = list(reader)
    if len(rows) != conformers:
        raise ValueError("incomplete conformer audit")
    for row in rows:
        if not all(math.isfinite(float(row[field])) for field in FIELDS):
            raise ValueError("non-finite conformer output")


def fingerprints(directory: Path) -> dict:
    return {
        "packed_sha256": sha(directory / "volumes.sigpack"),
        "conformers_sha256": sha(directory / "conformers.csv"),
        "stderr_sha256": sha(directory / "stderr.txt"),
        "stdout_normalized_sha256": hashlib.sha256(
            normalized_stdout((directory / "stdout.txt").read_bytes())
        ).hexdigest(),
    }


def assert_exact(expected: dict, actual: dict, name: str) -> None:
    if (
        actual["returncode"] != expected["returncode"]
        or actual["fingerprints"] != expected["fingerprints"]
    ):
        raise ValueError(f"complete numeric volume output changed: {name}")


def execute(binary: Path, case: dict, directory: Path, threads: int) -> dict:
    directory.mkdir(parents=True, exist_ok=False)
    argv = [
        str(binary),
        *case["argv"],
        "--output",
        str(directory / "volumes.sigpack"),
        "--per-conformer-output",
        str(directory / "conformers.csv"),
    ]
    environment = dict(
        LC_ALL="C", TZ="UTC", RAYON_NUM_THREADS=str(threads), OMP_NUM_THREADS="1"
    )
    env = dict(os.environ, **environment)
    env.pop("STERICX_PROFILE_PATH", None)
    with (
        (directory / "stdout.txt").open("wb") as out,
        (directory / "stderr.txt").open("wb") as err,
    ):
        result = subprocess.run(argv, cwd=REPO, env=env, stdout=out, stderr=err)
    if result.returncode:
        raise ValueError(f"volume oracle failed: {case['name']}; see {directory}")
    validate_artifacts(directory, case["records"], case["conformers"])
    record = dict(
        argv=argv,
        environment=environment,
        returncode=result.returncode,
        raw_stdout_sha256=sha(directory / "stdout.txt"),
        fingerprints=fingerprints(directory),
    )
    save(directory / "result.json", record)
    return record


def verify_inputs(manifest: dict) -> None:
    for entry in manifest["inputs"]:
        if sha(Path(entry["path"])) != entry["sha256"]:
            raise ValueError(f"frozen volume input changed: {entry['path']}")


def prepare(args: argparse.Namespace) -> None:
    root = args.root.resolve()
    if root.exists():
        raise ValueError(f"refusing to overwrite frozen oracle: {root}")
    root.mkdir(parents=True)
    source = REPO / ".stericx/profiling/workloads/sources/data"
    csv_path = source / "reactions_raw.csv"
    with csv_path.open(newline="") as stream:
        source_rows = list(csv.DictReader(stream))
    count = sum(len(row["Conformer_XYZ_Paths"].split(";")) for row in source_rows)
    if (len(source_rows), count) != (11, 56):
        raise ValueError("expected the frozen 11-reaction, 56-conformer corpus")
    base = ["buried-volume", "--csv", str(csv_path), "--xyz-dir", str(source)]
    cases = []
    for name, config in (
        ("default_56", []),
        ("descriptor_center_56", ["--center-distance", "2.28"]),
        ("dense_56", ["--density", "0.004"]),
        (
            "custom_coarse_56",
            [
                "--sphere-radius",
                "4.2",
                "--density",
                "0.07",
                "--center-distance",
                "3.1",
                "--radii-scale",
                "0.95",
            ],
        ),
    ):
        cases.append(
            dict(
                name=name,
                argv=base + config,
                records=len(source_rows),
                conformers=count,
            )
        )
    # A synthetic explicit-center case exercises a separate production API path.
    inputs = root / "inputs"
    inputs.mkdir()
    xyz = inputs / "phosphine.xyz"
    xyz.write_text(
        "5\nphosphine\nP 0 0 0\nC 1.4 0 0.45\nC -0.7 1.212 0.45\n"
        "C -0.7 -1.212 0.45\nC 2.8 0 0.7\n"
    )
    explicit = inputs / "explicit.csv"
    with explicit.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "Reaction_ID",
                "Ligand_XYZ_Path",
                "Attach_Atom_Idx",
                "Primary_Bond_Vector_Idx",
                "NBO_Charge",
                "IR_Frequency",
                "Temp_K",
                "Exp_ddG_kcal_mol",
                "Conformer_Coordination_Centers_Angstrom",
                "Coordination_Center_Method",
            ]
        )
        for index, center in enumerate(
            ["0,0,-2.1", "0.37,-0.19,-2.28", "-0.81,0.62,-3.1"]
        ):
            writer.writerow(
                [
                    f"EXPLICIT-{index}",
                    str(xyz),
                    0,
                    1,
                    0.8,
                    1650,
                    298.15,
                    1.2,
                    center,
                    "frozen_synthetic_center",
                ]
            )
    cases.append(
        dict(
            name="explicit_centers",
            argv=[
                "buried-volume",
                "--csv",
                str(explicit),
                "--xyz-dir",
                str(inputs),
                "--require-explicit-centers",
            ],
            records=3,
            conformers=3,
        )
    )
    paths = {csv_path, xyz, explicit}
    for row in source_rows:
        paths.update(source / p for p in row["Conformer_XYZ_Paths"].split(";"))
    baseline = args.baseline.resolve(strict=True)
    manifest = dict(
        schema_version=1,
        baseline_binary=str(baseline),
        baseline_binary_sha256=sha(baseline),
        script_sha256=sha(Path(__file__)),
        cases=cases,
        results={},
        normalized_stdout_exclusions=sorted(VOLATILE),
        exactness_scope=(
            "All bytes of all 32 f32 values per v2 record and all nine unrounded "
            "f32 per-conformer fields; no numeric tolerance."
        ),
        inputs=[dict(path=str(p), sha256=sha(p)) for p in sorted(paths)],
    )
    for case in cases:
        manifest["results"][case["name"]] = execute(
            baseline, case, root / "baseline" / case["name"], 1
        )
    save(root / "manifest.json", manifest)
    print(f"Frozen {len(cases)} complete-volume cases: {root}")


def compare(args: argparse.Namespace) -> None:
    root = args.root.resolve()
    manifest = json.loads((root / "manifest.json").read_text())
    verify_inputs(manifest)
    binary = args.candidate.resolve(strict=True)
    directory = root / "comparisons" / args.label
    directory.mkdir(parents=True, exist_ok=False)
    summary = dict(
        status="running",
        binary=str(binary),
        binary_sha256=sha(binary),
        manifest_sha256=sha(root / "manifest.json"),
        threads=args.threads,
        results={},
    )
    try:
        for case in manifest["cases"]:
            frozen = root / "baseline" / case["name"]
            validate_artifacts(frozen, case["records"], case["conformers"])
            if (
                fingerprints(frozen)
                != manifest["results"][case["name"]]["fingerprints"]
            ):
                raise ValueError(f"frozen volume output changed: {case['name']}")
            result = execute(binary, case, directory / case["name"], args.threads)
            summary["results"][case["name"]] = result
            assert_exact(manifest["results"][case["name"]], result, case["name"])
        summary["status"] = "complete"
    except BaseException as error:
        summary.update(status="failed", error=str(error))
        raise
    finally:
        save(directory / "summary.json", summary)
    print(f"PASS: complete v2 packed records and conformer volumes ({args.label})")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    actions = parser.add_subparsers(dest="action", required=True)
    preparation = actions.add_parser("prepare")
    preparation.add_argument("--baseline", type=Path, required=True)
    comparison = actions.add_parser("compare")
    comparison.add_argument("--candidate", type=Path, required=True)
    comparison.add_argument("--label", required=True)
    comparison.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    (prepare if args.action == "prepare" else compare)(args)


if __name__ == "__main__":
    main()
