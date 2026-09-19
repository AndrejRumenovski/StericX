#!/usr/bin/env python3
"""Freeze original-CLI screening behavior, including errors, then compare candidates.

This is a correctness oracle, not a benchmark. No output field is normalized or
discarded. Baseline and candidate use the identical immutable paths and models.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import struct
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULT_ROOT = REPO / ".stericx/profiling/optimization/error_oracles"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def verify_inputs(manifest: dict) -> None:
    for entry in manifest["inputs"]:
        path = Path(entry["path"])
        if not path.is_file() or sha(path) != entry["sha256"]:
            raise ValueError(f"frozen oracle input changed: {path}")


def execute(binary: Path, argv: list[str], directory: Path, threads: int) -> dict:
    directory.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ)
    env.update(LC_ALL="C", TZ="UTC", RAYON_NUM_THREADS=str(threads))
    env.pop("STERICX_PROFILE_PATH", None)
    with (
        (directory / "stdout.txt").open("wb") as stdout,
        (directory / "stderr.txt").open("wb") as stderr,
    ):
        result = subprocess.run(
            [str(binary), *argv], cwd=REPO, env=env, stdout=stdout, stderr=stderr
        )
    record = {
        "argv": [str(binary), *argv],
        "environment": {"LC_ALL": "C", "TZ": "UTC", "RAYON_NUM_THREADS": str(threads)},
        "returncode": result.returncode,
        "stdout_sha256": sha(directory / "stdout.txt"),
        "stderr_sha256": sha(directory / "stderr.txt"),
    }
    save(directory / "result.json", record)
    return record


def fingerprint(result: dict) -> dict:
    return {
        key: result[key] for key in ("returncode", "stdout_sha256", "stderr_sha256")
    }


def assert_exact(expected: dict, actual: dict, name: str) -> None:
    if fingerprint(expected) != fingerprint(actual):
        raise ValueError(f"screening behavior changed: {name}")


def validate_coverage(case: dict, result: dict, directory: Path) -> dict:
    """Check that fixtures exercised useful behavior, not merely two equal errors."""
    if result["returncode"] != case["expected_status"]:
        raise ValueError(f"unexpected baseline status for {case['name']}")
    stdout = (directory / "stdout.txt").read_text()
    stderr = (directory / "stderr.txt").read_text()
    if result["returncode"] != 0:
        if not stderr.startswith("error:"):
            raise ValueError(f"missing error message for {case['name']}")
        return {"error": stderr.strip()}
    report = json.loads(stdout)
    if report["screened"] < case.get("minimum_screened", 1):
        raise ValueError(
            f"fixture failed to screen intended valid geometry: {case['name']}"
        )
    if report["skipped"] < case.get("minimum_skipped", 0):
        raise ValueError(
            f"fixture did not exercise intended exclusions: {case['name']}"
        )
    details = "\n".join(exclusion["detail"] for exclusion in report["excluded"])
    for expected in case.get("exclusion_details", []):
        if expected not in details:
            raise ValueError(f"missing expected exclusion {expected!r}: {case['name']}")
    for hit in report["hits"]:
        for key in (
            "rank",
            "predicted_ddg_kcal_mol",
            "descriptors",
            "uncertainty",
            "prediction_interval_low",
            "prediction_interval_high",
            "leverage",
            "domain_verdict",
            "nearest_training_distance",
            "applicability",
        ):
            if key not in hit:
                raise ValueError(f"uncovered output field {key}: {case['name']}")
    if case.get("require_uncertainty", True) and not any(
        hit["uncertainty"] is not None and hit["prediction_interval_low"] is not None
        for hit in report["hits"]
    ):
        raise ValueError(f"fixture did not exercise uncertainty: {case['name']}")
    return {
        "screened": report["screened"],
        "skipped": report["skipped"],
        "returned": report["returned"],
        "exclusion_details": details,
        "ranked_ligands": [hit["ligand"] for hit in report["hits"]],
        "domain_verdicts": [hit["domain_verdict"] for hit in report["hits"]],
    }


def sdf_from_xyz(paths: list[Path]) -> str:
    blocks = []
    for path in paths:
        lines = path.read_text().splitlines()
        count = int(lines[0])
        atoms = []
        for line in lines[2 : 2 + count]:
            element, x, y, z = line.split()
            atoms.append(f"{x} {y} {z} {element} 0 0 0")
        blocks.append(
            "\n".join(
                [
                    path.stem,
                    "  frozen oracle",
                    "original coordinate decimal tokens",
                    f"{count:3d}  0  0  0  0  0            999 V2000",
                    *atoms,
                    "M  END",
                    "$$$$",
                    "",
                ]
            )
        )
    return "".join(blocks)


def make_model(baseline: Path, directory: Path) -> Path:
    """Train the same geometry-only synthetic fixture used in cli_screen.rs."""
    directory.mkdir(parents=True)
    packed = directory / "training.sigpack"
    metadata = directory / "labels.csv"
    with packed.open("wb") as matrix, metadata.open("w", newline="") as labels:
        writer = csv.writer(labels)
        writer.writerow(["Reaction_ID", "Dataset_Split", "Ligand_Group"])
        for index in range(24):
            length = 5.0 + 0.15 * index
            response = 0.5 + 0.6 * length + 0.01 * ((index * 11) % 5)
            matrix.write(
                struct.pack(
                    "=16f", length, 3, 7, -0.4, 1650, 298.15, response, *([0] * 9)
                )
            )
            writer.writerow(
                [
                    f"G{index:02}",
                    "train" if index < 22 else "blind",
                    f"group_{index % 4}",
                ]
            )
    model = directory / "geometry_model.json"
    argv = [
        "fit",
        "--data",
        str(packed),
        "--metadata",
        str(metadata),
        "--output",
        str(directory / "fit_report.json"),
        "--predictions",
        str(directory / "predictions.csv"),
        "--portable-model",
        str(model),
        "--model-id",
        "screening-semantics-geometry",
        "--optimize",
        "maximize",
        "--bootstrap",
        "30",
        "--permutations",
        "30",
        "--max-terms",
        "1",
        "--seed",
        "20260725",
    ]
    result = execute(baseline, argv, directory / "fit_execution", 1)
    if result["returncode"] != 0:
        raise ValueError(f"geometry model preparation failed: {directory}")
    if json.loads(model.read_text())["selected_features"] != ["L_boltz"]:
        raise ValueError("geometry fixture did not select only L_boltz")
    return model


def prepare(args: argparse.Namespace) -> None:
    root, baseline = args.root.resolve(), args.baseline.resolve(strict=True)
    if (root / "manifest.json").is_file():
        verify_inputs(json.loads((root / "manifest.json").read_text()))
        print(f"Existing frozen semantics oracle verified: {root}")
        return
    root.mkdir(parents=True, exist_ok=False)
    inputs = root / "inputs"
    inputs.mkdir()
    geometry_model = make_model(baseline, inputs / "geometry_model")
    electronic_model = inputs / "electronic_model.json"
    shutil.copyfile(
        REPO / "docs/study_001/stericx_portable_model.json", electronic_model
    )
    geometries = inputs / "geometries"
    geometries.mkdir()
    named = {}

    def xyz(name: str, body: str) -> Path:
        path = geometries / f"{name}.xyz"
        lines = body.strip().splitlines()
        path.write_text(f"{len(lines)}\n{name}\n" + "\n".join(lines) + "\n")
        named[name] = path
        return path

    valid = xyz(
        "valid_p", "P 0 0 0\nC 1.4 0 .45\nC -.7 1.212 .45\nC -.7 -1.212 .45\nC 2.8 0 .7"
    )
    xyz(
        "valid_n", "N 0 0 0\nC 1.4 0 .45\nC -.7 1.212 .45\nC -.7 -1.212 .45\nC 2.8 0 .7"
    )
    xyz("missing_donor", "C 0 0 0\nC 1.4 0 .45\nC -.7 1.212 .45\nC -.7 -1.212 .45")
    xyz(
        "ambiguous_donor",
        "P 0 0 0\nC 1.4 0 .45\nC -.7 1.212 .45\nC -.7 -1.212 .45\nP 8 8 8",
    )
    xyz("nontrivalent_donor", "P 0 0 0\nC 1.4 0 .45\nC -.7 1.212 .45")
    xyz("no_heavy_substituent", "P 0 0 0\nH 1.0 0 .3\nH -.5 .866 .3\nH -.5 -.866 .3")
    xyz(
        "planar_normal_fallback",
        "P 0 0 0\nC 1 0 0\nC -.5 .8660254 0\nC -.5 -.8660254 0",
    )
    xyz("collinear_frame", "P 0 0 0\nC 1 0 0\nC -1 0 0\nC 2 0 0")
    xyz("coincident_neighbor", "P 0 0 0\nC 0 0 0\nC 1.4 0 .45\nC -.7 1.212 .45")
    xyz(
        "nonfinite_coordinate",
        "P nan 0 0\nC 1.4 0 .45\nC -.7 1.212 .45\nC -.7 -1.212 .45",
    )
    truncated = geometries / "truncated.xyz"
    truncated.write_text("5\ntruncated\nP 0 0 0\n")
    named["truncated"] = truncated
    for name, original in (
        ("real_small", "data/xyz/SIG-NIHDA-401_9d42bff1.xyz"),
        ("real_large", "data/xyz/SIG-NIHDA-2062_d86f6436.xyz"),
    ):
        target = geometries / f"{name}.xyz"
        shutil.copyfile(REPO / original, target)
        named[name] = target
    for name, sources in (
        ("ensemble_valid", [valid, named["real_small"]]),
        ("ensemble_later_bad_topology", [valid, named["nontrivalent_donor"]]),
        ("ensemble_later_missing_donor", [valid, named["missing_donor"]]),
    ):
        path = geometries / f"{name}.sdf"
        path.write_text(sdf_from_xyz(sources))
        named[name] = path
    malformed_sdf = geometries / "ensemble_later_bad_parse.sdf"
    malformed_sdf.write_text(
        sdf_from_xyz([valid])
        + "bad\nprogram\ncomment\n  3  0  0  0  0  0 999 V2000\n0 0 0 P\n$$$$\n"
    )
    named["ensemble_later_bad_parse"] = malformed_sdf
    valid_names = [
        "valid_p",
        "real_small",
        "real_large",
        "ensemble_valid",
        "planar_normal_fallback",
    ]
    invalid_names = [name for name in named if name not in {*valid_names, "valid_n"}]

    def library(name: str, files: list[str], csv_fallback: bool = False) -> Path:
        if csv_fallback:
            path = inputs / f"{name}.csv"
            with path.open("w", newline="") as output:
                writer = csv.writer(output)
                writer.writerow(
                    ["Reaction_ID", "Ligand_XYZ_Path", "NBO_Charge", "IR_Frequency"]
                )
                for index, file in enumerate(files):
                    geometry = named.get(file, geometries / "does_not_exist.xyz")
                    writer.writerow(
                        [f"candidate_{index:02}_{file}", geometry, "0.74667798", "1650"]
                    )
            return path
        path = inputs / name
        path.mkdir()
        for file in files:
            shutil.copyfile(named[file], path / named[file].name)
        return path

    libraries = {}
    for loader in ("directory", "csv"):
        for collection, files in (
            ("valid", valid_names),
            ("mixed", valid_names + invalid_names),
            ("invalid", ["missing_donor"]),
            ("nonp", ["valid_n"]),
        ):
            if loader == "csv" and collection == "mixed":
                files = [*files, "missing_file"]
            libraries[loader, collection] = library(
                f"{loader}_{collection}", files, loader == "csv"
            )
    # Explicit descriptors bypass fallback, including invalid geometry/configuration.
    explicit = inputs / "explicit_descriptors.csv"
    explicit.write_text(
        "Reaction_ID,Ligand_XYZ_Path,L_boltz,B1_boltz,B5_boltz,NBO_Charge,IR_Frequency\n"
        f"explicit,{named['missing_donor']},5.1,3.0,7.0,0.75,1650\n"
        f"missing_charge,{valid},5.2,3.0,7.0,,1650\n"
    )
    cases = []

    def add(
        name: str,
        model: Path,
        library_path: Path,
        extra: list[str],
        status: int = 0,
        **checks: object,
    ) -> None:
        cases.append(
            {
                "name": name,
                "argv": [
                    "screen",
                    str(model),
                    str(library_path),
                    "--format",
                    "json",
                    *extra,
                ],
                "expected_status": status,
                **checks,
            }
        )

    config_failures = [
        ("radius_zero", "--sphere-radius=0"),
        ("radius_negative", "--sphere-radius=-1"),
        ("radius_nan", "--sphere-radius=NaN"),
        ("density_zero", "--density=0"),
        ("density_infinite", "--density=inf"),
        ("center_zero", "--center-distance=0"),
        ("center_negative", "--center-distance=-1"),
        ("radii_zero", "--radii-scale=0"),
        ("radii_nan", "--radii-scale=NaN"),
        ("empty_grid", "--density=1000000000"),
        ("symmetric_occupied", "--radii-scale=100"),
    ]
    for axis in ("bond", "coordination"):
        for loader in ("directory", "csv"):
            model = geometry_model if loader == "directory" else electronic_model
            common = ["--sterimol-axis", axis]
            for collection in ("valid", "mixed", "invalid", "nonp"):
                checks = {}
                if collection == "mixed":
                    checks["minimum_skipped"] = len(invalid_names) + (loader == "csv")
                    if loader == "directory":
                        checks["exclusion_details"] = [
                            "no P donor",
                            "found 2 P",
                            "not three-coordinate",
                            "no bonded heavy atom",
                            "conformer 1",
                            "collinear with the coordination axis",
                        ]
                add(
                    f"{loader}_{collection}_{axis}",
                    model,
                    libraries[loader, collection],
                    common + (["--donor-element", "N"] if collection == "nonp" else []),
                    2 if collection == "invalid" else 0,
                    **checks,
                )
            for condition, flag in config_failures:
                add(
                    f"{loader}_{condition}_{axis}",
                    model,
                    libraries[loader, "valid"],
                    [*common, flag],
                    2,
                )
            add(
                f"{loader}_custom_config_{axis}",
                model,
                libraries[loader, "valid"],
                [
                    *common,
                    "--sphere-radius",
                    "4.2",
                    "--density",
                    "0.07",
                    "--center-distance",
                    "3.1",
                    "--radii-scale",
                    "0.95",
                ],
            )
            add(
                f"{loader}_rank_top_{axis}",
                model,
                libraries[loader, "valid"],
                [*common, "--top", "2", "--ascending"],
            )
        add(
            f"csv_explicit_bypass_{axis}",
            electronic_model,
            explicit,
            ["--sterimol-axis", axis, "--density=0"],
            minimum_skipped=1,
        )
    results = {}
    for case in cases:
        directory = root / "baseline" / case["name"]
        result = execute(baseline, case["argv"], directory, 1)
        result["coverage"] = validate_coverage(case, result, directory)
        results[case["name"]] = result
        print(f"baseline {case['name']}: status {result['returncode']}", flush=True)
    manifest = {
        "schema_version": 1,
        "baseline_binary": str(baseline),
        "baseline_binary_sha256": sha(baseline),
        "script_sha256": sha(Path(__file__)),
        "comparison": (
            "Exact stdout bytes, stderr bytes, and process exit status; "
            "no normalization or tolerance."
        ),
        "cases": cases,
        "results": results,
        "inputs": [
            {"path": str(path), "sha256": sha(path)}
            for path in sorted(inputs.rglob("*"))
            if path.is_file()
        ],
    }
    save(root / "manifest.json", manifest)
    print(f"Frozen {len(cases)} screening cases: {root / 'manifest.json'}")


def compare(args: argparse.Namespace) -> None:
    root = args.root.resolve()
    manifest_path = root / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    verify_inputs(manifest)
    candidate = args.candidate.resolve(strict=True)
    destination = root / "comparisons" / args.label
    destination.mkdir(parents=True, exist_ok=False)
    summary = {
        "manifest_sha256": sha(manifest_path),
        "binary_sha256": sha(candidate),
        "binary": str(candidate),
        "threads": args.threads,
        "status": "running",
        "results": {},
    }
    try:
        for case in manifest["cases"]:
            result = execute(
                candidate, case["argv"], destination / case["name"], args.threads
            )
            summary["results"][case["name"]] = result
            assert_exact(manifest["results"][case["name"]], result, case["name"])
        summary["status"] = "complete"
    except BaseException as error:
        summary["status"] = "failed"
        summary["error"] = str(error)
        raise
    finally:
        save(destination / "summary.json", summary)
    print(
        f"PASS: all {len(summary['results'])} screening cases are byte-exact "
        f"({args.label})"
    )


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
