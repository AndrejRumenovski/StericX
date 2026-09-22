#!/usr/bin/env python3
"""Prepare new corrected-science profiling inputs; never run benchmarks.

Native parsing, fitting and a fresh database build are preparation. Historical
inputs/results and the old timing harness remain unchanged. Any native failure
or skipped selected geometry aborts, preserving evidence without a manifest.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

import profile_stericx as profile
import replay_scientific_remediation as replay

REPO = Path(__file__).resolve().parents[1]
AUDIT = REPO / "docs/scientific_accuracy_audit"
AUDIT_SHA256 = "c31b655641a352402b14dc2d4535261a9de9c9004979c5e571d78e45cac82f03"


def record(path: Path) -> dict:
    return {
        "path": str(path.resolve()),
        "bytes": path.stat().st_size,
        "sha256": profile.digest(path),
    }


def verify(item: dict) -> Path:
    path = Path(item["path"])
    if (
        not path.is_file()
        or path.stat().st_size != item["bytes"]
        or profile.digest(path) != item["sha256"]
    ):
        raise ValueError(f"bound input changed: {path}")
    return path


def file_records(value):
    if isinstance(value, dict):
        if {"path", "bytes", "sha256"} <= value.keys():
            yield value
        else:
            for child in value.values():
                yield from file_records(child)
    elif isinstance(value, list):
        for child in value:
            yield from file_records(child)


def same_content(left: dict, right: dict) -> bool:
    return all(left[key] == right[key] for key in ("bytes", "sha256"))


def copy_bound(item: dict, destination: Path) -> dict:
    source = verify(item)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("xb") as stream:
        stream.write(source.read_bytes())
    shutil.copymode(source, destination)
    copied = record(destination)
    if not same_content(item, copied):
        raise ValueError("copied input identity mismatch")
    return copied


def metadata_policy(fields, original, corrected_fields, corrected) -> None:
    """Check exact source tokens and explicit corrections, independent of paths."""
    additions = [
        "Historical_Recorded_Temp_K",
        "Conformer_Population_Source",
        "Target_Source_Note",
    ]
    if corrected_fields != [*fields, *additions] or len(original) != len(corrected):
        raise ValueError("corrected metadata schema/cardinality mismatch")
    for old, new in zip(original, corrected, strict=True):
        for key in fields:
            if key not in {"Temp_K", "Ligand_XYZ_Path", "Conformer_XYZ_Paths"}:
                if old[key] != new[key]:
                    raise ValueError(f"historical scientific value changed: {key}")
        expected = {
            "Temp_K": "353.15",
            "Historical_Recorded_Temp_K": "298.15",
            "Conformer_Population_Source": (
                "Supplied historical MMFF94 populations at 298.15 K; "
                "retained without reweighting"
            ),
            "Target_Source_Note": (
                "Published ddG_abs retained; source ee/DDG inconsistency "
                "remains unresolved"
                if old["Source_ID"] == "2064"
                else "Published ddG_abs retained"
            ),
        }
        if old["Temp_K"] != "298.15" or any(
            new[key] != value for key, value in expected.items()
        ):
            raise ValueError("response/population/source limitation policy changed")


def admit_metadata(input_manifest: Path, sealed) -> tuple[dict, dict, list[dict]]:
    data = json.loads(input_manifest.read_text())
    if (
        data.get("kind") != "Ni_hDA_response_temperature_metadata_correction"
        or data.get("temperature_k") != 353.15
        or data.get("conformer_population_source_temperature_k") != 298.15
        or data.get("conformer_populations_recomputed") is not False
        or data.get("not_a_new_experiment") is not True
        or data.get("pre_post_input_identities_verified") is not True
        or data.get("authority") != ["C21", "M67", "M68"]
        or data["audit_manifest"]["sha256"] != AUDIT_SHA256
    ):
        raise ValueError(
            "expected corrected response metadata and retained populations"
        )
    bindings = list(file_records(data))
    for item in bindings:
        verify(item)
    for key, path in (
        ("source", "frozen/repository/data/reactions_raw.csv"),
        ("primary_target_table", "frozen/repository/data/official/ni_hda_kraken.csv"),
        (
            "primary_target_provenance",
            "frozen/repository/data/official/provenance.json",
        ),
        ("audit_claims", "claims.json"),
    ):
        source = sealed(AUDIT / path)
        if not same_content(source, data[key]):
            raise ValueError(f"metadata {key} differs from sealed source")
        bindings.append(source)
    evidence = [
        sealed(AUDIT / name)
        for name in (
            "models/sources/correction_si.pdf",
            "models/sources/correction_si.txt",
            "models/inputs/corrected_ni_hda_tableS3_transcription_v2.json",
            "kinetics/results/ni_hda_target_temperature.csv",
        )
    ]
    if sorted((r["bytes"], r["sha256"]) for r in evidence) != sorted(
        (r["bytes"], r["sha256"]) for r in data["temperature_primary_evidence"]
    ):
        raise ValueError("response temperature evidence differs from sealed source")
    bindings.extend(evidence)
    fields, old = read_rows(Path(data["source"]["path"]))
    corrected_fields, corrected = read_rows(Path(data["output"]["path"]))
    metadata_policy(fields, old, corrected_fields, corrected)
    geometry = {}
    for before, after in zip(old, corrected, strict=True):
        for key in ("Ligand_XYZ_Path", "Conformer_XYZ_Paths"):
            old_paths, new_paths = before[key].split(";"), after[key].split(";")
            if len(old_paths) != len(new_paths):
                raise ValueError("source conformer count/order changed")
            for old_name, new_name in zip(old_paths, new_paths, strict=True):
                role = (REPO / "data" / old_name).resolve().relative_to(REPO)
                source = sealed(AUDIT / "frozen/repository" / role)
                supplied = record(Path(new_name))
                if not same_content(source, supplied):
                    raise ValueError(
                        "supplied geometry differs from exact source bytes"
                    )
                path = supplied["path"]
                if path in geometry and geometry[path]["role"] != str(role):
                    raise ValueError("ambiguous source geometry role")
                geometry[path] = {"record": supplied, "role": str(role)}
                bindings.extend([source, supplied])
    declared = {item["path"]: item for item in data["geometry_inputs"]}
    if declared != {path: value["record"] for path, value in geometry.items()}:
        raise ValueError("metadata geometry inventory differs from consumed inputs")
    return data, geometry, bindings


def eligible_ensembles(rows: list[dict]) -> tuple[dict[int, list[dict]], list[dict]]:
    """Whole-ensemble admission based on recorded source topology, before execution."""
    groups: dict[int, list[dict]] = defaultdict(list)
    for row in rows:
        groups[row["molecule_id"]].append(row)
    selected, excluded = {}, []
    for molecule, members in sorted(groups.items()):
        reasons = set()
        for row in members:
            if row.get("status") != "submitted":
                reasons.add(row.get("status", "unsubmitted"))
            donors = row.get("phosphorus_atoms", [])
            if len(donors) != 1:
                reasons.add("requires_exactly_one_phosphorus")
            elif donors[0].get("degree") != 3:
                reasons.add("requires_phosphorus_degree_three")
        if reasons:
            excluded.append(
                {
                    "molecule_id": molecule,
                    "reasons": sorted(reasons),
                    "inventory_rows": members,
                }
            )
        else:
            selected[molecule] = sorted(members, key=lambda row: row["conformer_id"])
    return selected, excluded


def read_rows(path: Path) -> tuple[list[str], list[dict]]:
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        return list(reader.fieldnames or []), list(reader)


def write_rows(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def prepare(args: argparse.Namespace) -> None:
    root = args.root.resolve()
    if root.is_relative_to(AUDIT):
        raise ValueError("historical audit is immutable")
    work = root / "workloads"
    if work.exists():
        raise FileExistsError(
            f"choose a new destination; existing evidence retained: {work}"
        )
    build = args.build.resolve()
    build_manifest = build / "manifest.json"
    build_data, _ = replay.verify_build(build_manifest)
    binary = verify(build_data["binaries"]["native"])
    input_manifest = args.inputs.resolve() / "manifest.json"
    seal = AUDIT / "manifest_final.json"
    if profile.digest(seal) != AUDIT_SHA256:
        raise ValueError("original audit seal changed")
    inventory = {item["path"]: item for item in json.loads(seal.read_text())["files"]}

    def sealed(path: Path) -> dict:
        item = inventory[str(path.relative_to(AUDIT))]
        bound = {**item, "path": str(path)}
        verify(bound)
        return bound

    input_data, geometry, metadata_bindings = admit_metadata(input_manifest, sealed)
    csv_source = verify(input_data["output"])
    inventory_path = AUDIT / "kraken/prepared_all/inventory.jsonl"
    sealed(inventory_path)
    corpus = [json.loads(line) for line in inventory_path.read_text().splitlines()]
    selected, excluded = eligible_ensembles(corpus)
    if not selected:
        raise ValueError("no source-topology eligible complete ensembles")
    work.mkdir(parents=True)
    sources = work / "sources"
    generated = work / "generated"
    commands_dir = work / "preparation_commands"
    for path in (sources, generated, commands_dir, work / "helpers"):
        path.mkdir()
    helper_paths = [
        Path(__file__),
        Path(profile.__file__),
        REPO / "scripts/profile_child.c",
        *replay.executing_helpers(),
    ]
    helper_bindings = [record(path) for path in helper_paths]
    bindings = [
        record(build_manifest),
        record(binary),
        record(input_manifest),
        record(csv_source),
        record(seal),
        record(inventory_path),
        record(build_manifest.with_suffix(".json.sha256")),
        *helper_bindings,
        *metadata_bindings,
    ]
    commands = []
    copies = []
    try:
        for item in helper_bindings:
            copied = copy_bound(item, work / "helpers" / Path(item["path"]).name)
            copies.append({"source": item, "copy": copied})
        native = work / "helpers/stericx"
        copies.append(
            {"source": record(binary), "copy": copy_bound(record(binary), native)}
        )
        frozen_paths = {}
        for path, value in geometry.items():
            target = sources / value["role"]
            copies.append(
                {"source": value["record"], "copy": copy_bound(value["record"], target)}
            )
            frozen_paths[path] = target
        fields, rows = read_rows(csv_source)
        if len(rows) != 11 or sum(r["Dataset_Split"] == "train" for r in rows) != 10:
            raise ValueError(
                "expected original 10-train/1-blind eleven-record partition"
            )
        for row in rows:
            for key in ("Ligand_XYZ_Path", "Conformer_XYZ_Paths"):
                row[key] = ";".join(
                    str(frozen_paths[str(Path(name).resolve())])
                    for name in row[key].split(";")
                )
        corrected_csv = generated / "reactions.csv"
        write_rows(corrected_csv, fields, rows)
        xyz_paths = sorted(p for p in frozen_paths.values() if p.parent.name == "xyz")
        conformers = sorted(p for p in frozen_paths.values() if "conformers" in p.parts)
        if len(xyz_paths) != 11 or len(conformers) != 56:
            raise ValueError("audited source geometry count differs from 11+56")
        by_size = sorted(
            xyz_paths, key=lambda p: (int(p.read_text().splitlines()[0]), str(p))
        )
        batch = generated / "batch"
        batch.mkdir()
        batch_paths = []
        for index in range(args.descriptor_count):
            target = batch / f"ligand_{index:05d}.xyz"
            shutil.copyfile(conformers[index % len(conformers)], target)
            batch_paths.append(target)
        csv_paths = {}
        for name, count in [
            ("parse_ensembles", args.parse_count),
            ("screen", args.screen_count),
        ]:
            expanded = []
            for index in range(count):
                row = dict(rows[index % len(rows)])
                row["Reaction_ID"] += f"-CORRECTED-PROFILE-{index:06d}"
                expanded.append(row)
            csv_paths[name] = generated / f"{name}.csv"
            write_rows(csv_paths[name], fields, expanded)

        # Preserve exact original SDF payloads and their explicit source bonds.
        kraken = sources / "kraken"
        source_sdfs = {}
        for molecule, members in selected.items():
            directory = kraken / str(molecule)
            directory.mkdir(parents=True)
            source_sdfs[molecule] = []
            for row in members:
                original = AUDIT / "kraken" / row["sdf_path"]
                bound = sealed(original)
                if bound["sha256"] != row["sdf_sha256"]:
                    raise ValueError("inventory/source SDF identity mismatch")
                target = directory / f"{row['conformer_id']}.sdf"
                target.symlink_to(original.resolve())
                source_sdfs[molecule].append(target)
                bindings.append(bound)
        selection = {
            "rule": (
                "all available complete ensembles with exactly one P of degree 3 "
                "in every conformer; no adaptive removal after native execution"
            ),
            "selected_molecule_ids": list(selected),
            "selected_molecules": len(selected),
            "selected_conformers": sum(map(len, selected.values())),
            "excluded": excluded,
            "topology_oracle_scope": (
                "separate full unfiltered oracle still includes multiple donors "
                "and unsupported cases"
            ),
        }
        profile.write_json(work / "corpus_selection.json", selection)
        ensemble_id = max(selected, key=lambda mid: (len(selected[mid]), -mid))
        ensemble = generated / "ensemble.sdf"
        with ensemble.open("xb") as stream:
            for path in source_sdfs[ensemble_id]:
                data = path.read_bytes()
                if data.count(b"$$$$") != 1 or not data.rstrip().endswith(b"$$$$"):
                    raise ValueError(
                        "expected exactly one delimited original SDF conformer"
                    )
                stream.write(data)
                if not data.endswith(b"\n"):
                    stream.write(b"\n")
        profile.write_json(
            work / "started.json",
            {
                "created_utc": profile.now(),
                "build": record(build_manifest),
                "input_manifest": record(input_manifest),
                "selection": selection,
                "executing_helpers": helper_bindings,
                "copied_inputs": copies,
                "preparation_environment": {
                    "LC_ALL": "C",
                    "TZ": "UTC",
                    "OMP_NUM_THREADS": "1",
                    "RAYON_NUM_THREADS": str(args.threads),
                },
                "purpose": "preparation, no timings or optimization acceptance",
            },
        )

        def native_command(argv: list[str]) -> bytes:
            for item in file_records(copies):
                verify(item)
            for item in helper_bindings:
                verify(item)
            result = subprocess.run(
                [str(native), *argv],
                capture_output=True,
                env={
                    **os.environ,
                    "LC_ALL": "C",
                    "TZ": "UTC",
                    "OMP_NUM_THREADS": "1",
                    "RAYON_NUM_THREADS": str(args.threads),
                },
            )
            index = len(commands)
            stdout = commands_dir / f"{index:02}.stdout"
            stderr = commands_dir / f"{index:02}.stderr"
            stdout.write_bytes(result.stdout)
            stderr.write_bytes(result.stderr)
            item = {
                "argv": [str(native), *argv],
                "returncode": result.returncode,
                "stdout": record(stdout),
                "stderr": record(stderr),
            }
            commands.append(item)
            profile.write_json(commands_dir / f"{index:02}.json", item)
            if result.returncode or b"skipped " in result.stderr:
                raise RuntimeError(
                    f"native preparation failed/skipped inputs: {argv[0]}"
                )
            return result.stdout

        packed = generated / "reactions.sigpack"
        native_command(
            [
                "parse",
                "--csv",
                str(corrected_csv),
                "--xyz-dir",
                str(sources / "data"),
                "--output",
                str(packed),
            ]
        )
        if packed.stat().st_size != 11 * 64:
            raise ValueError("fresh parse did not write eleven records")
        model = generated / "model.json"
        native_command(
            [
                "fit",
                "--data",
                str(packed),
                "--metadata",
                str(corrected_csv),
                "--output",
                str(generated / "fit-report.json"),
                "--predictions",
                str(generated / "frozen-predictions.csv"),
                "--portable-model",
                str(model),
                "--model-id",
                "corrected-profile-ni-hda",
                "--descriptor-aggregation",
                "supplied_weight_mean",
                "--response-temp-k",
                "353.15",
                "--response-sign-convention",
                "Published absolute ddG; no R/S assignment",
                "--optimize",
                "maximize",
                "--bootstrap",
                "1000",
                "--permutations",
                "500",
            ]
        )
        model_data = json.loads(model.read_text())
        if (
            model_data["training_count"] != 10
            or model_data["schema_version"] != 3
            or model_data["descriptor_aggregation"] != "supplied_weight_mean"
        ):
            raise ValueError("fresh model contract mismatch")
        native_command(["model", "validate", str(model), "--format", "json"])
        database = generated / "search_database.csv"
        native_command(
            [
                "db",
                "build",
                "--source",
                str(kraken),
                "--output",
                str(database),
                "--group-by-parent",
                "--label-from",
                "parent",
                "--sterimol-axis",
                "coordination",
            ]
        )
        db_manifest = json.loads(database.with_suffix(".manifest.json").read_text())
        if (
            db_manifest["ligands"] != len(selected)
            or db_manifest["geometries_skipped"] != 0
            or db_manifest["geometries_featurized"] != selection["selected_conformers"]
        ):
            raise ValueError("database omitted part of the predeclared corpus")
        inference = generated / "inference.sigpack"
        profile.repeat_sigpack(packed, inference, args.inference_count)
        weights = generated / "weights.json"
        profile.write_json(weights, model_data["weights"])
        workloads = {}

        def descriptors(name: str, paths: list[Path], **extra: object) -> None:
            workloads[name] = {
                "argv": [
                    "descriptors",
                    *map(str, paths),
                    "--sterimol-axis",
                    "coordination",
                    "--format",
                    "json",
                ],
                "kind": "descriptors",
                "expected_records": len(paths),
                **extra,
            }

        descriptors("single_small", [by_size[0]])
        descriptors("single_large", [by_size[-1]])
        descriptors(
            "ensemble_sdf", [ensemble], expected_conformers=len(selected[ensemble_id])
        )
        descriptors("conformers_56", conformers)
        descriptors("descriptors_10000", batch_paths)
        workloads["parse_ensembles_1000"] = {
            "argv": [
                "parse",
                "--csv",
                str(csv_paths["parse_ensembles"]),
                "--xyz-dir",
                str(sources / "data"),
                "--output",
                "{run_dir}/parsed.sigpack",
            ],
            "kind": "parse",
            "expected_records": args.parse_count,
        }
        workloads["predict_1000000"] = {
            "argv": ["predict", "--data", str(inference), "--weights", str(weights)],
            "kind": "predict",
            "expected_records": args.inference_count,
            "exactness_scope": (
                "CLI preview/aggregate only; full per-prediction independent "
                "gate remains required"
            ),
        }
        workloads["search_database"] = {
            "argv": [
                "search",
                "--similar-to",
                str(by_size[-1]),
                "--database",
                str(database),
                "--sterimol-axis",
                "coordination",
                "--top",
                str(len(selected)),
                "--format",
                "json",
            ],
            "kind": "search",
            "expected_records": len(selected),
        }
        workloads["screen_1000"] = {
            "argv": [
                "screen",
                str(model),
                "--library",
                str(csv_paths["screen"]),
                "--top",
                str(args.screen_count),
                "--temperature",
                "353.15",
                "--format",
                "json",
            ],
            "kind": "screen",
            "expected_records": args.screen_count,
        }
        conformer_root = sources / "data/conformers"
        workloads["db_build_ensembles"] = {
            "argv": [
                "db",
                "build",
                "--source",
                str(conformer_root),
                "--output",
                "{run_dir}/database.csv",
                "--group-by-parent",
                "--label-from",
                "parent",
                "--sterimol-axis",
                "coordination",
            ],
            "kind": "db",
            "expected_records": len({p.parent for p in conformers}),
            "expected_geometries": 56,
        }
        for item in bindings:
            verify(item)
        for item in file_records(copies):
            verify(item)
        final_build, _ = replay.verify_build(build_manifest)
        if final_build != build_data:
            raise ValueError("build receipt changed during preparation")
        # Preserve symlink *access paths* so the timing harness validates what it reads.
        files = [
            {
                "path": str(path),
                "bytes": path.stat().st_size,
                "sha256": profile.digest(path),
            }
            for path in sorted(work.rglob("*"))
            if path.is_file()
        ]
        manifest = {
            "schema_version": 1,
            "kind": "corrected_science_profile_inputs",
            "created_utc": profile.now(),
            "repo": str(REPO),
            "root": str(root),
            "files": files,
            "sources": bindings,
            "executing_helpers": helper_bindings,
            "copied_inputs": copies,
            "full_build_proof_verified_before_after": True,
            "metadata_admission": (
                "Exact sealed source scientific tokens and geometry bytes; only "
                "declared response metadata correction and path relocation allowed"
            ),
            "workloads": workloads,
            "preparation_commands": commands,
            "build_manifest": record(build_manifest),
            "input_manifest": record(input_manifest),
            "corpus_selection": record(work / "corpus_selection.json"),
            "sdf_sources": [record(p) for p in source_sdfs[ensemble_id]],
            "ensemble_molecule_id": ensemble_id,
            "diversity": {
                "distinct_ni_hda_ligands": 11,
                "distinct_ni_hda_conformers": 56,
                "search_molecules": len(selected),
                "search_conformers": selection["selected_conformers"],
                "scientific_scope": (
                    "Cyclic throughput repetition is not new chemistry; population "
                    "source 298.15 K retained, response 353.15 K. DB file means are "
                    "unweighted. No comparability claim to old invalid inputs."
                ),
            },
            "defaults": {
                "rayon_threads": 1,
                "affinity": [2],
                "warmups": 1,
                "repetitions": 7,
            },
            "stdout_keys_excluded_from_comparison": sorted(
                profile.VOLATILE_STDOUT_KEYS
            ),
            "db_manifest_keys_excluded_from_comparison": ["build_seconds"],
            "cache_policy": (
                "Sequential warm launches only when separately authorized; "
                "no benchmarks executed in preparation"
            ),
            "scientific_acceptance": (
                "pending independent final gates; input preparation does not "
                "establish scientific validity"
            ),
        }
        profile.write_json(work / "manifest.json", manifest)
        profile.verify_manifest(manifest)
        print(
            json.dumps(
                {
                    "manifest": str(work / "manifest.json"),
                    "workloads": len(workloads),
                    "search_molecules": len(selected),
                    "search_conformers": selection["selected_conformers"],
                }
            )
        )
    except Exception as error:
        profile.write_json(
            work / "preparation_failed.json",
            {
                "error": repr(error),
                "commands": commands,
                "bindings": bindings,
                "scientific_acceptance": False,
                "policy": "preserved; no adaptive exclusion or overwrite",
            },
        )
        raise


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--descriptor-count", type=int, default=10_000)
    parser.add_argument("--parse-count", type=int, default=1_000)
    parser.add_argument("--inference-count", type=int, default=1_000_000)
    parser.add_argument("--screen-count", type=int, default=1_000)
    args = parser.parse_args()
    if (
        min(
            args.threads,
            args.descriptor_count,
            args.parse_count,
            args.inference_count,
            args.screen_count,
        )
        <= 0
    ):
        parser.error("counts and threads must be positive")
    prepare(args)


if __name__ == "__main__":
    main()
