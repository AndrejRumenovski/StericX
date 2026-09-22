#!/usr/bin/env python3
"""Capture corrected SUT observations separately from immutable audit evidence.

This orchestrates measurements, not a scientific PASS. Intentional changes and
independent residuals require review. Every phase refuses an existing destination.
The original reference programs are copied unchanged; new conventions require
separately identified reference analysis, never edits to those programs.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
from itertools import zip_longest
from pathlib import Path

import check_current_scientific_equivalence as oracle
import recheck_current_scientific_references as refs

REPO = oracle.REPO
AUDIT = oracle.DEFAULT_AUDIT


def new_output(path: Path) -> Path:
    path = path.resolve()
    if path.is_relative_to(AUDIT.resolve()):
        raise ValueError("The original scientific audit is immutable")
    path.mkdir(parents=True, exist_ok=False)
    return path


def manifest(path: Path, value: dict) -> None:
    oracle.write_new(path, value)
    path.with_suffix(path.suffix + ".sha256").write_text(oracle.sha(path) + "\n")


def load_manifest(path: Path) -> dict:
    if (
        oracle.sha(path)
        != path.with_suffix(path.suffix + ".sha256").read_text().strip()
    ):
        raise ValueError(f"Manifest identity changed: {path}")
    return oracle.read(path)


def executing_helpers() -> list[Path]:
    return [Path(__file__), Path(oracle.__file__), Path(refs.__file__)]


def verify_build(path: Path) -> tuple[dict, dict]:
    built = load_manifest(path)
    if built.get("kind") != "scientific_remediation_build":
        raise ValueError("Expected a remediation build receipt")
    for item in (
        built["sources"]
        + built["python_sources"]
        + built["adapter_files"]
        + built["helpers"]
        + list(built["binaries"].values())
        + [built["oracle_manifest"], built["audit_manifest"]]
    ):
        oracle.verify(item)
    frozen = {Path(item["path"]).name: item for item in built["helpers"]}
    for path in executing_helpers():
        if oracle.sha(path) != frozen[path.name]["sha256"]:
            raise ValueError(f"Executing helper differs from build snapshot: {path}")
    prepared, _ = oracle.load_oracle(Path(built["oracle_manifest"]["path"]).parent)
    if prepared["evidence_lock"] != built["reference_evidence"]:
        raise ValueError("Build and oracle refer to different scientific evidence")
    if set(built["commands"]) != {"observer", "native"} or any(
        command["returncode"] != 0 for command in built["commands"].values()
    ):
        raise ValueError("Build did not complete both executables")
    return built, prepared


def topology_request(row: dict, dump: bool) -> dict:
    if len(row.get("neighbors", [])) != 3:
        raise ValueError("Topology replay requires three source neighbors")
    return {**row, "op": "geometry_topology", "dump": dump, "include_points": False}


def capture_counts(name: str, request: Path, response: Path) -> dict:
    errors = []
    count = private_frames = private_bins = 0
    for given, received in zip_longest(
        oracle.requests(request), oracle.requests(response)
    ):
        if given is None or received is None:
            raise RuntimeError(f"Missing or extra remediation observation: {name}")
        if given["id"] != received.get("id"):
            raise ValueError(f"Missing/misordered observation: {name}/{given['id']}")
        validation_request = dict(given)
        if validation_request.get("op") == "geometry_topology":
            validation_request["op"] = "geometry"
        oracle.validate_result(validation_request, received)
        count += 1
        if given.get("dump"):
            # Record errors for review, without treating their presence as PASS.
            allowed = {}
            if isinstance(received.get("dump"), dict) and "error" in received["dump"]:
                allowed[given["id"]] = received["dump"]
                errors.append({"id": given["id"], "dump": received["dump"]})
            oracle.validate_dump(
                validation_request,
                received,
                allowed,
                finite=name.startswith("kraken"),
            )
            for frame in received["dump"].get("orientations", []):
                private_frames += 1
                private_bins += len(frame["quadrants"]) + len(frame["octants"])
    return {
        "rows": count,
        "private_errors_for_review": errors,
        "private_frames": private_frames,
        "private_bins": private_bins,
    }


def verify_observations(path: Path) -> tuple[dict, dict]:
    observed = load_manifest(path)
    if observed.get(
        "kind"
    ) != "scientific_remediation_observations" or not observed.get("complete_capture"):
        raise ValueError("Complete corrected SUT observations are required")
    built, prepared = verify_build(oracle.verify(observed["build_manifest"]))
    if observed["binary"] != built["binaries"]["observer"]:
        raise ValueError("Observed executable differs from built observer")
    expected = set(prepared["lanes"]) | {"kraken_topology", "kraken_topology_all_bins"}
    if set(observed["lanes"]) != expected:
        raise ValueError("Incomplete or unexpected remediation lane set")
    for name, lane in observed["lanes"].items():
        request = oracle.verify(lane["requests"])
        response = oracle.verify(lane["stdout"])
        oracle.verify(lane["stderr"])
        if name.startswith("kraken_topology"):
            base = prepared["lanes"]["kraken"]
            if lane["source_requests"] != base["requests"]:
                raise ValueError("Derived topology lane has different source inputs")
            for original, derived in zip_longest(
                oracle.requests(oracle.verify(base["requests"])),
                oracle.requests(request),
            ):
                if original is None or derived != topology_request(
                    original, name.endswith("all_bins")
                ):
                    raise ValueError("Derived topology request changed")
        else:
            base = prepared["lanes"][name]
            if lane["requests"] != base["requests"]:
                raise ValueError(f"Historical request lane changed: {name}")
        counts = capture_counts(name, request, response)
        if lane["returncode"] != 0 or counts["rows"] != base["expected_rows"]:
            raise ValueError(f"Incomplete observation lane: {name}")
        if any(lane[key] != value for key, value in counts.items()):
            raise ValueError(f"Observation coverage receipt differs from rows: {name}")
    return observed, built


def execute(
    argv: list[str], output: Path, *, cwd: Path, env: dict | None = None
) -> dict:
    with (
        (output / "stdout").open("xb") as stdout,
        (output / "stderr").open("xb") as stderr,
    ):
        result = subprocess.run(argv, cwd=cwd, env=env, stdout=stdout, stderr=stderr)
    record = {
        "argv": argv,
        "cwd": str(cwd),
        "returncode": result.returncode,
        "stdout": oracle.file_record(output / "stdout"),
        "stderr": oracle.file_record(output / "stderr"),
    }
    manifest(output / "command.json", record)
    if result.returncode:
        raise RuntimeError(f"Command failed; retained logs: {output}")
    return record


def build(args) -> None:
    prepared, _ = oracle.load_oracle(args.oracle.resolve())
    output = new_output(args.output)
    source = output / "source"
    original = oracle.source_records(REPO)
    python_paths = [REPO / "pyproject.toml", REPO / "uv.lock"]
    for directory in ("scripts", "studies"):
        python_paths.extend(sorted((REPO / directory).glob("*.py")))
    python_original = [oracle.file_record(path) for path in python_paths]
    # The snapshot records all Rust scientific kernels and checked adapters.
    for item in original + python_original:
        relative = Path(item["path"]).relative_to(REPO)
        destination = source / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(oracle.verify(item), destination)
    adapter_source = args.adapter.resolve(strict=True)
    original_adapter = [
        oracle.file_record(path)
        for path in sorted(adapter_source.rglob("*"))
        if path.is_file()
    ]
    crate = output / "observer"
    shutil.copytree(adapter_source, crate)
    cargo = (AUDIT / "observer/Cargo.toml").read_text()
    old = 'steric_x = { path = "../frozen/repository" }'
    if cargo.count(old) != 1:
        raise ValueError("Unknown audited dependency declaration")
    (crate / "Cargo.toml").write_text(
        cargo.replace(old, "steric_x = { path = " + json.dumps(str(source)) + " }")
    )
    shutil.copyfile(AUDIT / "observer/Cargo.lock", crate / "Cargo.lock")
    oracle.check_dependency_locks(source / "Cargo.lock", crate / "Cargo.lock")
    prefix = (source / "src/geometry/buried_volume.rs").read_bytes()
    append = (AUDIT / "geometry/observer_append.rs").read_bytes()
    (crate / "src/observed_buried_volume.rs").write_bytes(prefix + append)
    helper = output / "helpers"
    helper.mkdir()
    for path in executing_helpers():
        shutil.copyfile(path, helper / path.name)
    started = {
        "kind": "scientific_remediation_build",
        "purpose": "Corrected SUT capture; no scientific pass inferred",
        "original_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPO, text=True
        ).strip(),
        "sources": oracle.source_records(source),
        "python_sources": [
            oracle.file_record(source / path.relative_to(REPO)) for path in python_paths
        ],
        "source": str(source),
        "adapter_files": [
            oracle.file_record(p) for p in sorted(crate.rglob("*")) if p.is_file()
        ],
        "helpers": [oracle.file_record(p) for p in sorted(helper.iterdir())],
        "oracle_manifest": oracle.file_record(args.oracle.resolve() / "manifest.json"),
        "audit_manifest": oracle.file_record(AUDIT / "manifest_final.json"),
        "reference_evidence": prepared["evidence_lock"],
        "bv_prefix_bytes": len(prefix),
        "bv_prefix_sha256": oracle.sha(source / "src/geometry/buried_volume.rs"),
        "unchanged_private_adapter": oracle.file_record(
            AUDIT / "geometry/observer_append.rs"
        ),
        "toolchain": {
            name: oracle.tool_version(argv)
            for name, argv in (("rustc", ["rustc", "-vV"]), ("cargo", ["cargo", "-V"]))
        },
        "environment": {
            name: os.environ.get(name)
            for name in (
                "RUSTFLAGS",
                "CARGO_ENCODED_RUSTFLAGS",
                "RUSTC",
                "RUSTC_WRAPPER",
                "CARGO_BUILD_TARGET",
            )
        },
    }
    manifest(output / "started.json", started)
    binaries = {}
    commands = {}
    for label, path, binary in (
        ("observer", crate, "stericx-audit-observer"),
        ("native", source, "stericx"),
    ):
        logs = output / (label + "-build")
        logs.mkdir()
        commands[label] = execute(
            [
                "cargo",
                "build",
                "--locked",
                "--offline",
                "--release",
                "--manifest-path",
                str(path / "Cargo.toml"),
                "--target-dir",
                str(output / "target"),
            ],
            logs,
            cwd=source,
        )
        destination = output / "bin" / binary
        destination.parent.mkdir(exist_ok=True)
        shutil.copyfile(output / "target/release" / binary, destination)
        destination.chmod(0o755)
        binaries[label] = oracle.file_record(destination)
    for item in (
        original
        + python_original
        + original_adapter
        + started["sources"]
        + started["python_sources"]
        + started["adapter_files"]
    ):
        oracle.verify(item)
    manifest(
        output / "manifest.json",
        {**started, "commands": commands, "binaries": binaries},
    )
    print(json.dumps({"built": str(output), "binaries": binaries}), flush=True)


def observe(args) -> None:
    build_path = args.build.resolve() / "manifest.json"
    build_identity = oracle.file_record(build_path)
    built, prepared = verify_build(build_path)
    binary = oracle.verify(built["binaries"]["observer"])
    output = new_output(args.output)
    plans = dict(prepared["lanes"])
    derived = output / "derived_requests"
    derived.mkdir()
    for name, dump in (("kraken_topology", False), ("kraken_topology_all_bins", True)):
        base = prepared["lanes"]["kraken"]
        source_request = oracle.verify(base["requests"])
        path = derived / (name + ".jsonl")
        with path.open("x") as stream:
            for row in oracle.requests(source_request):
                stream.write(json.dumps(topology_request(row, dump)) + "\n")
        plans[name] = {
            "requests": oracle.file_record(path),
            "expected_rows": base["expected_rows"],
            "source_requests": base["requests"],
            "derivation": (
                "Recorded explicit source neighbors; same atoms/configuration/IDs"
            ),
        }
    lanes = {}
    for name, plan in plans.items():
        directory = output / name
        directory.mkdir()
        request = oracle.verify(plan["requests"])
        print(f"Observing {name}: {plan['expected_rows']} rows", flush=True)
        with (
            request.open("rb") as stdin,
            (directory / "stdout.jsonl").open("xb") as stdout,
            (directory / "stderr").open("xb") as stderr,
        ):
            result = subprocess.run(
                [str(binary)], stdin=stdin, stdout=stdout, stderr=stderr
            )
        counts = capture_counts(name, request, directory / "stdout.jsonl")
        if result.returncode or counts["rows"] != plan["expected_rows"]:
            raise RuntimeError(f"Incomplete remediation observations: {name}")
        lanes[name] = {
            "requests": plan["requests"],
            "stdout": oracle.file_record(directory / "stdout.jsonl"),
            "stderr": oracle.file_record(directory / "stderr"),
            "returncode": result.returncode,
            **counts,
            "source_requests": plan.get("source_requests"),
            "derivation": plan.get("derivation"),
        }
        manifest(directory / "manifest.json", lanes[name])
    oracle.verify(build_identity)
    verify_build(build_path)
    for lane in lanes.values():
        for key in ("requests", "stdout", "stderr"):
            oracle.verify(lane[key])
    manifest(
        output / "manifest.json",
        {
            "kind": "scientific_remediation_observations",
            "complete_capture": True,
            "scientific_pass": False,
            "build_manifest": build_identity,
            "binary": built["binaries"]["observer"],
            "lanes": lanes,
        },
    )


def references(args) -> None:
    observation_path = args.observations.resolve() / "manifest.json"
    observation_identity = oracle.file_record(observation_path)
    observed, _built = verify_observations(observation_path)
    output = new_output(args.output)
    sealed = refs.SealedInputs(AUDIT)
    environment = refs.reference_environment(sealed)
    commands = []
    for component in args.components:
        stage = {"geometry": refs.stage_geometry, "kraken": refs.stage_kraken}[
            component
        ]
        selected = observed
        if component == "kraken":
            selected = {
                **observed,
                "lanes": {
                    **observed["lanes"],
                    "kraken": observed["lanes"][args.kraken_lane],
                },
            }
        commands.extend(stage(sealed, output, selected))
    staged = [
        oracle.file_record(path) for path in sorted(output.rglob("*")) if path.is_file()
    ]
    staged_links = {
        str(path): str(path.resolve())
        for path in sorted(output.rglob("*"))
        if path.is_symlink()
    }
    manifest(
        output / "started.json",
        {
            "kind": "unchanged_independent_references_after_scientific_remediation",
            "observations": observation_identity,
            "staged_files": staged,
            "staged_links": staged_links,
            "reference_environment": environment,
            "sealed_inputs": list(sealed.used.values()),
            "components": args.components,
            "kraken_observation_lane": args.kraken_lane,
            "commands": commands,
        },
    )
    started_identity = oracle.file_record(output / "started.json")
    runs = refs.run_commands(commands, output, args.timeout)
    for item in [
        *sealed.used.values(),
        *staged,
        started_identity,
        observation_identity,
    ]:
        oracle.verify(item)
    for path, target in staged_links.items():
        if not Path(path).is_symlink() or str(Path(path).resolve()) != target:
            raise ValueError(f"Staged reference input link changed: {path}")
    verify_observations(observation_path)
    refs.reference_environment(sealed)
    manifest(
        output / "result.json",
        {
            "commands_completed": True,
            "scientific_pass": False,
            "interpretation": (
                "Residuals and intentional convention changes require separate review"
            ),
            "started": oracle.file_record(output / "started.json"),
            "runs": runs,
            "artifacts": [
                oracle.file_record(p)
                for p in refs.scientific_files(output, args.components)
            ],
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="action", required=True)
    command = commands.add_parser("build")
    command.add_argument("--oracle", required=True, type=Path)
    command.add_argument("--adapter", required=True, type=Path)
    command.add_argument("--output", required=True, type=Path)
    command.set_defaults(run=build)
    command = commands.add_parser("observe")
    command.add_argument("--build", required=True, type=Path)
    command.add_argument("--output", required=True, type=Path)
    command.set_defaults(run=observe)
    command = commands.add_parser("references")
    command.add_argument("--observations", required=True, type=Path)
    command.add_argument("--output", required=True, type=Path)
    command.add_argument(
        "--components", nargs="+", choices=["geometry", "kraken"], required=True
    )
    command.add_argument("--timeout", type=float, default=7200.0)
    command.add_argument(
        "--kraken-lane", choices=["kraken", "kraken_topology"], default="kraken"
    )
    command.set_defaults(run=references)
    args = parser.parse_args()
    args.run(args)


if __name__ == "__main__":
    main()
