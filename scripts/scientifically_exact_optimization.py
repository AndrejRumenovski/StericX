#!/usr/bin/env python3
"""Freeze a reviewed corrected build and compare full optimization observations.

Scientific admission is a separately reviewed receipt, never inferred from a
successful subprocess or equality with an older implementation. Optimization
comparisons use exact complete observations; no tolerance is introduced here.
The historical audit and prior performance freezes are never modified.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
from itertools import zip_longest
from pathlib import Path

import check_current_scientific_equivalence as oracle
import performance_freeze as inventory
import profile_stericx as profile
import replay_scientific_remediation as replay

REPO = oracle.REPO


def bound_copy(source: Path, destination: Path) -> dict:
    before = oracle.file_record(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        raise ValueError(f"Refusing to overwrite evidence: {destination}")
    shutil.copy2(source, destination)
    oracle.verify(before)
    copied = oracle.file_record(destination)
    if (copied["bytes"], copied["sha256"]) != (before["bytes"], before["sha256"]):
        raise ValueError(f"Copy identity differs: {source}")
    return {"original": before, "snapshot": copied}


def source_set() -> list[Path]:
    """Include current added helpers/tests, not just already committed files."""
    selected = set(inventory.source_paths())
    for directory, suffixes in (
        ("src", {".rs"}),
        ("scripts", {".py", ".c", ".sh"}),
        ("studies", {".py"}),
        ("tests", {".py", ".rs"}),
        ("examples", {".rs"}),
    ):
        selected.update(
            path
            for path in (REPO / directory).rglob("*")
            if path.is_file() and path.suffix in suffixes
        )
    return sorted(selected)


def verify_frozen(path: Path) -> dict:
    frozen = replay.load_manifest(path)
    if frozen.get("kind") != "reviewed_corrected_optimization_baseline":
        raise ValueError("Expected corrected baseline")
    # Deliberately verify frozen source here, not live source: candidates change
    # live files, while the frozen executable/source must remain unchanged.
    for item in frozen["copies"]:
        oracle.verify(item["snapshot"])
    for item in frozen["evidence"]:
        oracle.verify(item)
    for item in frozen["raw_evidence"]:
        oracle.verify(item)
    oracle.verify(frozen["machine"])
    workload_path = oracle.verify(frozen["workload_manifest"])
    profile.verify_manifest(oracle.read(workload_path))
    return frozen


def corrected_workloads(path: Path, build_path: Path) -> dict:
    workloads = oracle.read(path)
    if (
        workloads.get("kind") != "corrected_science_profile_inputs"
        or workloads.get("build_manifest") != oracle.file_record(build_path)
        or workloads.get("full_build_proof_verified_before_after") is not True
        or any(item["returncode"] != 0 for item in workloads["preparation_commands"])
    ):
        raise ValueError(
            "Workloads must be prepared successfully by the admitted build"
        )
    profile.verify_manifest(workloads)
    return workloads


def freeze(args) -> None:
    build_path = args.build.resolve() / "manifest.json"
    built, _ = replay.verify_build(build_path)
    observation_path = args.observations.resolve() / "manifest.json"
    observed, _ = replay.verify_observations(observation_path)
    if observed["build_manifest"] != oracle.file_record(build_path):
        raise ValueError("Scientific observations belong to another build")
    admission_path = args.admission.resolve()
    admission = replay.load_manifest(admission_path)
    required = {"geometry", "thermodynamics", "models", "behavior", "engineering"}
    if (
        admission.get("kind") != "reviewed_scientific_remediation_admission"
        or admission.get("admitted_for_exact_optimization") is not True
        or set(admission.get("gates", {})) != required
        or any(item.get("passed") is not True for item in admission["gates"].values())
        or admission.get("build") != oracle.file_record(build_path)
        or admission.get("observations") != oracle.file_record(observation_path)
        or not admission.get("retained_limitations")
    ):
        raise ValueError("Explicit complete reviewed scientific admission required")
    evidence = [oracle.file_record(build_path), oracle.file_record(observation_path)]
    evidence.append(oracle.file_record(admission_path))
    raw_evidence = {}
    for gate in admission["gates"].values():
        if not gate.get("evidence") or not gate.get("artifact_roots"):
            raise ValueError("Gate lacks supporting evidence")
        for item in gate["evidence"]:
            oracle.verify(item)
            evidence.append(item)
        for root_name in gate["artifact_roots"]:
            root = Path(root_name).resolve(strict=True)
            if not root.is_dir():
                raise ValueError(f"Raw evidence directory required: {root}")
            for path in sorted(root.rglob("*")):
                if path.is_file() and "__pycache__" not in path.parts:
                    raw_evidence[str(path)] = oracle.file_record(path)
    if str(observation_path) not in raw_evidence:
        raise ValueError("Full observation directory must be included in gate evidence")
    # Build snapshot dependencies are immutable even after live candidate edits.
    for key in ("sources", "python_sources", "adapter_files", "helpers"):
        for item in built[key]:
            raw_evidence[item["path"]] = item
    for item in built["binaries"].values():
        raw_evidence[item["path"]] = item
    audit_manifest = oracle.verify(built["audit_manifest"])
    evidence.append(built["audit_manifest"])
    for item in oracle.read(audit_manifest)["files"]:
        absolute = {**item, "path": str(audit_manifest.parent / item["path"])}
        oracle.verify(absolute)
        raw_evidence[absolute["path"]] = absolute
    # Bind current production to the independently exercised source snapshot.
    source_root = Path(built["source"])
    for item in built["sources"]:
        relative = Path(item["path"]).relative_to(source_root)
        current = oracle.file_record(REPO / relative)
        if (current["bytes"], current["sha256"]) != (item["bytes"], item["sha256"]):
            raise ValueError(f"Unreviewed Rust source: {relative}")
    for item in built["python_sources"]:
        relative = Path(item["path"]).relative_to(source_root)
        if relative.parts[0] == "studies" or relative.as_posix() in {
            "scripts/prepare_data.py",
            "scripts/stericx_quantum.py",
            "scripts/prepare_quantum_data.py",
            "pyproject.toml",
            "uv.lock",
        }:
            if oracle.sha(REPO / relative) != item["sha256"]:
                raise ValueError(f"Unreviewed scientific Python source: {relative}")
    commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=REPO, text=True
    ).strip()
    if commit != args.expected_commit:
        raise ValueError("Unexpected baseline commit")
    workload_path = args.workloads.resolve() / "workloads/manifest.json"
    corrected_workloads(workload_path, build_path)
    output = replay.new_output(args.output)
    copies = [
        bound_copy(path, output / "source" / path.relative_to(REPO))
        for path in source_set()
    ]
    for name, item in built["binaries"].items():
        copies.append(bound_copy(oracle.verify(item), output / "bin" / name))
    for name, path in (
        ("admission.json", admission_path),
        ("build.json", build_path),
        ("observations.json", observation_path),
        ("workloads.json", workload_path),
    ):
        copies.append(bound_copy(path, output / "receipts" / name))
    for path in sorted((REPO / "docs/scientific_remediation").rglob("*")):
        if path.is_file() and path.suffix in {".md", ".json"}:
            copies.append(bound_copy(path, output / "reports" / path.relative_to(REPO)))
    machine_path = output / "machine.json"
    replay.manifest(machine_path, inventory.machine())
    for item in copies:
        oracle.verify(item["original"])
        oracle.verify(item["snapshot"])
    for item in evidence:
        oracle.verify(item)
    for item in raw_evidence.values():
        oracle.verify(item)
    profile.verify_manifest(oracle.read(workload_path))
    replay.manifest(
        output / "manifest.json",
        {
            "kind": "reviewed_corrected_optimization_baseline",
            "created_utc": profile.now(),
            "git_commit": commit,
            "scientific_scope": admission["scope"],
            "retained_limitations": admission["retained_limitations"],
            "optimization_equivalence": (
                "exact complete values, bits and errors; zero tolerance"
            ),
            "copies": copies,
            "evidence": evidence,
            "raw_evidence": list(raw_evidence.values()),
            "admitted_observations": oracle.file_record(observation_path),
            "machine": oracle.file_record(machine_path),
            "workload_manifest": oracle.file_record(workload_path),
            "thread_affinities": inventory.AFFINITIES,
            "native_repetitions": 7,
            "warmups": 1,
        },
    )
    verify_frozen(output / "manifest.json")
    print(f"Corrected scientific baseline frozen: {output}", flush=True)


def compare_observations(args) -> None:
    frozen = verify_frozen(args.snapshot.resolve() / "manifest.json")
    if (
        oracle.file_record(args.baseline.resolve() / "manifest.json")
        != frozen["admitted_observations"]
    ):
        raise ValueError("Comparison baseline is not the frozen admitted baseline")
    output = replay.new_output(args.output)
    manifests = [
        root.resolve() / "manifest.json" for root in (args.baseline, args.candidate)
    ]
    identities = [oracle.file_record(path) for path in manifests]
    observations = [replay.verify_observations(path)[0] for path in manifests]
    old, new = observations
    if set(old["lanes"]) != set(new["lanes"]):
        raise ValueError("Observation lane set changed")
    lanes = {}
    failures = []
    for name, left in old["lanes"].items():
        right = new["lanes"][name]
        if left["requests"]["sha256"] != right["requests"]["sha256"]:
            raise ValueError(f"Scientific requests differ: {name}")
        mismatches = []
        count = 0
        for count, (a, b) in enumerate(
            zip_longest(
                oracle.requests(oracle.verify(left["stdout"])),
                oracle.requests(oracle.verify(right["stdout"])),
            ),
            1,
        ):
            # Canonical bytes distinguish -0.0 from +0.0 and preserve every
            # reported IEEE bit, value, selected atom, frame and error string.
            if oracle.canonical(a) != oracle.canonical(b):
                mismatches.append({"row": count, "baseline": a, "candidate": b})
        stderr_equal = left["stderr"]["sha256"] == right["stderr"]["sha256"]
        same = (
            not mismatches
            and stderr_equal
            and left["returncode"] == right["returncode"]
        )
        if not same:
            failures.append(name)
        replay.manifest(output / f"{name}.json", {"mismatches": mismatches})
        lanes[name] = {
            "rows": count,
            "passed": same,
            "mismatches": len(mismatches),
            "stderr_equal": stderr_equal,
        }
    for item in identities:
        oracle.verify(item)
    replay.manifest(
        output / "comparison.json",
        {
            "kind": "exact_corrected_baseline_comparison",
            "baseline": identities[0],
            "candidate": identities[1],
            "executing_helper": oracle.file_record(Path(__file__)),
            "tolerance": 0,
            "lanes": lanes,
            "passed": not failures,
            "independent_reference_gate_required_separately": True,
        },
    )
    if failures:
        raise ValueError(f"Scientific equivalence FAILED: {failures}")
    total = sum(row["rows"] for row in lanes.values())
    print(f"Exact match: {total} complete observations")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_subparsers(dest="action", required=True)
    action = actions.add_parser("freeze")
    for name in ("build", "observations", "admission", "workloads", "output"):
        action.add_argument(f"--{name}", type=Path, required=True)
    action.add_argument("--expected-commit", required=True)
    action = actions.add_parser("compare-observations")
    for name in ("snapshot", "baseline", "candidate", "output"):
        action.add_argument(f"--{name}", type=Path, required=True)
    action = actions.add_parser("verify")
    action.add_argument("--snapshot", type=Path, required=True)
    args = parser.parse_args()
    if args.action == "freeze":
        freeze(args)
    elif args.action == "compare-observations":
        compare_observations(args)
    else:
        verify_frozen(args.snapshot.resolve() / "manifest.json")
        print("Frozen baseline identity verified")


if __name__ == "__main__":
    main()
