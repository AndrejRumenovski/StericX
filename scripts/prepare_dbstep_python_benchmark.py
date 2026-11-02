#!/usr/bin/env python3
"""Prepare source-bound matched Python environments; never run timing blocks."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import zipfile
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--candidate-install", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    old_receipt_path = (
        ROOT / "docs/paton_recommendations/python/distribution_receipt.json"
    )
    old = json.loads(old_receipt_path.read_text())
    install = json.loads((args.candidate_install / "receipt.json").read_text())
    previous = subprocess.check_output(
        ["git", "rev-parse", "f04b156f"], text=True
    ).strip()
    candidate_commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], text=True
    ).strip()
    source_checks = {}
    for name, expected in old["current_sources"].items():
        historical = hashlib.sha256(
            subprocess.check_output(["git", "show", f"{previous}:{name}"])
        ).hexdigest()
        current = sha(ROOT / name)
        source_checks[name] = {
            "historical_receipt": expected,
            "historical_git": historical,
            "current": current,
        }
        if historical != expected or current != expected:
            raise ValueError(f"Scientific/package source changed: {name}")
    for name in [
        "scripts/benchmark_paton_python.py",
        "examples/paton_ensemble_benchmark.rs",
        "python/examples/data/ligand_401_ensemble.xyz",
        "python/examples/data/ligand_401_options.json",
    ]:
        historical = hashlib.sha256(
            subprocess.check_output(["git", "show", f"{previous}:{name}"])
        ).hexdigest()
        source_checks[name] = {
            "historical_git": historical,
            "current": sha(ROOT / name),
        }
        if historical != source_checks[name]["current"]:
            raise ValueError(f"Benchmark controller/input changed: {name}")
    wheels = {
        "baseline": ROOT / old["wheel"]["path"],
        "candidate": Path(install["artifact"]["path"]),
    }
    expected_wheels = {
        "baseline": old["wheel"]["sha256"],
        "candidate": install["artifact"]["sha256"],
    }
    artifacts = {}
    for label, path in wheels.items():
        if sha(path) != expected_wheels[label]:
            raise ValueError(f"Historical or candidate wheel identity changed: {label}")
        with zipfile.ZipFile(path) as wheel:
            members = {
                name: hashlib.sha256(wheel.read(name)).hexdigest()
                for name in wheel.namelist()
                if name.startswith("stericx/")
            }
        artifacts[label] = {
            "path": str(path),
            "sha256": sha(path),
            "package_members": members,
        }
    if (
        artifacts["baseline"]["package_members"]
        != artifacts["candidate"]["package_members"]
    ):
        raise ValueError(
            "This preparation expects unchanged scientific package binaries"
        )
    candidate_environment = Path(
        next(
            c["argv"][-1]
            for c in install["commands"]
            if c["name"] == "create_environment"
        )
    )
    environment = output / "baseline_env"
    candidate_python = candidate_environment / "bin/python"
    receipt = {
        "started_utc": datetime.now(UTC).isoformat(),
        "command": sys.argv,
        "controller_sha256": sha(__file__),
        "baseline_commit": previous,
        "candidate_commit": candidate_commit,
        "source_checks": source_checks,
        "historical_distribution_receipt": {
            "path": str(old_receipt_path),
            "sha256": sha(old_receipt_path),
        },
        "candidate_install_receipt": {
            "path": str(args.candidate_install / "receipt.json"),
            "sha256": sha(args.candidate_install / "receipt.json"),
        },
        "artifacts": artifacts,
        "commands": [],
        "timings_run": False,
    }
    (output / "started.json").write_text(json.dumps(receipt, indent=2) + "\n")
    env = dict(
        os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", RAYON_NUM_THREADS="2"
    )

    def run(name, command):
        command = [str(p) for p in command]
        p = subprocess.run(
            ["taskset", "-c", "0,1", *command],
            env=env,
            cwd=ROOT,
            capture_output=True,
            text=True,
            timeout=900,
        )
        log = output / (name + ".log")
        log.write_text(p.stdout + p.stderr)
        receipt["commands"].append(
            {
                "name": name,
                "argv": ["taskset", "-c", "0,1", *command],
                "returncode": p.returncode,
                "log_sha256": sha(log),
            }
        )
        (output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
        if p.returncode:
            raise RuntimeError(f"{name} failed; evidence retained")
        return p.stdout

    run(
        "create_baseline",
        ["uv", "venv", "--seed", "--python", candidate_python, environment],
    )
    baseline_python = environment / "bin/python"
    run(
        "install_historical_wheel",
        [baseline_python, "-m", "pip", "install", "--no-deps", wheels["baseline"]],
    )
    packages = [
        f"{p['name']}=={p['version']}"
        for p in install["packages"]
        if p["name"].lower() != "stericx"
    ]
    run(
        "match_candidate_packages", [baseline_python, "-m", "pip", "install", *packages]
    )
    inspector = """import hashlib,importlib.metadata,json,platform,sys,pathlib,stericx
p=pathlib.Path(stericx.__file__).parent
packages=sorted([
    {'name':d.metadata['Name'],'version':d.version}
    for d in importlib.metadata.distributions()],key=lambda d:d['name'].lower())
files={str(f):hashlib.sha256(f.read_bytes()).hexdigest()
    for f in sorted(p.iterdir()) if f.is_file()}
print(json.dumps({'python':sys.version,'executable':sys.executable,
    'platform':platform.platform(),'packages':packages,'files':files},indent=2))
"""
    inspected = {}
    for label, python in [
        ("baseline", baseline_python),
        ("candidate", candidate_python),
    ]:
        inspected[label] = json.loads(
            run(label + "_identity", [python, "-I", "-c", inspector])
        )
        run(label + "_pip_check", [python, "-m", "pip", "check"])
        expected = artifacts[label]["package_members"]
        for path, digest in inspected[label]["files"].items():
            if expected["stericx/" + Path(path).name] != digest:
                raise ValueError(f"Installed artifact differs: {path}")
    if (
        inspected["baseline"]["python"] != inspected["candidate"]["python"]
        or inspected["baseline"]["packages"] != inspected["candidate"]["packages"]
    ):
        raise ValueError("Interpreter/package environments differ")
    receipt.update(
        environments=inspected,
        package_files_identical=True,
        source_changes={
            name: sha(ROOT / name)
            for name, row in source_checks.items()
            if sha(ROOT / name) != row["current"]
        },
        ready_for_quiet_benchmark=True,
        timings_run=False,
    )
    if receipt["source_changes"]:
        receipt["ready_for_quiet_benchmark"] = False
    (output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    if receipt["source_changes"]:
        raise ValueError("Bound sources changed during preparation")
    print(
        json.dumps(
            {
                "baseline_python": str(baseline_python),
                "candidate_python": str(candidate_python),
                "package_files_identical": True,
                "timings_run": False,
            }
        )
    )


if __name__ == "__main__":
    main()
