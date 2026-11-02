#!/usr/bin/env python3
"""Replay existing advanced scientific/package checks without changing their gates.

Execution completion and scientific acceptance are recorded separately. Historical
strict hard-sphere comparisons remain intact, including their expected failures.
Every phase uses a new destination and records sources, environment and logs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REFERENCE_PYTHON = ROOT / ".stericx/paton/reference/venv/bin/python"
BUILD_PYTHON = ROOT / ".stericx/paton/python-build/bin/python"
MATURIN = BUILD_PYTHON.parent / "maturin"
LIMITS = {
    "RAYON_NUM_THREADS": "2",
    "OPENBLAS_NUM_THREADS": "1",
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "CARGO_BUILD_JOBS": "2",
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def now():
    return datetime.now(UTC).isoformat()


def sources():
    paths = [ROOT / p for p in ["Cargo.toml", "Cargo.lock"]]
    for directory in ["src", "examples", "python"]:
        paths.extend(
            p
            for p in (ROOT / directory).rglob("*")
            if p.is_file()
            and not {"target", "__pycache__", ".pytest_cache", ".venv"}.intersection(
                p.relative_to(ROOT).parts
            )
            and p.suffix in {".rs", ".py", ".pyi", ".toml", ".lock", ".json", ".md"}
        )
    paths.extend((ROOT / "scripts").glob("validate_*.py"))
    paths.append(ROOT / "scripts/study_conformer_ensembles.py")
    paths.append(Path(__file__))
    paths.extend((ROOT / "docs/validation/dbstep/resolution_baseline/dbstep").glob("*"))
    return {str(p.relative_to(ROOT)): sha(p) for p in sorted(set(paths)) if p.is_file()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=["references", "distribution", "installs"])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--build", type=Path)
    parser.add_argument("--distribution", type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, **LIMITS)
    env["PATH"] = str(BUILD_PYTHON.parent) + os.pathsep + env["PATH"]
    env["VIRTUAL_ENV"] = str(BUILD_PYTHON.parent.parent)
    receipt = {
        "started_utc": now(),
        "command": sys.argv,
        "controller_sha256": sha(__file__),
        "source_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "git_status": subprocess.check_output(
            ["git", "status", "--short"], cwd=ROOT, text=True
        ),
        "sources_before": sources(),
        "resource_limits": {"cpu_affinity": [0, 1], **LIMITS},
        "commands": [],
        "acceptance_policy": "Original controllers unchanged; exit=0 is not admission.",
    }
    (output / "receipt_started.json").write_text(json.dumps(receipt, indent=2) + "\n")

    def run(name, command, timeout=7200):
        command = [str(s) for s in command]
        launched = ["taskset", "-c", "0,1", *command]
        print(json.dumps({"started": name, "argv": launched}), flush=True)
        start = time.monotonic()
        log = output / (name + ".log")
        with log.open("w") as stream:
            try:
                proc = subprocess.run(
                    launched,
                    cwd=ROOT,
                    env=env,
                    stdout=stream,
                    stderr=subprocess.STDOUT,
                    timeout=timeout,
                    check=False,
                )
                returncode = proc.returncode
            except subprocess.TimeoutExpired:
                returncode = "timeout"
        item = {
            "name": name,
            "argv": launched,
            "returncode": returncode,
            "elapsed_seconds": time.monotonic() - start,
            "log": str(log),
            "log_sha256": sha(log),
        }
        receipt["commands"].append(item)
        (output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
        print(json.dumps(item), flush=True)
        return returncode

    if args.phase == "references":
        if not args.build:
            parser.error("references requires --build")
        code = run(
            "build_examples",
            [
                "cargo",
                "build",
                "--locked",
                "--offline",
                "--example",
                "ensemble_report",
                "--example",
                "density_report",
                "--example",
                "volume_profile",
                "--example",
                "population_observer",
            ],
        )
        if code != 0:
            raise RuntimeError("Fresh example build failed; receipt retained")
        binaries = [
            "ensemble_report",
            "density_report",
            "volume_profile",
            "population_observer",
        ]
        receipt["binaries"] = {
            name: {
                "path": str(ROOT / "target/debug/examples" / name),
                "sha256": sha(ROOT / "target/debug/examples" / name),
            }
            for name in binaries
        }
        jobs = [
            (
                "sterimol",
                "validate_dbstep_sterimol.py",
                ["--engine", "target/debug/examples/ensemble_report"],
            ),
            ("volume", "validate_dbstep_volume.py", ["--build", args.build.resolve()]),
            (
                "density",
                "validate_density_volume.py",
                ["--engine", "target/debug/examples/density_report"],
            ),
            ("profiles", "validate_volume_profiles.py", []),
            (
                "populations",
                "validate_conformer_populations.py",
                ["--observer", "target/debug/examples/population_observer"],
            ),
            ("ensemble_reference", "study_conformer_ensembles.py", []),
        ]
        for name, script, extra in jobs:
            run(
                name,
                [
                    REFERENCE_PYTHON,
                    ROOT / "scripts" / script,
                    *extra,
                    "--output",
                    output / name,
                ],
            )
    elif args.phase == "distribution":
        run("maturin_version", [MATURIN, "--version"])
        run("build_packages", [BUILD_PYTHON, "-m", "pip", "list", "--format", "json"])
        run(
            "wheel_build",
            [
                MATURIN,
                "build",
                "--release",
                "--locked",
                "--manifest-path",
                "python/Cargo.toml",
                "--zig",
                "--compatibility",
                "manylinux_2_28",
                "--interpreter",
                BUILD_PYTHON,
                "--out",
                output / "wheels",
            ],
        )
        run(
            "sdist_build",
            [
                MATURIN,
                "sdist",
                "--manifest-path",
                "python/Cargo.toml",
                "--out",
                output / "sdist",
            ],
        )
        receipt["artifacts"] = {
            str(p.relative_to(output)): sha(p)
            for glob in ["wheels/*.whl", "sdist/*.tar.gz"]
            for p in output.glob(glob)
        }
    else:
        if not args.distribution:
            parser.error("installs requires --distribution")
        wheels = list(args.distribution.glob("wheels/*.whl"))
        archives = list(args.distribution.glob("sdist/*.tar.gz"))
        if len(wheels) != 1 or len(archives) != 1:
            raise ValueError("Exactly one newly built wheel and one sdist are required")
        jobs = [("wheel-" + v, v, wheels[0]) for v in ["3.11", "3.12", "3.13", "3.14"]]
        jobs.append(("sdist-3.12", "3.12", archives[0]))
        for name, version, artifact in jobs:
            run(
                name,
                [
                    BUILD_PYTHON,
                    ROOT / "scripts/validate_python_install.py",
                    "--python",
                    version,
                    "--artifact",
                    artifact.resolve(),
                    "--environment",
                    output / "environments" / name,
                    "--output",
                    output / name,
                ],
            )
    receipt["sources_after"] = sources()
    receipt["source_changes"] = {
        p: [digest, receipt["sources_after"].get(p)]
        for p, digest in receipt["sources_before"].items()
        if receipt["sources_after"].get(p) != digest
    }
    receipt["completed_utc"] = now()
    receipt["execution_completed"] = True
    receipt["all_commands_exit_zero"] = all(
        c["returncode"] == 0 for c in receipt["commands"]
    )
    (output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(
        json.dumps(
            {
                k: receipt[k]
                for k in [
                    "execution_completed",
                    "all_commands_exit_zero",
                    "source_changes",
                ]
            }
        ),
        flush=True,
    )
    return bool(receipt["source_changes"])


if __name__ == "__main__":
    raise SystemExit(main())
