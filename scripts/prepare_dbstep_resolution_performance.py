#!/usr/bin/env python3
"""Build separate frozen-source diagnostic binaries before quiet-window timing."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASELINE = "f04b156f7b8b31700e1c6ab978e337b00cc734e0"


def record(path):
    path = Path(path).resolve()
    return {
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    os.sched_setaffinity(0, {0, 1})
    env = {**os.environ, "CARGO_BUILD_JOBS": "2", "RAYON_NUM_THREADS": "2"}
    receipt = {
        "controller": record(__file__),
        "baseline_commit": BASELINE,
        "builds": {},
        "complete": False,
        "timing_started": False,
    }
    save(args.output / "receipt.json", receipt)
    for label, build_path in [
        ("baseline", args.baseline),
        ("candidate", args.candidate),
    ]:
        manifest_path = build_path.resolve() / "manifest.json"
        build = json.loads(manifest_path.read_text())
        for item in [*build["sources"], *build["binaries"].values()]:
            if record(item["path"]) != item:
                raise ValueError("Changed frozen build")
        directory = args.output / label
        source = directory / "source"
        source.mkdir(parents=True)
        for item in build["sources"]:
            rel = Path(item["path"]).relative_to(build["source"])
            target = source / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(item["path"], target)
            original = subprocess.check_output(
                ["git", "show", f"{BASELINE}:{rel.as_posix()}"], cwd=ROOT
            )
            if hashlib.sha256(original).hexdigest() != item["sha256"]:
                raise ValueError(
                    "Unexpected production source correction; review preparation"
                )
        for name in ["paton_ensemble_benchmark", "paton_ensemble_allocations"]:
            relative = f"examples/{name}.rs"
            data = subprocess.check_output(
                ["git", "show", f"{BASELINE}:{relative}"], cwd=ROOT
            )
            if data != (ROOT / relative).read_bytes():
                raise ValueError("Benchmark helper changed")
            (source / "examples").mkdir(exist_ok=True)
            (source / relative).write_bytes(data)
        bound = [record(p) for p in sorted(source.rglob("*")) if p.is_file()]
        entry = {
            "source_build": record(manifest_path),
            "sources": bound,
            "commands": [],
        }
        receipt["builds"][label] = entry
        for feature, targets in [
            (False, ["--example", "paton_ensemble_benchmark"]),
            (True, ["--bin", "stericx", "--example", "paton_ensemble_allocations"]),
        ]:
            target_dir = directory / (
                "instrumented-target" if feature else "native-target"
            )
            command = [
                "cargo",
                "build",
                "--release",
                "--locked",
                "--offline",
                "--manifest-path",
                str(source / "Cargo.toml"),
                "--target-dir",
                str(target_dir),
                *targets,
            ]
            if feature:
                command += ["--features", "profiling"]
            log = directory / (
                "instrumented-build.log" if feature else "native-build.log"
            )
            start = time.perf_counter()
            with log.open("w") as stream:
                result = subprocess.run(
                    command,
                    cwd=source,
                    env=env,
                    stdout=stream,
                    stderr=subprocess.STDOUT,
                )
            entry["commands"].append(
                {
                    "argv": command,
                    "returncode": result.returncode,
                    "elapsed_seconds": time.perf_counter() - start,
                    "log": record(log),
                }
            )
            save(args.output / "receipt.json", receipt)
            if result.returncode:
                raise RuntimeError(f"Build failed: {log}")
        if bound != [record(p["path"]) for p in bound]:
            raise ValueError("Diagnostic source changed during compilation")
        entry["binaries"] = {
            "allocation_cli": record(directory / "instrumented-target/release/stericx"),
            "ensemble_allocations": record(
                directory
                / "instrumented-target/release/examples/paton_ensemble_allocations"
            ),
            "ensemble_benchmark": record(
                directory / "native-target/release/examples/paton_ensemble_benchmark"
            ),
        }
        save(args.output / "receipt.json", receipt)
        print(f"Prepared {label} diagnostics", flush=True)
    receipt["complete"] = True
    save(args.output / "receipt.json", receipt)


if __name__ == "__main__":
    main()
