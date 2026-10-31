#!/usr/bin/env python3
"""Install one wheel/sdist in a clean interpreter, then test optional adapters."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python", required=True)
    parser.add_argument("--artifact", required=True, type=Path)
    parser.add_argument("--environment", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    environment = args.environment.resolve()
    if environment.exists():
        raise ValueError("clean environment path already exists")
    artifact = args.artifact.resolve()
    receipt = {
        "created_utc": datetime.now(UTC).isoformat(),
        "command": sys.argv,
        "controller_sha256": sha(__file__),
        "artifact": {"path": str(artifact), "sha256": sha(artifact)},
        "tests": {
            str(p.relative_to(ROOT / "python/tests")): sha(p)
            for p in sorted((ROOT / "python/tests").rglob("*"))
            if p.is_file() and "__pycache__" not in p.parts
        },
        "commands": [],
    }
    (output / "receipt_started.json").write_text(json.dumps(receipt, indent=2) + "\n")

    def run(name, command):
        print(name, flush=True)
        result = subprocess.run(command, text=True, capture_output=True, timeout=900)
        log = output / f"{name}.log"
        log.write_text(result.stdout + result.stderr)
        receipt["commands"].append(
            {
                "name": name,
                "argv": command,
                "returncode": result.returncode,
                "log": str(log),
                "sha256": sha(log),
            }
        )
        (output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
        if result.returncode:
            raise RuntimeError(f"{name} failed; complete output retained")
        return result.stdout

    run(
        "create_environment",
        ["uv", "venv", "--seed", "--python", args.python, str(environment)],
    )
    python = str(environment / "bin/python")
    tests = environment.parent / (environment.name + "-tests")
    if tests.exists():
        raise ValueError("test copy already exists")
    shutil.copytree(
        ROOT / "python/tests",
        tests,
        ignore=shutil.ignore_patterns("__pycache__", ".pytest_cache"),
    )
    run("install_base", [python, "-m", "pip", "install", "--no-deps", str(artifact)])
    smoke = """import importlib.util,json,sys,stericx
assert all(importlib.util.find_spec(n) is None for n in ['numpy','pandas','rdkit'])
m=stericx.Molecule.from_xyz_string('2\\naxis\\nC 0 0 0\\nC 0 0 1.5\\n')
r=m.descriptors(stericx.DescriptorOptions(donor_index=0,reference_index=1))
assert r.sterimol.require().l>3.19
print(json.dumps({'python':sys.version,'module':stericx.__file__,'version':stericx.__version__,'l':r.sterimol.value.l}))
"""
    run("base_smoke", [python, "-I", "-c", smoke])
    run(
        "install_test_adapters",
        [python, "-m", "pip", "install", "pytest", "numpy", "pandas", "rdkit"],
    )
    run("tests", [python, "-I", "-m", "pytest", "-q", str(tests)])
    run("pip_check", [python, "-m", "pip", "check"])
    packages = run("packages", [python, "-m", "pip", "list", "--format", "json"])
    receipt["packages"] = json.loads(packages)
    receipt["passed"] = True
    receipt["artifact_unchanged"] = sha(artifact) == receipt["artifact"]["sha256"]
    if not receipt["artifact_unchanged"]:
        raise ValueError("artifact changed during installation")
    (output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"python": args.python, "passed": True}))


if __name__ == "__main__":
    main()
