#!/usr/bin/env python3
"""Build and test the version-only RC using unchanged installation contracts."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import time
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASELINE = "497de873"
BUILD_PYTHON = ROOT / ".stericx/paton/python-build/bin/python"
MATURIN = BUILD_PYTHON.parent / "maturin"
VERSION_CHANGES = {
    "python/Cargo.toml": ('version = "0.4.0-dev.0"', 'version = "0.4.0-rc.1"'),
    "python/Cargo.lock": ('version = "0.4.0-dev.0"', 'version = "0.4.0-rc.1"'),
    "python/pyproject.toml": ('version = "0.4.0.dev0"', 'version = "0.4.0rc1"'),
    "python/stericx/__init__.py": (
        '__version__ = "0.4.0.dev0"',
        '__version__ = "0.4.0rc1"',
    ),
}
INVENTORY = """import importlib.metadata,json,platform,sys
print(json.dumps({'python':sys.version,'executable':sys.executable,
'platform':platform.platform(),'packages':sorted([
{'name':d.metadata['Name'],'version':d.version}
for d in importlib.metadata.distributions()],key=lambda d:d['name'].lower())}))
"""
VERSION_CHECK = """import importlib.metadata,json,stericx
assert stericx.__version__=='0.4.0rc1'
assert importlib.metadata.version('stericx')=='0.4.0rc1'
assert stericx._native.__version__=='0.4.0-rc.1'
print(json.dumps({'python_version':stericx.__version__,
'distribution_version':importlib.metadata.version('stericx'),
'native_binding_version':stericx._native.__version__,
'engine_version':'0.3.0','module':stericx.__file__,
'native_module':stericx._native.__file__}))
"""


def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def record(path):
    path = Path(path).resolve()
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": sha(path)}


def now():
    return datetime.now(UTC).isoformat()


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT)


def source_paths():
    return sorted(
        git(
            "ls-files",
            "Cargo.toml",
            "Cargo.lock",
            "src",
            "examples",
            "python",
            "scripts/validate_python_install.py",
            "docs/validation/dbstep/ACCEPTANCE_POLICY.md",
            "docs/validation/dbstep/POLICY_APPROVAL.md",
        )
        .decode()
        .splitlines()
    )


def version_only(paths):
    checks = {}
    for name in paths:
        original = git("show", f"{BASELINE}:{name}")
        expected = original
        if name in VERSION_CHANGES:
            before, after = (value.encode() for value in VERSION_CHANGES[name])
            if original.count(before) != 1:
                raise ValueError(
                    f"Version substitution must occur exactly once: {name}"
                )
            expected = original.replace(before, after)
        current = (ROOT / name).read_bytes()
        if current != expected:
            raise ValueError(
                f"Unexpected source change since validated baseline: {name}"
            )
        checks[name] = {
            "baseline_sha256": hashlib.sha256(original).hexdigest(),
            "current_sha256": hashlib.sha256(current).hexdigest(),
            "version_substitution": name in VERSION_CHANGES,
        }
    if not set(VERSION_CHANGES).issubset(checks):
        raise ValueError("Version metadata coverage incomplete")
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--archive", required=True, type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    paths = source_paths()
    checks = version_only(paths)
    receipt = {
        "kind": "dbstep_resolution_local_release_candidate_validation",
        "started_utc": now(),
        "source_commit": git("rev-parse", "HEAD").decode().strip(),
        "validated_pre_version_commit": git("rev-parse", BASELINE).decode().strip(),
        "controller": record(__file__),
        "commands": [],
        "source_checks": checks,
        "sources_before": {p: sha(ROOT / p) for p in paths},
        "version_only_verified": True,
        "source_changes": None,
        "artifacts": {},
        "installations": {},
        "complete": False,
        "passed": False,
        "publication_performed": False,
    }
    receipt_path = output / "receipt.json"

    def save():
        receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")

    save()
    env = dict(
        os.environ,
        RAYON_NUM_THREADS="2",
        OPENBLAS_NUM_THREADS="1",
        OMP_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
        NUMEXPR_NUM_THREADS="1",
        CARGO_BUILD_JOBS="2",
        CARGO_TARGET_DIR=str(output / "build_target"),
        VIRTUAL_ENV=str(BUILD_PYTHON.parent.parent),
    )
    env["PATH"] = str(BUILD_PYTHON.parent) + os.pathsep + env["PATH"]
    receipt["resource_configuration"] = {
        "cpu_affinity": [0, 1],
        **{k: env[k] for k in env if k.endswith("NUM_THREADS")},
        "CARGO_BUILD_JOBS": env["CARGO_BUILD_JOBS"],
        "CARGO_TARGET_DIR": env["CARGO_TARGET_DIR"],
    }

    def run(name, argv, timeout=1800):
        argv = ["taskset", "-c", "0,1", *[str(x) for x in argv]]
        log = output / f"{name}.log"
        command = {"name": name, "argv": argv, "started_utc": now()}
        receipt["commands"].append(command)
        save()
        print(json.dumps({"started": name}), flush=True)
        started = time.monotonic()
        with log.open("w") as stream:
            try:
                proc = subprocess.run(
                    argv,
                    cwd=ROOT,
                    env=env,
                    stdout=stream,
                    stderr=subprocess.STDOUT,
                    timeout=timeout,
                    check=False,
                )
                code = proc.returncode
            except subprocess.TimeoutExpired:
                code = "timeout"
        command.update(
            returncode=code,
            elapsed_seconds=time.monotonic() - started,
            log=record(log),
        )
        save()
        if code != 0:
            raise RuntimeError(f"RC validation failed: {name}; evidence preserved")
        return log.read_text()

    run("rust_version", ["rustc", "--version", "--verbose"])
    run("cargo_version", ["cargo", "--version"])
    run("maturin_version", [MATURIN, "--version"])
    receipt["build_environment"] = json.loads(
        run("build_inventory", [BUILD_PYTHON, "-c", INVENTORY])
    )
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
    wheels = list((output / "wheels").glob("*.whl"))
    sdists = list((output / "sdist").glob("*.tar.gz"))
    if len(wheels) != 1 or len(sdists) != 1:
        raise ValueError("Exactly one wheel and source distribution are required")
    if wheels[0].name != "stericx-0.4.0rc1-cp311-abi3-manylinux_2_28_x86_64.whl":
        raise ValueError("Unexpected RC wheel identity")
    if sdists[0].name != "stericx-0.4.0rc1.tar.gz":
        raise ValueError("Unexpected RC source distribution identity")
    receipt["artifacts"] = {"wheel": record(wheels[0]), "sdist": record(sdists[0])}
    save()
    jobs = [("wheel-" + v, v, wheels[0]) for v in ["3.11", "3.12", "3.13", "3.14"]]
    jobs.append(("sdist-3.12", "3.12", sdists[0]))
    for name, version, artifact in jobs:
        environment = output / "environments" / name
        installation = output / "installations" / name
        run(
            name,
            [
                BUILD_PYTHON,
                ROOT / "scripts/validate_python_install.py",
                "--python",
                version,
                "--artifact",
                artifact,
                "--environment",
                environment,
                "--output",
                installation,
            ],
        )
        installed = json.loads((installation / "receipt.json").read_text())
        if not installed["passed"] or not installed["artifact_unchanged"]:
            raise ValueError(f"Installation not admitted: {name}")
        if any(c["returncode"] for c in installed["commands"]):
            raise ValueError(f"Installation contains failed command: {name}")
        test_text = (installation / "tests.log").read_text()
        if "71 passed" not in test_text or "skipped" in test_text:
            raise ValueError(f"Original 71-test coverage not complete: {name}")
        versions = json.loads(
            run(
                name + "_versions",
                [environment / "bin/python", "-I", "-c", VERSION_CHECK],
            )
        )
        receipt["installations"][name] = {
            "receipt": record(installation / "receipt.json"),
            "tests": record(installation / "tests.log"),
            "passed": True,
            "tests_passed": 71,
            "versions": versions,
        }
        save()
    after = {p: sha(ROOT / p) for p in paths}
    receipt["sources_after"] = after
    receipt["source_changes"] = {
        p: [receipt["sources_before"][p], after.get(p)]
        for p in paths
        if receipt["sources_before"][p] != after.get(p)
    }
    if source_paths() != paths or receipt["source_changes"]:
        raise ValueError("Source changed during RC validation")
    if version_only(paths) != checks:
        raise ValueError("Version-only source comparison changed")
    for item in receipt["artifacts"].values():
        if record(item["path"]) != item:
            raise ValueError("Built artifact changed during validation")
    if record(__file__) != receipt["controller"]:
        raise ValueError("RC controller changed during validation")
    receipt.update(completed_utc=now(), complete=True, passed=True)
    save()
    archive = args.archive.resolve()
    archive.mkdir(parents=True, exist_ok=False)
    included = [receipt_path, Path(__file__), *output.glob("*.log")]
    included += [
        p
        for p in (output / "installations").rglob("*")
        if p.is_file() and p.suffix in {".json", ".log"}
    ]
    archived = []
    for path in sorted(set(included)):
        relative = (
            Path("controller.py")
            if path == Path(__file__)
            else path.relative_to(output)
        )
        dest = archive / relative
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, dest)
        if sha(dest) != sha(path):
            raise ValueError("Archive copy differs")
        archived.append({"source": record(path), "archive_relative": str(relative)})
    (archive / "archive_manifest.json").write_text(
        json.dumps(archived, indent=2) + "\n"
    )
    print(
        json.dumps({"complete": True, "passed": True, "installations": 5}), flush=True
    )


if __name__ == "__main__":
    main()
