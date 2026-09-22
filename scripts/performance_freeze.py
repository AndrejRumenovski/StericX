"""Freeze and profile the current audited implementation without claiming validation.

The immutable historical workload/audit manifests remain authoritative at their
original paths. This helper verifies and binds them, copies current source and
executables, and runs the existing resource/oracle harness sequentially.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tomllib
from pathlib import Path

import profile_stericx as harness

REPO = Path(__file__).resolve().parents[1]
AUDIT = REPO / "docs/scientific_accuracy_audit"
AFFINITIES = {1: "2", 2: "2,3", 4: "0,1,2,3", 6: "0,1,2,3,4,5"}
ENV_KEYS = {
    "PATH",
    "LANG",
    "LC_ALL",
    "LC_CTYPE",
    "TZ",
    "RUSTFLAGS",
    "CARGO_ENCODED_RUSTFLAGS",
    "RUSTDOCFLAGS",
    "CARGO_BUILD_TARGET",
    "CARGO_BUILD_RUSTFLAGS",
    "CARGO_PROFILE_RELEASE_DEBUG",
    "CARGO_PROFILE_RELEASE_LTO",
    "CARGO_PROFILE_RELEASE_CODEGEN_UNITS",
    "RAYON_NUM_THREADS",
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "MALLOC_ARENA_MAX",
    "LD_PRELOAD",
    "LD_LIBRARY_PATH",
    "STERICX_PROFILE_PATH",
}


def save(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def capture(argv: list[str]) -> dict:
    try:
        result = subprocess.run(argv, cwd=REPO, capture_output=True, text=True)
        return dict(
            argv=argv,
            returncode=result.returncode,
            stdout=result.stdout,
            stderr=result.stderr,
        )
    except OSError as error:
        return dict(argv=argv, error=str(error))


def record(path: Path, relative_to: Path | None = None) -> dict:
    return {
        "path": str(path.relative_to(relative_to) if relative_to else path),
        "bytes": path.stat().st_size,
        "sha256": harness.digest(path),
    }


def verify_records(base: Path, rows: list[dict]) -> dict:
    failures = []
    for row in rows:
        path = base / row["path"]
        if (
            not path.is_file()
            or path.stat().st_size != row["bytes"]
            or harness.digest(path) != row["sha256"]
        ):
            failures.append(row["path"])
    return {
        "verified_utc": harness.now(),
        "files_checked": len(rows),
        "bytes_checked": sum(row["bytes"] for row in rows),
        "failures": failures,
        "passed": not failures,
    }


def require_passed(result: dict) -> None:
    if not result["passed"]:
        raise RuntimeError(f"frozen evidence changed: {result['failures'][:10]}")


def source_paths() -> list[Path]:
    tracked = subprocess.check_output(["git", "ls-files", "-z"], cwd=REPO)
    names = [name for name in tracked.decode().split("\0") if name]
    roots = {"src", "scripts", "tests", "examples", "studies", ".cargo", ".github"}
    selected = {
        REPO / name
        for name in names
        if Path(name).parts[0] in roots
        or (
            len(Path(name).parts) == 1
            and (
                Path(name).suffix in {".toml", ".lock", ".rs", ".py", ".sh"}
                or name == "README.md"
            )
        )
    }
    selected.add(Path(__file__).resolve())
    return sorted(selected)


def machine() -> dict:
    files = [
        Path(name)
        for name in (
            "/proc/cpuinfo",
            "/proc/meminfo",
            "/proc/loadavg",
            "/etc/os-release",
            "/proc/sys/kernel/perf_event_paranoid",
            "/proc/sys/kernel/randomize_va_space",
            "/sys/devices/system/cpu/cpufreq/boost",
        )
    ]
    for cpu in sorted(os.sched_getaffinity(0)):
        root = Path(f"/sys/devices/system/cpu/cpu{cpu}")
        files.extend(
            root / suffix
            for suffix in (
                "topology/core_id",
                "topology/thread_siblings_list",
                "cpufreq/scaling_driver",
                "cpufreq/scaling_governor",
                "cpufreq/energy_performance_preference",
            )
        )
        for cache in sorted((root / "cache").glob("index*")):
            files.extend(
                cache / field
                for field in (
                    "level",
                    "type",
                    "size",
                    "ways_of_associativity",
                    "coherency_line_size",
                    "shared_cpu_list",
                )
            )
    package_code = (
        "import importlib.metadata,json,sys; "
        "print(json.dumps({'executable':sys.executable,'version':sys.version,"
        "'packages':dict(sorted((d.metadata['Name'],d.version) "
        "for d in importlib.metadata.distributions()))},sort_keys=True))"
    )
    return {
        "captured_utc": harness.now(),
        "allowed_cpus": sorted(os.sched_getaffinity(0)),
        "environment": {key: os.environ.get(key) for key in sorted(ENV_KEYS)},
        "system_files": {
            str(path): path.read_text() for path in files if path.is_file()
        },
        "commands": [
            capture(argv)
            for argv in (
                ["uname", "-a"],
                ["lscpu"],
                ["rustc", "-Vv"],
                ["cargo", "-V"],
                ["cc", "--version"],
                ["ld", "--version"],
                ["perf", "--version"],
                [
                    "perf",
                    "stat",
                    "-e",
                    "cycles,instructions,cache-misses,branch-misses",
                    "--",
                    "/usr/bin/true",
                ],
                ["perf", "stat", "-e", "task-clock", "--", "/usr/bin/true"],
                [sys.executable, "-c", package_code],
                [str(REPO / ".venv/bin/python"), "-c", package_code],
                ["git", "status", "--short"],
                ["git", "diff", "--binary", "HEAD"],
            )
        ],
        "cargo_locked_packages": tomllib.loads((REPO / "Cargo.lock").read_text())[
            "package"
        ],
    }


def freeze(args: argparse.Namespace) -> None:
    output = args.snapshot.resolve()
    output.mkdir(parents=True, exist_ok=False)
    commit = (
        subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO).decode().strip()
    )
    if commit != args.expected_commit:
        raise RuntimeError(f"unexpected commit: {commit}")
    source = [record(path, REPO) for path in source_paths()]
    for row in source:
        destination = output / "source" / row["path"]
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPO / row["path"], destination)
    require_passed(verify_records(output / "source", source))
    manifest = args.root.resolve() / "workloads/manifest.json"
    workloads = json.loads(manifest.read_text())
    harness.verify_manifest(workloads)
    shutil.copyfile(manifest, output / "workload_manifest.json")
    audit_manifest = AUDIT / "manifest_final.json"
    audit = json.loads(audit_manifest.read_text())
    print(f"Verifying {len(audit['files'])} sealed audit files", flush=True)
    audit_check = verify_records(AUDIT, audit["files"])
    save(output / "audit_verification.json", audit_check)
    require_passed(audit_check)
    shutil.copyfile(audit_manifest, output / "audit_manifest_final.json")
    shutil.copyfile(
        AUDIT / "manifest_initial.json", output / "audit_manifest_initial.json"
    )
    initial = json.loads((AUDIT / "manifest_initial.json").read_text())
    audited_source = [
        dict(path=name, bytes=item["bytes"], sha256=item["sha256"])
        for name, item in initial["repository_files"].items()
        if name.startswith("src/") or name in {"Cargo.toml", "Cargo.lock"}
    ]
    audited_source_check = verify_records(REPO, audited_source)
    save(output / "audited_production_source_verification.json", audited_source_check)
    require_passed(audited_source_check)
    reports = []
    for relative in (
        "SCIENTIFIC_ACCURACY_AUDIT.md",
        "CLAIMS.md",
        "claims.json",
        "COMPLETION_CHECKLIST.md",
        "REPRODUCE.md",
        "geometry/REPORT.md",
        "kinetics/KINETICS_CONFORMERS.md",
        "kraken/REPORT.md",
        "models/REPORT.md",
    ):
        destination = output / "audit_reports" / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(AUDIT / relative, destination)
        reports.append(record(destination, output))
    save(output / "environment.json", machine())
    identity = {
        "schema_version": 1,
        "created_utc": harness.now(),
        "git_commit": commit,
        "classification": (
            "current audited implementation; scientific limitations retained"
        ),
        "scientifically_corrected_or_validated_baseline": False,
        "repository": str(REPO),
        "run_root": str(args.root.resolve()),
        "run_prefix": args.run_prefix,
        "source": source,
        "workload_manifest": record(manifest),
        "audit_manifest": record(audit_manifest),
        "audit_root": str(AUDIT),
        "audit_reports": reports,
        "snapshot_files": [
            record(path, output)
            for path in sorted(output.rglob("*"))
            if path.is_file() and "source" not in path.relative_to(output).parts
        ],
        "audit_binding_policy": (
            "All final-manifest payloads freshly hashed; "
            "original evidence remains in place."
        ),
        "thread_affinities": AFFINITIES,
        "native_repetitions": 7,
        "diagnostic_repetitions": 5,
        "warmups": 1,
        "timing_policy": (
            "Sequential warm-cache process launches; no concurrent builds "
            "or profilers; no sample exclusion."
        ),
    }
    save(output / "identity.json", identity)
    (output / "identity.sha256").write_text(
        harness.digest(output / "identity.json") + "\n"
    )
    print(f"Current-audited snapshot frozen: {output}", flush=True)


def verify(snapshot: Path, audit: bool = False) -> dict:
    identity_file = snapshot / "identity.json"
    if (
        harness.digest(identity_file)
        != (snapshot / "identity.sha256").read_text().strip()
    ):
        raise RuntimeError("snapshot identity changed")
    identity = json.loads(identity_file.read_text())
    checks = {
        "source_current": verify_records(REPO, identity["source"]),
        "source_frozen": verify_records(snapshot / "source", identity["source"]),
        "snapshot_files": verify_records(snapshot, identity["snapshot_files"]),
        "bound_manifests": verify_records(
            Path("/"), [identity["workload_manifest"], identity["audit_manifest"]]
        ),
    }
    harness.verify_manifest(
        json.loads((snapshot / "workload_manifest.json").read_text())
    )
    if audit:
        checks["audit"] = verify_records(
            Path(identity["audit_root"]),
            json.loads((snapshot / "audit_manifest_final.json").read_text())["files"],
        )
    if (snapshot / "builds.json").exists():
        if (
            harness.digest(snapshot / "builds.json")
            != (snapshot / "builds.sha256").read_text().strip()
        ):
            raise RuntimeError("build identity changed")
        builds = json.loads((snapshot / "builds.json").read_text())
        checks["binaries"] = verify_records(snapshot, builds["binaries"])
    for result in checks.values():
        require_passed(result)
    return {"verified_utc": harness.now(), "passed": True, "checks": checks}


def build(args: argparse.Namespace) -> None:
    snapshot = args.snapshot.resolve()
    verify(snapshot)
    if (snapshot / "builds.json").exists():
        raise RuntimeError("build evidence already exists")
    builds = {
        "created_utc": harness.now(),
        "status": "running",
        "commands": [],
        "binaries": [],
    }
    try:
        for kind in ("native", "diagnostic"):
            target = args.target_root.resolve() / kind
            if target == snapshot or snapshot in target.parents:
                raise RuntimeError("build target must be outside snapshot")
            argv = [
                "cargo",
                "build",
                "--release",
                "--locked",
                "--bin",
                "stericx",
                "--example",
                "profile_predictions",
                "--target-dir",
                str(target),
            ]
            env = dict(os.environ)
            env.pop("STERICX_PROFILE_PATH", None)
            if kind == "diagnostic":
                argv.extend(["--features", "profiling"])
                env["CARGO_PROFILE_RELEASE_DEBUG"] = "1"
            log = snapshot / f"build_{kind}.log"
            item = {
                "kind": kind,
                "argv": argv,
                "started_utc": harness.now(),
                "environment": {key: env.get(key) for key in sorted(ENV_KEYS)},
            }
            builds["commands"].append(item)
            print(f"Building {kind}: {target}", flush=True)
            with log.open("w") as output:
                result = subprocess.run(
                    argv, cwd=REPO, env=env, stdout=output, stderr=subprocess.STDOUT
                )
            item.update(
                returncode=result.returncode,
                finished_utc=harness.now(),
                log=record(log, snapshot),
            )
            if result.returncode:
                raise RuntimeError(f"{kind} build failed; see {log}")
            verify(snapshot)
            for relative in ("stericx", "examples/profile_predictions"):
                destination = snapshot / "bin" / kind / Path(relative).name
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(target / "release" / relative, destination)
                builds["binaries"].append(record(destination, snapshot))
        builds["status"] = "complete"
    except BaseException as error:
        builds.update(status="failed", error=str(error))
        raise
    finally:
        builds["finished_utc"] = harness.now()
        save(snapshot / "builds.json", builds)
        (snapshot / "builds.sha256").write_text(
            harness.digest(snapshot / "builds.json") + "\n"
        )


def run(args: argparse.Namespace) -> None:
    snapshot = args.snapshot.resolve()
    verify(snapshot)
    identity = json.loads((snapshot / "identity.json").read_text())
    builds = json.loads((snapshot / "builds.json").read_text())
    if builds["status"] != "complete":
        raise RuntimeError("complete frozen builds required")
    root = Path(identity["run_root"])
    prefix = identity["run_prefix"]
    kind = args.kind
    phase_dir = snapshot / "measurements" / kind
    phase_dir.mkdir(parents=True, exist_ok=False)
    state = {
        "status": "running",
        "created_utc": harness.now(),
        "kind": kind,
        "runs": [],
    }
    try:
        for threads, affinity in AFFINITIES.items():
            label = f"{prefix}_{kind}_t{threads}"
            argv = [
                sys.executable,
                str(REPO / "scripts/profile_stericx.py"),
                "--root",
                str(root),
                "run",
                "--binary",
                str(snapshot / "bin" / kind / "stericx"),
                "--label",
                label,
                "--threads",
                str(threads),
                "--affinity",
                affinity,
                "--warmups",
                "1",
                "--reps",
                "7" if kind == "native" else "5",
            ]
            if kind == "diagnostic":
                argv.extend(
                    [
                        "--profile-env",
                        "STERICX_PROFILE_PATH",
                        "--compare-to",
                        f"{prefix}_native_t{threads}",
                    ]
                )
            log = phase_dir / f"t{threads}.log"
            item = {
                "label": label,
                "argv": argv,
                "started_utc": harness.now(),
                "loadavg_before": Path("/proc/loadavg").read_text(),
            }
            state["runs"].append(item)
            save(phase_dir / "status.json", state)
            print(f"Starting {label}: all ten complete workloads", flush=True)
            with log.open("w") as output:
                result = subprocess.run(
                    argv, cwd=REPO, stdout=output, stderr=subprocess.STDOUT
                )
            item.update(
                returncode=result.returncode,
                finished_utc=harness.now(),
                loadavg_after=Path("/proc/loadavg").read_text(),
            )
            directory = root / "runs" / label
            if directory.exists():
                item["artifacts"] = [
                    record(path)
                    for path in sorted(directory.rglob("*"))
                    if path.is_file()
                ]
                shutil.copyfile(
                    directory / "summary.json", phase_dir / f"t{threads}_summary.json"
                )
            if result.returncode:
                raise RuntimeError(
                    f"{label} failed; original outputs retained at {directory}"
                )
            verify(snapshot)
            print(
                f"Completed {label}; frozen source and binaries still match", flush=True
            )
        state["status"] = "complete"
    except BaseException as error:
        state.update(status="failed", error=str(error))
        raise
    finally:
        state["finished_utc"] = harness.now()
        save(phase_dir / "status.json", state)
        (phase_dir / "status.sha256").write_text(
            harness.digest(phase_dir / "status.json") + "\n"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", type=Path, required=True)
    actions = parser.add_subparsers(dest="action", required=True)
    freezer = actions.add_parser("freeze")
    freezer.add_argument("--root", type=Path, default=harness.DEFAULT_ROOT)
    freezer.add_argument("--expected-commit", required=True)
    freezer.add_argument("--run-prefix", required=True)
    builder = actions.add_parser("build")
    builder.add_argument("--target-root", type=Path, required=True)
    runner = actions.add_parser("run")
    runner.add_argument("--kind", choices=("native", "diagnostic"), required=True)
    checker = actions.add_parser("verify")
    checker.add_argument("--audit", action="store_true")
    args = parser.parse_args()
    if args.action == "verify":
        print(json.dumps(verify(args.snapshot.resolve(), args.audit), indent=2))
    else:
        {"freeze": freeze, "build": build, "run": run}[args.action](args)


if __name__ == "__main__":
    main()
