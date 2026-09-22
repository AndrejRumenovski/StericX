#!/usr/bin/env python3
"""Freeze tracked scientific inputs and measure one StericX child at a time.

Preparation never executes StericX. Timings include process startup and output to
regular files. Linux wait4 supplies *per-child*, not cumulative, resource usage.
All generated inputs and raw outputs remain under .stericx/profiling. No source
optimization, cache flushing, network access, or numerical tolerance is used.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import shutil
import signal
import statistics
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULT_ROOT = REPO / ".stericx/profiling"
# Deliberately exact keys: never strip an arbitrary number or an entire JSON tree.
VOLATILE_STDOUT_KEYS = frozenset(
    {
        "weights_load_ms",
        "memory_map_ms",
        "prediction_latency_ms",
        "total_ms",
        "throughput_records_per_second",
        "rss_start_bytes",
        "rss_end_bytes",
        "rss_delta_bytes",
        "csv_and_geometry_ms",
        "geometry_compute_ms",
        "binary_export_ms",
        "sigpack_output",
        "database",
        "manifest",
        "build_seconds",
    }
)
CONTROLLED_ENV = {"LC_ALL": "C", "TZ": "UTC", "OMP_NUM_THREADS": "1"}
LAUNCHER_SOURCE = (REPO / "scripts/profile_child.c").read_text()


def build_launcher(directory: Path) -> Path:
    """A tiny native parent avoids inheriting Python's large RSS high-water mark.

    Linux records memory used immediately before exec in child rusage. Direct
    wait4 from a Python controller therefore overstates RSS of tiny targets even
    with posix_spawn. The native parent forks only after its small image starts.
    """
    source_hash = hashlib.sha256(LAUNCHER_SOURCE.encode()).hexdigest()
    directory.mkdir(parents=True, exist_ok=True)
    source = directory / f"wait4-{source_hash}.c"
    binary = source.with_suffix("")
    if not binary.is_file():
        source.write_text(LAUNCHER_SOURCE)
        command = [
            "cc",
            "-O2",
            "-std=c11",
            "-Wall",
            "-Wextra",
            str(source),
            "-o",
            str(binary),
        ]
        subprocess.run(command, check=True)
        write_json(
            binary.with_suffix(".identity.json"),
            {
                "source_sha256": source_hash,
                "binary_sha256": digest(binary),
                "compiler_identity": subprocess.check_output(
                    ["cc", "--version"]
                ).decode(),
                "compile_argv": command,
            },
        )
    return binary


def now() -> str:
    return datetime.now(UTC).isoformat()


def digest(path: Path) -> str:
    sha = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            sha.update(chunk)
    return sha.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def xyz_to_sdf(paths: list[Path]) -> str:
    """Keep original decimal coordinate tokens; four-decimal conversion is lossy.

    StericX's V2000 parser accepts whitespace-delimited atom fields. These blocks
    use that supported representation with zero explicit bonds; topology is
    inferred by exactly the same kernel used on XYZ files.
    """
    blocks = []
    for path in paths:
        lines = path.read_text().splitlines()
        count = int(lines[0])
        atoms = [line.split() for line in lines[2 : 2 + count]]
        if len(atoms) != count or any(len(atom) != 4 for atom in atoms):
            raise ValueError(f"malformed XYZ: {path}")
        body = [f"{x} {y} {z} {element} 0 0 0" for element, x, y, z in atoms]
        blocks.append(
            "\n".join(
                [
                    path.name,
                    "  StericX profiling",
                    "lossless original XYZ decimals",
                    f"{count:3d}  0  0  0  0  0            999 V2000",
                    *body,
                    "M  END",
                    "$$$$",
                    "",
                ]
            )
        )
    return "".join(blocks)


def repeat_sigpack(source: Path, target: Path, count: int) -> int:
    data = source.read_bytes()
    if not data or len(data) % 64:
        raise ValueError("expected nonempty native-endian sigpack v1, 64-byte records")
    source_count = len(data) // 64
    cycles, tail = divmod(count, source_count)
    with target.open("wb") as output:
        while cycles:
            chunk = min(cycles, 4096)
            output.write(data * chunk)
            cycles -= chunk
        output.write(data[: tail * 64])
    return source_count


def prepare(args: argparse.Namespace) -> None:
    root = args.root.resolve()
    work = root / "workloads"
    manifest_path = work / "manifest.json"
    if manifest_path.exists():
        verify_manifest(json.loads(manifest_path.read_text()))
        print(f"Existing immutable inputs verified: {manifest_path}")
        return
    if work.exists():
        raise RuntimeError(
            f"incomplete preparation exists; inspect before removing: {work}"
        )
    work.mkdir(parents=True)
    tracked = set(
        subprocess.check_output(["git", "ls-files", "-z"], cwd=REPO)
        .decode()
        .split("\0")
    )
    xyz_names = sorted(
        p for p in tracked if p.startswith("data/xyz/") and p.endswith(".xyz")
    )
    conformer_names = sorted(
        p for p in tracked if p.startswith("data/conformers/") and p.endswith(".xyz")
    )
    required = [
        "data/reactions_raw.csv",
        "data/reactions.sigpack",
        "data/ligand_db/kraken_phosphines.csv",
        "docs/study_001/stericx_portable_model.json",
        *xyz_names,
        *conformer_names,
    ]
    if not xyz_names or not conformer_names:
        raise RuntimeError("tracked XYZ and conformer inputs are required")
    sources = []
    for name in required:
        if name not in tracked:
            raise RuntimeError(f"input is not tracked: {name}")
        source, target = REPO / name, work / "sources" / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        sources.append(
            {"original": name, "path": str(target), "sha256": digest(target)}
        )
    generated = work / "generated"
    generated.mkdir()
    xyz_paths = [work / "sources" / name for name in xyz_names]
    conformers = [work / "sources" / name for name in conformer_names]
    atom_counts = {
        str(path): int(path.read_text().splitlines()[0]) for path in xyz_paths
    }
    by_size = sorted(xyz_paths, key=lambda path: (atom_counts[str(path)], str(path)))
    groups: dict[Path, list[Path]] = {}
    for path in conformers:
        groups.setdefault(path.parent, []).append(path)
    ensemble = max(groups.values(), key=lambda paths: (len(paths), str(paths[0])))
    sdf = generated / "ensemble.sdf"
    sdf.write_text(xyz_to_sdf(ensemble))

    batch = generated / "batch"
    batch.mkdir()
    batch_paths = []
    for index in range(args.descriptor_count):
        target = batch / f"ligand_{index:05d}.xyz"
        shutil.copyfile(conformers[index % len(conformers)], target)
        batch_paths.append(target)

    source_csv = work / "sources/data/reactions_raw.csv"
    with source_csv.open(newline="") as source:
        reader = csv.DictReader(source)
        rows = list(reader)
        fields = reader.fieldnames
    # Absolute frozen paths keep CSV resolution identical across result folders.
    for row in rows:
        for key in ("Ligand_XYZ_Path", "Conformer_XYZ_Paths"):
            row[key] = ";".join(
                str(work / "sources/data" / name)
                for name in row[key].split(";")
                if name
            )
            for name in row[key].split(";"):
                if name and not Path(name).is_file():
                    raise RuntimeError(f"unfrozen CSV coordinate reference: {name}")
    csv_paths = {}
    for name, count in (
        ("parse_ensembles", args.parse_count),
        ("screen", args.screen_count),
    ):
        target = generated / f"{name}.csv"
        with target.open("w", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=fields)
            writer.writeheader()
            for index in range(count):
                row = dict(rows[index % len(rows)])
                row["Reaction_ID"] += f"-PROFILE-{index:06d}"
                writer.writerow(row)
        csv_paths[name] = target
    packed = generated / "inference.sigpack"
    packed_source_count = repeat_sigpack(
        work / "sources/data/reactions.sigpack", packed, args.inference_count
    )
    model = work / "sources/docs/study_001/stericx_portable_model.json"
    weights = generated / "weights.json"
    write_json(weights, json.loads(model.read_text())["weights"])
    database = work / "sources/data/ligand_db/kraken_phosphines.csv"
    with database.open(newline="") as source:
        database_count = sum(1 for _ in csv.DictReader(source))
    workloads = {}

    def descriptors(name: str, paths: list[Path], **extra: object) -> None:
        workloads[name] = {
            "argv": ["descriptors", *map(str, paths), "--format", "json"],
            "kind": "descriptors",
            "expected_records": len(paths),
            **extra,
        }

    descriptors("single_small", [by_size[0]], atoms=atom_counts[str(by_size[0])])
    descriptors("single_large", [by_size[-1]], atoms=atom_counts[str(by_size[-1])])
    descriptors("ensemble_sdf", [sdf], expected_conformers=len(ensemble))
    descriptors("conformers_56", conformers)
    descriptors("descriptors_10000", batch_paths)
    workloads["parse_ensembles_1000"] = {
        "argv": [
            "parse",
            "--csv",
            str(csv_paths["parse_ensembles"]),
            "--xyz-dir",
            str(work / "sources/data"),
            "--output",
            "{run_dir}/parsed.sigpack",
        ],
        "kind": "parse",
        "expected_records": args.parse_count,
    }
    workloads["predict_1000000"] = {
        "argv": ["predict", "--data", str(packed), "--weights", str(weights)],
        "kind": "predict",
        "expected_records": args.inference_count,
        "exactness_scope": (
            "CLI preview and aggregate only; complete f32 prediction audit "
            "is required separately"
        ),
    }
    workloads["search_database"] = {
        "argv": [
            "search",
            "--similar-to",
            str(by_size[-1]),
            "--database",
            str(database),
            "--top",
            str(database_count),
            "--format",
            "json",
        ],
        "kind": "search",
        "expected_records": database_count,
    }
    workloads["screen_1000"] = {
        "argv": [
            "screen",
            str(model),
            str(csv_paths["screen"]),
            "--top",
            str(args.screen_count),
            "--format",
            "json",
        ],
        "kind": "screen",
        "expected_records": args.screen_count,
    }
    workloads["db_build_ensembles"] = {
        "argv": [
            "db",
            "build",
            "--source",
            str(work / "sources/data/conformers"),
            "--output",
            "{run_dir}/database.csv",
            "--group-by-parent",
            "--label-from",
            "parent",
        ],
        "kind": "db",
        "expected_records": len(groups),
        "expected_geometries": len(conformers),
    }
    files = [
        {"path": str(path), "bytes": path.stat().st_size, "sha256": digest(path)}
        for path in sorted(work.rglob("*"))
        if path.is_file()
    ]
    manifest = {
        "schema_version": 1,
        "created_utc": now(),
        "repo": str(REPO),
        "root": str(root),
        "git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO)
        .decode()
        .strip(),
        "sources": sources,
        "files": files,
        "workloads": workloads,
        "diversity": {
            "distinct_ligands": len(xyz_paths),
            "distinct_conformers": len(conformers),
            "packed_source_records": packed_source_count,
            "descriptor_batch_repetition": (
                "cyclic copies of sorted real conformer files"
            ),
            "reaction_batch_repetition": (
                "cyclic real reaction rows; unique IDs, original weights "
                "and decimal strings"
            ),
            "scientific_scope": (
                "Repeated throughput workloads do not represent "
                "10000 distinct chemistries."
            ),
        },
        "sdf_sources": list(map(str, ensemble)),
        "xyz_atom_counts": atom_counts,
        "defaults": {
            "rayon_threads": 1,
            "affinity": [2],
            "warmups": 1,
            "repetitions": 7,
        },
        "stdout_keys_excluded_from_comparison": sorted(VOLATILE_STDOUT_KEYS),
        "db_manifest_keys_excluded_from_comparison": ["build_seconds"],
        "cache_policy": (
            "Inputs hash-checked before run; first observed run retained separately; "
            "warm repeated process launches, no forced cache dropping or "
            "cold-cache claim."
        ),
    }
    write_json(manifest_path, manifest)
    print(
        f"Frozen {len(files)} input files and {len(workloads)} workloads: "
        f"{manifest_path}"
    )


def verify_manifest(manifest: dict) -> None:
    for item in manifest["files"]:
        path = Path(item["path"])
        if (
            not path.is_file()
            or path.stat().st_size != item["bytes"]
            or digest(path) != item["sha256"]
        ):
            raise RuntimeError(f"frozen input changed or missing: {path}")


def normalized_stdout(data: bytes, kind: str) -> bytes:
    if kind not in ("parse", "predict", "db"):
        return data
    return b"".join(
        line
        for line in data.splitlines(keepends=True)
        if line.decode().partition("=")[0] not in VOLATILE_STDOUT_KEYS
    )


def check_result(workload: dict, directory: Path) -> dict:
    stdout = (directory / "stdout.txt").read_bytes()
    stderr = (directory / "stderr.txt").read_text()
    if "skipped " in stderr or " skipped)" in stderr:
        raise RuntimeError(f"workload skipped inputs: {directory}")
    kind, expected = workload["kind"], workload["expected_records"]
    if kind in ("descriptors", "search", "screen"):
        parsed = json.loads(stdout)
        if kind == "descriptors":
            if len(parsed) != expected:
                raise RuntimeError(f"descriptor count mismatch: {directory}")
            if (
                "expected_conformers" in workload
                and parsed[0]["conformers"] != workload["expected_conformers"]
            ):
                raise RuntimeError(f"conformer count mismatch: {directory}")
        elif kind == "screen":
            if (
                parsed["screened"] != expected
                or parsed["returned"] != expected
                or parsed["skipped"] != 0
                or parsed["excluded"]
            ):
                raise RuntimeError(f"incomplete screen: {directory}")
        elif parsed["library_size"] != expected or len(parsed["hits"]) != expected:
            raise RuntimeError(f"incomplete database search: {directory}")
    else:
        values = dict(
            line.split("=", 1) for line in stdout.decode().splitlines() if "=" in line
        )
        key = {
            "parse": "records_processed",
            "predict": "records_predicted",
            "db": "ligands",
        }[kind]
        if int(values[key]) != expected:
            raise RuntimeError(f"record count mismatch: {directory}")
        if kind == "db" and (
            int(values["geometries_skipped"]) != 0
            or int(values["geometries_featurized"]) != workload["expected_geometries"]
        ):
            raise RuntimeError(f"incomplete database build: {directory}")
    output = {
        "stdout_sha256": hashlib.sha256(normalized_stdout(stdout, kind)).hexdigest()
    }
    if kind == "parse":
        packed = directory / "parsed.sigpack"
        if packed.stat().st_size != expected * 64:
            raise RuntimeError(f"sigpack byte count mismatch: {directory}")
        output["sigpack_sha256"] = digest(packed)
    if kind == "db":
        output["table_sha256"] = digest(directory / "database.csv")
        manifest = json.loads((directory / "database.manifest.json").read_text())
        del manifest["build_seconds"]
        output["manifest_sha256"] = hashlib.sha256(
            json.dumps(manifest, sort_keys=True).encode()
        ).hexdigest()
    return output


def measure(
    argv: list[str],
    env: dict[str, str],
    affinity: set[int],
    directory: Path,
    launcher: Path | None = None,
) -> dict:
    """Native wait4 measures that exact child, after a small clean parent forks."""
    launcher = launcher or build_launcher(directory.parent / "launcher")
    directory.mkdir(parents=True)
    usage_path = directory / "native_usage.json"
    with (
        (directory / "stdout.txt").open("wb") as stdout,
        (directory / "stderr.txt").open("wb") as stderr,
    ):
        previous = os.sched_getaffinity(0)
        os.sched_setaffinity(0, affinity)
        started = time.perf_counter_ns()
        try:
            process = subprocess.Popen(
                [str(launcher), str(usage_path), *argv],
                env=env,
                stdin=subprocess.DEVNULL,
                stdout=stdout,
                stderr=stderr,
                start_new_session=True,
            )
        finally:
            os.sched_setaffinity(0, previous)
        try:
            process.wait()
        except BaseException:
            os.killpg(process.pid, signal.SIGTERM)
            process.wait()
            raise
        controller_elapsed = time.perf_counter_ns() - started
    if not usage_path.is_file():
        raise RuntimeError(
            f"native resource launcher failed ({process.returncode}): {directory}"
        )
    metrics = json.loads(usage_path.read_text())
    metrics.update(
        {
            "argv": argv,
            "launcher_pid": process.pid,
            "controller_wall_ns": controller_elapsed,
            "cpu_percent_one_core": (metrics["cpu_user_s"] + metrics["cpu_system_s"])
            / (metrics["wall_ns"] / 1e9)
            * 100,
            "raw_stdout_sha256": digest(directory / "stdout.txt"),
            "raw_stderr_sha256": digest(directory / "stderr.txt"),
        }
    )
    write_json(directory / "metrics.json", metrics)
    if process.returncode != 0 or metrics["returncode"] != 0:
        raise RuntimeError(
            f"child failed ({process.returncode}); see {directory / 'stderr.txt'}"
        )
    return metrics


def summarize(records: list[dict]) -> dict:
    def distribution(key: str) -> dict:
        values = sorted(record[key] for record in records)
        median = statistics.median(values)
        return {
            "median": median,
            "mad": statistics.median(abs(v - median) for v in values),
            "minimum": min(values),
            "maximum": max(values),
            "p90_nearest_rank": values[math.ceil(0.9 * len(values)) - 1],
        }

    return {
        key: distribution(key)
        for key in (
            "wall_ns",
            "cpu_user_s",
            "cpu_system_s",
            "cpu_percent_one_core",
            "peak_rss_bytes",
            "minor_faults",
            "major_faults",
            "voluntary_context_switches",
            "involuntary_context_switches",
        )
    }


def run(args: argparse.Namespace) -> None:
    if sys.platform != "linux":
        raise RuntimeError(
            "this harness requires Linux wait4 RSS units and sched_setaffinity"
        )
    manifest_path = args.root.resolve() / "workloads/manifest.json"
    manifest = json.loads(manifest_path.read_text())
    verify_manifest(manifest)
    selected = (
        args.workloads.split(",") if args.workloads else list(manifest["workloads"])
    )
    if len(selected) != len(set(selected)) or any(
        name not in manifest["workloads"] for name in selected
    ):
        raise RuntimeError("unknown or duplicate workload name")
    if args.reps < 1 or args.warmups < 1 or args.threads < 1:
        raise RuntimeError("reps, warmups, and threads must be positive")
    affinity = set(map(int, args.affinity.split(",")))
    if not affinity or not affinity <= os.sched_getaffinity(0):
        raise RuntimeError("requested CPUs are unavailable")
    binary = args.binary.resolve(strict=True)
    launcher = build_launcher(args.root.resolve() / "launcher")
    target = args.root.resolve() / "runs" / args.label
    target.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ)
    env.update(CONTROLLED_ENV, RAYON_NUM_THREADS=str(args.threads))
    env.pop("STERICX_PROFILE_PATH", None)
    settings = {**CONTROLLED_ENV, "RAYON_NUM_THREADS": str(args.threads)}
    info = {
        "created_utc": now(),
        "binary": str(binary),
        "binary_sha256": digest(binary),
        "manifest_sha256": digest(manifest_path),
        "manifest": str(manifest_path),
        "harness_sha256": digest(Path(__file__)),
        "platform": platform.platform(),
        "launcher_sha256": digest(launcher),
        "launcher_source_sha256": hashlib.sha256(LAUNCHER_SOURCE.encode()).hexdigest(),
        "launcher_identity": json.loads(
            launcher.with_suffix(".identity.json").read_text()
        ),
        "measurement": (
            "Native parent fork/exec through exact-child wait4; process startup "
            "and output included. Python/launcher startup excluded from wall_ns "
            "and retained as controller_wall_ns."
        ),
        "python": sys.version,
        "workloads": selected,
        "environment": settings,
        "affinity": sorted(affinity),
        "repetitions": args.reps,
        "warmups": args.warmups,
        "profile_env": args.profile_env,
        "status": "running",
        "results": {},
        "cache_policy": manifest["cache_policy"],
    }
    comparison = None
    if args.compare_to:
        comparison = json.loads(
            (
                args.root.resolve() / "runs" / args.compare_to / "summary.json"
            ).read_text()
        )
        if comparison["manifest_sha256"] != info["manifest_sha256"]:
            raise RuntimeError("comparison used different frozen inputs")
    write_json(target / "summary.json", info)
    os.chdir(REPO)
    try:
        for name in selected:
            workload = manifest["workloads"][name]
            records, expected = [], None
            for index in range(args.warmups + args.reps):
                phase = "warmup" if index < args.warmups else "measured"
                repetition = index if phase == "warmup" else index - args.warmups
                directory = target / name / f"{phase}_{repetition:02d}"
                argv = [
                    str(binary),
                    *(
                        value.replace("{run_dir}", str(directory))
                        for value in workload["argv"]
                    ),
                ]
                child_env = dict(env)
                if args.profile_env:
                    child_env[args.profile_env] = str(directory / "stages.json")
                metrics = measure(argv, child_env, affinity, directory, launcher)
                if args.profile_env:
                    report_path = directory / "stages.json"
                    report = json.loads(report_path.read_text())
                    if (
                        report.get("schema_version") != 1
                        or report.get("dropped_scopes") != 0
                    ):
                        raise RuntimeError(
                            f"incomplete or unsupported stage profile: {directory}"
                        )
                    metrics["stages_sha256"] = digest(report_path)
                fingerprints = check_result(workload, directory)
                if expected is None:
                    expected = fingerprints
                if fingerprints != expected:
                    raise RuntimeError(f"numerical/output nondeterminism: {directory}")
                if (
                    comparison
                    and fingerprints != comparison["results"][name]["fingerprints"]
                ):
                    raise RuntimeError(
                        f"exact result differs from {args.compare_to}: {directory}"
                    )
                metrics["fingerprints"] = fingerprints
                metrics["phase"] = phase
                write_json(directory / "metrics.json", metrics)
                if phase == "measured":
                    records.append(metrics)
                print(
                    f"{args.label} {name} {phase} {repetition + 1}: "
                    f"{metrics['wall_ns'] / 1e6:.3f} ms; "
                    f"RSS {metrics['peak_rss_bytes']} bytes",
                    flush=True,
                )
            info["results"][name] = {
                "summary": summarize(records),
                "fingerprints": expected,
                "exactness_scope": workload.get(
                    "exactness_scope", "all emitted values byte-for-byte"
                ),
                "measured_runs": len(records),
            }
            write_json(target / "summary.json", info)
        info["status"] = "complete"
    except BaseException as error:
        info["status"] = "failed"
        info["error"] = str(error)
        raise
    finally:
        info["finished_utc"] = now()
        write_json(target / "summary.json", info)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    actions = parser.add_subparsers(dest="action", required=True)
    preparation = actions.add_parser(
        "prepare", help="freeze tracked inputs; never execute StericX"
    )
    preparation.add_argument("--descriptor-count", type=int, default=10_000)
    preparation.add_argument("--parse-count", type=int, default=1_000)
    preparation.add_argument("--inference-count", type=int, default=1_000_000)
    preparation.add_argument("--screen-count", type=int, default=1_000)
    runner = actions.add_parser("run", help="measure sequential warm process launches")
    runner.add_argument("--binary", type=Path, required=True)
    runner.add_argument("--label", required=True)
    runner.add_argument(
        "--workloads", help="comma-separated manifest names; default all"
    )
    runner.add_argument("--reps", type=int, default=7)
    runner.add_argument("--warmups", type=int, default=1)
    runner.add_argument("--threads", type=int, default=1)
    runner.add_argument("--affinity", default="2", help="comma-separated CPU IDs")
    runner.add_argument(
        "--compare-to", help="prior run label for strict output hash comparison"
    )
    runner.add_argument(
        "--profile-env",
        help="environment variable accepting a per-run stages.json path",
    )
    args = parser.parse_args()
    (prepare if args.action == "prepare" else run)(args)


if __name__ == "__main__":
    main()
