#!/usr/bin/env python3
"""Run admitted, matched native/Python/allocation measurements in one quiet window."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import benchmark_paton_extensions as benchmark
import profile_stericx as profile

ROOT = Path(__file__).resolve().parents[1]

PYTHON_INSPECTOR = """import hashlib,importlib.metadata,json,platform,sys,pathlib
import stericx
p=pathlib.Path(stericx.__file__).parent
packages=sorted([
    {'name':d.metadata['Name'],'version':d.version}
    for d in importlib.metadata.distributions()],key=lambda d:d['name'].lower())
files={str(f):hashlib.sha256(f.read_bytes()).hexdigest()
    for f in sorted(p.iterdir()) if f.is_file()}
print(json.dumps({'python':sys.version,'executable':sys.executable,
    'platform':platform.platform(),'packages':packages,'files':files}))
"""


def verify_python_preparation(receipt, executables):
    if not receipt["ready_for_quiet_benchmark"] or receipt["source_changes"]:
        raise ValueError("Python performance preparation is not admitted")
    for side, executable in executables.items():
        expected = receipt["environments"][side]
        if executable.absolute() != Path(expected["executable"]):
            raise ValueError("A different Python environment was supplied")
        observed = json.loads(
            subprocess.check_output(
                [str(executable.absolute()), "-c", PYTHON_INSPECTOR], cwd=ROOT
            )
        )
        if observed != expected:
            raise ValueError(f"Python interpreter/package/files changed: {side}")
        artifact = receipt["artifacts"][side]
        if benchmark.record(artifact["path"])["sha256"] != artifact["sha256"]:
            raise ValueError("Python wheel artifact changed")
    for name, checks in receipt["source_checks"].items():
        if benchmark.record(ROOT / name)["sha256"] != checks["current"]:
            raise ValueError(f"Python benchmark source changed: {name}")


def read(path):
    return Path(path).read_text()


def host_state():
    return {
        "utc": profile.now(),
        "processes": subprocess.check_output(
            ["ps", "-eo", "pid,ppid,comm,pcpu,pmem,stat", "--sort=-pcpu"], text=True
        ),
        "proc": {
            name: read(Path("/proc") / name)
            for name in [
                "stat",
                "meminfo",
                "diskstats",
                "loadavg",
                "pressure/cpu",
                "pressure/io",
                "pressure/memory",
            ]
        },
    }


def main():
    parser = argparse.ArgumentParser()
    for name in [
        "baseline",
        "candidate",
        "preparation",
        "python-baseline",
        "python-candidate",
        "python-preparation",
        "review",
        "engineering",
        "quiet-window-receipt",
        "output",
    ]:
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    quiet = json.loads(args.quiet_window_receipt.read_text())
    if not quiet.get("user_confirmed") or not quiet.get("agent_heavy_work_complete"):
        raise ValueError("Quiet window has not been confirmed")
    review = json.loads(args.review.read_text())
    if not review.get("descriptor_performance_admitted"):
        raise ValueError("Scientific timing admission is incomplete")
    prepared = json.loads((args.preparation / "receipt.json").read_text())
    if not prepared["complete"]:
        raise ValueError("Diagnostic builds are incomplete")
    for build in prepared["builds"].values():
        for item in [*build["sources"], *build["binaries"].values()]:
            benchmark.verify(item)
    python_preparation = json.loads(args.python_preparation.read_text())
    python_executables = {
        "baseline": args.python_baseline,
        "candidate": args.python_candidate,
    }
    verify_python_preparation(python_preparation, python_executables)
    protocol = (
        ROOT / "docs/validation/dbstep/resolution_analysis/PERFORMANCE_PROTOCOL.md"
    )
    workload = (
        ROOT
        / ".stericx/profiling/scientifically_validated_optimization"
        / "corrected_workloads_v2/workloads/manifest.json"
    )
    input_path = ROOT / "python/examples/data/ligand_401_ensemble.xyz"
    options_path = ROOT / "python/examples/data/ligand_401_options.json"
    bindings = [
        benchmark.record(p)
        for p in [
            Path(__file__),
            protocol,
            args.review,
            args.engineering,
            args.quiet_window_receipt,
            args.preparation / "receipt.json",
            args.python_preparation,
            workload,
            input_path,
            options_path,
            ROOT / "scripts/benchmark_paton_python.py",
            ROOT / "scripts/benchmark_paton_extensions.py",
        ]
    ]
    journal = {
        "bindings": bindings,
        "quiet_window": quiet,
        "host_before": host_state(),
        "commands": [],
        "complete": False,
        "python_campaign_order": ["baseline", "candidate", "candidate", "baseline"],
    }
    path = args.output / "journal.json"
    profile.write_json(path, journal)

    def run(name, command):
        log = args.output / (name + ".log")
        row = {"name": name, "argv": command, "started_utc": profile.now()}
        journal["commands"].append(row)
        profile.write_json(path, journal)
        start = time.perf_counter()
        with log.open("w") as stream:
            result = subprocess.run(
                command,
                cwd=ROOT,
                stdout=stream,
                stderr=subprocess.STDOUT,
                env={**os.environ, **profile.CONTROLLED_ENV},
            )
        row.update(
            returncode=result.returncode,
            elapsed_seconds=time.perf_counter() - start,
            log=benchmark.record(log),
        )
        profile.write_json(path, journal)
        if result.returncode:
            raise RuntimeError(f"Measurement failed: {name}; all output retained")
        print(f"Completed {name}", flush=True)

    binaries = {
        side: prepared["builds"][side]["binaries"] for side in ["baseline", "candidate"]
    }
    run(
        "native",
        [
            sys.executable,
            str(ROOT / "scripts/benchmark_paton_extensions.py"),
            "--baseline",
            str(args.baseline),
            "--candidate",
            str(args.candidate),
            "--workloads",
            str(workload),
            "--review",
            str(args.review),
            "--engineering",
            str(args.engineering),
            "--output",
            str(args.output / "native"),
            "--allocation-baseline",
            binaries["baseline"]["allocation_cli"]["path"],
            "--allocation-candidate",
            binaries["candidate"]["allocation_cli"]["path"],
        ],
    )
    for index, side in enumerate(journal["python_campaign_order"]):
        python = args.python_baseline if side == "baseline" else args.python_candidate
        run(
            f"python_{index}_{side}",
            [
                str(python.absolute()),
                str(ROOT / "scripts/benchmark_paton_python.py"),
                "--native",
                binaries[side]["ensemble_benchmark"]["path"],
                "--input",
                str(input_path),
                "--options",
                str(options_path),
                "--iterations",
                "100",
                "--output",
                str(args.output / f"python_{index}_{side}"),
            ],
        )
        journal["commands"][-1]["produced_records"] = [
            benchmark.record(p)
            for p in sorted((args.output / f"python_{index}_{side}").glob("*.json"))
        ]
        profile.write_json(path, journal)
    launcher = profile.build_launcher(args.output / "ensemble_launcher")
    journal["ensemble_allocation_runs"] = []
    for side in ["baseline", "candidate"]:
        directory = args.output / "ensemble_allocations" / side
        command = [
            binaries[side]["ensemble_allocations"]["path"],
            str(input_path),
            str(options_path),
        ]
        metrics = profile.measure(
            command,
            {**os.environ, **profile.CONTROLLED_ENV, "RAYON_NUM_THREADS": "1"},
            {2},
            directory,
            launcher,
        )
        profile.write_json(directory / "metrics.json", metrics)
        value = json.loads((directory / "stdout.txt").read_text())
        journal["ensemble_allocation_runs"].append(
            {
                "side": side,
                "metrics": benchmark.record(directory / "metrics.json"),
                "stdout": benchmark.record(directory / "stdout.txt"),
                "values": value,
            }
        )
        profile.write_json(path, journal)
    for item in bindings:
        benchmark.verify(item)
    for build in prepared["builds"].values():
        for item in [*build["sources"], *build["binaries"].values()]:
            benchmark.verify(item)
    verify_python_preparation(python_preparation, python_executables)
    journal.update(complete=True, host_after=host_state())
    profile.write_json(path, journal)


if __name__ == "__main__":
    main()
