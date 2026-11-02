#!/usr/bin/env python3
"""Run unchanged scientific controllers with new, source-bound evidence paths.

This is orchestration only. It does not change reference equations, tolerances,
scientific inputs or acceptance policy. Successful execution is not scientific
admission. All subprocess logs, failures, commands and durations are retained.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

import replay_scientific_remediation as replay

ROOT = Path(__file__).resolve().parents[1]


def now():
    return datetime.now(UTC).isoformat()


def commands(phase, build, output):
    python = str(ROOT / ".venv/bin/python")
    ref_python = str(ROOT / ".stericx/paton/reference/venv/bin/python")
    observe = output / "observations"
    cli = output / "cli"
    base = [python, "scripts/replay_scientific_remediation.py"]
    if phase == "engineering":
        return [
            ("rust_format", ["cargo", "fmt", "--all", "--check"], ROOT),
            (
                "clippy",
                [
                    "cargo",
                    "clippy",
                    "--all-targets",
                    "--all-features",
                    "--",
                    "-D",
                    "warnings",
                ],
                ROOT,
            ),
            ("rustdoc", ["cargo", "doc", "--no-deps", "--all-features"], ROOT),
            ("rust_tests", ["cargo", "test", "--all-targets", "--all-features"], ROOT),
            (
                "screening_regression",
                [
                    "cargo",
                    "test",
                    "--test",
                    "screening_regression",
                    "--",
                    "--nocapture",
                ],
                ROOT,
            ),
            (
                "python_tests",
                [
                    python,
                    "-m",
                    "unittest",
                    "discover",
                    "-s",
                    "tests",
                    "-p",
                    "test_*.py",
                ],
                ROOT,
            ),
            (
                "dbstep_upstream_tests",
                [ref_python, "-m", "pytest", "-q", "tests"],
                ROOT / ".stericx/paton/reference/DBSTEP",
            ),
            (
                "archived_morfeus_evidence",
                [python, "docs/benchmarks/morfeus_current/verify.py"],
                ROOT,
            ),
            (
                "retrospective_screening",
                [python, "studies/study_011_ligand_ranking.py", "--verify-only"],
                ROOT,
            ),
        ]
    if phase == "references":
        return [
            (
                "references_inferred",
                [
                    "taskset",
                    "-c",
                    "0,1",
                    *base,
                    "references",
                    "--observations",
                    str(observe),
                    "--components",
                    "geometry",
                    "kraken",
                    "--output",
                    str(output / "references_inferred"),
                ],
                ROOT,
            ),
            (
                "references_topology",
                [
                    "taskset",
                    "-c",
                    "0,1",
                    *base,
                    "references",
                    "--observations",
                    str(observe),
                    "--components",
                    "kraken",
                    "--kraken-lane",
                    "kraken_topology",
                    "--output",
                    str(output / "references_topology"),
                ],
                ROOT,
            ),
        ]
    if phase == "volume":
        prefix = "docs/scientific_remediation/geometry/"
        return [
            (
                "volume_reference",
                [
                    "taskset",
                    "-c",
                    "0,1",
                    python,
                    prefix + "recompute_volume_means.py",
                    "--requests",
                    str(observe / "derived_requests/kraken_topology_all_bins.jsonl"),
                    "--observations",
                    str(observe / "kraken_topology_all_bins/stdout.jsonl"),
                    "--expected-rows",
                    "31721",
                    "--workers",
                    "2",
                    "--output",
                    str(output / "volume_reference"),
                ],
                ROOT,
            ),
            (
                "all_bins",
                [
                    python,
                    prefix + "compare_all_volume_bins.py",
                    "--references",
                    str(output / "volume_reference"),
                    "--observations",
                    str(observe),
                    "--output",
                    str(output / "all_bins"),
                ],
                ROOT,
            ),
            (
                "fixed_union",
                [
                    "taskset",
                    "-c",
                    "0,1",
                    python,
                    prefix + "replay_fixed_variant_union.py",
                    "--observations",
                    str(observe),
                    "--inferred",
                    str(output / "references_inferred/kraken"),
                    "--topology",
                    str(output / "references_topology/kraken"),
                    "--workers",
                    "2",
                    "--output",
                    str(output / "fixed_union"),
                ],
                ROOT,
            ),
        ]
    if phase in ("cli", "models_thermo"):
        capture = [python, "scripts/capture_remediation_cli.py", "--root", str(cli)]
        result = [
            (
                "cli_prepare",
                [
                    *capture,
                    "prepare",
                    "--historical",
                    str(
                        ROOT
                        / ".stericx/profiling/scientifically_validated_optimization"
                        / "cli_oracle_v5"
                    ),
                ],
                ROOT,
            )
        ]
        for threads in (1, 2, 4, 6):
            result.append(
                (
                    f"cli_{threads}",
                    [
                        *capture,
                        "observe",
                        "--build",
                        str(build),
                        "--label",
                        f"current_{threads}",
                        "--threads",
                        str(threads),
                    ],
                    ROOT,
                )
            )
        for threads in (2, 4, 6):
            result.append(
                (
                    f"cli_compare_{threads}",
                    [
                        python,
                        "scripts/compare_corrected_cli.py",
                        "--root",
                        str(cli),
                        "--baseline",
                        "current_1",
                        "--candidate",
                        f"current_{threads}",
                        "--mode",
                        "threads",
                        "--output",
                        str(output / f"cli_comparison_{threads}"),
                    ],
                    ROOT,
                )
            )
        if phase == "cli":
            return result
        result = []
        demo = output / "model_demo"
        models = output / "models"
        result.extend(
            [
                (
                    "model_demo",
                    [
                        python,
                        "scripts/replay_paton_model_demo.py",
                        "--build",
                        str(build),
                        "--output",
                        str(demo),
                    ],
                    ROOT,
                ),
                (
                    "models",
                    [
                        python,
                        "docs/scientific_remediation/models/recheck.py",
                        "--build",
                        str(build),
                        "--observations",
                        str(observe),
                        "--cli-run",
                        str(cli / "runs/current_1"),
                        "--demo",
                        str(demo / "demo"),
                        "--input-bundle",
                        str(
                            ROOT
                            / ".stericx/scientific_remediation"
                            / "ni_hda_response_metadata_v2"
                        ),
                        "--output",
                        str(models),
                    ],
                    ROOT,
                ),
                (
                    "model_supplement",
                    [
                        python,
                        "scripts/replay_paton_model_supplement.py",
                        "--reference-run",
                        str(models),
                        "--build",
                        str(build),
                        "--cli-run",
                        str(cli / "runs/current_1"),
                        "--cli-plan",
                        str(cli),
                        "--demo",
                        str(demo / "demo"),
                        "--output",
                        str(output / "model_supplement"),
                    ],
                    ROOT,
                ),
                (
                    "thermodynamics",
                    [
                        python,
                        "docs/scientific_remediation/thermodynamics/recheck.py",
                        "--binary",
                        str(build / "bin/stericx"),
                        "--observer",
                        str(build / "bin/stericx-audit-observer"),
                        "--out",
                        str(output / "thermodynamics"),
                    ],
                    ROOT,
                ),
            ]
        )
        return result
    raise ValueError(phase)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--logs", type=Path, required=True)
    parser.add_argument(
        "--phase",
        choices=("engineering", "references", "volume", "cli", "models_thermo"),
        required=True,
    )
    args = parser.parse_args()
    args.build, args.output, args.logs = (
        p.resolve() for p in (args.build, args.output, args.logs)
    )
    phase_dir = args.output / "execution" / args.phase
    phase_dir.mkdir(parents=True, exist_ok=False)
    args.logs.mkdir(parents=True, exist_ok=True)
    built, _ = replay.verify_build(args.build / "manifest.json")
    for item in built["sources"]:
        relative = Path(item["path"]).relative_to(built["source"])
        if replay.oracle.sha(ROOT / relative) != item["sha256"]:
            raise ValueError(f"Current scientific source differs: {relative}")
    controller = phase_dir / "controller.py"
    with controller.open("xb") as archived:
        archived.write(Path(__file__).read_bytes())
    bindings = [
        replay.oracle.file_record(args.build / "manifest.json"),
        replay.oracle.file_record(controller),
    ]
    env = dict(
        os.environ,
        OPENBLAS_NUM_THREADS="1",
        OMP_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
        RAYON_NUM_THREADS="2",
        CARGO_BUILD_JOBS="2",
        RUSTDOCFLAGS="-D warnings",
        LC_ALL="C",
        TZ="UTC",
    )
    journal = dict(
        phase=args.phase,
        started_utc=now(),
        command=sys.argv,
        build=bindings[0],
        controller=bindings[1],
        complete=False,
        scientific_global_pass=False,
        runs=[],
        interpretation=(
            "Execution receipt only; independent artifact review is required."
        ),
    )
    path = phase_dir / "journal.json"

    def save():
        path.write_text(json.dumps(journal, indent=2) + "\n")

    save()
    for name, argv, cwd in commands(args.phase, args.build, args.output):
        stdout, stderr = (
            args.logs / f"{name}.stdout.log",
            args.logs / f"{name}.stderr.log",
        )
        row = dict(
            name=name, argv=argv, cwd=str(cwd), started_utc=now(), status="running"
        )
        journal["runs"].append(row)
        save()
        started = time.monotonic()
        with stdout.open("xb") as out, stderr.open("xb") as err:
            result = subprocess.run(argv, cwd=cwd, env=env, stdout=out, stderr=err)
        row.update(
            returncode=result.returncode,
            seconds=time.monotonic() - started,
            completed_utc=now(),
            status="completed" if result.returncode == 0 else "failed",
            stdout=replay.oracle.file_record(stdout),
            stderr=replay.oracle.file_record(stderr),
        )
        save()
        print(
            json.dumps({k: row[k] for k in ("name", "returncode", "seconds")}),
            flush=True,
        )
        if result.returncode:
            raise SystemExit(result.returncode)
    replay.verify_build(args.build / "manifest.json")
    for item in bindings:
        replay.oracle.verify(item)
    journal.update(complete=True, completed_utc=now())
    save()


if __name__ == "__main__":
    main()
