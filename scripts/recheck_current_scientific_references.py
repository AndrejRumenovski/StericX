#!/usr/bin/env python3
r"""Recalculate independent audit references against NEW live-source captures.

Run in the pinned science environment, after the observation helpers finish:
  .venv/bin/python3 scripts/recheck_current_scientific_references.py \
      --observations /tmp/candidate-observations --output /tmp/candidate-references \
      --components geometry kinetics kraken --cli-run /tmp/cli/runs/candidate \
      --cli-build-snapshot /tmp/candidate-build-snapshot

Every output root must be new and outside the sealed audit. Components may be
run separately. Kraken reruns its full-corpus analysis, Morfeus/analytic kernels,
and interpretation; it can take hours. This helper never builds or observes a
SUT. It refuses incomplete captures and never substitutes historical SUT output.
Exact per-value historical comparisons are diagnostics, not scientific PASS
claims. Known failures and source/experimental limitations remain in force.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import importlib.util
import json
import os
import shutil
import signal
import subprocess
import sys
import time
from datetime import UTC, datetime
from decimal import Decimal, InvalidOperation
from itertools import zip_longest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
AUDIT = REPO / "docs/scientific_accuracy_audit"
KINETICS_SCRIPT = "scripts/kinetics_audit.py"
GEOMETRY_SCRIPTS = (
    "scripts/geometry_reference.py",
    "scripts/geometry_run.py",
    "scripts/geometry_summarize.py",
)
KRAKEN_SCRIPTS = (
    "kraken_analyze.py",
    "kraken_morfeus_outliers.py",
    "kraken_full_analytic_reference.py",
    "kraken_interpret.py",
    "kraken_matched_percent.py",
    "kraken_delta_outliers.py",
    "kraken_finish_delta.py",
)
GEOMETRY_RESULTS = (
    "reference_results.jsonl",
    "descriptor_comparisons.csv",
    "bin_comparisons.csv",
    "failures.csv",
    "metrics.json",
    "invariance.csv",
    "connectivity.json",
)
KINETICS_RESULTS = (
    "kinetics.csv",
    "mmff_weights.csv",
    "crest_weights.csv",
    "aggregation.csv",
    "ee_to_ddg.csv",
    "native_packed_ensembles.csv",
    "ni_hda_target_temperature.csv",
    "same_selectivity_different_rates.json",
    "summary.json",
)
KRAKEN_RESULT_DIRS = (
    "analysis",
    "morfeus_outliers",
    "full_analytic_reference",
    "interpretation",
    "delta_outliers",
    "delta_interpretation",
)
# These are provenance-only files, never rows or scientific values.
PROVENANCE_FILES = frozenset(
    {
        "manifest.json",
        "manifest_frozen.json",
        "plan_frozen.json",
        "selection_frozen.json",
        "analytic_raw_frozen.json",
    }
)
LIMITATIONS = [
    "This is new independent numerical evaluation, not experimental validation.",
    "Known INCORRECT/UNCERTAIN audit classifications are not promoted by equality.",
    "Original source access, ensemble provenance, donor chemistry and "
    "predictive-validity limits remain.",
    "Focused/alignment/historical geometry supplements, model references, and "
    "Kraken focused rounding diagnostics are retained historical evidence, "
    "not recomputed by this helper.",
    "The additional kraken_all_bins observation lane is retained by the observer "
    "oracle; this helper does not add independent full-corpus regional quadrature "
    "beyond the original reference campaigns.",
    "Unchanged reference scripts contain historical explanatory prose. "
    "New numerical rows and provenance establish what was rerun.",
    "No observed maximum is used as a global tolerance. Every historical "
    "scientific value is compared exactly; changed values and missing records "
    "are retained.",
]


def sha(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read(path: Path):
    return json.loads(path.read_text())


def write_new(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def record(path: Path) -> dict:
    return {
        "path": str(path.resolve()),
        "bytes": path.stat().st_size,
        "sha256": sha(path),
    }


def verify(item: dict) -> Path:
    path = Path(item["path"]).resolve(strict=True)
    if not path.is_file() or sha(path) != item["sha256"]:
        raise ValueError(f"Changed captured file: {path}")
    if "bytes" in item and path.stat().st_size != item["bytes"]:
        raise ValueError(f"Changed captured size: {path}")
    return path


class SealedInputs:
    """Verify bytes against the original seal before linking/copying any input."""

    def __init__(self, audit: Path):
        self.audit = audit
        self.manifest = record(audit / "manifest_final.json")
        self.files = {
            r["path"]: r for r in read(audit / "manifest_final.json")["files"]
        }
        self.used: dict[str, dict] = {}

    def check(self, relative: str) -> Path:
        if relative not in self.files:
            raise ValueError(f"Input is not in the sealed audit: {relative}")
        item = self.files[relative]
        path = self.audit / relative
        actual = verify({**item, "path": str(path)})
        self.used[relative] = {**item, "path": str(actual)}
        return actual

    def copy(self, relative: str, destination: Path) -> None:
        source = self.check(relative)
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists() or destination.is_symlink():
            raise ValueError(f"Refusing to replace staged file: {destination}")
        shutil.copyfile(source, destination)

    def link_input_tree(self, relative: str, destination: Path) -> None:
        # Only scientific input trees are allowed, never prior SUT/derived trees.
        allowed = {"kraken/primary", "kraken/prepared_all", "kraken/raw_sources"}
        if relative not in allowed:
            raise ValueError(f"Not a permitted immutable input tree: {relative}")
        names = [p for p in self.files if p.startswith(relative + "/")]
        if not names:
            raise ValueError(f"Missing immutable tree: {relative}")
        for name in names:
            self.check(name)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.symlink_to(self.audit / relative, target_is_directory=True)


def reference_environment(sealed: SealedInputs) -> dict:
    initial = read(sealed.check("manifest_initial.json"))
    expected = initial["python"]["packages"]
    packages = {
        name: importlib.metadata.version(name)
        for name in (
            "numpy",
            "scipy",
            "pandas",
            "matplotlib",
            "rdkit",
            "morfeus-ml",
        )
    }
    if sys.version != initial["python"]["version"] or any(
        version != expected[name] for name, version in packages.items()
    ):
        raise ValueError(
            "Reference environment differs from the pinned audit environment"
        )
    spec = importlib.util.find_spec("morfeus")
    if spec is None or spec.origin is None:
        raise ValueError("Morfeus reference package is unavailable")
    root = Path(spec.origin).parent
    sources = {}
    for name, item in initial["reference_tools"]["morfeus-ml"][
        "python_source_files"
    ].items():
        sources[name] = record(root / name)
        if sources[name]["sha256"] != item["sha256"]:
            raise ValueError(f"Independent Morfeus source changed: {name}")
    return {
        "python": sys.version,
        "executable": sys.executable,
        "packages": packages,
        "morfeus_sources": sources,
    }


def observations(
    root: Path, audit: Path, components: list[str]
) -> tuple[dict, dict, dict]:
    manifest = read(root / "manifest.json")
    if manifest.get("kind") != "live_scientific_observations" or not manifest.get(
        "complete"
    ):
        raise ValueError("A complete NEW live-source observer capture is required")
    if root.is_relative_to(audit):
        raise ValueError(
            "Historical observations cannot be used as a candidate capture"
        )
    prepared = read(verify(manifest["oracle_manifest"]))
    if (
        prepared.get("kind") != "live_scientific_oracle_inputs"
        or Path(prepared["audit"]).resolve() != audit
    ):
        raise ValueError("Observer oracle uses a different audit")
    verify(manifest["helper"])
    verify(prepared["helper"])
    for item in prepared["locked_files"]:
        verify(item)
    evidence = read(verify(manifest["evidence_lock"]))
    if manifest["evidence_lock"] != prepared["evidence_lock"]:
        raise ValueError("Observer evidence lock differs from its prepared oracle")
    verify(evidence["sealed_manifest"])
    for item in evidence["files"]:
        verify(item)
    built = read(verify(manifest["build_manifest"]))
    if (
        built.get("kind") != "live_scientific_observer_build"
        or built.get("returncode") != 0
    ):
        raise ValueError(
            "Candidate capture lacks successful live-source build provenance"
        )
    if (
        built["oracle_manifest"] != manifest["oracle_manifest"]
        or built["binary"] != manifest["binary"]
    ):
        raise ValueError("Candidate capture and build provenance do not match")
    verify(built["binary"])
    for item in built["sources"] + built["adapter_files"]:
        verify(item)
    required = []
    for component in components:
        required.extend(
            {
                "geometry": ["geometry"],
                "kinetics": ["kinetics", "aggregation"],
                "kraken": ["kraken"],
            }[component]
        )
    for name in required:
        lane = manifest["lanes"].get(name)
        if lane is None or lane["returncode"] != 0 or lane.get("timed_out"):
            raise ValueError(f"Missing successful candidate observer lane: {name}")
        plan = prepared["lanes"][name]
        if (
            lane["requests"] != plan["requests"]
            or lane["coverage"]["rows"] != plan["expected_rows"]
        ):
            raise ValueError(
                f"Candidate input/coverage does not match locked plan: {name}"
            )
        for field in ("requests", "stdout", "stderr"):
            path = verify(lane[field])
            if field != "requests" and not path.is_relative_to(root):
                raise ValueError(
                    f"Candidate stream is not in its fresh capture root: {path}"
                )
    return manifest, built, evidence


def candidate_copy(item: dict, target: Path) -> None:
    source = verify(item)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        raise ValueError(f"Refusing to replace candidate evidence: {target}")
    shutil.copyfile(source, target)


def stage_geometry(sealed: SealedInputs, work: Path, observed: dict) -> list[list[str]]:
    for relative in GEOMETRY_SCRIPTS:
        sealed.copy(relative, work / relative)
    for name in ("inputs.json", "supplement_inputs.json"):
        sealed.copy("geometry/" + name, work / "geometry" / name)
    lane = observed["lanes"]["geometry"]
    out = work / "geometry"
    candidate_copy(lane["stdout"], out / "sut_outputs.jsonl")
    candidate_copy(lane["stderr"], out / "sut_stderr.txt")
    candidate_copy(lane["requests"], out / "requests.jsonl")
    write_new(
        out / "sut_freeze.json",
        {
            "kind": "new_candidate_capture_for_unchanged_independent_reference",
            "output_sha256": sha(out / "sut_outputs.jsonl"),
            "input_hashes": {
                name: sha(out / name)
                for name in ("inputs.json", "supplement_inputs.json")
            },
            "candidate_lane": lane,
        },
    )
    return [
        [str(work / "scripts/geometry_run.py"), "compare"],
        [str(work / "scripts/geometry_summarize.py")],
    ]


def stage_kinetics(
    sealed: SealedInputs,
    work: Path,
    observed: dict,
    built: dict,
    cli_run: Path | None,
    cli_snapshot: Path | None,
) -> tuple[list[list[str]], dict]:
    if cli_run is None:
        raise ValueError(
            "Kinetics requires --cli-run with new CLI/Python/packed captures"
        )
    if cli_run.is_relative_to(sealed.audit):
        raise ValueError(
            "Historical CLI captures cannot substitute for new candidate observations"
        )
    helper = REPO / "scripts/check_scientific_cli_equivalence.py"
    spec = importlib.util.spec_from_file_location("candidate_cli_oracle", helper)
    if spec is None or spec.loader is None:
        raise ValueError("CLI observation verifier is unavailable")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    cli_root = cli_run.parent.parent
    identity = read(cli_run / "identity.json")
    completed = read(cli_run / "complete.json")
    cli_plan = read(cli_root / "manifest.json")
    if identity["manifest_sha256"] != sha(cli_root / "manifest.json") or identity[
        "harness_sha256"
    ] != sha(helper):
        raise ValueError("CLI capture helper/input identity changed")
    module.verify(cli_plan)
    if sha(Path(identity["binary"])) != identity["binary_sha256"]:
        raise ValueError("Candidate CLI executable changed after capture")
    build_receipt = bind_cli_build(cli_snapshot, identity, built, sealed)
    for name, expected in identity["python_sut"].items():
        if sha(Path(built["source"]) / "scripts" / name) != expected:
            raise ValueError(
                f"Candidate Python capture does not match live source: {name}"
            )
    actual = {
        case["id"]: module.fingerprint(cli_run / case["id"])
        for case in cli_plan["cases"]
    }
    actual["python"] = sha(cli_run / "python/python.json")
    if actual != completed["results"]:
        raise ValueError("CLI/Python observations changed after capture")
    sealed.copy(KINETICS_SCRIPT, work / KINETICS_SCRIPT)
    for relative in sealed.files:
        if relative.startswith("kinetics/inputs/"):
            sealed.copy(relative, work / relative)
    sealed.copy("kinetics/inputs_manifest.json", work / "kinetics/inputs_manifest.json")
    # analyze() uses this one frozen primary target table, never Python SUT code.
    relative = "frozen/repository/data/official/ni_hda_kraken.csv"
    sealed.copy(relative, work / relative)
    out = work / "kinetics/frozen_outputs"
    out.mkdir()
    for name in ("kinetics", "aggregation"):
        candidate_copy(observed["lanes"][name]["stdout"], out / (name + ".jsonl"))
        candidate_copy(observed["lanes"][name]["stderr"], out / (name + ".stderr"))
    candidate_copy(record(cli_run / "python/python.json"), out / "python.json")
    runs = []
    for case in cli_plan["cases"]:
        if not case["id"].startswith("kinetics_"):
            continue
        directory = cli_run / case["id"]
        runs.append(
            {
                "id": case["id"].removeprefix("kinetics_"),
                "argv": [identity["binary"], *case["argv"]],
                **read(directory / "result.json"),
                "stdout": (directory / "stdout").read_text(),
                "stderr": (directory / "stderr").read_text(),
            }
        )
        for packed in (directory / "artifacts").glob("*.sigpack"):
            candidate_copy(record(packed), out / packed.name)
    expected_ids = {
        r["id"] for r in read(sealed.check("kinetics/frozen_outputs/cli_runs.json"))
    }
    if {r["id"] for r in runs} != expected_ids:
        raise ValueError("Incomplete new kinetics CLI cases")
    write_new(out / "cli_runs.json", runs)
    write_new(
        work / "kinetics/sut_manifest.json",
        {
            "description": "New candidate observations frozen before analysis",
            "files": [
                {**record(p), "path": str(p.relative_to(work))}
                for p in sorted(out.iterdir())
            ],
        },
    )
    provenance = {
        "run": str(cli_run),
        "identity": identity,
        "capture_files": [
            record(cli_run / name) for name in ("identity.json", "complete.json")
        ],
        "helper": record(helper),
        "build_receipt": build_receipt,
    }
    return [[str(work / KINETICS_SCRIPT), "analyze"]], provenance


def bind_cli_build(
    snapshot: Path | None,
    cli_identity: dict,
    observer_build: dict,
    sealed: SealedInputs,
) -> dict:
    """Bind the CLI executable to the same source as the live observer library."""
    if snapshot is None:
        raise ValueError(
            "Kinetics requires --cli-build-snapshot for executable/source provenance"
        )
    for name in ("identity", "builds"):
        if (
            sha(snapshot / (name + ".json"))
            != (snapshot / (name + ".sha256")).read_text().strip()
        ):
            raise ValueError(f"CLI build receipt changed: {name}")
    identity = read(snapshot / "identity.json")
    builds = read(snapshot / "builds.json")
    if builds["status"] != "complete" or any(
        r["returncode"] != 0 for r in builds["commands"]
    ):
        raise ValueError("CLI build receipt is incomplete")
    if identity["audit_manifest"]["sha256"] != sealed.manifest["sha256"]:
        raise ValueError("CLI build receipt is bound to a different scientific audit")
    binaries = [
        r for r in builds["binaries"] if r["sha256"] == cli_identity["binary_sha256"]
    ]
    if len(binaries) != 1:
        raise ValueError("CLI executable does not uniquely match its build receipt")
    verify({**binaries[0], "path": str(snapshot / binaries[0]["path"])})
    source = {r["path"]: r for r in identity["source"]}
    source_root = Path(observer_build["source"])
    checked = []
    for item in observer_build["sources"]:
        relative = str(Path(item["path"]).relative_to(source_root))
        expected = source.get(relative)
        if expected is None or (item["sha256"], item["bytes"]) != (
            expected["sha256"],
            expected["bytes"],
        ):
            raise ValueError(
                f"CLI and observer were built from different source: {relative}"
            )
        verify({**expected, "path": str(snapshot / "source" / relative)})
        checked.append(relative)
    for name, digest in cli_identity["python_sut"].items():
        relative = "scripts/" + name
        if source[relative]["sha256"] != digest:
            raise ValueError(f"Python SUT differs from build snapshot: {name}")
        verify({**source[relative], "path": str(snapshot / "source" / relative)})
    return {
        "snapshot": str(snapshot),
        "matched_binary": binaries[0],
        "matching_observer_source_files": checked,
        "receipts": [
            record(snapshot / name)
            for name in (
                "identity.json",
                "identity.sha256",
                "builds.json",
                "builds.sha256",
            )
        ],
        "scope": "CLI binary and observer source match the same build receipt",
    }


def stage_kraken(sealed: SealedInputs, work: Path, observed: dict) -> list[list[str]]:
    for name in KRAKEN_SCRIPTS:
        relative = "kraken/scripts/" + name
        sealed.copy(relative, work / relative)
    for relative in ("kraken/primary", "kraken/prepared_all", "kraken/raw_sources"):
        sealed.link_input_tree(relative, work / relative)
    out = work / "kraken/sut"
    out.mkdir()
    lane = observed["lanes"]["kraken"]
    candidate_copy(lane["stdout"], out / "stdout.jsonl")
    candidate_copy(lane["stderr"], out / "stderr.txt")
    write_new(
        out / "manifest_frozen.json",
        {
            "returncode": lane["returncode"],
            "frozen_before_reference_comparison": True,
            "scope": "New candidate SUT; independent references are recalculated",
            "raw_output_files": {
                name: record(out / name) for name in ("stdout.jsonl", "stderr.txt")
            },
            "candidate_lane": lane,
        },
    )
    return [[str(work / "kraken/scripts" / name)] for name in KRAKEN_SCRIPTS]


def normalize(value, old_root: Path, new_root: Path):
    # Traceback path spelling is not numerical evidence. No values are rounded.
    if isinstance(value, str):
        return value.replace(str(old_root), "<AUDIT_ROOT>").replace(
            str(new_root), "<AUDIT_ROOT>"
        )
    if isinstance(value, list):
        return [normalize(v, old_root, new_root) for v in value]
    if isinstance(value, dict):
        return {k: normalize(v, old_root, new_root) for k, v in value.items()}
    return value


def differences(before, after, path=""):
    if isinstance(before, dict) and isinstance(after, dict):
        for key in sorted(before.keys() | after.keys()):
            child = path + "/" + key.replace("~", "~0").replace("/", "~1")
            if key not in before or key not in after:
                yield {
                    "pointer": child,
                    "historical_present": key in before,
                    "candidate_present": key in after,
                    "historical": before.get(key),
                    "candidate": after.get(key),
                }
            else:
                yield from differences(before[key], after[key], child)
    elif isinstance(before, list) and isinstance(after, list):
        for index, (a, b) in enumerate(zip_longest(before, after)):
            if index >= len(before) or index >= len(after):
                yield {
                    "pointer": path + f"/{index}",
                    "historical_present": index < len(before),
                    "candidate_present": index < len(after),
                    "historical": a,
                    "candidate": b,
                }
            else:
                yield from differences(a, b, path + f"/{index}")
    elif type(before) is not type(after) or before != after:
        # int/float JSON representations may denote exactly the same number.
        numeric = type(before) in (int, float) and type(after) in (int, float)
        if not numeric or before != after:
            yield {"pointer": path, "historical": before, "candidate": after}


def scientific_files(work: Path, components: list[str]) -> list[Path]:
    result = []
    if "geometry" in components:
        result.extend(work / "geometry" / name for name in GEOMETRY_RESULTS)
    if "kinetics" in components:
        result.extend(work / "kinetics/results" / name for name in KINETICS_RESULTS)
    if "kraken" in components:
        for directory in KRAKEN_RESULT_DIRS:
            result.extend(
                p
                for p in sorted((work / "kraken" / directory).iterdir())
                if p.suffix in {".json", ".jsonl", ".csv"}
                and p.name not in PROVENANCE_FILES
            )
    return result


def compare_artifacts(
    sealed: SealedInputs, work: Path, components: list[str], out: Path
) -> dict:
    summaries = []
    with (out / "historical_value_differences.jsonl").open("x") as dst:
        for candidate in scientific_files(work, components):
            relative = str(candidate.relative_to(work))
            historical = sealed.check(relative)
            changes = rows = provenance_changes = 0
            with historical.open() as left, candidate.open() as right:
                if candidate.suffix == ".json":
                    pairs = [(json.load(left), json.load(right))]
                elif candidate.suffix == ".jsonl":
                    pairs = zip_longest(map(json.loads, left), map(json.loads, right))
                else:
                    pairs = zip_longest(csv.DictReader(left), csv.DictReader(right))
                for rows, (before, after) in enumerate(pairs, 1):
                    before = normalize(before, sealed.audit, work)
                    after = normalize(after, sealed.audit, work)
                    for delta in differences(before, after):
                        if (
                            relative
                            == "kraken/interpretation/matched_default_percent.json"
                            and delta["pointer"] == "/source_manifest_sha256"
                        ):
                            provenance_changes += 1
                            dst.write(
                                json.dumps(
                                    {
                                        "file": relative,
                                        "record": rows,
                                        "kind": "new_reference_manifest_identity",
                                        **delta,
                                    }
                                )
                                + "\n"
                            )
                            continue
                        # Compare exact decimal values, allowing CSV formatting.
                        if (
                            candidate.suffix == ".csv"
                            and "historical_present" not in delta
                        ):
                            a, b = delta["historical"], delta["candidate"]
                            try:
                                if (
                                    isinstance(a, str)
                                    and isinstance(b, str)
                                    and Decimal(a).is_finite()
                                    and Decimal(b).is_finite()
                                    and Decimal(a) == Decimal(b)
                                ):
                                    continue
                            except InvalidOperation:
                                pass
                        changes += 1
                        dst.write(
                            json.dumps(
                                {"file": relative, "record": rows, **delta},
                                allow_nan=False,
                            )
                            + "\n"
                        )
            summaries.append(
                {
                    "file": relative,
                    "records_compared": rows,
                    "changed_values_or_structure": changes,
                    "provenance_only_changes": provenance_changes,
                    "historical": record(historical),
                    "candidate": record(candidate),
                }
            )
    return {
        "comparison": "Exact scientific values/statuses; no tolerance or retuning",
        "normalization": [
            "audit root in text/tracebacks",
            "exactly equal decimal CSV formatting",
            "matched_default_percent.json source_manifest_sha256: provenance only",
        ],
        "files": summaries,
        "changed_values_or_structure": sum(
            r["changed_values_or_structure"] for r in summaries
        ),
    }


def run_commands(commands: list[list[str]], output: Path, timeout: float) -> list[dict]:
    log = output / "commands"
    log.mkdir()
    env = dict(
        os.environ,
        OPENBLAS_NUM_THREADS="1",
        OMP_NUM_THREADS="1",
        LC_ALL="C",
        TZ="UTC",
        PYTHONDONTWRITEBYTECODE="1",
        MPLCONFIGDIR=str(output / "matplotlib-cache"),
    )
    records = []
    for index, arguments in enumerate(commands):
        command = [sys.executable, *arguments]
        print(
            f"Reference {index + 1}/{len(commands)}: {Path(arguments[0]).name}",
            flush=True,
        )
        before = time.monotonic()
        with (
            (log / f"{index:02}.stdout").open("xb") as stdout,
            (log / f"{index:02}.stderr").open("xb") as stderr,
        ):
            try:
                process = subprocess.Popen(
                    command,
                    env=env,
                    stdout=stdout,
                    stderr=stderr,
                    start_new_session=True,
                )
                status = {
                    "returncode": process.wait(timeout=timeout),
                    "timed_out": False,
                }
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
                status = {"returncode": None, "timed_out": True}
            except BaseException:
                if "process" in locals() and process.poll() is None:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                raise
        row = {
            "argv": command,
            **status,
            "seconds": time.monotonic() - before,
            "stdout": record(log / f"{index:02}.stdout"),
            "stderr": record(log / f"{index:02}.stderr"),
        }
        records.append(row)
        write_new(log / f"{index:02}.json", row)
        if status["returncode"] != 0:
            raise RuntimeError(
                f"Reference failed; all evidence retained: {arguments[0]}"
            )
    return records


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--observations", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--audit", type=Path, default=AUDIT)
    parser.add_argument(
        "--components",
        nargs="+",
        choices=["geometry", "kinetics", "kraken"],
        required=True,
    )
    parser.add_argument("--cli-run", type=Path)
    parser.add_argument(
        "--cli-build-snapshot",
        type=Path,
        help="performance_freeze identity/builds receipts; required for kinetics",
    )
    parser.add_argument(
        "--command-timeout",
        type=float,
        default=21600,
        help="Per-reference subprocess timeout in seconds; Kraken is expensive",
    )
    args = parser.parse_args()
    if len(set(args.components)) != len(args.components) or args.command_timeout <= 0:
        parser.error("Components must be unique and timeout positive")
    audit = args.audit.resolve(strict=True)
    candidate = args.observations.resolve(strict=True)
    output = args.output.resolve()
    cli_run = args.cli_run.resolve(strict=True) if args.cli_run else None
    cli_snapshot = (
        args.cli_build_snapshot.resolve(strict=True)
        if args.cli_build_snapshot
        else None
    )
    protected = [audit, candidate] + ([cli_run.parent.parent] if cli_run else [])
    if any(output.is_relative_to(p) or p.is_relative_to(output) for p in protected):
        raise ValueError("Output must be separate from sealed audit and capture roots")
    output.mkdir(parents=True, exist_ok=False)
    start = {
        "kind": "new_independent_scientific_reference_recheck",
        "created_utc": datetime.now(UTC).isoformat(),
        "helper": record(Path(__file__)),
        "components": args.components,
        "candidate_observation_manifest": record(candidate / "manifest.json"),
        "limitations": LIMITATIONS,
    }
    write_new(output / "started.json", start)
    sealed = SealedInputs(audit)
    work = output / "recalculation"
    work.mkdir()
    commands = []
    result = {**start, "complete": False, "scientific_pass_claimed": False}
    try:
        observed, built, evidence = observations(candidate, audit, args.components)
        result["reference_environment"] = reference_environment(sealed)
        result["candidate_build_manifest"] = observed["build_manifest"]
        result["unchanged_claim_classifications"] = evidence["classification_counts"]
        if "geometry" in args.components:
            commands.extend(stage_geometry(sealed, work, observed))
        if "kinetics" in args.components:
            calls, provenance = stage_kinetics(
                sealed, work, observed, built, cli_run, cli_snapshot
            )
            commands.extend(calls)
            result["candidate_cli_provenance"] = provenance
        if "kraken" in args.components:
            commands.extend(stage_kraken(sealed, work, observed))
        write_new(
            output / "reference_plan.json",
            {
                "commands": commands,
                "input_files": list(sealed.used.values()),
                "original_sealed_manifest": sealed.manifest,
            },
        )
        result["commands"] = run_commands(commands, output, args.command_timeout)
        result["historical_comparison"] = compare_artifacts(
            sealed, work, args.components, output
        )
        # Detect any altered sources/inputs while reference processes were executing.
        for item in sealed.used.values():
            verify(item)
        verify(sealed.manifest)
        observations(candidate, audit, args.components)
        result["complete"] = True
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
    finally:
        result["original_sealed_manifest"] = sealed.manifest
        result["used_sealed_files"] = list(sealed.used.values())
        write_new(output / "result.json", result)
        # Linked immutable trees are locked separately; do not traverse them.
        files = []
        for directory, subdirs, names in os.walk(output):
            subdirs[:] = [d for d in subdirs if not (Path(directory) / d).is_symlink()]
            for name in sorted(names):
                path = Path(directory) / name
                if not path.is_symlink():
                    files.append(record(path))
        write_new(
            output / "manifest.json",
            {"kind": "independent_reference_recheck_files", "files": files},
        )
    print(
        json.dumps(
            {
                "output": str(output),
                "complete": result["complete"],
                "changed_values_or_structure": result.get(
                    "historical_comparison", {}
                ).get("changed_values_or_structure"),
                "scientific_pass_claimed": False,
                "error": result.get("error"),
            },
            indent=2,
        )
    )
    if not result["complete"]:
        raise SystemExit(1)
    if result["historical_comparison"]["changed_values_or_structure"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
