#!/usr/bin/env python3
"""Replay audited CLI/Python behavior without modifying the sealed audit.

Passing this oracle means unchanged behavior, including known scientific failures.
It does not turn an audited INCORRECT or UNCERTAIN claim into a scientific pass.
Only explicitly named process metrics and portable-model creation time are ignored.
All raw observations are retained. Output paths are identical between builds.
"""

from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import importlib.util
import json
import os
import shutil
import struct
import subprocess
import sys
from datetime import datetime
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
AUDIT = REPO / "docs/scientific_accuracy_audit"
METRICS = frozenset(
    {
        "fit_ms",
        "total_ms",
        "total_microseconds",
        "rss_start_bytes",
        "rss_end_bytes",
        "rss_delta_bytes",
        "csv_and_geometry_ms",
        "geometry_compute_ms",
        "binary_export_ms",
        "throughput_records_per_second",
    }
)


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def save(path: Path, value: object) -> None:
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def scientific_stdout(raw: bytes) -> bytes:
    return b"".join(
        line
        for line in raw.splitlines(keepends=True)
        if line.partition(b"=")[0].decode(errors="replace") not in METRICS
    )


def artifact_bytes(path: Path) -> bytes:
    raw = path.read_bytes()
    if path.name == "portable.json":
        value = json.loads(raw)
        # This exact field comes from SystemTime::now(), not the fitted science.
        datetime.fromisoformat(value["created"]["created_utc"].replace("Z", "+00:00"))
        value["created"]["created_utc"] = "<process creation time>"
        return json.dumps(value, sort_keys=True, allow_nan=False).encode()
    return raw


def fingerprint(directory: Path) -> dict:
    result = json.loads((directory / "result.json").read_text())
    return {
        "returncode": result["returncode"],
        "stderr": digest(directory / "stderr"),
        "stdout": hashlib.sha256(
            scientific_stdout((directory / "stdout").read_bytes())
        ).hexdigest(),
        "artifacts": {
            str(path.relative_to(directory / "artifacts")): hashlib.sha256(
                artifact_bytes(path)
            ).hexdigest()
            for path in sorted((directory / "artifacts").rglob("*"))
            if path.is_file()
        },
    }


def validate_outputs(case: dict, directory: Path) -> None:
    result = json.loads((directory / "result.json").read_text())
    if result["returncode"] != case["expected_status"]:
        raise ValueError(f"audited status changed: {case['id']}")
    required = set(case["expected_artifacts"])
    if case["expected_status"] == 0:
        required.update(Path(path).name for path in case["outputs"])
    for name in required:
        path = directory / "artifacts" / name
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"missing or empty required output: {case['id']}/{name}")
    if case["id"].startswith("deck_"):
        validate_deck(directory)
    if "batch_contract" in case:
        for stream in ("stdout", "stderr"):
            if (directory / stream).read_bytes() != case[
                f"historical_{stream}"
            ].encode():
                raise ValueError(
                    f"descriptor batch {stream} contract changed: {case['id']}"
                )


def descriptor_batch_stdout(rows: list[dict], output_format: str) -> str:
    """Render the reviewed CLI contract from audited values, never a new executable.

    JSON retains the audit's exact numeric spellings. Text/CSV formatting first
    restores each audited shortest-roundtrip number to its original f32 value.
    This avoids rounding a decimal approximation instead of the stored value.
    """
    if not rows:
        return ""
    if output_format == "json":
        records = []
        for row in rows:
            raw = row["historical_stdout"]
            if not raw.startswith("[\n  {\n") or not raw.endswith("\n  }\n]\n"):
                raise ValueError("audited single-file JSON layout changed")
            records.append(raw[len("[\n") : -len("\n]\n")])
        return "[\n" + ",\n".join(records) + "\n]\n"

    def number(record: dict, key: str, places: int) -> str:
        exact = struct.unpack("=f", struct.pack("=f", record[key]))[0]
        return f"{exact:.{places}f}"

    records = [json.loads(row["historical_stdout"])[0] for row in rows]
    if any(record["conformers"] != 1 for record in records):
        raise ValueError("batch oracle expects reviewed single-conformer fixtures")
    if output_format == "text":
        rendered = []
        for r in records:
            rendered.append(
                f"{r['file']}\n"
                f"  donor          {r['donor_element']} (atom {r['donor_index']})\n"
                f"  substituents   {', '.join(r['substituents'])}\n"
                "  conformers     1\n"
                f"  Sterimol      L {number(r, 'sterimol_l', 2)}   "
                f"B1 {number(r, 'sterimol_b1', 2)}   "
                f"B5 {number(r, 'sterimol_b5', 2)}   Å\n"
                f"  buried volume  Vbur {number(r, 'percent_buried_volume', 1)}%   "
                f"({number(r, 'buried_volume', 1)} Å³)\n"
                f"                 qvbur_min {number(r, 'qvbur_min', 2)}   "
                f"qvbur_max {number(r, 'qvbur_max', 2)}   "
                f"max_delta_qvbur {number(r, 'max_delta_qvbur', 2)}   Å³\n"
                f"  pyramidalization pyr_P {number(r, 'pyr_p', 3)}   "
                f"pyr_alpha {number(r, 'pyr_alpha', 2)}°\n"
            )
        return "\n".join(rendered)
    if output_format != "csv":
        raise ValueError(f"unsupported descriptor oracle format: {output_format}")
    columns = (
        "file",
        "conformers",
        "donor_element",
        "donor_index",
        "substituents",
        "sterimol_l",
        "sterimol_b1",
        "sterimol_b5",
        "percent_buried_volume",
        "buried_volume",
        "qvbur_min",
        "qvbur_max",
        "max_delta_qvbur",
        "max_delta_qvbur_min",
        "pyr_p",
        "pyr_alpha",
    )

    def csv_field(value: str) -> str:
        return (
            '"' + value.replace('"', '""') + '"'
            if any(char in value for char in ',"\n')
            else value
        )

    lines = [",".join(columns)]
    for r in records:
        fields = [
            csv_field(r["file"]),
            str(r["conformers"]),
            r["donor_element"],
            str(r["donor_index"]),
            csv_field(" ".join(r["substituents"])),
            *(number(r, key, 4) for key in columns[5:]),
        ]
        lines.append(",".join(fields))
    return "\n".join(lines) + "\n"


def descriptor_batch_cases(audited_cases: list[dict]) -> list[dict]:
    """Fixed order/duplicates/errors derived from sealed single-file observations."""
    by_id = {case["id"]: case for case in audited_cases}
    selections = {
        "mixed": (
            "nearly_collinear",
            "triphenylphosphine",
            "multiple_phosphorus",
            "methylphosphine",
            "triphenylphosphine",
            "multiple_phosphorus",
            "asymmetric_phosphine",
        ),
        "all_failures": ("nearly_collinear", "nearly_collinear", "multiple_phosphorus"),
        "donor_index_precheck": ("multiple_phosphorus", "triphenylphosphine"),
    }
    terminal_error = "error: no ligand files could be featurized\n"
    cases = []
    for kind, names in selections.items():
        rows = [by_id[f"geometry_{name}"] for name in names]
        successful = [row for row in rows if row["expected_status"] == 0]
        skips = []
        for row in rows:
            if row["argv"][0] != "descriptors" or row["argv"][2:] != [
                "--format",
                "json",
                "--donor-element",
                "P",
            ]:
                raise ValueError("audited descriptor options changed")
            if row["expected_status"] == 0:
                if (
                    row["historical_stderr"]
                    or len(json.loads(row["historical_stdout"])) != 1
                ):
                    raise ValueError("audited descriptor success contract changed")
            elif (
                row["expected_status"] == 2
                and not row["historical_stdout"]
                and row["historical_stderr"].startswith(f"skipped {row['argv'][1]}: ")
                and row["historical_stderr"].endswith(terminal_error)
            ):
                skips.append(row["historical_stderr"][: -len(terminal_error)])
            else:
                raise ValueError("audited descriptor failure contract changed")
        stderr = "".join(skips)
        if successful:
            stderr += (
                f"featurized {len(successful)} of {len(rows)} files "
                f"({len(skips)} skipped)\n"
            )
        else:
            stderr += terminal_error
        status = 0 if successful else 2
        if kind == "donor_index_precheck":
            successful = []
            status = 2
            stderr = (
                "error: --donor-index applies to a single file; "
                "omit it for batch runs\n"
            )
        for output_format in ("json", "text", "csv"):
            cases.append(
                {
                    "id": f"batch_{kind}_{output_format}",
                    "argv": [
                        "descriptors",
                        *(row["argv"][1] for row in rows),
                        "--format",
                        output_format,
                        "--donor-element",
                        "P",
                        *(
                            ["--donor-index", "0"]
                            if kind == "donor_index_precheck"
                            else []
                        ),
                    ],
                    "expected_status": status,
                    "historical_stdout": descriptor_batch_stdout(
                        successful, output_format
                    ),
                    "historical_stderr": stderr,
                    "batch_contract": {
                        "basis": (
                            "Sealed single-file audit values/errors and reviewed "
                            "frozen CLI emitters"
                        ),
                        "source_case_ids": [row["id"] for row in rows],
                        "successful_records": len(successful),
                        "skipped_records": 0
                        if kind == "donor_index_precheck"
                        else len(skips),
                        "precheck_before_file_processing": kind
                        == "donor_index_precheck",
                        "output_format": output_format,
                    },
                }
            )
    return cases


def validate_deck(directory: Path) -> None:
    deck = directory / "artifacts/deck.csv"
    metadata = json.loads((directory / "artifacts/deck.meta.json").read_text())
    report = json.loads((directory / "stdout").read_text())
    with deck.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if metadata["deck_sha256"] != digest(deck):
        raise ValueError("deck checksum does not match its metadata")
    if len(rows) != report["returned"] or len(rows) != len(report["hits"]):
        raise ValueError("deck does not contain every returned candidate")
    if metadata["provenance"] != report["provenance"]:
        raise ValueError("deck provenance differs from the screening report")
    for row, hit in zip(rows, report["hits"], strict=True):
        if int(row["rank"]) != hit["rank"] or row["ligand_id"] != hit["ligand"]:
            raise ValueError("deck ranking differs from the screening report")
        for key in (
            "predicted_ddg_kcal_mol",
            "predicted_ee_percent",
            "prediction_interval_low",
            "prediction_interval_high",
            "maximum_extrapolation",
            "nearest_training_distance",
            "nearest_training_threshold",
            "leverage",
        ):
            value = None if row[key] == "" else float(row[key])
            if value != hit[key]:
                raise ValueError(f"deck scientific value differs: {key}")
        for key in ("domain_verdict", "trust", "nearest_training_ligand"):
            if row[key] != (hit[key] or ""):
                raise ValueError(f"deck classification differs: {key}")
        uncertainty = hit["uncertainty"] or {}
        for key in ("lower", "upper", "level", "replicates"):
            value = row[f"uncertainty_{key}"]
            if (None if value == "" else float(value)) != uncertainty.get(key):
                raise ValueError(f"deck uncertainty differs: {key}")
        if row["uncertainty_method"] != (uncertainty.get("method") or ""):
            raise ValueError("deck uncertainty method differs")


def raw_hashes(run: Path, manifest: dict) -> dict:
    paths = [run / "python/python.json"]
    for case in manifest["cases"]:
        directory = run / case["id"]
        paths.extend(directory / name for name in ("stdout", "stderr", "result.json"))
        paths.extend(p for p in (directory / "artifacts").rglob("*") if p.is_file())
    return {str(p.relative_to(run)): digest(p) for p in sorted(paths)}


def verify_run(root: Path, manifest: dict, run: Path) -> tuple[dict, dict]:
    identity = json.loads((run / "identity.json").read_text())
    if identity["manifest_sha256"] != digest(root / "manifest.json"):
        raise ValueError(f"manifest changed after capture: {run}")
    if not (
        identity["harness_sha256"]
        == digest(root / "harness.py")
        == manifest["harness_sha256"]
    ):
        raise ValueError(f"captured harness identity differs: {run}")
    if identity["binary_sha256"] != digest(run / "bin/stericx"):
        raise ValueError(f"captured executable changed: {run}")
    for name, expected in identity["python_sut"].items():
        if digest(run / "python_source" / name) != expected:
            raise ValueError(f"captured Python source changed: {run}/{name}")
    complete = json.loads((run / "complete.json").read_text())
    observed = {}
    for case in manifest["cases"]:
        validate_outputs(case, run / case["id"])
        observed[case["id"]] = fingerprint(run / case["id"])
    observed["python"] = digest(run / "python/python.json")
    if (
        observed != complete["results"]
        or raw_hashes(run, manifest) != complete["raw_hashes"]
    ):
        raise ValueError(f"observation changed after capture: {run}")
    return identity, observed


def prepare(root: Path) -> None:
    root.mkdir(parents=True, exist_ok=False)
    (root / "scratch").mkdir()
    shutil.copyfile(Path(__file__), root / "harness.py")
    cases = []
    sources = [
        AUDIT / "geometry/focused/cli_outputs.json",
        AUDIT / "kinetics/frozen_outputs/cli_runs.json",
    ]
    for prefix, source in zip(("geometry", "kinetics"), sources, strict=True):
        for row in json.loads(source.read_text()):
            cases.append(
                {
                    "id": f"{prefix}_{row['id']}",
                    "argv": row.get("argv", row.get("command"))[1:],
                    "expected_status": row["returncode"],
                    "historical_stdout": row["stdout"],
                    "historical_stderr": row["stderr"],
                }
            )
    commands = sorted((AUDIT / "models/raw").rglob("command.json"))
    if len(commands) != 33 or len(cases) != 49:
        raise ValueError("audited CLI coverage changed; inspect before proceeding")
    cases.extend(descriptor_batch_cases(cases))
    sources.extend(
        AUDIT / "frozen/repository/src" / name for name in ("descriptors.rs", "main.rs")
    )
    sources.extend(commands)
    for command in commands:
        row = json.loads(command.read_text())
        cases.append(
            {
                "id": f"models_{command.parent.name}",
                "argv": row["argv"][1:],
                "expected_status": row["exit_code"],
                "historical_stdout": (command.parent / "stdout").read_text(),
                "historical_stderr": (command.parent / "stderr").read_text(),
            }
        )
        sources.extend((command.parent / "stdout", command.parent / "stderr"))
    for name in ("screen_v2_default", "screen_v2_diverse"):
        source = next(case for case in cases if case["id"] == f"models_{name}")
        cases.append(
            {
                "id": f"deck_{name}",
                "argv": [*source["argv"], "--export-deck", "deck.csv"],
                "expected_status": source["expected_status"],
            }
        )
    inputs = set(sources)
    for case in cases:
        argv = case["argv"]
        flags = {"--output", "--export-deck"}
        if argv[0] == "fit":
            flags.update(("--predictions", "--portable-model"))
        replacements = {}
        outputs = []
        expected_artifacts = {}
        for index, token in enumerate(argv):
            if index and argv[index - 1] in flags:
                target = root / "scratch" / case["id"] / Path(token).name
                replacements[token] = str(target)
                outputs.append(str(target))
                old = Path(token)
                if old.is_file():
                    inputs.add(old.resolve())
                    expected_artifacts[old.name] = str(old.resolve())
                continue
            path = Path(token)
            if path.is_file():
                inputs.add(path.resolve())
            elif path.is_dir():
                inputs.update(p.resolve() for p in path.rglob("*") if p.is_file())
        case["argv"] = [replacements.get(token, token) for token in argv]
        case["outputs"] = outputs
        case["expected_artifacts"] = expected_artifacts
        case["output_path_replacements"] = replacements
    # Inputs also include indirect XYZ files named within CSVs and Python fixtures.
    for directory in ("geometry/focused/xyz", "kinetics/inputs", "models/inputs"):
        inputs.update(p for p in (AUDIT / directory).rglob("*") if p.is_file())
    inputs.update(
        (
            AUDIT / "scripts/kinetics_audit.py",
            AUDIT / "kinetics/frozen_outputs/python.json",
            AUDIT / "CLAIMS.md",
            AUDIT / "SCIENTIFIC_ACCURACY_AUDIT.md",
        )
    )
    save(
        root / "manifest.json",
        {
            "scope": "Exact audited behavior, including documented failures",
            "ignored_process_metric_keys": sorted(METRICS),
            "ignored_json_fields": ["portable.json:/created/created_utc"],
            "cases": cases,
            "inputs": {str(path): digest(path) for path in sorted(inputs)},
            "harness_sha256": digest(Path(__file__)),
        },
    )
    print(f"Prepared {len(cases)} CLI cases and live Python fixture replay: {root}")


def verify(manifest: dict) -> None:
    if digest(Path(__file__)) != manifest["harness_sha256"]:
        raise ValueError(
            "harness changed; prepare a new oracle without overwriting the old one"
        )
    for name, expected in manifest["inputs"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"immutable input changed: {name}")


def python_observation(repo: Path, destination: Path) -> None:
    source = AUDIT / "scripts/kinetics_audit.py"
    spec = importlib.util.spec_from_file_location("audit_python_observation", source)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    # Reuse only the audited input driver. Its scientific SUT imports now use repo.
    module.FROZEN = repo
    module.OUT = destination
    module.observe_python()


def observe(root: Path, binary: Path, repo: Path, label: str, threads: int) -> None:
    manifest = json.loads((root / "manifest.json").read_text())
    verify(manifest)
    run = root / "runs" / label
    run.mkdir(parents=True, exist_ok=False)
    if digest(root / "harness.py") != manifest["harness_sha256"]:
        raise ValueError("frozen harness changed")
    (run / "bin").mkdir()
    shutil.copy2(binary, run / "bin/stericx")
    (run / "python_source").mkdir()
    for name in ("prepare_data.py", "stericx_quantum.py"):
        shutil.copy2(repo / "scripts" / name, run / "python_source" / name)
    identity = {
        "binary": str(binary),
        "binary_sha256": digest(binary),
        "manifest_sha256": digest(root / "manifest.json"),
        "harness_sha256": digest(Path(__file__)),
        "python_executable": sys.executable,
        "python_version": sys.version,
        "threads": threads,
        "python_sut": {
            name: digest(repo / "scripts" / name)
            for name in ("prepare_data.py", "stericx_quantum.py")
        },
    }
    save(run / "identity.json", identity)
    env = dict(os.environ, LC_ALL="C", TZ="UTC", RAYON_NUM_THREADS=str(threads))
    env.pop("STERICX_PROFILE_PATH", None)
    results = {}
    for case in manifest["cases"]:
        directory = run / case["id"]
        directory.mkdir()
        scratch = root / "scratch" / case["id"]
        scratch.mkdir(exist_ok=False)
        try:
            with (
                (directory / "stdout").open("xb") as stdout,
                (directory / "stderr").open("xb") as stderr,
            ):
                result = subprocess.run(
                    [str(binary), *case["argv"]],
                    cwd=REPO,
                    env=env,
                    stdout=stdout,
                    stderr=stderr,
                    check=False,
                )
            shutil.copytree(scratch, directory / "artifacts")
            save(directory / "result.json", {"returncode": result.returncode})
            validate_outputs(case, directory)
            results[case["id"]] = fingerprint(directory)
        finally:
            # Only the fresh, harness-owned case directory is removed.
            shutil.rmtree(scratch)
    python_dir = run / "python"
    python_dir.mkdir()
    python_observation(repo, python_dir)
    results["python"] = digest(python_dir / "python.json")
    verify(manifest)
    if digest(binary) != identity["binary_sha256"]:
        raise ValueError("executable changed while observing")
    for name, expected in identity["python_sut"].items():
        if digest(repo / "scripts" / name) != expected:
            raise ValueError("Python SUT changed while observing")
    save(
        run / "complete.json",
        {
            "results": results,
            "raw_hashes": raw_hashes(run, manifest),
            "scope": manifest["scope"],
        },
    )
    verify_run(root, manifest, run)
    print(f"Observed {len(results) - 1} CLI cases plus Python: {run}")


def compare(root: Path, baseline: str, candidate: str, output: Path) -> None:
    manifest = json.loads((root / "manifest.json").read_text())
    verify(manifest)
    before = root / "runs" / baseline
    after = root / "runs" / candidate
    checked = [verify_run(root, manifest, p) for p in (before, after)]
    identities = [item[0] for item in checked]
    for key in ("manifest_sha256", "threads", "python_version"):
        if identities[0][key] != identities[1][key]:
            raise ValueError(f"comparison environment differs: {key}")
    rows = [item[1] for item in checked]
    differences = [
        name
        for name in sorted(rows[0].keys() | rows[1].keys())
        if rows[0].get(name) != rows[1].get(name)
    ]
    save(
        output,
        {
            "baseline": baseline,
            "candidate": candidate,
            "differences": differences,
            "manifest_sha256": digest(root / "manifest.json"),
            "identities": identities,
            "rust_source_linkage": (
                "Requires separate immutable build receipts; this oracle checks "
                "executable identities and outputs."
            ),
        },
    )
    if differences:
        raise ValueError(f"scientific/behavioral differences: {differences}")
    print(
        f"Exact audited behavior retained for {len(rows[0])} lanes; "
        "known failures remain"
    )


def fidelity(root: Path, label: str, output: Path) -> None:
    """Check the new baseline against original audited observations, never retune it."""
    manifest = json.loads((root / "manifest.json").read_text())
    verify(manifest)
    run = root / "runs" / label
    verify_run(root, manifest, run)
    complete = json.loads((run / "complete.json").read_text())
    differences = []
    count = 0
    derived_count = 0
    for case in manifest["cases"]:
        if "historical_stdout" not in case:
            continue
        if "batch_contract" in case:
            derived_count += 1
        else:
            count += 1
        directory = run / case["id"]
        expected_stdout = case["historical_stdout"]
        expected_stderr = case["historical_stderr"]
        for old, new in case["output_path_replacements"].items():
            expected_stdout = expected_stdout.replace(old, new)
            expected_stderr = expected_stderr.replace(old, new)
        observed = fingerprint(directory)
        expected_artifacts = {
            name: hashlib.sha256(artifact_bytes(Path(source))).hexdigest()
            for name, source in case["expected_artifacts"].items()
        }
        expected = {
            "returncode": case["expected_status"],
            "stdout": hashlib.sha256(
                scientific_stdout(expected_stdout.encode())
            ).hexdigest(),
            "stderr": hashlib.sha256(expected_stderr.encode()).hexdigest(),
            "artifacts": expected_artifacts,
        }
        if observed != complete["results"][case["id"]] or observed != expected:
            differences.append(
                {"id": case["id"], "expected": expected, "observed": observed}
            )
    python = digest(run / "python/python.json")
    if python != digest(AUDIT / "kinetics/frozen_outputs/python.json"):
        differences.append({"id": "python", "observed": python})
    save(
        output,
        {
            "historical_cli_cases": count,
            "derived_descriptor_batch_cases": derived_count,
            "differences": differences,
            "manifest_sha256": digest(root / "manifest.json"),
            "capture_identity_sha256": digest(run / "identity.json"),
            "capture_complete_sha256": digest(run / "complete.json"),
        },
    )
    if differences:
        raise ValueError(
            f"historical fidelity differences: {[r['id'] for r in differences]}"
        )
    print(
        f"Historical observation fidelity: {count} audited + {derived_count} derived "
        "CLI cases plus Python; exact"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    sub = parser.add_subparsers(dest="action", required=True)
    sub.add_parser("prepare")
    observation = sub.add_parser("observe")
    observation.add_argument("--binary", type=Path, required=True)
    observation.add_argument("--repo", type=Path, default=REPO)
    observation.add_argument("--label", required=True)
    observation.add_argument("--threads", type=int, default=1)
    comparison = sub.add_parser("compare")
    comparison.add_argument("--baseline", required=True)
    comparison.add_argument("--candidate", required=True)
    comparison.add_argument("--output", type=Path, required=True)
    historical = sub.add_parser("fidelity")
    historical.add_argument("--label", required=True)
    historical.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    if args.action == "prepare":
        prepare(root)
        return
    with (root / "lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if args.action == "observe":
            if args.threads < 1 or Path(args.label).name != args.label:
                raise ValueError("invalid threads or run label")
            observe(
                root,
                args.binary.resolve(),
                args.repo.resolve(),
                args.label,
                args.threads,
            )
        elif args.action == "compare":
            compare(root, args.baseline, args.candidate, args.output)
        else:
            fidelity(root, args.label, args.output)


if __name__ == "__main__":
    main()
