#!/usr/bin/env python3
"""Review downloaded CI evidence, wheel RECORD hashes, and source archive identity."""

from __future__ import annotations

import argparse
import ast
import base64
import csv
import gzip
import hashlib
import io
import json
import re
import shutil
import tarfile
import zipfile
from datetime import UTC, datetime
from pathlib import Path


def sha_bytes(data):
    return hashlib.sha256(data).hexdigest()


def sha(path):
    return sha_bytes(path.read_bytes())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    raw, output = args.input.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    started = json.loads((raw / "started.json").read_text())
    commit, source_hashes = started["source_commit"], started["source_hashes"]
    supplement = json.loads((raw / "source_bindings.json").read_text())
    if supplement["source_commit"] != commit:
        raise ValueError("Source binding commits disagree")
    source_hashes.update(supplement["source_hashes"])
    failures = []
    runs = {}
    for label, expected_jobs in [("wheels", 13), ("rust", 3)]:
        run = json.loads((raw / (label + "_final.json")).read_text())
        if run["headSha"] != commit or len(run["jobs"]) != expected_jobs:
            failures.append({"run": label, "error": "source/job count mismatch"})
        if run["conclusion"] != "success" or any(
            j["conclusion"] != "success" for j in run["jobs"]
        ):
            failures.append({"run": label, "error": "run/job conclusion not success"})
        lines = (raw / (label + "_full.log")).read_text().splitlines()
        outcomes = [line for line in lines if re.search(r"\b\d+ passed in ", line)]
        if label == "wheels" and (
            len(outcomes) != expected_jobs
            or any(not re.search(r"\b71 passed in ", line) for line in outcomes)
        ):
            failures.append(
                {"run": label, "error": "installed test outcomes incomplete"}
            )
        runs[label] = {"metadata": run, "python_test_outcomes": outcomes}
    artifacts = []
    paths = sorted((raw / "artifacts").rglob("*.whl"))
    if len(paths) != 12:
        failures.append(
            {"error": "Expected twelve wheel artifacts", "actual": len(paths)}
        )
    for path in paths:
        entry = {
            "path": str(path.relative_to(raw)),
            "bytes": path.stat().st_size,
            "sha256": sha(path),
            "files": {},
            "python_source_matches": {},
            "python_source_serialization": {},
        }
        if "-cp311-abi3-" not in path.name:
            failures.append({"path": entry["path"], "error": "wheel ABI tag mismatch"})
        with zipfile.ZipFile(path) as wheel:
            names = wheel.namelist()
            records = [n for n in names if n.endswith(".dist-info/RECORD")]
            if len(records) != 1:
                raise ValueError("Expected one wheel RECORD")
            for name, digest, size in csv.reader(
                io.StringIO(wheel.read(records[0]).decode())
            ):
                data = wheel.read(name)
                entry["files"][name] = {"sha256": sha_bytes(data), "bytes": len(data)}
                if digest:
                    algorithm, encoded = digest.split("=", 1)
                    actual = (
                        base64.urlsafe_b64encode(hashlib.new(algorithm, data).digest())
                        .rstrip(b"=")
                        .decode()
                    )
                    if actual != encoded or str(len(data)) != size:
                        failures.append(
                            {
                                "path": entry["path"],
                                "member": name,
                                "error": "wheel RECORD mismatch",
                            }
                        )
                source = "python/" + name
                if source in source_hashes:
                    match = entry["files"][name]["sha256"] == source_hashes[source]
                    normalized = data.replace(b"\r\n", b"\n")
                    crlf_equivalent = (
                        not match
                        and name.endswith((".py", ".pyi"))
                        and b"\r\n" in data
                        and sha_bytes(normalized) == source_hashes[source]
                        and ast.dump(ast.parse(data), include_attributes=False)
                        == ast.dump(ast.parse(normalized), include_attributes=False)
                    )
                    entry["python_source_serialization"][source] = {
                        "exact_bytes": match,
                        "exact_crlf_to_lf_equivalence": crlf_equivalent,
                        "crlf_pairs": data.count(b"\r\n"),
                        "normalized_sha256": sha_bytes(normalized),
                        "git_source_sha256": source_hashes[source],
                        "raw_sha256": entry["files"][name]["sha256"],
                    }
                    match = match or crlf_equivalent
                    entry["python_source_matches"][source] = match
                    if not match:
                        failures.append(
                            {
                                "path": entry["path"],
                                "member": name,
                                "error": "Python source differs from CI commit",
                            }
                        )
            if not entry["python_source_matches"]:
                failures.append(
                    {"path": entry["path"], "error": "No Python source matches checked"}
                )
        artifacts.append(entry)
    sdists = list((raw / "artifacts").rglob("*.tar.gz"))
    if len(sdists) != 1:
        failures.append({"error": "Expected one source archive", "actual": len(sdists)})
    for path in sdists:
        entry = {
            "path": str(path.relative_to(raw)),
            "sha256": sha(path),
            "bytes": path.stat().st_size,
            "scientific_source_matches": {},
        }
        with tarfile.open(path) as tar:
            files = {m.name.split("/", 1)[1]: m for m in tar.getmembers() if m.isfile()}
            expected = {
                name: digest
                for name, digest in source_hashes.items()
                if name.startswith(("src/", "python/rust/", "python/tests/"))
            }
            for name, digest in expected.items():
                actual = (
                    sha_bytes(tar.extractfile(files[name]).read())
                    if name in files
                    else None
                )
                entry["scientific_source_matches"][name] = actual == digest
                if actual != digest:
                    failures.append(
                        {
                            "path": entry["path"],
                            "member": name,
                            "error": "source archive differs from CI commit",
                        }
                    )
        artifacts.append(entry)
    summary = {
        "created_utc": datetime.now(UTC).isoformat(),
        "reviewer_sha256": sha(Path(__file__)),
        "source_commit": commit,
        "runs": runs,
        "artifacts": artifacts,
        "failures": failures,
        "passed": not failures,
        "scope": "Actual wheel/source installation and Rust CI; no publication.",
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    review = {
        "source_commit": commit,
        "passed": not failures,
        "source_files_bound": len(source_hashes),
        "artifact_count": len(artifacts),
        "wheel_count": len(paths),
        "source_archive_count": len(sdists),
        "wheel_record_hashes_checked": sum(len(a.get("files", {})) for a in artifacts),
        "wheel_python_sources_checked": sum(
            len(a.get("python_source_matches", {})) for a in artifacts
        ),
        "source_archive_scientific_files_checked": sum(
            len(a.get("scientific_source_matches", {})) for a in artifacts
        ),
        "rust_run": {
            k: runs["rust"]["metadata"][k]
            for k in ["databaseId", "headSha", "status", "conclusion", "url", "jobs"]
        },
        "wheel_run": {
            k: runs["wheels"]["metadata"][k]
            for k in ["databaseId", "headSha", "status", "conclusion", "url", "jobs"]
        },
        "installed_python_test_outcomes": runs["wheels"]["python_test_outcomes"],
        "failures": failures,
        "release_claim": False,
    }
    (output / "ci_review.json").write_text(json.dumps(review, indent=2) + "\n")
    for path in raw.glob("*.json"):
        shutil.copyfile(path, output / path.name)
    for path in raw.glob("*.txt"):
        shutil.copyfile(path, output / path.name)
    for path in raw.glob("*.log"):
        with (output / (path.name + ".gz")).open("wb") as stream:
            with gzip.GzipFile(fileobj=stream, mode="wb", filename="", mtime=0) as z:
                z.write(path.read_bytes())
    manifest = {
        str(p.relative_to(output)): sha(p)
        for p in sorted(output.iterdir())
        if p.is_file()
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(
        json.dumps(
            {
                "source_commit": commit,
                "passed": not failures,
                "artifacts": len(artifacts),
                "failures": failures,
            },
            indent=2,
        )
    )
    return bool(failures)


if __name__ == "__main__":
    raise SystemExit(main())
