#!/usr/bin/env python3
"""Verify raw receipts and require exact equality of legacy scientific streams."""

import argparse
import hashlib
import json
import sys
from datetime import UTC, datetime
from pathlib import Path


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    manifests = [json.loads(p.read_text()) for p in (args.baseline, args.candidate)]
    failures, lanes = [], {}
    checked = {}
    for label, manifest in zip(["baseline", "candidate"], manifests, strict=True):
        if not manifest.get("complete_capture"):
            failures.append(f"{label}: incomplete capture")
        for name, lane in manifest["lanes"].items():
            if lane["returncode"] != 0:
                failures.append(f"{label}/{name}: nonzero return code")
            for field, record in lane.items():
                if isinstance(record, dict) and "sha256" in record and "path" in record:
                    path = record["path"]
                    if path not in checked:
                        checked[path] = digest(path)
                    if checked[path] != record["sha256"]:
                        failures.append(
                            f"{label}/{name}/{field}: recorded hash invalid"
                        )
    baseline, candidate = manifests
    if baseline["lanes"].keys() != candidate["lanes"].keys():
        failures.append("lane membership changed")
    for name, old in baseline["lanes"].items():
        new = candidate["lanes"].get(name)
        if new is None:
            continue
        same = old["rows"] == new["rows"]
        for field in ["requests", "stdout", "stderr"]:
            same = same and old[field]["sha256"] == new[field]["sha256"]
        if not same:
            failures.append(f"{name}: input or output stream changed")
        lanes[name] = {
            "rows": new["rows"],
            "exact_match": same,
            "stdout_sha256": new["stdout"]["sha256"],
        }
    report = {
        "created_utc": datetime.now(UTC).isoformat(),
        "command": sys.argv,
        "script_sha256": digest(__file__),
        "baseline_manifest_sha256": digest(args.baseline),
        "candidate_manifest_sha256": digest(args.candidate),
        "total_rows": sum(lane["rows"] for lane in lanes.values()),
        "lanes": lanes,
        "verified_files": checked,
        "failures": failures,
        "exact_legacy_replay_pass": not failures,
        "interpretation": "Raw numerical preservation only; not independent physical or predictive validation.",  # noqa: E501
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps(
            {
                k: report[k]
                for k in ["total_rows", "failures", "exact_legacy_replay_pass"]
            }
        )
    )
    return bool(failures)


if __name__ == "__main__":
    raise SystemExit(main())
