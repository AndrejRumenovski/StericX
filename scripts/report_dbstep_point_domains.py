#!/usr/bin/env python3
"""Verify raw point-study evidence and export compact reproducible records."""

import argparse
import json
import shutil
from pathlib import Path

from analyze_dbstep_point_domains import record, save


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads((args.input / "manifest.json").read_text())
    for item in manifest:
        if record(item["path"]) != item:
            raise ValueError(f"Evidence changed: {item['path']}")
    args.output.mkdir(parents=True, exist_ok=False)
    rows = [
        json.loads(path.read_text())
        for path in sorted((args.input / "cases").rglob("plane_*.json"))
    ]
    save(args.output / "frames.json", rows)
    for name in ["summary.json", "pre_execution.json", "manifest.json"]:
        shutil.copyfile(args.input / name, args.output / name)
    save(
        args.output / "receipt.json",
        {
            "raw_manifest": record(args.input / "manifest.json"),
            "raw_artifacts_verified": len(manifest),
            "raw_bytes_verified": sum(item["bytes"] for item in manifest),
            "frames": len(rows),
            "controller": record(__file__),
            "copied_outputs": [record(p) for p in sorted(args.output.iterdir())],
            "release_admission_scored": False,
        },
    )
    print(json.dumps({"verified": len(manifest), "frames": len(rows)}))


if __name__ == "__main__":
    main()
