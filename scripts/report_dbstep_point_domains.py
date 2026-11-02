#!/usr/bin/env python3
"""Verify raw point-study evidence and export compact reproducible records."""

import argparse
import json
import math
import shutil
from pathlib import Path

import numpy as np
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
    summaries = []
    with (args.input / "native/responses.jsonl").open() as stream:
        for line in stream:
            if '"dump":' not in line:
                continue
            captured = json.loads(line)
            frames = [row for row in rows if row["case_id"] == captured["id"]]
            if len(frames) != 3:
                raise ValueError("Public reduction requires all three native frames")
            radius = float(np.float32(frames[0]["config"]["sphere_radius"]))
            sphere_volume = 4 * math.pi * radius**3 / 3
            fraction = (
                math.fsum(
                    frame["dbstep_count"] / frame["point_count"] for frame in frames
                )
                / 3
            )
            for field, expected, limit in [
                (
                    "buried_volume",
                    sphere_volume * fraction,
                    32 * 2.0**-23 * sphere_volume,
                ),
                ("percent_buried_volume", 100 * fraction, 32 * 2.0**-23 * 100),
            ]:
                observed = float(np.float32(captured["buried_volume"][field]))
                summaries.append(
                    {
                        "case_id": captured["id"],
                        "field": field,
                        "native": observed,
                        "same_domain_reference": expected,
                        "residual": observed - expected,
                        "unchanged_arithmetic_limit": limit,
                        "within_old_arithmetic_limit": abs(observed - expected)
                        <= limit,
                    }
                )
    if len(summaries) != 2 * len(rows) // 3:
        raise ValueError("Missing public summaries")
    save(args.output / "public_summaries.json", summaries)
    for name in ["summary.json", "pre_execution.json", "manifest.json"]:
        shutil.copyfile(args.input / name, args.output / name)
    save(
        args.output / "receipt.json",
        {
            "raw_manifest": record(args.input / "manifest.json"),
            "raw_artifacts_verified": len(manifest),
            "raw_bytes_verified": sum(item["bytes"] for item in manifest),
            "frames": len(rows),
            "public_summary_comparisons": len(summaries),
            "public_summary_arithmetic_failures": sum(
                not row["within_old_arithmetic_limit"] for row in summaries
            ),
            "controller": record(__file__),
            "copied_outputs": [record(p) for p in sorted(args.output.iterdir())],
            "release_admission_scored": False,
        },
    )
    print(json.dumps({"verified": len(manifest), "frames": len(rows)}))


if __name__ == "__main__":
    main()
