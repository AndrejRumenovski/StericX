#!/usr/bin/env python3
"""Recombine independent per-plane references under the corrected mean convention.

Accepts geometry reference JSONL, preserving old per-plane values and reporting
new-mean and old-first-frame residuals separately. This is evidence, not a gate.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sealed_reference import AUDIT, MANIFEST, refuse_audit_output, verify


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(args):
    refuse_audit_output(args.output)
    if args.references.resolve().is_relative_to(AUDIT.resolve()):
        verify(args.references)
    observations = (
        {
            row["id"]: row
            for row in map(json.loads, args.observations.read_text().splitlines())
        }
        if args.observations
        else {}
    )
    records = []
    with args.references.open() as handle:
        for line in handle:
            reference = json.loads(line)
            variants = {
                key: value
                for key, value in reference.items()
                if isinstance(value, dict) and "orientations" in value
            }
            if not variants:
                records.append(
                    {
                        "id": reference["id"],
                        "status": "missing independent three-frame volume",
                    }
                )
            for variant, volume in variants.items():
                frames = volume["orientations"]
                if len(frames) != 3:
                    records.append(
                        {
                            "id": reference["id"],
                            "variant": variant,
                            "status": "not three independent frames",
                        }
                    )
                    continue
                values = [frame.get("independent", frame) for frame in frames]
                keys = [
                    key
                    for key in [
                        "buried_volume",
                        "percent_buried_volume",
                        "near_vbur",
                        "far_vbur",
                    ]
                    if all(key in value for value in values)
                ]
                mean = {key: sum(value[key] for value in values) / 3 for key in keys}
                sut = observations.get(reference["id"], {}).get("buried_volume", {})
                records.append(
                    {
                        "id": reference["id"],
                        "variant": variant,
                        "status": "compared"
                        if sut and "error" not in sut
                        else "independent recombination only",
                        "independent_mean": mean,
                        "independent_first": values[0],
                        "original_variant_summary": {
                            key: volume.get(key) for key in keys
                        },
                        "mean_residual": {
                            key: sut[key] - value for key, value in mean.items()
                        }
                        if sut and "error" not in sut
                        else None,
                        "planes": [
                            frame.get("plane", frame.get("plane_atom"))
                            for frame in frames
                        ],
                        "reference_center": reference.get("sut_center")
                        if variant.startswith("matched_sut")
                        else reference.get("primary_center", volume.get("center")),
                    }
                )
    result = {
        "kind": "independent_three_frame_recombination",
        "global_scientific_pass": False,
        "inputs": {
            str(path): sha(path)
            for path in [
                args.references,
                args.observations,
                MANIFEST,
                Path(__file__),
                Path(__file__).with_name("sealed_reference.py"),
            ]
            if path
        },
        "records": records,
    }
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["references", "output"]:
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--observations", type=Path)
    main(parser.parse_args())
