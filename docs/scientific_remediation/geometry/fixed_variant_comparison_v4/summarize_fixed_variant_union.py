#!/usr/bin/env python3
"""Summarize independently recombined three-plane means at every fixed variant.

No residual becomes a tolerance. Raw historical variant semantics remain in the
union evidence; incompatible plane sets and unavailable SUT values stay explicit.
"""

import argparse
import json
import math
import shutil
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from sealed_reference import refuse_audit_output, sha


def rows(path):
    with path.open() as handle:
        for line in handle:
            yield json.loads(line)


def main(args):
    refuse_audit_output(args.output)
    args.output.mkdir(parents=True, exist_ok=False)
    shutil.copy2(__file__, args.output / Path(__file__).name)
    receipt_path = args.references / "complete.json"
    receipt = json.loads(receipt_path.read_text())
    if not receipt["complete"] or receipt["union_count"] != 2517:
        raise ValueError("Incomplete fixed historical union")
    inputs = {
        str(receipt_path): sha(receipt_path),
        str(Path(__file__)): sha(Path(__file__)),
    }
    stats, status, exceptions = defaultdict(list), Counter(), []
    with (args.output / "per_id_residuals.jsonl").open("x") as dst:
        for mode in ["inferred", "topology"]:
            directory = args.references / mode
            for name, expected in receipt["lanes"][mode]["artifacts"].items():
                path = directory / name
                if sha(path) != expected:
                    raise ValueError("Changed fixed reference evidence")
                inputs[str(path)] = expected
            pairs = {
                x["request"]["id"]: x for x in rows(directory / "selected_inputs.jsonl")
            }
            means = {x["id"]: x for x in rows(directory / "mean_recombination.jsonl")}
            count = 0
            for row in rows(directory / "reference_raw.jsonl"):
                count += 1
                id_ = row["id"]
                sut = pairs[id_]["sut"]["buried_volume"]
                for variant in [
                    "matched_sut_001",
                    "matched_sut_0001",
                    "primary_001",
                    "primary_0001",
                ]:
                    raw, reduced = row.get(variant), means[id_]["variants"].get(variant)
                    key = mode + "/" + variant
                    if (
                        raw is None
                        or "error" in raw
                        or reduced is None
                        or "error" in sut
                    ):
                        reason = (
                            "missing_reference"
                            if raw is None
                            else "reference_error"
                            if "error" in raw
                            else "sut_error"
                        )
                        status[key + "/" + reason] += 1
                        exceptions.append(
                            {
                                "mode": mode,
                                "id": id_,
                                "variant": variant,
                                "reason": reason,
                            }
                        )
                        continue
                    if not reduced["exact_three_topological_planes"]:
                        status[key + "/different_plane_set"] += 1
                        exceptions.append(
                            {
                                "mode": mode,
                                "id": id_,
                                "variant": variant,
                                "reason": "different_plane_set",
                                "plane_atoms": reduced["plane_atoms"],
                            }
                        )
                        continue
                    expected = dict(reduced["independent_three_plane_mean"])
                    expected["percent_buried_volume"] = (
                        expected["buried_volume"] * 100 / (4 * math.pi * 3.5**3 / 3)
                    )
                    expected.update(
                        {
                            field: raw[field]
                            for field in [
                                "qvbur_min",
                                "qvbur_max",
                                "ovbur_min",
                                "ovbur_max",
                                "max_delta_qvbur",
                            ]
                        }
                    )
                    residual = {
                        field: sut[field] - value for field, value in expected.items()
                    }
                    status[key + "/compared"] += 1
                    dst.write(
                        json.dumps(
                            {
                                "mode": mode,
                                "id": id_,
                                "variant": variant,
                                "independent": expected,
                                "residual": residual,
                            },
                            allow_nan=False,
                        )
                        + "\n"
                    )
                    for field, delta in residual.items():
                        stats[key + "/" + field].append(
                            {
                                "id": id_,
                                "residual": delta,
                                "sut": sut[field],
                                "independent": expected[field],
                            }
                        )
            if count != receipt["union_count"]:
                raise ValueError("Incomplete union rows")
    metrics, top = {}, {}
    for key, values in stats.items():
        errors = np.array([row["residual"] for row in values])
        metrics[key] = {
            "n": len(values),
            "max_abs": float(np.abs(errors).max()),
            "rms": float(np.sqrt(np.mean(errors**2))),
            "mean_abs": float(np.abs(errors).mean()),
            "median_abs": float(np.median(np.abs(errors))),
        }
        top[key] = sorted(values, key=lambda row: abs(row["residual"]), reverse=True)[
            :20
        ]
    for name, value in [
        ("metrics.json", metrics),
        ("top20.json", top),
        ("exceptions.json", exceptions),
    ]:
        (args.output / name).write_text(
            json.dumps(value, indent=2, allow_nan=False) + "\n"
        )
    complete = {
        "kind": "fixed_variant_mean_comparison",
        "global_scientific_pass": False,
        "thresholds_added": False,
        "inputs": inputs,
        "status_counts": dict(status),
        "artifacts": {p.name: sha(p) for p in args.output.iterdir() if p.is_file()},
    }
    (args.output / "complete.json").write_text(json.dumps(complete, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--references", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    main(parser.parse_args())
