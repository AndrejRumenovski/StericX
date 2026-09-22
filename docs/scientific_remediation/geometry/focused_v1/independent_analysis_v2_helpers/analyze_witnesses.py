#!/usr/bin/env python3
"""Independent equation replay for focused corrections; not a global pass gate.

Uses unchanged sealed analytic support, determinant, and per-plane cubature
references. Only the explicitly new symmetric mean is added here. The supplied
center is controlled to separate center convention from occupancy/projection.
No measured residual becomes a new tolerance.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from sealed_reference import MANIFEST, REFERENCE, load_reference, refuse_audit_output

reference = load_reference()


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return {row["id"]: row for row in map(json.loads, path.read_text().splitlines())}


def per_frame_mean(volume):
    """New scientific reduction applied to unchanged independent frame values."""
    frames = volume["orientations"]
    return {
        key: sum(frame["independent"][key] for frame in frames) / len(frames)
        for key in ["buried_volume", "percent_buried_volume", "near_vbur", "far_vbur"]
    }


def analyze(args):
    refuse_audit_output(args.output)
    requests, before, after = [
        read(path) for path in (args.requests, args.before, args.after)
    ]
    assert set(requests) == set(after) == set(before)
    output = {
        "kind": "focused_independent_geometry_equations",
        "global_scientific_pass": False,
        "inputs": {
            str(p): sha(p)
            for p in (
                args.requests,
                args.before,
                args.after,
                REFERENCE,
                MANIFEST,
                Path(__file__),
                Path(__file__).with_name("sealed_reference.py"),
            )
        },
        "coordinate_contract": (
            "Identical input f32 coordinates and SUT-declared radii for supports; "
            "sealed independent Bondi radii and float64 grids for volume reference. "
            "Controlled current center."
        ),
        "cases": [],
        "ordinary_change_distribution": {},
    }
    ordinary = [
        "methylphosphine",
        "dimethylphosphine",
        "trimethylphosphine",
        "triethylphosphine",
    ]
    changes = {}
    for id_, request in requests.items():
        old, new = before[id_], after[id_]
        result = {
            "id": id_,
            "explicit_errors": {
                key: value
                for key, value in new.items()
                if isinstance(value, dict) and "error" in value
            },
        }
        xyz = np.array(
            [atom["position"] for atom in request["atoms"]], dtype=np.float32
        ).astype(np.float64)
        radii = np.array(
            [atom["vdw_radius"] for atom in new["atoms"]], dtype=np.float64
        )
        if np.all(np.isfinite(xyz)) and np.all(np.isfinite(radii)):
            for key, dummy, base, excluded in [
                (
                    "sterimol_bond",
                    xyz[request.get("attach", request["donor"])],
                    request.get("sterimol_neighbor", request["reference"]),
                    request.get("attach", request["donor"]),
                ),
                ("sterimol_dummy_raw", new.get("center"), request["donor"], None),
            ]:
                if (
                    not isinstance(dummy, (list, np.ndarray))
                    or not isinstance(new.get(key), dict)
                    or "error" in new[key]
                ):
                    continue
                exact = reference.sterimol_exact(
                    xyz, radii, np.array(dummy), base, excluded
                )
                result[key] = {
                    "independent": exact,
                    "before_residual": {
                        k: old.get(key, {}).get(k) - exact[k]
                        if isinstance(old.get(key, {}).get(k), (int, float))
                        else None
                        for k in exact
                    },
                    "after_residual": {k: new[key][k] - exact[k] for k in exact},
                }
            if "error" not in new["pyramidalization"]:
                result["pyramidalization"] = reference.pyramidalization_ref(
                    request, xyz
                )
                exact = result["pyramidalization"]["independent"]
                result["pyramidalization"]["after_residual"] = {
                    key: new["pyramidalization"][key] - value
                    for key, value in exact.items()
                }
        # Bound the focused cubature work: ordinary molecules, PH3, and its
        # density series (exclude the 1e-4 lane already retained in the audit).
        if id_ in [
            *ordinary,
            "phosphine_PH3",
            "ph3_lens_0.1",
            "ph3_lens_0.01",
            "ph3_lens_0.001",
        ] and "error" not in new.get("buried_volume", {"error": True}):
            independent = reference.bv_refs(request, xyz, center=new["center"])
            mean = per_frame_mean(independent)
            result["buried_volume"] = {
                "per_plane_independent": independent["orientations"],
                "new_independent_mean": mean,
                "new_mean_residual": {
                    key: new["buried_volume"][key] - value
                    for key, value in mean.items()
                },
                "old_first_frame_reference": independent["independent"],
            }
        if id_ in ordinary:
            for key in [
                "sterimol_bond",
                "sterimol_dummy_raw",
                "pyramidalization",
                "buried_volume",
            ]:
                for field, value in new.get(key, {}).items():
                    if isinstance(value, (int, float)) and isinstance(
                        old.get(key, {}).get(field), (int, float)
                    ):
                        changes.setdefault(f"{key}.{field}", []).append(
                            value - old[key][field]
                        )
        output["cases"].append(result)
    output["ordinary_change_distribution"] = {
        key: {
            "count": len(values),
            "deltas": values,
            "max_abs": max(map(abs, values)),
            "median_abs": float(np.median(np.abs(values))),
        }
        for key, values in changes.items()
    }
    with args.output.open("x") as handle:
        json.dump(output, handle, indent=2, allow_nan=False)
        handle.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["requests", "before", "after", "output"]:
        parser.add_argument("--" + name, type=Path, required=True)
    analyze(parser.parse_args())
