#!/usr/bin/env python3
"""Stream fresh independent three-frame BV references; never calls SUT kernels.

Retains the sealed float64 grid, radius table, strict-open octant convention and
union-of-balls reference. Inputs use identical f32 coordinates and the captured
SUT center, isolating projection/occupancy from center-model differences. Output
is new reference evidence, not an acceptance verdict or a new residual limit.
"""

import argparse
import functools
import hashlib
import importlib.metadata
import itertools
import json
import os
import platform
import sys
from collections import deque
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from sealed_reference import (
    AUDIT,
    MANIFEST,
    MANIFEST_SHA256,
    REFERENCE,
    ROOT,
    load_reference,
    refuse_audit_output,
)

reference = load_reference()


def sha(path):
    result = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


@functools.lru_cache(maxsize=8)
def points(radius, density):
    # Cache only independent immutable grids, not any occupancy observations.
    return reference.grid(radius, density)


def calculate(pair):
    request, observed = pair
    row = {
        "id": request["id"],
        "coordinate_contract": (
            "f32 input coordinates converted to float64; controlled captured center"
        ),
    }
    try:
        xyz = np.array(
            [atom["position"] for atom in request["atoms"]], dtype=np.float32
        ).astype(np.float64)
        center = np.array(observed["center"], dtype=np.float64)
        if (
            center.shape != (3,)
            or not np.all(np.isfinite(xyz))
            or not np.all(np.isfinite(center))
        ):
            raise ValueError("nonfinite coordinates or unavailable center")
        cfg = request["config"]
        radii = reference.radii_for(request) * cfg["radii_scale"]
        include = np.array(
            [
                cfg["include_hydrogens"] or atom["element"].upper() != "H"
                for atom in request["atoms"]
            ]
        )
        cloud = points(cfg["sphere_radius"], cfg["density"])
        orientations = []
        z = reference.unit(center - xyz[request["donor"]])
        for neighbor in request["neighbors"]:
            plane = xyz[neighbor] - center
            x = reference.unit(plane - np.dot(plane, z) * z)
            # Exactly the sealed bv_refs convention: z and x are normalized;
            # y is their cross product, without an additional normalization.
            y = np.cross(z, x)
            aligned = (xyz - center) @ np.array([x, y, z]).T
            volume = reference.volume_on_grid(
                cloud, aligned[include], radii[include], cfg["sphere_radius"]
            )
            orientations.append({"plane": neighbor, "independent": volume})
        if len(orientations) != 3:
            raise ValueError("reference requires exactly three explicit neighbors")
        keys = ["buried_volume", "percent_buried_volume", "near_vbur", "far_vbur"]
        mean = {
            key: sum(frame["independent"][key] for frame in orientations) / 3
            for key in keys
        }
        q = np.array([frame["independent"]["quadrants"] for frame in orientations])
        octants = np.array([frame["independent"]["octants"] for frame in orientations])
        summary = mean | {
            "qvbur_min": float(q.min()),
            "qvbur_max": float(q.max()),
            "ovbur_min": float(octants.min()),
            "ovbur_max": float(octants.max()),
            "max_delta_qvbur": float(np.abs(q - np.roll(q, 1, axis=1)).max()),
        }
        row["buried_volume"] = {
            "center": center.tolist(),
            "orientations": orientations,
            "independent": summary,
            "original_first_frame_scalars": {
                key: orientations[0]["independent"][key] for key in keys
            },
        }
        public = observed.get("buried_volume", {})
        row["sut_error"] = public.get("error")
        row["residual"] = (
            {key: public[key] - value for key, value in summary.items()}
            if "error" not in public
            else None
        )
    except (ValueError, KeyError, TypeError, FloatingPointError) as error:
        row["reference_error"] = str(error)
    return row


def requests(args, coverage):
    with args.requests.open() as req, args.observations.open() as obs:
        count = 0
        seen = set()
        for left, right in itertools.zip_longest(req, obs):
            if left is None or right is None:
                raise ValueError("request/observation lengths differ")
            request, observed = json.loads(left), json.loads(right)
            if request["id"] != observed["id"]:
                raise ValueError("request/observation ID order differs")
            if request["id"] in seen:
                raise ValueError("duplicate source request ID")
            seen.add(request["id"])
            count += 1
            if args.limit is None or count <= args.limit:
                yield request, observed
        if count != args.expected_rows:
            raise ValueError(
                f"Expected {args.expected_rows} source rows, observed {count}"
            )
        coverage["source_rows"] = count
        coverage["complete_source_stream_validated"] = True


def main(args):
    refuse_audit_output(args.output)
    args.output.mkdir(parents=True, exist_ok=False)
    identities = {
        str(path): sha(path)
        for path in [
            args.requests,
            args.observations,
            REFERENCE,
            MANIFEST,
            Path(__file__),
            Path(__file__).with_name("sealed_reference.py"),
            ROOT / "uv.lock",
            ROOT / "pyproject.toml",
        ]
    }
    manifest = {
        "kind": "new_independent_three_frame_volume",
        "global_scientific_pass": False,
        "inputs": identities,
        "argv": sys.argv,
        "workers": args.workers,
        "limit": args.limit,
        "expected_rows": args.expected_rows,
        "audit_seal_sha256": MANIFEST_SHA256,
        "audit": str(AUDIT),
        "numpy": np.__version__,
        "python": sys.version,
        "platform": platform.platform(),
        "packages": {
            name: importlib.metadata.version(name)
            for name in ["numpy", "scipy", "morfeus-ml"]
        },
        "environment": {
            name: os.environ.get(name)
            for name in [
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS",
                "NUMEXPR_NUM_THREADS",
                "PYTHONHASHSEED",
            ]
        },
        "independent_reference_scope": "full supplied request stream"
        if not args.limit
        else "explicitly limited prefix",
    }
    (args.output / "started.json").write_text(json.dumps(manifest, indent=2) + "\n")
    count = errors = compared = 0
    seen = set()
    coverage = {}
    with (
        (args.output / "reference_results.jsonl").open("x") as output,
        ProcessPoolExecutor(max_workers=args.workers) as executor,
    ):
        pending = deque()
        for pair in requests(args, coverage):
            pending.append(executor.submit(calculate, pair))
            if len(pending) >= 2 * args.workers:
                row = pending.popleft().result()
                if row["id"] in seen:
                    raise ValueError("duplicate request ID")
                seen.add(row["id"])
                output.write(json.dumps(row, allow_nan=False) + "\n")
                count += 1
                errors += "reference_error" in row
                compared += row.get("residual") is not None
        while pending:
            row = pending.popleft().result()
            if row["id"] in seen:
                raise ValueError("duplicate request ID")
            seen.add(row["id"])
            output.write(json.dumps(row, allow_nan=False) + "\n")
            count += 1
            errors += "reference_error" in row
            compared += row.get("residual") is not None
    if identities != {path: sha(Path(path)) for path in identities}:
        raise ValueError("input/reference bytes changed during replay")
    expected_computations = (
        min(args.limit, args.expected_rows) if args.limit else args.expected_rows
    )
    if count != expected_computations or not coverage.get(
        "complete_source_stream_validated"
    ):
        raise ValueError("Incomplete selected reference scope")
    manifest.update(
        {
            "complete": True,
            "requests": count,
            "reference_errors": errors,
            "compared": compared,
            "reference_results_sha256": sha(args.output / "reference_results.jsonl"),
            **coverage,
        }
    )
    (args.output / "complete.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["requests", "observations", "output"]:
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--expected-rows", type=int, required=True)
    args = parser.parse_args()
    if (
        args.workers < 1
        or args.expected_rows < 1
        or (args.limit is not None and not 1 <= args.limit <= args.expected_rows)
    ):
        parser.error("workers and optional limit must be positive")
    main(args)
