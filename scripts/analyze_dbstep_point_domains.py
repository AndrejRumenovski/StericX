#!/usr/bin/env python3
"""Observe Rust decisions; compare independent predicates on identical nodes.

This research controller does not change the historical gate or implement a
production descriptor. New release admission is deliberately not scored here.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import subprocess
import time
from contextlib import redirect_stdout
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from dbstep import sterics
from dbstep_resolution_reference import (
    f32_predicate_certificate,
    union_mask_exact,
    union_mask_f64,
)
from validate_dbstep_volume import load_cases, request

ROOT = Path(__file__).resolve().parents[1]


def record(path):
    path = Path(path).resolve()
    with path.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest}


def save(path, value):
    Path(path).write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )


class Observer:
    def __init__(self, executable, directory):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True)
        self.requests = (self.directory / "requests.jsonl").open("w")
        self.responses = (self.directory / "responses.jsonl").open("w")
        self.stderr = (self.directory / "stderr").open("w")
        self.process = subprocess.Popen(
            [str(executable)],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=self.stderr,
            text=True,
            bufsize=1,
        )

    def call(self, payload):
        line = json.dumps(payload, allow_nan=False) + "\n"
        self.requests.write(line)
        self.requests.flush()
        self.process.stdin.write(line)
        self.process.stdin.flush()
        response = self.process.stdout.readline()
        self.responses.write(response)
        self.responses.flush()
        if not response:
            raise RuntimeError(f"Observer stopped: {self.process.poll()}")
        parsed = json.loads(response)
        if "error" in parsed:
            raise RuntimeError(parsed)
        return parsed

    def close(self):
        self.process.stdin.close()
        code = self.process.wait()
        self.requests.close()
        self.responses.close()
        self.stderr.close()
        if code:
            raise RuntimeError(f"Observer exit {code}")


def trace_f32(point, center, radius):
    """Diagnostic IEEE graph; never the independent scientific oracle."""
    f = np.float32
    differences = f(f(point) - f(center))
    squares = f(differences * differences)
    partial = f(squares[0] + squares[1])
    total = f(partial + squares[2])
    r2 = f(f(radius) * f(radius))
    return {
        "differences": differences.tolist(),
        "squares": squares.tolist(),
        "xy_sum": float(partial),
        "squared_distance": float(total),
        "radius_squared": float(r2),
        "inside": bool(total <= r2),
        "squared_distance_bits": int(total.view(np.uint32)),
        "radius_squared_bits": int(r2.view(np.uint32)),
    }


def dbstep_mask(points, centers, radii):
    options = SimpleNamespace(grid=1.0, verbose=False, debug=False)
    with redirect_stdout(io.StringIO()):
        occupied, tree, _ = sterics.occupied(
            points, centers, radii, np.zeros(3), options
        )
    # DBSTEP returns unique rows of the supplied array. Map every row exactly;
    # counts alone could hide compensated errors. Duplicate input points rejected.
    mapping = {tuple(point): i for i, point in enumerate(points)}
    if len(mapping) != len(points):
        raise ValueError("Duplicate input sample nodes")
    mask = np.zeros(len(points), dtype=bool)
    for point in occupied:
        mask[mapping[tuple(point)]] = True
    return mask, tree


def analyze_mask(points, centers, radii, observed):
    native = np.asarray(observed["occupied_mask"], dtype=bool)
    reference = union_mask_f64(points, centers, radii)
    dbstep, tree = dbstep_mask(points, centers, radii)
    disputed = np.flatnonzero((native != reference) | (native != dbstep))
    witnesses, unexplained = [], []
    for index in disputed:
        point = points[index]
        exact = bool(union_mask_exact([point], centers, radii)[0])
        atoms = []
        for center, radius in zip(centers, radii, strict=True):
            certificate = f32_predicate_certificate(point, center, radius)
            atoms.append(
                {
                    "center": center.tolist(),
                    "radius": float(radius),
                    "certificate": certificate,
                    "diagnostic_ieee_trace": trace_f32(point, center, radius),
                }
            )
        stable_inside = any(
            a["certificate"]["classification"] == "stable_inside" for a in atoms
        )
        stable_outside = all(
            a["certificate"]["classification"] == "stable_outside" for a in atoms
        )
        traced = any(a["diagnostic_ieee_trace"]["inside"] for a in atoms)
        stable_disagreement = (stable_inside and not native[index]) or (
            stable_outside and native[index]
        )
        explained = (
            not stable_disagreement
            and traced == native[index]
            and dbstep[index] == exact
        )
        witness = {
            "index": int(index),
            "point": point.tolist(),
            "native": bool(native[index]),
            "direct_f64": bool(reference[index]),
            "dbstep": bool(dbstep[index]),
            "exact_rational": exact,
            "atoms": atoms,
            "stable_disagreement": bool(stable_disagreement),
            "classified_binary32_roundoff": bool(explained),
        }
        witnesses.append(witness)
        if not explained:
            unexplained.append(witness)
    return (
        {
            "native_count": int(native.sum()),
            "dbstep_count": int(dbstep.sum()),
            "direct_f64_count": int(reference.sum()),
            "point_count": len(points),
            "native_dbstep_mask_differences": int(np.count_nonzero(native != dbstep)),
            "native_direct_mask_differences": int(
                np.count_nonzero(native != reference)
            ),
            "disputed_points": witnesses,
            "unexplained_count": len(unexplained),
        },
        tree,
        dbstep,
    )


def sphere_witnesses(points, radius, tree):
    tree_inside = np.zeros(len(points), dtype=bool)
    tree_inside[tree.query_ball_point(np.zeros(3), radius)] = True
    # Native grid admitted all supplied nodes. Inspect every removed point and
    # every point near the independent a priori roundoff band.
    squared = np.einsum("ij,ij->i", points, points)
    u = 2.0**-24
    bound = 5 * u / (1 - 5 * u) * squared + u * radius**2
    near = np.flatnonzero((abs(squared - radius**2) <= bound * 1.000001) | ~tree_inside)
    witnesses = []
    for index in near:
        point = points[index]
        certificate = f32_predicate_certificate(point, [0, 0, 0], radius)
        trace = trace_f32(point, [0, 0, 0], radius)
        witnesses.append(
            {
                "index": int(index),
                "point": point.tolist(),
                "native_grid_included": True,
                "dbstep_tree_included": bool(tree_inside[index]),
                "certificate": certificate,
                "diagnostic_ieee_trace": trace,
                "explained": trace["inside"]
                and certificate["classification"] != "stable_outside",
            }
        )
    return tree_inside, witnesses


def verify_counts(observed):
    mask = np.asarray(observed["occupied_mask"], dtype=bool)
    bins = np.asarray(observed["point_bins"])
    f = np.float32
    volume = f(observed["sphere_volume"])
    reconstructed = int(f(f(mask.sum()) / f(len(mask)) * volume).view(np.uint32))
    errors = []
    if reconstructed != observed["aggregate"]["_bits"]["buried_volume"]:
        errors.append("aggregate count arithmetic")
    for column, name, count in [(0, "quadrant", 4), (1, "octant", 8)]:
        expected_totals = np.bincount(bins[:, column], minlength=count)
        expected_occupied = np.bincount(bins[mask, column], minlength=count)
        if not np.array_equal(expected_totals, observed[f"{name}_total"]):
            errors.append(f"{name} totals")
        if not np.array_equal(expected_occupied, observed[f"{name}_occupied"]):
            errors.append(f"{name} occupation")
        for index in range(count):
            expected = (
                f(0)
                if expected_totals[index] == 0
                else f(
                    f(expected_occupied[index])
                    / f(expected_totals[index])
                    * f(volume / f(count))
                )
            )
            if (
                int(expected.view(np.uint32))
                != observed["aggregate"]["_bits"][name + "s"][index]
            ):
                errors.append(f"{name} {index} count arithmetic")
    return errors


def run_case(case, observer, folder):
    folder.mkdir(parents=True)
    captured = observer.call(request(case, True))
    config = captured["config"] if "config" in captured else case["options"]["volume"]
    points = np.asarray(captured["dump"]["points"], dtype=np.float32).astype(float)
    radius = float(np.float32(config["sphere_radius"]))
    retained = [
        a
        for a in captured["atoms"]
        if config["include_hydrogens"] or a["element"].upper() != "H"
    ]
    radii = np.asarray(
        [
            np.float32(a["vdw_radius"]) * np.float32(config["radii_scale"])
            for a in retained
        ],
        dtype=np.float32,
    ).astype(float)
    results = []
    for frame in captured["dump"]["orientations"]:
        centers = np.asarray(
            [a["position"] for a in frame["aligned_atoms"]], dtype=np.float32
        ).astype(float)
        expected_r2 = np.float32(np.float32(radii) * np.float32(radii))
        actual_r2 = np.float32([a["radius_squared"] for a in frame["aligned_atoms"]])
        if not np.array_equal(expected_r2.view(np.uint32), actual_r2.view(np.uint32)):
            raise ValueError(
                "Captured radius square differs from accepted scaled radius"
            )
        observed = observer.call(
            {
                "id": f"{case['id']}/mask/{frame['plane']}",
                "op": "raw_point_decisions",
                "config": config,
                "points": points.tolist(),
                "atoms": frame["aligned_atoms"],
                "include_points": False,
            }
        )
        analysis, tree, mask = analyze_mask(points, centers, radii, observed)
        sphere, sphere_records = sphere_witnesses(points, radius, tree)
        volume = 4 * math.pi * radius**3 / 3
        native_volume = float(np.float32(frame["buried_volume"]))
        matched_volume = float(mask.sum()) / len(points) * volume
        reclipped_volume = (
            float(np.count_nonzero(mask & sphere)) / int(sphere.sum()) * volume
        )
        legacy_limit = 32 * 2.0**-23 * volume
        arithmetic_errors = verify_counts(observed)
        if observed["aggregate"]["_bits"]["buried_volume"] != int(
            np.float32(native_volume).view(np.uint32)
        ):
            arithmetic_errors.append("per-point adapter vs public frame")
        result = {
            "case_id": case["id"],
            "plane": frame["plane"],
            "config": config,
            **analysis,
            "sphere_witnesses": sphere_records,
            "native_volume": native_volume,
            "same_domain_dbstep_volume": matched_volume,
            "legacy_reclipped_dbstep_volume": reclipped_volume,
            "legacy_sphere_count": int(sphere.sum()),
            "legacy_occupied_count": int(np.count_nonzero(mask & sphere)),
            "same_domain_residual": native_volume - matched_volume,
            "legacy_residual": native_volume - reclipped_volume,
            "unchanged_arithmetic_limit": legacy_limit,
            "same_domain_within_old_arithmetic_limit": abs(
                native_volume - matched_volume
            )
            <= legacy_limit,
            "legacy_gate_pass": len(points) == int(sphere.sum())
            and abs(native_volume - reclipped_volume) <= legacy_limit,
            "arithmetic_errors": arithmetic_errors,
            "unexplained_sphere_count": sum(
                not row["explained"] for row in sphere_records
            ),
            "accepted_radii": radii.tolist(),
            "accepted_centers": centers.tolist(),
        }
        save(folder / f"plane_{frame['plane']}.json", result)
        results.append(result)
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--observer", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    cases = load_cases()
    save(args.output / "inputs.json", cases)
    sources = [
        Path(__file__),
        ROOT / "scripts/dbstep_resolution_reference.py",
        ROOT / "scripts/validate_dbstep_volume.py",
        ROOT / "docs/validation/dbstep/resolution_analysis/ANALYSIS_PROTOCOL.md",
        args.observer,
    ]
    before = [record(path) for path in sources]
    save(
        args.output / "pre_execution.json",
        {
            "sources": before,
            "case_count": len(cases),
            "release_admission_scored": False,
        },
    )
    observer = Observer(args.observer.resolve(), args.output / "native")
    rows, failures = [], []
    start = time.perf_counter()
    try:
        for index, case in enumerate(cases):
            try:
                rows.extend(
                    run_case(case, observer, args.output / "cases" / f"{index:03d}")
                )
            except Exception as error:
                failures.append({"case_id": case["id"], "exception": repr(error)})
            print(
                json.dumps(
                    {
                        "completed": index + 1,
                        "total": len(cases),
                        "execution_failures": len(failures),
                    }
                ),
                flush=True,
            )
    finally:
        observer.close()
    summary = {
        "case_count": len(cases),
        "frame_count": len(rows),
        "execution_failures": failures,
        "legacy_failed_frames": [
            {
                "case_id": r["case_id"],
                "plane": r["plane"],
                "residual": r["legacy_residual"],
            }
            for r in rows
            if not r["legacy_gate_pass"]
        ],
        "same_domain_arithmetic_failed_frames": sum(
            not r["same_domain_within_old_arithmetic_limit"] for r in rows
        ),
        "disputed_atomic_points": sum(len(r["disputed_points"]) for r in rows),
        "unexplained_atomic_points": sum(r["unexplained_count"] for r in rows),
        "unexplained_sphere_points": sum(r["unexplained_sphere_count"] for r in rows),
        "arithmetic_error_count": sum(len(r["arithmetic_errors"]) for r in rows),
        "elapsed_seconds": time.perf_counter() - start,
        "release_admission_scored": False,
        "sources_unchanged": before == [record(path) for path in sources],
    }
    save(args.output / "summary.json", summary)
    save(
        args.output / "manifest.json",
        [record(path) for path in sorted(args.output.rglob("*")) if path.is_file()],
    )
    print(json.dumps(summary, indent=2))
    return int(
        bool(failures)
        or not summary["sources_unchanged"]
        or bool(
            summary["unexplained_atomic_points"]
            + summary["unexplained_sphere_points"]
            + summary["arithmetic_error_count"]
        )
    )


if __name__ == "__main__":
    raise SystemExit(main())
