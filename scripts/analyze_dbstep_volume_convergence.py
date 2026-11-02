#!/usr/bin/env python3
"""Fixed-resolution independent volume study; observes SUTs, never scores release."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import itertools
import json
import math
import os
import shutil
import subprocess
import sys
import time
from contextlib import redirect_stderr, redirect_stdout
from fractions import Fraction
from importlib.metadata import version
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from dbstep import sterics
from dbstep_resolution_reference import (
    certified_union_volume,
    f32_predicate_certificate,
    pi_interval,
    sphere_intersection_interval,
    union_mask_exact,
    union_mask_f64,
)
from morfeus import BuriedVolume

ROOT = Path(__file__).resolve().parents[1]
DENSITIES = (0.08, 0.01, 0.00125, 0.00015625)
SPACINGS = (0.4, 0.2, 0.1, 0.05)
PIN = "4365353e0a2f24d05336e619759e350bc630a3cb"


def fraction(value):
    return value if isinstance(value, Fraction) else Fraction(float(value))


def upward(value):
    value = fraction(value)
    result = float(value)
    return math.nextafter(result, math.inf) if fraction(result) < value else result


def downward(value):
    return -upward(-fraction(value))


def record(path):
    path = Path(path).resolve()
    with path.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest}


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def float_bits(value):
    return int(np.asarray(value, dtype=np.float32).view(np.uint32))


def tensor_points(axes):
    return np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)


def point_indices(points, axes):
    points = np.asarray(points, dtype=float).reshape((-1, 3))
    indices = []
    for dim, axis in enumerate(axes):
        index = np.searchsorted(axis, points[:, dim])
        if np.any(index == len(axis)) or not np.array_equal(
            axis[index], points[:, dim]
        ):
            raise ValueError("Captured point does not belong to declared actual axis")
        indices.append(index)
    return np.ravel_multi_index(tuple(indices), tuple(map(len, axes)))


def cell_geometry(axes, radius):
    """Exact dyadic Voronoi width extrema and an outward cell-radius bound."""
    minima, maxima, steps, half_spans = [], [], [], []
    coverage = []
    for axis in axes:
        axis = [fraction(x) for x in axis]
        if len(axis) < 2 or any(a >= b for a, b in itertools.pairwise(axis)):
            raise ValueError("Strictly increasing tensor axes are required")
        boundaries = [axis[0] - (axis[1] - axis[0]) / 2]
        boundaries += [(a + b) / 2 for a, b in itertools.pairwise(axis)]
        boundaries.append(axis[-1] + (axis[-1] - axis[-2]) / 2)
        if boundaries[0] > -fraction(radius) or boundaries[-1] < fraction(radius):
            raise ValueError("Finite tensor cells do not cover the integration sphere")
        widths = [b - a for a, b in itertools.pairwise(boundaries)]
        minima.append(min(widths))
        maxima.append(max(widths))
        steps.append((axis[-1] - axis[0]) / (len(axis) - 1))
        half_spans.append(
            max(
                max(x - a, b - x)
                for x, a, b in zip(axis, boundaries[:-1], boundaries[1:], strict=True)
            )
        )
        coverage.append([str(boundaries[0]), str(boundaries[-1])])
    radius_squared = sum(x * x for x in half_spans)
    delta = math.sqrt(upward(radius_squared))
    if fraction(delta) ** 2 < radius_squared:
        delta = math.nextafter(delta, math.inf)
    minimum, maximum, nominal = math.prod(minima), math.prod(maxima), math.prod(steps)
    return {
        "cell_radius_upper": fraction(delta),
        "minimum_cell_volume": minimum,
        "maximum_cell_volume": maximum,
        "count_cell_volume": nominal,
        "weight_error_per_selected_node": max(
            abs(minimum - nominal), abs(maximum - nominal)
        ),
        "coverage_exact": coverage,
    }


def shell_upper(radius, cell_radius):
    """Shell bound includes a priori binary64 distance/radius roundoff."""
    r = fraction(radius)
    u = Fraction(1, 2**53)
    gamma = 5 * u / (1 - 5 * u)
    # The direct reference evaluates three subtracts/squares and two additions.
    # This algebra proves its radial uncertainty is at most 8*u*r.
    assert (1 + u) / (1 - gamma) <= (1 + 8 * u) ** 2
    assert (1 - u) / (1 + gamma) >= (1 - 8 * u) ** 2
    delta = fraction(cell_radius) + 8 * u * r
    coefficient = Fraction(4, 3) * ((r + delta) ** 3 - max(Fraction(0), r - delta) ** 3)
    return coefficient * pi_interval()[1]


def estimator_enclosure(axes, centers, radii, radius, denominator, numerator, reported):
    """Bound an actual finite-count estimator using actual tensor-cell geometry.

    Actual and reference numerator/denominator masks may differ. Every such
    decision is charged one maximum cell volume. Reference predicates have an
    a priori roundoff shell, so they are not treated as exact classifications.
    """
    points = tensor_points(axes)
    denominator = np.asarray(denominator, dtype=bool)
    numerator = np.asarray(numerator, dtype=bool)
    if denominator.shape != (len(points),) or numerator.shape != denominator.shape:
        raise ValueError("One explicit decision per tensor node is required")
    if (
        not np.isfinite(points).all()
        or max(
            abs(points).max(initial=0),
            np.max(np.abs(centers), initial=0),
            max(radii, default=0),
            radius,
        )
        > 1e100
    ):
        raise ValueError("Outside the finite molecular-scale certificate domain")
    # All campaign coordinates/radii and nonzero differences are normal f64.
    smallest = np.finfo(float).tiny
    for center, r in zip([np.zeros(3), *centers], [radius, *radii], strict=True):
        offsets = np.abs(points - center)
        if np.any((offsets != 0) & (offsets < math.sqrt(smallest) * 2)) or (
            r and r < math.sqrt(smallest) * 2
        ):
            raise ValueError("Subnormal reference arithmetic is not certified")
    reference_sphere = union_mask_f64(points, [[0, 0, 0]], [radius])
    reference_occupied = reference_sphere & union_mask_f64(points, centers, radii)
    ds = int(np.count_nonzero(reference_sphere != denominator))
    do = int(np.count_nonzero(reference_occupied != numerator))
    ns, no = int(denominator.sum()), int(numerator.sum())
    if ns == 0:
        raise ValueError("A positive denominator is required")
    cells = cell_geometry(axes, radius)
    sphere_shell = shell_upper(radius, cells["cell_radius_upper"])
    atom_shells = sum(shell_upper(r, cells["cell_radius_upper"]) for r in radii if r)
    e_sphere = (
        sphere_shell
        + ds * cells["maximum_cell_volume"]
        + ns * cells["weight_error_per_selected_node"]
    )
    e_occupied = (
        (sphere_shell + atom_shells if len(radii) else Fraction(0))
        + do * cells["maximum_cell_volume"]
        + no * cells["weight_error_per_selected_node"]
    )
    ratio_bound = e_occupied + Fraction(no, ns) * e_sphere
    pi_low, pi_high = pi_interval()
    ratio_coefficient = Fraction(4, 3) * fraction(radius) ** 3 * Fraction(no, ns)
    ratio_low, ratio_high = ratio_coefficient * pi_low, ratio_coefficient * pi_high
    arithmetic_or_clamp = max(
        abs(fraction(reported) - ratio_low), abs(fraction(reported) - ratio_high)
    )
    error_bound = ratio_bound + arithmetic_or_clamp
    result = {
        "method": "actual tensor-cell shells with measured decision differences",
        "sphere_nodes": ns,
        "occupied_nodes": no,
        "reference_sphere_nodes": int(reference_sphere.sum()),
        "reference_occupied_nodes": int(reference_occupied.sum()),
        "denominator_reference_mask_disagreements": ds,
        "numerator_reference_mask_disagreements": do,
        "cell_geometry": {
            key: str(value) if isinstance(value, Fraction) else value
            for key, value in cells.items()
        },
        "sphere_quadrature_error_upper": upward(e_sphere),
        "occupied_quadrature_error_upper": upward(e_occupied),
        "ratio_discretization_error_upper": upward(ratio_bound),
        "normalization_arithmetic_or_clamp_error_upper": upward(arithmetic_or_clamp),
        "total_absolute_error_upper": upward(error_bound),
        "certified_volume_lower": downward(fraction(reported) - error_bound),
        "certified_volume_upper": upward(fraction(reported) + error_bound),
        "reference_roundoff_radial_enlargement": "8 * 2^-53 * radius",
        "release_gate_scored": False,
    }
    differences = {
        "denominator_disagreement_indices": np.flatnonzero(
            reference_sphere != denominator
        ),
        "numerator_disagreement_indices": np.flatnonzero(
            reference_occupied != numerator
        ),
    }
    return result, differences


def native_ball_trace(point, center, radius):
    """Explicit IEEE arithmetic witness, not a replacement geometry kernel."""
    certificate = f32_predicate_certificate(point, center, radius)
    p, c, r = (
        np.asarray(point, dtype=np.float32),
        np.asarray(center, dtype=np.float32),
        np.float32(radius),
    )
    dx, dy, dz = np.asarray(p - c, dtype=np.float32)
    x2, y2, z2 = np.float32(dx * dx), np.float32(dy * dy), np.float32(dz * dz)
    xy = np.float32(x2 + y2)
    distance_squared = np.float32(xy + z2)
    radius_squared = np.float32(r * r)
    inside = bool(distance_squared <= radius_squared)
    trace = {
        name: {"value": float(value), "bits": float_bits(value)}
        for name, value in (
            ("dx", dx),
            ("dy", dy),
            ("dz", dz),
            ("x2", x2),
            ("y2", y2),
            ("z2", z2),
            ("xy", xy),
            ("distance_squared", distance_squared),
            ("radius_squared", radius_squared),
        )
    }
    return {
        "center": list(center),
        "radius": float(radius),
        "certificate": certificate,
        "trace": trace,
        "trace_inside": inside,
        "trace_exact_dispute_within_bound": inside == certificate["exact_inside"]
        or certificate["classification"] == "ambiguous",
    }


def native_dispute_witnesses(case, axes, denominator, numerator, differences):
    indices = np.union1d(
        differences["denominator_disagreement_indices"],
        differences["numerator_disagreement_indices"],
    )
    rows = []
    for index in indices:
        xyz = np.unravel_index(index, tuple(map(len, axes)))
        point = [float(axis[i]) for axis, i in zip(axes, xyz, strict=True)]
        sphere = native_ball_trace(point, [0, 0, 0], case["radius"])
        atoms = [
            native_ball_trace(point, c, r)
            for c, r in zip(case["centers"], case["radii"], strict=True)
        ]
        trace_sphere = sphere["trace_inside"]
        trace_occupied = trace_sphere and any(atom["trace_inside"] for atom in atoms)
        exact_sphere = sphere["certificate"]["exact_inside"]
        exact_occupied = exact_sphere and any(
            atom["certificate"]["exact_inside"] for atom in atoms
        )
        matches = trace_sphere == bool(denominator[index]) and trace_occupied == bool(
            numerator[index]
        )
        explained = (
            matches
            and sphere["trace_exact_dispute_within_bound"]
            and all(atom["trace_exact_dispute_within_bound"] for atom in atoms)
        )
        rows.append(
            {
                "tensor_index": int(index),
                "point": point,
                "observed_sphere": bool(denominator[index]),
                "observed_numerator": bool(numerator[index]),
                "exact_sphere": exact_sphere,
                "exact_numerator": exact_occupied,
                "sphere_witness": sphere,
                "atomic_witnesses": atoms,
                "traces_reproduce_observed_decisions": matches,
                "explained_within_a_priori_bound": explained,
            }
        )
    return {
        "rows": rows,
        "disputed_node_count": len(rows),
        "unexplained_count": sum(
            not row["explained_within_a_priori_bound"] for row in rows
        ),
    }


def primitive_cases():
    def case(name, centers, radii, truth=None):
        return {
            "id": name,
            "kind": "mathematical_primitive",
            "radius": 2.5,
            "centers": np.asarray(centers, dtype=np.float32)
            .reshape((-1, 3))
            .astype(float)
            .tolist(),
            "radii": np.asarray(radii, dtype=np.float32).astype(float).tolist(),
            "truth": truth,
        }

    cases = [case("empty", [], [], {"kind": "zero"})]
    singles = [
        ("concentric", 0, 1),
        ("enclosed", 0.5, 1),
        ("crossing", 2.2, 1.1),
        ("disjoint", 4, 0.5),
        ("covering", 0.25, 4),
        ("equal_spheres", 0, 2.5),
        ("internal_tangent", 1.5, 1),
        ("external_tangent", 3.5, 1),
    ]
    for label, d in (("internal", 1.5), ("external", 3.5)):
        for direction, target in (("below", -np.inf), ("above", np.inf)):
            singles.append(
                (
                    f"{label}_tangent_{direction}",
                    float(np.nextafter(np.float32(d), np.float32(target))),
                    1,
                )
            )
    for name, distance, r in singles:
        cases.append(case(name, [[distance, 0, 0]], [r], {"kind": "single"}))
    cases += [
        case(
            "separated",
            [[-0.75, 0, 0], [0.75, 0, 0]],
            [0.5, 0.5],
            {"kind": "two_inside"},
        ),
        case(
            "overlapping", [[-0.5, 0, 0], [0.5, 0, 0]], [1, 1], {"kind": "two_inside"}
        ),
        case("coincident", [[0, 0, 0], [0, 0, 0]], [1, 1], {"kind": "two_inside"}),
        case(
            "asymmetric",
            [[0.4, -0.3, 0.2], [1.9, 0.5, -0.7], [-1.3, 0.4, 0.9]],
            [1.3, 0.8, 1.1],
        ),
    ]
    return cases


def analytical_interval(case):
    truth = case["truth"]
    if not truth:
        return None
    if truth["kind"] == "zero":
        return [0.0, 0.0]
    centers, radii, radius = case["centers"], case["radii"], case["radius"]
    if truth["kind"] == "single":
        assert centers[0][1:] == [0, 0]
        return list(sphere_intersection_interval(radius, radii[0], abs(centers[0][0])))
    assert truth["kind"] == "two_inside"
    assert all(
        c[1:] == [0, 0] and abs(c[0]) + r <= radius
        for c, r in zip(centers, radii, strict=True)
    )
    first = sphere_intersection_interval(radius, radii[0], abs(centers[0][0]))
    second = sphere_intersection_interval(radius, radii[1], abs(centers[1][0]))
    lens = sphere_intersection_interval(
        radii[0], radii[1], abs(centers[0][0] - centers[1][0])
    )
    return [
        downward(fraction(first[0]) + fraction(second[0]) - fraction(lens[1])),
        upward(fraction(first[1]) + fraction(second[1]) - fraction(lens[0])),
    ]


def molecular_cases():
    source = ROOT / "docs/validation/dbstep/volume_v2"
    identity = "SIG-NIHDA-401/0/sphere_radius=2.5"
    request = next(
        json.loads(line)
        for line in (source / "requests.jsonl").read_text().splitlines()
        if json.loads(line)["id"] == identity
    )
    result = next(
        json.loads(line)
        for line in (source / "native_stdout.jsonl").read_text().splitlines()
        if json.loads(line)["id"] == identity
    )
    retained = [
        a
        for a in result["atoms"]
        if request["config"]["include_hydrogens"] or a["element"].upper() != "H"
    ]
    radii = [
        float(
            np.float32(a["vdw_radius"]) * np.float32(request["config"]["radii_scale"])
        )
        for a in retained
    ]
    cases = []
    for frame in result["dump"]["orientations"]:
        assert [float(np.float32(r) * np.float32(r)) for r in radii] == [
            a["radius_squared"] for a in frame["aligned_atoms"]
        ]
        cases.append(
            {
                "id": f"original401_plane{frame['plane']}",
                "kind": "captured_molecular_frame",
                "source_case_id": identity,
                "source_records": [
                    record(source / p)
                    for p in ("requests.jsonl", "native_stdout.jsonl")
                ],
                "radius": 2.5,
                "centers": [a["position"] for a in frame["aligned_atoms"]],
                "radii": radii,
                "truth": None,
                "source_frame": frame,
                "retained_source_atom_indices": [
                    i
                    for i, a in enumerate(result["atoms"])
                    if request["config"]["include_hydrogens"]
                    or a["element"].upper() != "H"
                ],
            }
        )
    return cases


class Observer:
    def __init__(self, binary, output):
        self.stderr = (output / "observer_stderr.txt").open("w")
        self.process = subprocess.Popen(
            [str(binary.resolve())],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=self.stderr,
            text=True,
        )

    def call(self, request, destination):
        payload = json.dumps(request, separators=(",", ":")) + "\n"
        with gzip.open(
            destination.with_suffix(".request.json.gz"), "wt", compresslevel=1
        ) as stream:
            stream.write(payload)
        start = time.perf_counter()
        self.process.stdin.write(payload)
        self.process.stdin.flush()
        line = self.process.stdout.readline()
        elapsed = time.perf_counter() - start
        with gzip.open(
            destination.with_suffix(".response.json.gz"), "wt", compresslevel=1
        ) as stream:
            stream.write(line)
        result = json.loads(line)
        if "error" in result:
            raise RuntimeError(result)
        return result, elapsed

    def close(self):
        self.process.stdin.close()
        code = self.process.wait()
        self.stderr.close()
        if code:
            raise RuntimeError(f"Observation process exited {code}")


def native_axes(points, radius):
    """Recover grid metadata, never the scientific occupancy predicate.

    The true side count is inferred by exact membership of all observed axes;
    source operation order is used only to supply unselected endpoint metadata.
    """
    observed = np.unique(points[:, 0])
    matches = []
    for count in range(max(2, len(observed)), len(observed) + 4):
        step = np.float32(np.float32(2) * np.float32(radius) / np.float32(count - 1))
        axis = (-np.float32(radius) + np.arange(count, dtype=np.float32) * step).astype(
            float
        )
        if np.isin(observed, axis).all():
            try:
                point_indices(points, [axis] * 3)
            except ValueError:
                continue
            matches.append(axis)
    if len(matches) != 1:
        raise ValueError(f"Native axis reconstruction is ambiguous: {len(matches)}")
    return [matches[0]] * 3


def native_run(case, density, observer, directory):
    request = {
        "id": case["id"],
        "op": "grid_occupancy",
        "include_points": True,
        "config": {"sphere_radius": case["radius"], "density": density},
        "atoms": [
            {"position": c, "radius_squared": float(np.float32(r) * np.float32(r))}
            for c, r in zip(case["centers"], case["radii"], strict=True)
        ],
    }
    captured, elapsed = observer.call(request, directory / "native")
    points = np.asarray(captured["points"], dtype=np.float32).astype(float)
    axes = native_axes(points, case["radius"])
    denominator = np.zeros(math.prod(map(len, axes)), dtype=bool)
    ids = point_indices(points, axes)
    assert len(np.unique(ids)) == len(ids)
    denominator[ids] = True
    numerator = np.zeros_like(denominator)
    numerator[ids] = captured["occupied_mask"]
    result = {
        "method": "stericx_native",
        "control": density,
        "elapsed_seconds": elapsed,
        "cost_scope": (
            "observation, per-point calls and JSON transfer; not production latency"
        ),
        "volume": float(np.float32(captured["aggregate"]["buried_volume"])),
        "percent": 100
        * captured["aggregate"]["buried_volume"]
        / captured["sphere_volume"],
        "sphere_nodes": int(denominator.sum()),
        "occupied_nodes": int(numerator.sum()),
        "accepted_config": captured["config"],
        "aggregate": captured["aggregate"],
        "percent_definition": "accepted native aggregate / native sphere volume * 100",
    }
    return [(result, axes, denominator, numerator)]


def dbstep_run(case, spacing, directory):
    # Explicit common cube covers the entire sampling sphere. Occupation outside
    # it cannot affect a buried-volume calculation. This is the documented
    # gridsize/rounded-linspace convention, not molecule-dependent box fitting.
    bound = 1.2 * case["radius"]
    axis = np.round(np.linspace(-bound, bound, 1 + round(2 * bound / spacing)), 8)
    axes = [axis] * 3
    points = tensor_points(axes)
    centers = np.asarray(case["centers"]).reshape((-1, 3))
    radii = np.asarray(case["radii"])
    options = SimpleNamespace(
        grid=spacing, radius=case["radius"], verbose=False, debug=False
    )
    stdout, stderr = io.StringIO(), io.StringIO()
    start = time.perf_counter()
    with redirect_stdout(stdout), redirect_stderr(stderr):
        occupied, tree, molecular_volume = sterics.occupied(
            points, centers, radii, np.zeros(3), options
        )
        percent, _ = sterics.buried_vol(
            occupied, tree, np.zeros(3), case["radius"], 0, options
        )
    elapsed = time.perf_counter() - start
    denominator = np.zeros(len(points), dtype=bool)
    denominator[tree.query_ball_point(np.zeros(3), case["radius"])] = True
    occupied_ids = point_indices(occupied, axes)
    numerator = np.zeros_like(denominator)
    numerator[occupied_ids[np.sum(occupied**2, axis=1) <= case["radius"] ** 2]] = True
    volume = percent / 100 * (4 * math.pi * case["radius"] ** 3 / 3)
    common = {
        "control": spacing,
        "axis_bounds": [-bound, bound],
        "axis_rounding_decimals": 8,
        "public_api": False,
        "api": "unmodified dbstep.sterics.occupied + buried_vol",
        "sphere_nodes": int(denominator.sum()),
        "occupied_nodes": int(numerator.sum()),
        "whole_grid_occupied_volume_not_buried_volume": float(molecular_volume),
    }
    tree_result = {
        **common,
        "method": "dbstep_tree_common_radii",
        "elapsed_seconds": elapsed,
        "cost_scope": "occupied + buried_vol, excluding grid construction",
        "volume": volume,
        "percent": float(percent),
    }
    start = time.perf_counter()
    with redirect_stdout(stdout), redirect_stderr(stderr):
        direct, direct_volume = sterics.occupied_direct(
            centers, radii, np.zeros(3), *axes, options
        )
        direct_percent, _ = sterics.buried_vol(
            direct, None, np.zeros(3), case["radius"], 0, options, grid_axes=axes
        )
    direct_elapsed = time.perf_counter() - start
    # Actual fast denominator's documented subtractive operation order.
    # Validate its full mask against the unmodified count routine exactly.
    x2, y2, z2 = (a * a for a in axes)
    fast_den = (
        z2[None, None, :]
        <= (case["radius"] ** 2 - x2[:, None, None]) - y2[None, :, None]
    ).reshape(-1)
    assert fast_den.sum() == sterics.count_grid_points_in_sphere(
        *axes, np.zeros(3), case["radius"]
    )
    direct_num = np.zeros_like(fast_den)
    direct_ids = point_indices(direct, axes)
    direct_num[direct_ids[np.sum(direct**2, axis=1) <= case["radius"] ** 2]] = True
    raw_percent = 100 * int(direct_num.sum()) / int(fast_den.sum())
    direct_result = {
        **common,
        "method": "dbstep_fast_common_radii",
        "elapsed_seconds": direct_elapsed,
        "cost_scope": "occupied_direct + buried_vol, excluding grid construction",
        "api": "unmodified dbstep.sterics.occupied_direct + buried_vol",
        "sphere_nodes": int(fast_den.sum()),
        "occupied_nodes": int(direct_num.sum()),
        "whole_grid_occupied_volume_not_buried_volume": float(direct_volume),
        "volume": direct_percent / 100 * (4 * math.pi * case["radius"] ** 3 / 3),
        "percent": float(direct_percent),
        "unclamped_percent": raw_percent,
        "clamped": raw_percent > 100,
        "tree_direct_occupancy_mask_differences": int(
            np.count_nonzero(direct_num != numerator)
        ),
        "denominator_mask_provenance": (
            "reconstructed subtractive mask; count checked against unmodified DBSTEP; "
            "not an independently observed per-node mask"
        ),
    }
    (directory / "dbstep_stdout.txt").write_text(stdout.getvalue())
    (directory / "dbstep_stderr.txt").write_text(stderr.getvalue())
    np.savez_compressed(
        directory / "dbstep_returned_points.npz",
        occupied_tree=occupied,
        occupied_direct=direct,
    )
    return [
        (tree_result, axes, denominator, numerator),
        (direct_result, axes, fast_den, direct_num),
    ]


def morfeus_run(case, density, directory):
    stdout, stderr = io.StringIO(), io.StringIO()
    start = time.perf_counter()
    with redirect_stdout(stdout), redirect_stderr(stderr):
        result = BuriedVolume(
            [0] + [6] * len(case["radii"]),
            [[0, 0, 0], *case["centers"]],
            metal_index=1,
            radii=[0, *case["radii"]],
            include_hs=True,
            radius=case["radius"],
            radii_scale=1,
            density=density,
        )
    elapsed = time.perf_counter() - start
    count = round((result._sphere.volume / density * 6 / math.pi) ** (1 / 3))
    axes = [np.linspace(-case["radius"], case["radius"], count)] * 3
    denominator = np.zeros(count**3, dtype=bool)
    numerator = np.zeros_like(denominator)
    denominator[point_indices(result._sphere.points, axes)] = True
    numerator[point_indices(result._buried_points, axes)] = True
    (directory / "morfeus_stdout.txt").write_text(stdout.getvalue())
    (directory / "morfeus_stderr.txt").write_text(stderr.getvalue())
    np.savez_compressed(
        directory / "morfeus_returned_points.npz",
        sphere=result._sphere.points,
        buried=result._buried_points,
    )
    row = {
        "method": "morfeus_public_common_radii",
        "control": density,
        "public_api": True,
        "api": "morfeus.BuriedVolume",
        "center_index": 1,
        "center_element": 0,
        "radii_scale": 1,
        "include_hs": True,
        "elapsed_seconds": elapsed,
        "cost_scope": "public BuriedVolume construction",
        "volume": result.buried_volume,
        "percent": 100 * result.fraction_buried_volume,
        "sphere_nodes": int(denominator.sum()),
        "occupied_nodes": int(numerator.sum()),
    }
    return [(row, axes, denominator, numerator)]


def error_metrics(volume, reference):
    value = fraction(volume)
    lower, upper = map(fraction, reference)
    error = [value - upper, value - lower]
    absolute = [
        Fraction(0) if lower <= value <= upper else min(map(abs, error)),
        max(map(abs, error)),
    ]
    relative = (
        None
        if lower <= 0
        else [downward(absolute[0] / upper), upward(absolute[1] / lower)]
    )
    return {
        "signed_error_interval": [downward(error[0]), upward(error[1])],
        "absolute_error_interval": [downward(absolute[0]), upward(absolute[1])],
        "relative_error_interval": relative,
        "relative_error_reason": "undefined at a zero or zero-containing reference"
        if relative is None
        else None,
    }


def execute_case(case, observer, output):
    directory = output / case["id"]
    directory.mkdir()
    save(directory / "input.json", case)
    start = time.perf_counter()
    depth = 11 if case["kind"] == "captured_molecular_frame" else 8
    bounds = certified_union_volume(
        case["centers"], case["radii"], case["radius"], max_depth=depth
    )
    bounds["elapsed_seconds"] = time.perf_counter() - start
    analytic = analytical_interval(case)
    bounds["analytical_interval"] = analytic
    nested = all(
        a["lower"] <= b["lower"] <= b["upper"] <= a["upper"]
        for a, b in itertools.pairwise(bounds["history"])
    )
    bounds["nested"] = nested
    bounds["analytical_enclosed"] = (
        None
        if analytic is None
        else bounds["lower"] <= analytic[0] <= analytic[1] <= bounds["upper"]
    )
    save(directory / "continuous_reference.json", bounds)
    rows = []
    previous = {}
    for index, (density, spacing) in enumerate(zip(DENSITIES, SPACINGS, strict=True)):
        level = directory / f"level{index}"
        level.mkdir()
        for operation, control in (
            (native_run, density),
            (morfeus_run, density),
            (dbstep_run, spacing),
        ):
            try:
                runs = (
                    operation(case, control, observer, level)
                    if operation is native_run
                    else operation(case, control, level)
                )
                for row, axes, denominator, numerator in runs:
                    start = time.perf_counter()
                    enclosure, differences = estimator_enclosure(
                        axes,
                        case["centers"],
                        case["radii"],
                        case["radius"],
                        denominator,
                        numerator,
                        row["volume"],
                    )
                    enclosure["construction_elapsed_seconds"] = (
                        time.perf_counter() - start
                    )
                    row["cell_enclosure"] = enclosure
                    if row["method"] == "stericx_native":
                        witnesses = native_dispute_witnesses(
                            case, axes, denominator, numerator, differences
                        )
                        with gzip.open(
                            level / "native_dispute_witnesses.json.gz",
                            "wt",
                            compresslevel=1,
                        ) as stream:
                            json.dump(witnesses, stream, allow_nan=False)
                        row["native_disputed_node_count"] = witnesses[
                            "disputed_node_count"
                        ]
                        row["native_unexplained_disputed_node_count"] = witnesses[
                            "unexplained_count"
                        ]
                    row.update(
                        error_metrics(
                            row["volume"],
                            analytic or [bounds["lower"], bounds["upper"]],
                        )
                    )
                    row["analytical_truth_inside_cell_enclosure"] = (
                        None
                        if analytic is None
                        else enclosure["certified_volume_lower"]
                        <= analytic[0]
                        <= analytic[1]
                        <= enclosure["certified_volume_upper"]
                    )
                    row["successive_volume_difference"] = (
                        None
                        if row["method"] not in previous
                        else row["volume"] - previous[row["method"]]
                    )
                    previous[row["method"]] = row["volume"]
                    row["case_id"] = case["id"]
                    row["level"] = index
                    row["actual_axis_count"] = list(map(len, axes))
                    np.savez_compressed(
                        level / (row["method"] + "_decisions.npz"),
                        x=axes[0],
                        y=axes[1],
                        z=axes[2],
                        denominator=denominator,
                        numerator=numerator,
                        **differences,
                    )
                    save(level / (row["method"] + ".json"), row)
                    rows.append(row)
            except Exception as error:
                failure = {
                    "case_id": case["id"],
                    "level": index,
                    "method": operation.__name__,
                    "control": control,
                    "status": "unsuccessful_comparison",
                    "error_type": type(error).__name__,
                    "error": str(error),
                }
                save(level / (operation.__name__ + "_failure.json"), failure)
                rows.append(failure)
        print(
            json.dumps(
                {"case": case["id"], "level_completed": index, "comparisons": len(rows)}
            ),
            flush=True,
        )
    save(directory / "results.json", rows)
    return {"case": case["id"], "continuous_reference": bounds, "rows": rows}


def boundary_witnesses(observer, output):
    """Actual one-point decisions plus exact rational margins and IEEE traces."""
    directory = output / "exact_point_boundaries"
    directory.mkdir()
    one = np.float32(1)
    below = float(np.nextafter(one, np.float32(-np.inf)))
    above = float(np.nextafter(one, np.float32(np.inf)))
    points = np.asarray(
        [
            [1, 0, 0],
            [below, 0, 0],
            [above, 0, 0],
            [-1, -0.0, 0],
            [-below, 0, 0],
            [-above, 0, 0],
            [0.6, 0.8, 0],
            [0.8, 0.6, 0],
            [1, 2**-15, 0],
            [0, 0, 0],
        ],
        dtype=np.float32,
    ).astype(float)
    request = {
        "id": "exact_point_boundaries",
        "op": "raw_point_decisions",
        "points": points.tolist(),
        "atoms": [{"position": [0, 0, 0], "radius_squared": 1}],
        "config": {"sphere_radius": 2.5},
        "include_points": True,
    }
    actual, elapsed = observer.call(request, directory / "native")
    options = SimpleNamespace(grid=0.1, radius=2.5, verbose=False, debug=False)
    occupied, _, _ = sterics.occupied(
        points, np.zeros((1, 3)), np.ones(1), np.zeros(3), options
    )
    db_mask = np.asarray([any(np.array_equal(p, q) for q in occupied) for p in points])
    exact = union_mask_exact(points, [[0, 0, 0]], [1])
    direct = union_mask_f64(points, [[0, 0, 0]], [1])
    rows = []
    for index, point in enumerate(points):
        dx, dy, dz = point.astype(np.float32)
        x2, y2, z2 = np.float32(dx * dx), np.float32(dy * dy), np.float32(dz * dz)
        xy, squared = np.float32(x2 + y2), np.float32(np.float32(x2 + y2) + z2)
        decision = bool(squared <= np.float32(1))
        certificate = f32_predicate_certificate(point, [0, 0, 0], 1)
        disputed = bool(actual["occupied_mask"][index]) != bool(exact[index])
        row = {
            "point": point.tolist(),
            "point_bits": [float_bits(v) for v in point],
            "native": actual["occupied_mask"][index],
            "dbstep": bool(db_mask[index]),
            "independent_f64": bool(direct[index]),
            "independent_exact": bool(exact[index]),
            "native_differs_from_exact": disputed,
            "certificate": certificate,
            "ieee_trace": {
                name: {"value": float(value), "bits": float_bits(value)}
                for name, value in (
                    ("dx", dx),
                    ("dy", dy),
                    ("dz", dz),
                    ("dx_squared", x2),
                    ("dy_squared", y2),
                    ("dz_squared", z2),
                    ("xy_squared", xy),
                    ("distance_squared", squared),
                )
            },
            "trace_matches_observed": decision == actual["occupied_mask"][index],
            "disagreement_explained": not disputed
            or (
                decision == actual["occupied_mask"][index]
                and certificate["classification"] == "ambiguous"
            ),
        }
        rows.append(row)
    result = {
        "rows": rows,
        "elapsed_native_observation_seconds": elapsed,
        "morfeus": "not applicable: public API chooses its grid; no point-mask input",
        "native_exact_disagreements": sum(
            row["native_differs_from_exact"] for row in rows
        ),
        "unexplained": [
            row
            for row in rows
            if not row["trace_matches_observed"] or not row["disagreement_explained"]
        ],
        "release_gate_scored": False,
    }
    save(directory / "results.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--observer", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--case", action="append")
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    os.sched_setaffinity(0, {0, 1})
    reference = ROOT / ".stericx/paton/reference/DBSTEP"
    commit = subprocess.check_output(
        ["git", "-C", str(reference), "rev-parse", "HEAD"], text=True
    ).strip()
    assert commit == PIN
    assert not subprocess.check_output(
        ["git", "-C", str(reference), "status", "--porcelain"], text=True
    )
    assert (
        record(sterics.__file__)["sha256"]
        == record(reference / "dbstep/sterics.py")["sha256"]
    )
    build_manifest_path = args.observer.resolve().parents[1] / "manifest.json"
    build_manifest = json.loads(build_manifest_path.read_text())
    assert record(args.observer) == build_manifest["binaries"]["observer"]
    assert (
        record(ROOT / "src/geometry/buried_volume.rs")["sha256"]
        == build_manifest["bv_prefix_sha256"]
    )
    bindings = [
        record(__file__),
        record(ROOT / "scripts/dbstep_resolution_reference.py"),
        record(
            ROOT / "docs/validation/dbstep/resolution_analysis/ANALYSIS_PROTOCOL.md"
        ),
        record(args.observer),
        record(reference / "dbstep/sterics.py"),
        record(sterics.__file__),
        record(sys.modules["morfeus.buried_volume"].__file__),
        record(sys.modules["morfeus.geometry"].__file__),
        record(build_manifest_path),
        record(ROOT / "tests/test_dbstep_volume_convergence.py"),
        record(
            ROOT
            / "docs/validation/dbstep/resolution_analysis/CONVERGENCE_CELL_BOUNDS.md"
        ),
    ]
    started = {
        "kind": "fixed_volume_convergence",
        "release_gate_scored": False,
        "bindings": bindings,
        "argv": sys.argv,
        "python": sys.version,
        "dbstep_commit": commit,
        "dependencies": {
            name: version(name) for name in ("numpy", "scipy", "dbstep", "morfeus-ml")
        },
        "affinity": sorted(os.sched_getaffinity(0)),
        "thread_environment": {
            name: os.environ.get(name)
            for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
        },
        "densities": DENSITIES,
        "spacings": SPACINGS,
    }
    save(output / "started.json", started)
    snapshots = output / "source_snapshots"
    snapshots.mkdir()
    for index, binding in enumerate(bindings):
        shutil.copyfile(
            binding["path"], snapshots / f"{index}_{Path(binding['path']).name}"
        )
    cases = primitive_cases() + molecular_cases()
    if args.case:
        cases = [case for case in cases if case["id"] in args.case]
        if len(cases) != len(set(args.case)):
            raise ValueError("Unknown or duplicated requested case")
    save(output / "cases.json", cases)
    observer = Observer(args.observer, output)
    results = []
    try:
        for case in cases:
            results.append(execute_case(case, observer, output))
        boundaries = boundary_witnesses(observer, output)
    finally:
        observer.close()
    for binding in bindings:
        if record(binding["path"]) != binding:
            raise ValueError(f"Immutable source changed: {binding['path']}")
    files = [record(path) for path in sorted(output.rglob("*")) if path.is_file()]
    summary = {
        "case_count": len(results),
        "comparison_count": sum(len(case["rows"]) for case in results),
        "unsuccessful_comparisons": [
            row for case in results for row in case["rows"] if row.get("status")
        ],
        "analytical_box_failures": [
            case["case"]
            for case in results
            if case["continuous_reference"]["analytical_enclosed"] is False
        ],
        "non_nested_box_failures": [
            case["case"]
            for case in results
            if not case["continuous_reference"]["nested"]
        ],
        "analytical_cell_bound_failures": [
            row
            for case in results
            for row in case["rows"]
            if row.get("analytical_truth_inside_cell_enclosure") is False
        ],
        "release_gate_scored": False,
        "exact_boundary_native_disagreements": boundaries["native_exact_disagreements"],
        "exact_boundary_unexplained": boundaries["unexplained"],
        "native_unexplained_disputed_node_count": sum(
            row.get("native_unexplained_disputed_node_count", 0)
            for case in results
            for row in case["rows"]
        ),
    }
    save(output / "summary.json", summary)
    save(
        output / "manifest.json",
        {**started, "files": files, "summary": record(output / "summary.json")},
    )
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
