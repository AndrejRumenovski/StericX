#!/usr/bin/env python3
"""Diagnose DBSTEP sphere-count consistency without modifying either engine.

Exact rational predicates refer to the represented float64 coordinates/radii.
They are a diagnostic reference, not a replacement StericX chemistry kernel.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import itertools
import json
import subprocess
import sys
from contextlib import redirect_stderr, redirect_stdout
from fractions import Fraction
from pathlib import Path

import numpy as np
from dbstep import Dbstep, constants, sterics
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[4]
REFERENCE = ROOT / ".stericx/paton/reference/DBSTEP"
PIN = "4365353e0a2f24d05336e619759e350bc630a3cb"


def record(path):
    path = Path(path).resolve()
    with path.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest}


def rational(value):
    value = Fraction(value)
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
        "float": float(value),
    }


def lens_fraction(radius, atom_radius, distance):
    """Analytic two-sphere intersection divided by the evaluation-sphere volume.

    Pi cancels, permitting an exact rational result for these axis-aligned cases.
    """
    r, a, d = map(Fraction, (radius, atom_radius, distance))
    if d >= r + a:
        return Fraction(0)
    if d <= abs(r - a):
        return min(r, a) ** 3 / r**3
    return (
        (r + a - d) ** 2 * (d**2 + 2 * d * (r + a) - 3 * (r - a) ** 2) / (16 * d * r**3)
    )


def minimal_boundary():
    point = np.array([[0.6, 0.8, 0.0]])
    axes = [point[:, index] for index in range(3)]
    margin = sum(Fraction(float(x)) ** 2 for x in point[0]) - 1
    return {
        "point": point[0].tolist(),
        "point_hex": [float(x).hex() for x in point[0]],
        "radius": 1.0,
        "exact_squared_distance_minus_radius_squared": rational(margin),
        "exact_included": margin <= 0,
        "fast_denominator_included": bool(
            sterics.count_grid_points_in_sphere(*axes, np.zeros(3), 1.0)
        ),
        "summed_numerator_included": bool(np.sum(point**2) <= 1.0),
        "tree_denominator_included": bool(
            cKDTree(point).query_ball_point(np.zeros(3), 1.0)
        ),
        "subtractive_remainder": float(1.0 - 0.6**2 - 0.8**2),
        "summed_squared_distance": float(np.sum(point**2)),
    }


def public_cases(spacing):
    radius = 1.0
    # Explicit covering box fixes the actual lattice for every atom position.
    axis = np.round(np.linspace(-3.0, 3.0, 1 + round(6.0 / spacing)), 8)
    grid = np.array(np.meshgrid(axis, axis, axis, indexing="ij")).reshape(3, -1).T
    sphere = np.sum(grid**2, axis=1) <= radius**2
    sphere_grid = grid[sphere]
    fast_denominator = sterics.count_grid_points_in_sphere(
        axis, axis, axis, np.zeros(3), radius
    )
    tree_denominator = len(cKDTree(grid).query_ball_point(np.zeros(3), radius))
    # Exact sphere membership can only occur within the coordinate-wise box.
    small_axis = [Fraction(float(x)) for x in axis if abs(x) <= radius]
    exact_sphere = [
        point
        for point in itertools.product(small_axis, repeat=3)
        if sum(x * x for x in point) <= Fraction(radius) ** 2
    ]
    rows = []
    for position in (
        [0.0, 0.0, 0.0],
        [0.4, 0.0, 0.0],
        [0.8, 0.0, 0.0],
        [0.0, 0.0, 1.0],
    ):
        settings = dict(
            atom1=1,
            volume=True,
            sterimol=False,
            r=radius,
            grid=spacing,
            scalevdw=1.0,
            radii="bondi",
            quiet=True,
            gridsize="-3,3,-3,3,-3,3",
        )
        stdout, stderr = io.StringIO(), io.StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            result = Dbstep.from_coords(
                ["Bq", "C"], np.array([[0.0, 0.0, 0.0], position]), **settings
            )
        atom_radius = float(constants.bondi["C"])
        atom = np.sum((sphere_grid - position) ** 2, axis=1) <= atom_radius**2
        occupied = int(np.count_nonzero(atom))
        exact_center = list(map(Fraction, position))
        exact_occupied = sum(
            sum((x - c) ** 2 for x, c in zip(point, exact_center, strict=True))
            <= Fraction(atom_radius) ** 2
            for point in exact_sphere
        )
        raw_percent = occupied * spacing**3 / (fast_denominator * spacing**3) * 100
        expected_public = min(100.0, raw_percent)
        assert (
            abs(float(result.bur_vol) - expected_public)
            <= 64 * np.finfo(float).eps * 100
        )
        analytical = lens_fraction(radius, atom_radius, max(map(abs, position)))
        rows.append(
            {
                "geometry": {
                    "elements": ["Bq", "C"],
                    "coordinates": [[0, 0, 0], position],
                },
                "settings": settings,
                "axis": axis.tolist(),
                "atom_radius_angstrom": atom_radius,
                "sphere_count_fast_denominator": int(fast_denominator),
                "sphere_count_summed_numerator_predicate": len(sphere_grid),
                "sphere_count_tree": tree_denominator,
                "sphere_count_exact_rational": len(exact_sphere),
                "occupied_count_summed": occupied,
                "occupied_count_exact_rational": exact_occupied,
                "summed_sphere_points_unoccupied": len(sphere_grid) - occupied,
                "raw_mixed_predicate_percent": raw_percent,
                "public_percent": float(result.bur_vol),
                "public_clamped": raw_percent > 100,
                "consistent_summed_predicate_percent": 100
                * occupied
                / len(sphere_grid),
                "consistent_exact_predicate_percent": 100
                * exact_occupied
                / len(exact_sphere),
                "analytical_continuum_percent": float(100 * analytical),
                "analytical_continuum_fraction": rational(analytical),
                "whole_grid_occupied_volume": float(result.occ_vol),
                "stdout": stdout.getvalue(),
                "stderr": stderr.getvalue(),
            }
        )
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    commit = subprocess.check_output(
        ["git", "-C", str(REFERENCE), "rev-parse", "HEAD"], text=True
    ).strip()
    assert commit == PIN
    assert not subprocess.check_output(
        ["git", "-C", str(REFERENCE), "status", "--porcelain"], text=True
    ).strip()
    files = [Path(__file__), Path(sterics.__file__), Path(Dbstep.__file__)]
    for installed in files[1:]:
        assert (
            record(installed)["sha256"]
            == record(REFERENCE / "dbstep" / installed.name)["sha256"]
        )
    result = {
        "kind": "unmodified_dbstep_sphere_predicate_diagnostic",
        "reference_commit": commit,
        "command": sys.argv,
        "python": sys.version,
        "bindings": [record(path) for path in files],
        "minimal_boundary": minimal_boundary(),
        "public_cases": [],
        "scientific_scope": (
            "DBSTEP internal predicate consistency and exact binary-coordinate point "
            "classification; continuum values are separate from finite-grid values. "
            "No StericX or DBSTEP implementation, radius or acceptance gate is changed."
        ),
    }
    for spacing in (0.2, 0.1, 0.05):
        result["public_cases"].extend(public_cases(spacing))
    result["complete"] = True
    destination = args.output / "result.json"
    destination.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {"cases": len(result["public_cases"]), "result": record(destination)}
        )
    )


if __name__ == "__main__":
    main()
