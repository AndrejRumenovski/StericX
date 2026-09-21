"""Analytic self-checks of the independent reference, not tests of StericX."""

import json
from pathlib import Path

import numpy as np
from geometry_reference import pyramidalization_ref, sterimol_exact, volume_on_grid
from scipy.spatial.transform import Rotation


def main():
    checks = []
    xyz = np.array([[0.0, 0, 0], [0, 0, 1], [3, 0, 1], [-3, 0, 1]])
    radii = np.array([1.2, 1.7, 1.7, 1.7])
    p = sterimol_exact(xyz, radii, xyz[0], 1, 0)
    assert np.allclose([p["l"], p["b1"], p["b5"]], [2.7, 1.7, 4.7], rtol=0, atol=1e-13)
    checks.append("Closed-form four-sphere Sterimol L=2.7 B1=1.7 B5=4.7")
    for rotation in Rotation.random(25, random_state=1964).as_matrix():
        transformed = xyz @ rotation.T + [7.0, -8.0, 3.0]
        q = sterimol_exact(transformed, radii, transformed[0], 1, 0)
        assert np.allclose(list(p.values()), list(q.values()), rtol=0, atol=1e-12)
    checks.append(
        "Analytic reference support minimum invariant over25 independent rotations"
    )
    pyramid = {"donor": 0, "neighbors": [1, 2, 3]}
    tetra = np.array([[0.0, 0, 0], [1, 1, 1], [1, -1, -1], [-1, 1, -1]])
    q = pyramidalization_ref(pyramid, tetra)["independent"]
    assert abs(q["pyr_p"] - 4 / (3 * np.sqrt(3))) < 1e-14
    assert abs(q["pyr_alpha"] - np.rad2deg(np.arccos(np.sqrt(2 / 3)))) < 1e-12
    checks.append("Tetrahedral determinant4/(3sqrt3), alpha35.26438968 degrees")
    points = np.array(
        [[x, y, z] for x in [-0.25, 0.25] for y in [-0.25, 0.25] for z in [-0.25, 0.25]]
    )
    full = volume_on_grid(points, np.zeros((1, 3)), np.array([1.0]), 1.0)
    assert abs(full["buried_volume"] - 4 * np.pi / 3) < 1e-14
    assert np.allclose(full["quadrants"], np.pi / 3)
    assert np.allclose(full["octants"], np.pi / 6)
    empty = volume_on_grid(points, np.array([[9.0, 9.0, 9.0]]), np.array([1.0]), 1.0)
    assert empty["buried_volume"] == 0
    checks.append("Analytic fully occupied/empty sphere and all regional fractions")
    root = Path(__file__).resolve().parents[1] / "geometry"
    rows = [
        json.loads(line)
        for line in (root / "reference_results.jsonl").read_text().splitlines()
    ]
    matched = 0
    max_error = 0.0
    for row in rows:
        data = row.get("buried_volume", {})
        if "independent" not in data:
            continue
        for key, value in data["independent"].items():
            error = abs(value - data["morfeus"][key])
            max_error = max(max_error, error)
        matched += 1
    assert max_error < 1e-10
    checks.append(
        f"Independent volume reference vs actual Morfeus: {matched} cases, maxdiff{max_error}"
    )
    (root / "reference_selfcheck.json").write_text(
        json.dumps(
            {
                "checks": checks,
                "passed": len(checks),
                "meaning": "Analytic sanity checks of audit reference; not evidence of SUT correctness.",
            },
            indent=2,
        )
        + "\n"
    )
    print("\n".join(checks))


if __name__ == "__main__":
    main()
