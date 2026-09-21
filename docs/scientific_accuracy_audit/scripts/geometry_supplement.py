# ruff: noqa: E501
"""Additional prespecified numerical witnesses; freezes before SUT evaluation."""

import copy
import hashlib
import json
from pathlib import Path

import numpy as np

R = Path(__file__).resolve().parents[1] / "geometry"
base = next(
    c
    for c in json.loads((R / "inputs.json").read_text())["cases"]
    if c["id"] == "asymmetric_toy"
)
cases = []
for element in ["Xe", "He", "Fe", "Na", "Se", "B", "Si", "Cl", "Br", "I"]:
    c = copy.deepcopy(base)
    c.update(id="radius_" + element, category="radii")
    c["atoms"].append({"element": element, "position": [3.2, 2.2, -1.0]})
    cases.append(c)
for epsilon in [0.0, 1e-8, 1e-6, 1e-4, 0.001]:
    c = {
        "id": f"axis_antiparallel_{epsilon}",
        "category": "numerical",
        "atoms": [
            {"element": "H", "position": [0, 0, 0]},
            {"element": "C", "position": [epsilon, 0, -1.0]},
            {"element": "O", "position": [2.0, 1.0, -1.0]},
        ],
        "donor": 0,
        "reference": 1,
        "neighbors": [1, 2, 1],
        "sterimol_only": True,
        "attach": 0,
        "sterimol_neighbor": 1,
        "config": base["config"],
    }
    cases.append(c)
for length in [
    float(np.nextafter(np.float32(np.sqrt(np.finfo(np.float32).eps)), np.float32(0))),
    float(np.sqrt(np.finfo(np.float32).eps)),
    float(
        np.nextafter(np.float32(np.sqrt(np.finfo(np.float32).eps)), np.float32(np.inf))
    ),
]:
    c = copy.deepcopy(cases[-1])
    c.update(id=f"short_axis_{length}", category="numerical")
    c["atoms"][1]["position"] = [0, 0, length]
    cases.append(c)
cases.append(
    {
        "id": "phosphine_PH3",
        "category": "chemical_edge",
        "atoms": [
            {"element": "P", "position": [0, 0, 0]},
            {"element": "H", "position": [1.2, 0, 0.72]},
            {"element": "H", "position": [-0.6, 1.0392304845413263, 0.72]},
            {"element": "H", "position": [-0.6, -1.0392304845413263, 0.72]},
        ],
        "donor": 0,
        "reference": 1,
        "neighbors": [1, 2, 3],
        "expected_bonded_neighbors": [1, 2, 3],
        "config": copy.deepcopy(base["config"]),
    }
)
c = copy.deepcopy(base)
c.update(id="small_sphere_fully_inside_donor", category="chemical_edge")
c["config"]["center_distance"] = 1.0
c["config"]["sphere_radius"] = 0.5
c["config"]["density"] = 0.0001
cases.append(c)
raw = []
for side, distance in [
    ("inside", float(np.nextafter(np.float32(1), np.float32(0)))),
    ("on", 1.0),
    ("outside", float(np.nextafter(np.float32(1), np.float32(np.inf)))),
]:
    raw.append(
        {
            "id": "atom_boundary_" + side,
            "points": [[distance, 0, 0]],
            "atoms": [{"position": [0, 0, 0], "radius_squared": 1.0}],
            "sphere_radius": 3.5,
            "expected_occupied": side != "outside",
        }
    )
for axis in range(3):
    for coordinate in [
        -float(np.nextafter(np.float32(0), np.float32(1))),
        -0.0,
        0.0,
        float(np.nextafter(np.float32(0), np.float32(1))),
    ]:
        p = [0.25, 0.25, 0.25]
        p[axis] = coordinate
        raw.append(
            {
                "id": f"partition_axis{axis}_{float(coordinate).hex()}",
                "points": [p],
                "atoms": [{"position": [0, 0, 0], "radius_squared": 1.0}],
                "sphere_radius": 3.5,
            }
        )
for radius in [
    3.5,
    float(np.nextafter(np.float32(3.5), np.float32(0))),
    float(np.nextafter(np.float32(3.5), np.float32(np.inf))),
]:
    c = copy.deepcopy(base)
    c.update(id=f"sphere_boundary_{radius}", category="boundary", include_points=True)
    c["config"]["sphere_radius"] = radius
    c["config"]["density"] = float(8 * radius**3 / 5**3)
    cases.append(c)
path = R / "supplement_inputs.json"
encoded = json.dumps(
        {
            "cases": cases,
            "raw": raw,
            "rationale": "Prespecified radius-fallback, antiparallel/small-axis and exactfloat boundary witnesses before SUT/reference evaluation.",
        },
        indent=2,
    ) + "\n"
if path.exists():
    assert path.read_text() == encoded, "Frozen supplemental inputs changed"
else:
    path.write_text(encoded)
(R / "supplement_input_freeze.json").write_text(
    json.dumps(
        {
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "n_cases": len(cases),
            "n_raw": len(raw),
        },
        indent=2,
    )
    + "\n"
)
