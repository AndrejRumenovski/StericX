"""Check preserved near-axis failures after an exact coordinate permutation.

This is an audit experiment, not a production implementation. The independent
reference uses direct axial/radial projections and analytic support envelopes.
A 90-degree X rotation only permutes/signs coordinates, avoiding coordinate
rounding changes from trigonometric rotation of the input.
"""

import copy
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
from geometry_reference import sterimol_exact

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "geometry/alignment/rotation"
OUT.mkdir(exist_ok=True)
SOURCE = ROOT / "geometry/alignment/requests.jsonl"
BINARY = ROOT / "frozen/bin/stericx-audit-observer"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rotate(vector):
    x, y, z = vector
    return [x, -z, y]


requests = []
for original in map(json.loads, SOURCE.read_text().splitlines()):
    transformed = copy.deepcopy(original)
    transformed["id"] += "__rot90x"
    for atom in transformed["atoms"]:
        atom["position"] = rotate(atom["position"])
    transformed["center"] = rotate(transformed["center"])
    requests.append(transformed)
input_path = OUT / "requests.jsonl"
serialized = "".join(json.dumps(request) + "\n" for request in requests)
if input_path.exists():
    assert input_path.read_text() == serialized
else:
    input_path.write_text(serialized)
input_freeze = {
    "source_sha256": sha(SOURCE),
    "requests_sha256": sha(input_path),
    "observer_sha256": sha(BINARY),
    "transform": "[x,y,z] -> [x,-z,y], proper rigid 90-degree X rotation",
}
input_manifest = OUT / "input_freeze.json"
if input_manifest.exists():
    assert json.loads(input_manifest.read_text()) == input_freeze
else:
    input_manifest.write_text(json.dumps(input_freeze, indent=2) + "\n")

output_path = OUT / "sut_outputs.jsonl"
output_manifest = OUT / "sut_freeze.json"
if not output_path.exists():
    with input_path.open("rb") as src, output_path.open("wb") as dst:
        with (OUT / "stderr.txt").open("wb") as err:
            subprocess.run([str(BINARY)], stdin=src, stdout=dst, stderr=err, check=True)
if output_manifest.exists():
    assert json.loads(output_manifest.read_text())["sut_sha256"] == sha(output_path)
else:
    output_manifest.write_text(
        json.dumps({**input_freeze, "sut_sha256": sha(output_path)}, indent=2) + "\n"
    )

original_observations = {
    row["id"]: row
    for row in map(json.loads, (ROOT / "geometry/alignment/sut_outputs.jsonl").read_text().splitlines())
}
observations = {
    row["id"]: row for row in map(json.loads, output_path.read_text().splitlines())
}
results = []
for request in requests:
    observed = observations[request["id"]]
    original = original_observations[request["id"].removesuffix("__rot90x")]
    xyz = np.array([atom["position"] for atom in observed["atoms"]])
    radii = np.array([atom["vdw_radius"] for atom in observed["atoms"]])
    reference = sterimol_exact(xyz, radii, np.array(observed["center"]), request["donor"])
    original_xyz = np.array([atom["position"] for atom in original["atoms"]])
    original_reference = sterimol_exact(
        original_xyz, radii, np.array(original["center"]), request["donor"]
    )
    results.append({
        "id": request["id"],
        "reference": reference,
        "sut": observed["sterimol_dummy_raw"],
        "errors": {key: observed["sterimol_dummy_raw"][key] - reference[key] for key in reference},
        "sut_rigid_rotation_change": {
            key: observed["sterimol_dummy_raw"][key] - original["sterimol_dummy_raw"][key]
            for key in reference
        },
        "independent_rigid_rotation_change": {
            key: reference[key] - original_reference[key] for key in reference
        },
    })
(OUT / "reference_results.json").write_text(json.dumps(results, indent=2) + "\n")
print(json.dumps(results, indent=2))
