"""Frozen invalid-number and boundary diagnostics; these never adjust production."""

import copy
import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "geometry" / "numerical"
OUT.mkdir(exist_ok=True)
base = next(
    c
    for c in json.loads((ROOT / "geometry/inputs.json").read_text())["cases"]
    if c["id"] == "asymmetric_toy"
)
cases = []
for label, value in [
    ("nan", "NaN"),
    ("positive_inf", "Infinity"),
    ("negative_inf", "-Infinity"),
    ("huge_finite", 1e38),
]:
    for target in ["coordinate", "radius"]:
        c = copy.deepcopy(base)
        c["id"] = label + "_" + target
        if target == "coordinate":
            c["atoms"][-1]["position"][0] = value
        else:
            c["atoms"][-1]["radius"] = value
        cases.append(c)
# Atom centers immediately inside/on/outside the integration-sphere boundary.
for label, x in [
    ("inside", 3.499999761581421),
    ("on", 3.5),
    ("outside", 3.500000238418579),
]:
    c = copy.deepcopy(base)
    c["id"] = "atom_center_sphere_boundary_" + label
    c["center"] = [0, 0, -2.28]
    c["atoms"].append({"element": "C", "position": [x, 0, -2.28]})
    c["dump"] = True
    cases.append(c)
path = OUT / "inputs.json"
data = json.dumps(cases, indent=2) + "\n"
if path.exists():
    assert path.read_text() == data
else:
    path.write_text(data)
req = OUT / "requests.jsonl"
text = "".join(json.dumps(c) + "\n" for c in cases)
if req.exists():
    assert req.read_text() == text
else:
    req.write_text(text)
output = OUT / "sut_outputs.jsonl"
binary = ROOT / "frozen/bin/stericx-audit-observer"
if not output.exists():
    with (
        req.open("rb") as src,
        output.open("wb") as dst,
        (OUT / "stderr.txt").open("wb") as err,
    ):
        subprocess.run([str(binary)], stdin=src, stdout=dst, stderr=err, check=True)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


(OUT / "freeze.json").write_text(
    json.dumps({p.name: sha(p) for p in [path, req, output, binary]}, indent=2) + "\n"
)
print(output)
