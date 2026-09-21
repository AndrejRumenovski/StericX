# ruff: noqa: E501
"""Minimized analytic B1, chemistry/frame witnesses and CLI auto-selection audit."""

import copy
import hashlib
import itertools
import json
import math
import subprocess
from pathlib import Path

import numpy as np
from geometry_reference import bv_refs, sterimol_refs, unit

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "geometry" / "focused"
OUT.mkdir(exist_ok=True)
OB = ROOT / "frozen/bin/stericx-audit-observer"
CLI = ROOT / "frozen/bin/stericx"


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def atom(e, p):
    return {"element": e, "position": list(p)}


base = json.loads((ROOT / "geometry/inputs.json").read_text())["cases"]
cfg = base[0]["config"]
cases = []
angle = math.radians(0.5)
x = 3 * math.cos(angle)
y = 3 * math.sin(angle)
cases.append(
    {
        "id": "minimal_b1_half_degree",
        "atoms": [
            atom("H", [0, 0, 0]),
            atom("C", [0, 0, 1]),
            atom("C", [x, y, 1]),
            atom("C", [-x, -y, 1]),
        ],
        "donor": 0,
        "reference": 1,
        "neighbors": [1, 2, 3],
        "sterimol_only": True,
        "config": cfg,
    }
)
tie = {
    "id": "tied_axis",
    "atoms": [
        atom("P", [0, 0, 0]),
        atom("C", [1.8, 0, 0]),
        atom("C", [0, 1.8, 0]),
        atom("C", [0, 0, 1.8]),
        atom("C", [3.3, 0, 0]),
        atom("F", [4.6, 0, 0]),
    ],
    "donor": 0,
    "reference": 1,
    "neighbors": [1, 2, 3],
    "config": cfg,
}
cases.append(tie)
for perm in itertools.permutations([1, 2, 3]):
    c = copy.deepcopy(tie)
    c["id"] = "tied_axis_permutation_" + "".join(map(str, perm))
    c["atoms"] = [tie["atoms"][0]] + [tie["atoms"][i] for i in perm] + tie["atoms"][4:]
    cases.append(c)
# Analytic perfect tetrahedron, orthogonal, acute, and coincident vectors.
for name, xyz in [
    ("tetrahedron", [[1, 1, 1], [1, -1, -1], [-1, 1, -1]]),
    ("orthogonal", [[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
    ("acute", [[1, 0, 0.1], [1, 0.1, 0], [1, -0.1, 0]]),
]:
    cases.append(
        {
            "id": name,
            "atoms": [atom("P", [0, 0, 0])] + [atom("C", p) for p in xyz],
            "donor": 0,
            "reference": 1,
            "neighbors": [1, 2, 3],
            "config": cfg,
        }
    )
# For scalar sphere intersection check analytic lens at all densities, expose raw bins even if error.
ph3 = next(
    c
    for c in json.loads((ROOT / "geometry/supplement_inputs.json").read_text())["cases"]
    if c["id"] == "phosphine_PH3"
)
for density in [0.1, 0.01, 0.001, 0.0001]:
    c = copy.deepcopy(ph3)
    c["id"] = "ph3_lens_" + str(density)
    c["config"]["density"] = density
    c["dump"] = True
    cases.append(c)
# Freeze before evaluation. Refuse to replace existing captured inputs/output.
inputs = OUT / "inputs.json"
blob = json.dumps(cases, indent=2) + "\n"
if inputs.exists():
    assert inputs.read_text() == blob
else:
    inputs.write_text(blob)
request = OUT / "requests.jsonl"
text = "".join(json.dumps(c) + "\n" for c in cases)
if request.exists():
    assert request.read_text() == text
else:
    request.write_text(text)
output = OUT / "sut_outputs.jsonl"
if not output.exists():
    with (
        request.open("rb") as src,
        output.open("wb") as dst,
        (OUT / "stderr.txt").open("wb") as err,
    ):
        subprocess.run([str(OB)], stdin=src, stdout=dst, stderr=err, check=True)
(OUT / "freeze.json").write_text(
    json.dumps(
        {
            "input_sha256": sha(inputs),
            "request_sha256": sha(request),
            "output_sha256": sha(output),
            "observer_sha256": sha(OB),
        },
        indent=2,
    )
    + "\n"
)
sut = {x["id"]: x for x in map(json.loads, output.read_text().splitlines())}
# Independently check all analytic witness values after SUT capture.
references = {}
for c in cases[:8]:
    s = sut[c["id"]]
    xyz = np.array([a["position"] for a in s["atoms"]])
    references[c["id"]] = sterimol_refs(
        c,
        xyz,
        controlled_center=s.get("center")
        if isinstance(s.get("center"), list)
        else None,
    )
R = 3.5
r = 1.8 * 1.17
d = 2.28
lens = (
    math.pi * (R + r - d) ** 2 * (d * d + 2 * d * (R + r) - 3 * (R - r) ** 2) / (12 * d)
)
references["analytic_lens"] = {
    "sphere_radius": R,
    "atom_radius": r,
    "separation": d,
    "volume": lens,
    "percent": lens / (4 * math.pi * R**3 / 3) * 100,
    "equation": "pi*(R+r-d)^2*(d^2+2d(R+r)-3(R-r)^2)/(12d)",
    "errors": {
        c["id"]: sut[c["id"]]["dump"]["orientations"][0]["buried_volume"] - lens
        for c in cases
        if c["id"].startswith("ph3_lens")
    },
}
# Historical selector reconstruction with independent graph truth and SAME mathematical integration.
history = []
for name in [
    "trimethylphosphine",
    "history_methylphosphine_nearby_nonbonded",
    "history_dimethylphosphine_nearby_nonbonded",
]:
    c = next(x for x in base if x["id"] == name)
    xyz = np.array([a["position"] for a in c["atoms"]], dtype=np.float32).astype(float)
    donor = c["donor"]
    heavy = [i for i, a in enumerate(c["atoms"]) if i != donor and a["element"] != "H"]
    old = sorted(heavy, key=lambda i: np.linalg.norm(xyz[i] - xyz[donor]))[:3]
    new = c["neighbors"]

    def center(indices, xyz=xyz, donor=donor):
        return xyz[donor] + cfg["center_distance"] * unit(
            -sum(unit(xyz[j] - xyz[donor]) for j in indices)
        )

    correct = center(new)
    wrong = center(old)
    history.append(
        {
            "id": name,
            "chemical_graph_neighbors": new,
            "historical_nearest_heavy": old,
            "correct_center": correct.tolist(),
            "historical_center": wrong.tolist(),
            "center_distance_error": float(np.linalg.norm(wrong - correct)),
            "correct_reference": bv_refs(c, xyz, correct)["independent"],
            "historical_center_reference": bv_refs(c, xyz, wrong)["independent"],
            "limitation": "Chemical graph determines bonded substituents; neither geometric center estimates a measured metal position.",
        }
    )
references["historical_frame"] = history
(OUT / "reference_results.json").write_text(json.dumps(references, indent=2) + "\n")
# Freeze XYZ input identities and raw CLI responses before interpreting donor/axis results.
cli_cases = (
    [
        c
        for c in base
        if c["category"] == "base"
        and not c.get("sterimol_only")
        and not any("radius" in a for a in c["atoms"])
    ]
    + [c for c in cases if c["id"].startswith("tied_axis")]
    + [ph3]
)
xyzdir = OUT / "xyz"
xyzdir.mkdir(exist_ok=True)
cmds = []
for c in cli_cases:
    path = xyzdir / (c["id"] + ".xyz")
    blob = (
        str(len(c["atoms"]))
        + "\n"
        + c["id"]
        + "\n"
        + "".join(
            a["element"]
            + " "
            + " ".join(format(v, ".17g") for v in a["position"])
            + "\n"
            for a in c["atoms"]
        )
    )
    if path.exists():
        assert path.read_text() == blob
    else:
        path.write_text(blob)
    cmd = [
        str(CLI),
        "descriptors",
        str(path),
        "--format",
        "json",
        "--donor-element",
        c["atoms"][c["donor"]]["element"],
    ]
    # Only tied-axis family uses chemically permuted order as target; no explicit donor/reference override.
    cmds.append({"id": c["id"], "input_sha256": sha(path), "argv": cmd})
(OUT / "cli_input_manifest.json").write_text(json.dumps(cmds, indent=2) + "\n")
cliout = OUT / "cli_outputs.json"
if not cliout.exists():
    rows = []
    for spec in cmds:
        p = subprocess.run(spec["argv"], capture_output=True, text=True)
        rows.append(
            spec | {"returncode": p.returncode, "stdout": p.stdout, "stderr": p.stderr}
        )
    cliout.write_text(json.dumps(rows, indent=2) + "\n")
(OUT / "cli_freeze.json").write_text(
    json.dumps(
        {
            "executable_sha256": sha(CLI),
            "input_manifest_sha256": sha(OUT / "cli_input_manifest.json"),
            "output_sha256": sha(cliout),
        },
        indent=2,
    )
    + "\n"
)
print(
    "focused witness inputs", len(cases), "CLI cases", len(cmds), "analytic lens", lens
)
