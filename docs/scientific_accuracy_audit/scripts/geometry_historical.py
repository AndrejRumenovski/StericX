"""Reconstruct the six actual historical zero-asymmetry ligands independently."""

import copy
import csv
import hashlib
import io
import json
import subprocess
from pathlib import Path

import numpy as np
from geometry_reference import bv_refs, unit
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "geometry/historical"
OUT.mkdir(exist_ok=True)
if (OUT / "freeze.json").exists():
    frozen = json.loads((OUT / "freeze.json").read_text())
    revision = frozen["historical_revision"]
    oldcsv = (OUT / "historical_comparison.csv").read_bytes()
    oldsrc = (OUT / "historical_buried_volume.rs").read_bytes()
    for name, content in [("historical_comparison.csv", oldcsv), ("historical_buried_volume.rs", oldsrc)]:
        assert hashlib.sha256(content).hexdigest() == frozen["files"][name]
else:
    revision = subprocess.check_output(["git", "rev-parse", "f9a8ac8^"], text=True).strip()
    oldcsv = subprocess.check_output(
        ["git", "show", revision + ":docs/study_004/kraken_dft_scaled_comparison.csv"]
    )
    oldsrc = subprocess.check_output(
        ["git", "show", revision + ":src/geometry/buried_volume.rs"]
    )
    (OUT / "historical_comparison.csv").write_bytes(oldcsv)
    (OUT / "historical_buried_volume.rs").write_bytes(oldsrc)
ids = [
    int(r["Source_ID"])
    for r in csv.DictReader(io.StringIO(oldcsv.decode()))
    if float(r["stericx_on_dft"]) == 0
]
inventory = list(
    map(
        json.loads,
        (ROOT / "kraken/prepared_all/inventory.jsonl").read_text().splitlines(),
    )
)
cases = []
for item in inventory:
    if item["molecule_id"] not in ids:
        continue
    path = ROOT / "kraken" / item["sdf_path"]
    mol = next(
        iter(
            Chem.ForwardSDMolSupplier(
                io.BytesIO(path.read_bytes()), removeHs=False, sanitize=False
            )
        )
    )
    xyz = mol.GetConformer().GetPositions()
    d = item["donor"]
    neighbors = [a.GetIdx() for a in mol.GetAtomWithIdx(d).GetNeighbors()]
    heavy = [i for i, a in enumerate(mol.GetAtoms()) if i != d and a.GetSymbol() != "H"]
    old = sorted(heavy, key=lambda i: float(np.linalg.norm(xyz[i] - xyz[d])))[:3]
    cases.append(
        {
            "id": item["id"],
            "source_sdf": str(path.relative_to(ROOT)),
            "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "atoms": [
                {"element": a.GetSymbol(), "position": xyz[a.GetIdx()].tolist()}
                for a in mol.GetAtoms()
            ],
            "donor": d,
            "reference": next(
                i for i in neighbors if mol.GetAtomWithIdx(i).GetSymbol() != "H"
            ),
            "neighbors": neighbors,
            "historical_heavy_neighbors": old,
            "config": {
                "sphere_radius": 3.5,
                "density": 0.01,
                "center_distance": 2.28,
                "radii_scale": 1.17,
                "include_hydrogens": False,
            },
        }
    )
requests = OUT / "requests.jsonl"
blob = "".join(json.dumps(c) + "\n" for c in cases)
if requests.exists():
    assert requests.read_text() == blob
else:
    requests.write_text(blob)
output = OUT / "sut_outputs.jsonl"
binary = ROOT / "frozen/bin/stericx-audit-observer"
if not output.exists():
    with (
        requests.open("rb") as src,
        output.open("wb") as dst,
        (OUT / "stderr.txt").open("wb") as err,
    ):
        subprocess.run([str(binary)], stdin=src, stdout=dst, stderr=err, check=True)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


(OUT / "freeze.json").write_text(
    json.dumps(
        {
            "historical_revision": revision,
            "zero_ligand_ids": ids,
            "n_conformers": len(cases),
            "files": {
                p.name: sha(p)
                for p in [
                    requests,
                    output,
                    binary,
                    OUT / "historical_comparison.csv",
                    OUT / "historical_buried_volume.rs",
                ]
            },
        },
        indent=2,
    )
    + "\n"
)
sut = {x["id"]: x for x in map(json.loads, output.read_text().splitlines())}
results = []
for c in cases:
    s = sut[c["id"]]
    xyz = np.array([a["position"] for a in s["atoms"]])
    d = c["donor"]
    old = c["historical_heavy_neighbors"]
    new = c["neighbors"]

    def center(indices, xyz=xyz, d=d):
        return xyz[d] + 2.28 * unit(-sum(unit(xyz[i] - xyz[d]) for i in indices))

    correct = center(new)
    historical = center(old)
    legacy = copy.deepcopy(c)
    legacy["neighbors"] = old
    results.append(
        {
            "id": c["id"],
            "chemical_bond_neighbors": new,
            "neighbor_elements": [c["atoms"][i]["element"] for i in new],
            "historical_heavy_neighbors": old,
            "sut_bonded_neighbors": [x["index"] for x in s["bonded_neighbors"]],
            "center_displacement": float(np.linalg.norm(correct - historical)),
            "current_center_reference": bv_refs(c, xyz, correct)["independent"],
            "historical_frame_reference": bv_refs(legacy, xyz, historical)[
                "independent"
            ],
            "current_sut": s["buried_volume"],
        }
    )
(OUT / "reference_results.json").write_text(json.dumps(results, indent=2) + "\n")
print("frozen", len(cases), "conformers", ids)
for row in results:
    print(
        row["id"],
        row["center_displacement"],
        row["historical_frame_reference"]["max_delta_qvbur"],
        row["current_center_reference"]["max_delta_qvbur"],
        row["current_sut"].get("max_delta_qvbur"),
    )
