# ruff: noqa: E501
"""Freeze adversarial inputs before SUT/reference evaluation; independent RDKit graphs."""

import copy
import hashlib
import json
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "geometry"


def atom(element, position, radius=None):
    a = {"element": element, "position": [float(x) for x in position]}
    if radius is not None:
        a["radius"] = radius
    return a


def case(name, elements, xyz, donor=0, neighbors=None, **kw):
    return {
        "id": name,
        "atoms": [atom(e, p) for e, p in zip(elements, xyz, strict=True)],
        "donor": donor,
        "neighbors": neighbors or [1, 2, 3],
        "reference": (neighbors or [1, 2, 3])[0],
        "config": {
            "sphere_radius": 3.5,
            "density": 0.01,
            "center_distance": 2.28,
            "radii_scale": 1.17,
            "include_hydrogens": False,
        },
        "category": "base",
        **kw,
    }


def main():
    rng = np.random.default_rng(20260919)
    bases = []
    specs = [
        ("methylphosphine", "CP", "P"),
        ("dimethylphosphine", "CPC", "P"),
        ("trimethylphosphine", "CP(C)C", "P"),
        ("triethylphosphine", "CCP(CC)CC", "P"),
        ("triisopropylphosphine", "CC(C)P(C(C)C)C(C)C", "P"),
        ("triphenylphosphine", "c1ccccc1P(c1ccccc1)c1ccccc1", "P"),
        ("asymmetric_phosphine", "CP(CC)c1ccccc1", "P"),
        ("fluorophosphine", "FP(F)F", "P"),
        ("methylamine", "CN", "N"),
        ("dimethylamine", "CNC", "N"),
        ("trimethylamine", "CN(C)C", "N"),
        ("aniline", "Nc1ccccc1", "N"),
    ]
    for i, (name, smiles, donor_element) in enumerate(specs):
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        params = AllChem.ETKDGv3()
        params.randomSeed = 20260919 + i
        assert AllChem.EmbedMolecule(mol, params) == 0
        status = AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
        donor = next(
            a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == donor_element
        )
        neighbors = [a.GetIdx() for a in mol.GetAtomWithIdx(donor).GetNeighbors()]
        heavy = [j for j in neighbors if mol.GetAtomWithIdx(j).GetSymbol() != "H"]
        c = case(
            name,
            [a.GetSymbol() for a in mol.GetAtoms()],
            mol.GetConformer().GetPositions(),
            donor,
            neighbors,
            smiles=smiles,
            chemical_bonds=[
                [b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondTypeAsDouble()]
                for b in mol.GetBonds()
            ],
            mmff_status=status,
            expected_donor=donor,
            expected_bonded_neighbors=neighbors,
        )
        c["reference"] = heavy[0]
        bases.append(c)
    # Analytic toy structures: graph labels intentionally explicit, never inferred from SUT.
    xyz = np.array(
        [
            [0, 0, 0],
            [1.8, 0, 0.3],
            [-0.9, 1.5588457268, 0.3],
            [-0.9, -1.5588457268, 0.3],
        ]
    )
    bases.append(
        case(
            "symmetric_tertiary",
            ["P", "C", "C", "C"],
            xyz,
            expected_bonded_neighbors=[1, 2, 3],
        )
    )
    t = copy.deepcopy(bases[-1])
    t["id"] = "asymmetric_toy"
    t["atoms"] += [atom("C", [3.1, 0.7, 0.8]), atom("O", [-2.9, -0.3, 0.5])]
    bases.append(t)
    bases.append(
        case(
            "planar_donor",
            ["P", "C", "C", "C"],
            [[0, 0, 0], [1.8, 0, 0], [-0.9, 1.5588457268, 0], [-0.9, -1.5588457268, 0]],
            expected_bonded_neighbors=[1, 2, 3],
        )
    )
    for z in [1e-7, 1e-4, 0.01, 2.0]:
        c = copy.deepcopy(next(x for x in bases if x["id"] == "planar_donor"))
        c["id"] = f"pyramid_height_{z}"
        for a in c["atoms"][1:]:
            a["position"][2] = z
        bases.append(c)
    bases.append(
        case(
            "nearly_collinear",
            ["P", "C", "C", "C"],
            [[0, 0, 0], [1.8, 0, 0], [1.8, 1e-5, 0], [-1.8, 0, 1e-5]],
            expected_bonded_neighbors=[1, 2, 3],
        )
    )
    bases.append(
        case(
            "coincident_neighbor",
            ["P", "C", "C", "C"],
            [[0, 0, 0], [0, 0, 0], [-0.9, 1.5, 0.3], [-0.9, -1.5, 0.3]],
            expected_bonded_neighbors=[1, 2, 3],
        )
    )
    c = copy.deepcopy(next(x for x in bases if x["id"] == "asymmetric_toy"))
    c["id"] = "full_sphere_cover"
    for a in c["atoms"]:
        a["radius"] = 10.0
    bases.append(c)
    # Connected P-H graph with distant nonbonded heavy contacts: recreates historical nearest-heavy failure.
    for name in ["methylphosphine", "dimethylphosphine"]:
        c = copy.deepcopy(next(x for x in bases if x["id"] == name))
        c["id"] = f"history_{name}_nearby_nonbonded"
        d = np.array(c["atoms"][c["donor"]]["position"])
        c["atoms"] += [
            atom("C", d + np.array([2.7, 0, 0])),
            atom("C", d + np.array([0, 2.8, 0])),
        ]
        bases.append(c)
    # Threshold perturbations define geometry-only ambiguous contacts, not experimental bond truth.
    c = copy.deepcopy(next(x for x in bases if x["id"] == "asymmetric_toy"))
    for side, value in [
        ("below", np.nextafter(np.float32(1.3 * (1.07 + 0.76)), np.float32(0))),
        ("at", np.float32(1.3 * (1.07 + 0.76))),
        ("above", np.nextafter(np.float32(1.3 * (1.07 + 0.76)), np.float32(np.inf))),
    ]:
        d = copy.deepcopy(c)
        d["id"] = "bond_threshold_" + side
        d["atoms"][1]["position"] = [float(value), 0, 0]
        d["expected_bonded_neighbors"] = None
        d["chemical_expectation"] = (
            "No bond identity can be established from one cutoff perturbation alone."
        )
        bases.append(d)
    c = copy.deepcopy(next(x for x in bases if x["id"] == "asymmetric_toy"))
    c["id"] = "multiple_phosphorus"
    c["atoms"].append(atom("P", [9, 0, 0]))
    bases.append(c)
    # Independent small substituent Sterimol cases (BV trivalence errors are expected).
    for name, smiles in [
        ("methane", "C"),
        ("ethane", "CC"),
        ("isobutane", "CC(C)C"),
        ("tertbutylbenzene", "CC(C)(C)c1ccccc1"),
    ]:
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        AllChem.EmbedMolecule(mol, randomSeed=20260919)
        AllChem.MMFFOptimizeMolecule(mol)
        coords = mol.GetConformer().GetPositions()
        h = (
            next(
                a.GetIdx()
                for a in mol.GetAtomWithIdx(0).GetNeighbors()
                if a.GetSymbol() == "H"
            )
            if name != "tertbutylbenzene"
            else 4
        )
        c = case(
            name,
            [a.GetSymbol() for a in mol.GetAtoms()],
            coords,
            donor=0,
            neighbors=[h],
            sterimol_only=True,
        )
        c["attach"] = h
        c["sterimol_neighbor"] = 0
        bases.append(c)
    cases = copy.deepcopy(bases)
    chosen = [
        "methylphosphine",
        "trimethylphosphine",
        "triphenylphosphine",
        "asymmetric_phosphine",
        "trimethylamine",
        "planar_donor",
    ]
    for name in chosen:
        b = next(c for c in bases if c["id"] == name)
        xyz = np.array([a["position"] for a in b["atoms"]])
        matrices = list(
            Rotation.from_euler(
                "xyz", [[0.7, 0, 0], [0, 0.7, 0], [0, 0, 0.7]]
            ).as_matrix()
        ) + list(Rotation.random(100, random_state=rng).as_matrix())
        for i, R in enumerate(matrices):
            c = copy.deepcopy(b)
            c.update(
                id=f"{name}__rigid_{i:03}",
                category="rigid",
                parent=name,
                transform_index=i,
            )
            shift = rng.uniform(-50, 50, 3)
            for a, p in zip(c["atoms"], xyz @ R.T + shift, strict=True):
                a["position"] = p.tolist()
            c["rotation"] = R.tolist()
            c["translation"] = shift.tolist()
            cases.append(c)
        for i in range(20):
            c = copy.deepcopy(b)
            c.update(
                id=f"{name}__permutation_{i:02}", category="permutation", parent=name
            )
            order = rng.permutation(len(b["atoms"]))
            inverse = np.argsort(order)
            c["atoms"] = [copy.deepcopy(b["atoms"][j]) for j in order]
            c["donor"] = int(inverse[b["donor"]])
            c["neighbors"] = [int(inverse[j]) for j in b["neighbors"]]
            c["reference"] = int(inverse[b["reference"]])
            c["permutation"] = order.tolist()
            cases.append(c)
    for name in [
        "methylphosphine",
        "trimethylphosphine",
        "triphenylphosphine",
        "asymmetric_toy",
    ]:
        b = next(c for c in bases if c["id"] == name)
        for density in [0.1, 0.01, 0.001, 0.0001, 0.008]:
            c = copy.deepcopy(b)
            c.update(
                id=f"{name}__density_{density}", parent=name, category="convergence"
            )
            c["config"]["density"] = density
            cases.append(c)
    for name in ["asymmetric_toy", "methylphosphine"]:
        for magnitude in [1e3, 1e5, 1e7, 1e10, 1e20]:
            c = copy.deepcopy(next(c for c in bases if c["id"] == name))
            c.update(
                id=f"{name}__translation_{magnitude}", category="precision", parent=name
            )
            for a in c["atoms"]:
                a["position"] = (
                    np.array(a["position"])
                    + np.array([magnitude, -magnitude, magnitude])
                ).tolist()
            cases.append(c)
    # Separate protocol cases probe radii choice, H occupancy, and explicit center.
    for name in ["methylphosphine", "trimethylamine", "asymmetric_toy"]:
        for flag in [True, False]:
            c = copy.deepcopy(next(c for c in bases if c["id"] == name))
            c.update(id=f"{name}__include_h_{flag}", category="protocol", parent=name)
            c["config"]["include_hydrogens"] = flag
            cases.append(c)
    data = {
        "seed": 20260919,
        "rationale": "Graphs from explicit SMILES; synthetic cases defined by coordinates; no SUT outputs used for selection. All rotations incl100 random rotations/representative. No outlier exclusion.",
        "cases": cases,
    }
    path = OUT / "inputs.json"
    encoded = json.dumps(data, indent=2) + "\n"
    if path.exists():
        assert path.read_text() == encoded, (
            "Frozen inputs changed; use a new audit directory"
        )
    else:
        path.write_text(encoded)
    manifest = {
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "n_cases": len(cases),
        "n_base": len(bases),
        "file": "geometry/inputs.json",
    }
    (OUT / "input_freeze.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(manifest)


if __name__ == "__main__":
    main()
