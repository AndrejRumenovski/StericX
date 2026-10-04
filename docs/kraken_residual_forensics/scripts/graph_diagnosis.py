"""Distinguish graph representation discrepancies from descriptor input changes."""
import collections
import json
import tarfile

import pandas as pd
from rdkit import Chem

from analyze import read_jsonl, references
from freeze import OUT, save, sha


def adjacency_code(mol, strip_metals=False):
    graph = Chem.RWMol()
    mapping = {}
    for atom in mol.GetAtoms():
        z = atom.GetAtomicNum()
        if z == 1 or (strip_metals and z in {26, 27, 28, 29, 30, 44, 45, 46, 47, 48, 76, 77, 78, 79, 80}):
            continue
        a = Chem.Atom(z)
        a.SetNoImplicit(True)
        mapping[atom.GetIdx()] = graph.AddAtom(a)
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in mapping and b in mapping:
            graph.AddBond(mapping[a], mapping[b], Chem.BondType.SINGLE)
    return Chem.MolToSmiles(graph.GetMol(), canonical=True, isomericSmiles=False)


def main():
    out = OUT / "provenance_analysis/graph_diagnosis"
    out.mkdir(exist_ok=False)
    table = pd.read_csv(OUT / "provenance_analysis/conformer_provenance.csv.gz")
    ids = set(table.loc[table.matches_API_connectivity == False, "id"])
    inventory = {r["sdf_path"]: r for r in read_jsonl(OUT / "frozen/inventory.jsonl.gz") if r["id"] in ids}
    _, metadata = references()
    results = []
    with tarfile.open(OUT / "frozen/inputs_and_reference.tar.gz") as archive:
        for member in archive:
            if member.name not in inventory:
                continue
            inv = inventory[member.name]
            mol = Chem.MolFromMolBlock(archive.extractfile(member).read().decode(), removeHs=False)
            ref = Chem.MolFromSmiles(metadata[inv["molecule_id"]]["smiles"])
            same = adjacency_code(mol) == adjacency_code(ref)
            nonmetal_same = adjacency_code(mol, True) == adjacency_code(ref, True)
            typ = "bond_order_or_formal_valence_representation" if same else ("metal_coordination_representation" if nonmetal_same else "nonmetal_connectivity_discrepancy")
            results.append({"id": inv["id"], "molecule_id": inv["molecule_id"], "classification": typ,
                "heavy_adjacency_equal": same, "metal_stripped_adjacency_equal": nonmetal_same,
                "descriptor_effect": "No direct effect in this experiment: donor neighbors, atomic identities/radii and positions are the descriptor inputs; SDF and primary donor neighbor sets already verified equal.",
                "limitation": "Atom identity and historical chemical provenance not proved by a graph representation match."})
    pd.DataFrame(results).to_csv(out / "cases.csv", index=False)
    save(out / "complete.json", {"N": len(results), "counts": dict(collections.Counter(r["classification"] for r in results)),
        "script_sha256": sha(__file__), "cases_sha256": sha(out / "cases.csv")})
    print((out / "complete.json").read_text())


if __name__ == "__main__":
    main()
