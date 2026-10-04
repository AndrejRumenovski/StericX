"""Input, graph, frame, conformer and independent reference-snapshot provenance."""
import collections
import datetime
import gzip
import json
import math
import tarfile

import numpy as np
import pandas as pd
from rdkit import Chem

from analyze import FIELDS, read_jsonl, references
from freeze import AUDIT, OUT, save, sha


def rmsd(a, b):
    aa, bb = a - a.mean(axis=0), b - b.mean(axis=0)
    u, _, vt = np.linalg.svd(aa.T @ bb)
    fix = np.eye(3)
    fix[2, 2] = np.linalg.det(u @ vt)
    rotation = u @ fix @ vt
    return float(np.sqrt(np.mean(np.sum((aa @ rotation - bb)**2, axis=1))))


def main():
    out = OUT / "provenance_analysis"
    out.mkdir(exist_ok=False)
    inventory = {r["sdf_path"]: r for r in read_jsonl(OUT / "frozen/inventory.jsonl.gz")}
    requests = {r["id"]: r for r in read_jsonl(OUT / "frozen/requests.jsonl.gz")}
    original = {r["id"]: r for r in read_jsonl(OUT / "frozen/prior_morfeus_analytic.jsonl.gz")}
    native = {r["id"]: r for r in read_jsonl(OUT / "frozen/current_predictions.jsonl.gz")}
    refs, molecules = references()
    expected_graph = {}
    for mid, m in molecules.items():
        mol = Chem.MolFromSmiles(m["smiles"])
        expected_graph[mid] = Chem.MolToSmiles(mol, isomericSmiles=False) if mol is not None else None
    save(out / "plan.json", {"created_utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "script_sha256": sha(__file__), "frozen_manifest_sha256": sha(OUT / "frozen/manifest.json"),
        "scope": "All 1541 ligands and all 31611 exported conformers, without selection by residual",
        "historical_identity_standard": "Class A requires affirmative mapping to actual historical calculation coordinates/retention records, not a DFT label or low descriptor error.",
        "RMSD": "Input float64 SDF to actual native f32 coordinates; intralibrary first-conformer heavy-atom RMSD uses one exact graph isomorphism (upper bound to symmetry-minimized RMSD). Historical RMSD is unavailable, never filled with zero.",
        "frames": "Actual native center plus independently reconstructed raw-vector source frame; original native private frame stream is independently hashed in its frozen validation receipt. Source frame xz convention is derived without fitting."})
    rows, first_mol = [], {}
    with tarfile.open(OUT / "frozen/inputs_and_reference.tar.gz") as archive, gzip.open(out / "frames.jsonl.gz", "wt") as stream:
        for member in archive:
            if member.name not in inventory:
                continue
            inv = inventory[member.name]
            key, mid = inv["id"], inv["molecule_id"]
            req = requests[key]
            data = archive.extractfile(member).read().decode()
            mol = Chem.MolFromMolBlock(data, sanitize=False, removeHs=False)
            sanitized = Chem.Mol(mol)
            flag = str(Chem.SanitizeMol(sanitized, catchErrors=True))
            graph = None
            heavy = None
            if flag == "SANITIZE_NONE":
                heavy = Chem.RemoveHs(sanitized)
                graph = Chem.MolToSmiles(heavy, isomericSmiles=False)
            mapped_rmsd = None
            mapping = None
            if heavy is not None:
                if mid not in first_mol:
                    first_mol[mid] = (key, heavy)
                refkey, first = first_mol[mid]
                if heavy.GetNumAtoms() == first.GetNumAtoms() and heavy.GetNumBonds() == first.GetNumBonds():
                    match = heavy.GetSubstructMatch(first, useChirality=False)
                    if len(match) == first.GetNumAtoms():
                        aa = heavy.GetConformer().GetPositions()[list(match)]
                        bb = first.GetConformer().GetPositions()
                        mapped_rmsd = rmsd(aa, bb)
                        mapping = list(match)
            else:
                refkey = None
            xyz = np.array([a["position"] for a in req["atoms"]])
            f32xyz = np.array([a["position"] for a in native[key]["atoms"]])
            d = req["donor"]
            native_center = np.array(native[key]["center"])
            source_center = np.array(original[key]["primary_center"])
            native_z = native_center - f32xyz[d]
            source_z = source_center - xyz[d]
            native_z /= np.linalg.norm(native_z)
            source_z /= np.linalg.norm(source_z)
            angle = float(np.rad2deg(np.arccos(np.clip(native_z @ source_z, -1, 1))))
            frames = []
            for n in original[key]["primary_neighbors"]:
                w = xyz[n] - source_center
                x = w - (w @ source_z) * source_z
                x /= np.linalg.norm(x)
                y = np.cross(source_z, x)
                frames.append({"plane_atom": n, "source_basis_rows_xyz": [x.tolist(), y.tolist(), source_z.tolist()]})
            stream.write(json.dumps({"id": key, "donor": d, "native_neighbors": req["neighbors"],
                "source_neighbors": original[key]["primary_neighbors"], "neighbor_elements": inv["neighbor_elements"],
                "native_center": native_center.tolist(), "source_center": source_center.tolist(),
                "native_BV_z": native_z.tolist(), "source_BV_z": source_z.tolist(),
                "native_Sterimol_axis": (-native_z).tolist(), "source_Sterimol_axis": (-source_z).tolist(),
                "frames": frames, "source_last_plane": original[key]["primary_neighbors"][-1],
                "native_reduction": "permutation-symmetric mean for total/near/far; extrema over all three frames",
                "source_reduction": "last plane for total/near/far; extrema over all three frames",
                "quadrants_xy": ["++", "-+", "--", "+-"],
                "native_octants": "same cyclic quadrants z>=0 then z<0; zero planes positive",
                "source_octants": "strict-open Morfeus octants; key/insertion ordering separately captured in source and Morfeus outputs",
                "heavy_atom_mapping_to_first": mapping}, separators=(",", ":")) + "\n")
            rows.append({"id": key, "molecule_id": mid, "conformer_id": inv["conformer_id"],
                "sdf_sha256": inv["sdf_sha256"], "source_atom_count": len(xyz),
                "native_atom_identity_equal": all(a["element"] == b["element"] for a, b in zip(req["atoms"], native[key]["atoms"], strict=True)),
                "canonical_connectivity": graph, "matches_API_connectivity": graph == expected_graph[mid] if graph is not None else None,
                "sanitization": flag, "source_native_neighbor_sets_equal": set(req["neighbors"]) == set(original[key]["primary_neighbors"]),
                "export_to_native_RMSD_A": rmsd(f32xyz, xyz),
                "export_to_native_unaligned_RMSD_A": float(np.sqrt(np.mean(np.sum((f32xyz-xyz)**2, axis=1)))),
                "first_conformer_id_for_RMSD": refkey, "within_ensemble_mapped_RMSD_A": mapped_rmsd,
                "frame_center_shift_A": float(np.linalg.norm(native_center-source_center)), "frame_axis_angle_deg": angle,
                "native_donor": d, "source_donor": d, "donor_H": inv["P_H_count"],
                "historical_geometry_RMSD_A": None, "historical_conformer_energy": None,
                "historical_ensemble_membership": "unverified", "historical_weight": None})
    df = pd.DataFrame(rows)
    df.to_csv(out / "conformer_provenance.csv.gz", index=False)
    ligrows = []
    for mid, g in df.groupby("molecule_id"):
        api_ids = molecules[mid]["conformers_id"]
        assert set(api_ids) == set(g.conformer_id)
        ligrows.append({"molecule_id": int(mid), "provenance_class": "F", "historical_exact_geometry_verified": False,
            "geometry_origin": "Public Kraken DFT SDF exports; not regenerated by StericX",
            "classification_reason": "Exact historical coordinate precision and retained conformer IDs not supplied; no residual-based classification as exact",
            "status": "UNRESOLVABLE FROM AVAILABLE DATA", "available_conformer_count": len(g),
            "source_API_membership_equal": True, "historical_conformer_count": None,
            "connectivity_mismatch_count": int((g.matches_API_connectivity == False).sum()),
            "unsanitized_count": int((g.sanitization != "SANITIZE_NONE").sum()),
            "unique_connectivity_count": int(g.canonical_connectivity.nunique()),
            "maximum_intracohort_RMSD_A": g.within_ensemble_mapped_RMSD_A.max(),
            "maximum_native_coordinate_RMSD_A": g.export_to_native_RMSD_A.max(),
            "maximum_frame_axis_angle_deg": g.frame_axis_angle_deg.max(),
            "suspected_ensemble_mismatch_evidence": "Earlier independently frozen two-ID-block dossier; not proof of historical inclusion" if mid == 369 else None})
    lf = pd.DataFrame(ligrows)
    lf.to_csv(out / "ligand_provenance.csv", index=False)
    hist = pd.read_csv(OUT / "frozen/historical_reference.csv", float_precision="round_trip").set_index("Unnamed: 0")
    si_path = AUDIT / "primary_si/DFT_data.csv"
    si = pd.read_csv(si_path, float_precision="round_trip").set_index("ID")
    comparisons = []
    for mid in sorted(molecules):
        for desc, (_, _, prop, _) in FIELDS.items():
            for red in ("min", "max", "delta", "vburminconf"):
                key = prop + "_" + red
                api = refs[mid][prop][red]
                old = hist.loc[mid, key] if key in hist else None
                original_si = si.loc[mid, key] if mid in si.index and key in si else None
                comparisons.append({"molecule_id": mid, "metric": desc + "_" + red, "API": api,
                    "historical_CSV": old, "original_SI": original_si,
                    "API_minus_historical": api-old if old is not None else None,
                    "API_minus_SI": api-original_si if original_si is not None else None})
    pd.DataFrame(comparisons).to_csv(out / "reference_snapshot_comparison.csv.gz", index=False)
    save(out / "complete.json", {"ligands": len(lf), "conformers": len(df),
        "provenance_counts": dict(collections.Counter(lf.provenance_class)),
        "graph_mismatches": int((df.matches_API_connectivity == False).sum()),
        "sanitization_failures": int((df.sanitization != "SANITIZE_NONE").sum()),
        "all_neighbor_sets_equal": bool(df.source_native_neighbor_sets_equal.all()),
        "maximum_native_coordinate_RMSD_A": float(df.export_to_native_RMSD_A.max()),
        "maximum_frame_axis_angle_deg": float(df.frame_axis_angle_deg.max()),
        "SI_path": str(si_path), "SI_sha256": sha(si_path),
        "original_SI_is_not_available_for_ids": sorted(set(molecules)-set(si.index)),
        "raw_records_sha256": sha(out / "conformer_provenance.csv.gz")})
    print((out / "complete.json").read_text())


if __name__ == "__main__":
    main()
