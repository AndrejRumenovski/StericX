"""Record actual native private frames and independently constructed source frames."""
import gzip
import json
import math

import numpy as np
from morfeus.geometry import kabsch_rotation_matrix

from analyze import read_jsonl
from freeze import OUT, save, sha


def transverse(axis):
    z = axis / np.linalg.norm(axis)
    t = math.hypot(z[0], z[1])
    if t == 0:
        return [[1. if z[2] >= 0 else -1., 0., 0.], [0., 1., 0.], z.tolist()]
    rotaxis = np.array([-z[1]/t, z[0]/t, 0.])
    if z[2] >= 0:
        c = math.sqrt((1+z[2])*.5)
        s = t/(2*c)
    else:
        s = math.sqrt((1-z[2])*.5)
        c = t/(2*s)
    qx, qy, _ = rotaxis*s
    x = [1-2*qy*qy, 2*qx*qy, -2*qy*c]
    y = [2*qx*qy, 1-2*qx*qx, 2*qx*c]
    return [x, y, z.tolist()]


def main():
    out = OUT / "frame_forensics"
    out.mkdir(exist_ok=False)
    manifest = json.loads((OUT / "frozen/current_observations_manifest.json").read_text())
    source = manifest["lanes"]["kraken_topology_all_bins"]["stdout"]
    assert sha(source["path"]) == source["sha256"]
    reqs = {r["id"]: r for r in read_jsonl(OUT / "campaigns/recovered_precision/requests.jsonl.gz")}
    refs = {r["id"]: r for r in read_jsonl(OUT / "campaigns/recovered_precision/morfeus_observations.jsonl.gz")}
    save(out / "plan.json", {"script_sha256": sha(__file__), "native_private_source": source,
        "source_convention_sha256": sha(OUT / "conventions/plan.json"),
        "native_frame": "Actual captured native basis/origin and per-frame bins from frozen validated baseline; baseline XYZ precision and radius table",
        "source_frame": "Primary raw-center/donor-neighbor convention independently reconstructed on the maximum precision API export; Morfeus Kabsch Sterimol frame explicitly recorded",
        "native_sterimol_frame": "Shortest-arc basis reconstructed from actual captured f32 center/coordinates promoted to f64; mathematical convention, not a claim of Python/Rust intermediate bit identity"})
    n = 0
    with gzip.open(out / "all_frames.jsonl.gz", "wt", compresslevel=3) as dst:
        for r in read_jsonl(source["path"]):
            if r["id"] not in reqs:
                continue
            req, ref = reqs[r["id"]], refs[r["id"]]
            xyz = np.array([a["position"] for a in req["atoms"]])
            c = np.array(ref["center"])
            z = c-xyz[req["donor"]]
            z /= np.linalg.norm(z)
            frames = []
            for neighbor in ref["neighbors"]:
                w = xyz[neighbor]-c
                x = w-(w@z)*z
                x /= np.linalg.norm(x)
                y = np.cross(z, x)
                frames.append({"plane": neighbor, "basis_rows_xyz": [x.tolist(), y.tolist(), z.tolist()]})
            nc = np.array(r["center"])
            donor32 = np.array(r["atoms"][req["donor"]]["position"])
            axis = -z
            R = kabsch_rotation_matrix(axis.reshape(1, 3), np.array([[1., 0., 0.]]), center=False)
            dst.write(json.dumps({"id": req["id"], "donor": req["donor"], "donor_neighbor_elements": [req["atoms"][i]["element"] for i in req["neighbors"]],
                "native_atom_input_source": "frozen/current_predictions.jsonl.gz (same ID)",
                "source_atom_input_source": "campaigns/recovered_precision/requests.jsonl.gz (same ID)",
                "native_center": r["center"], "native_neighbors": r["dump"]["neighbors"],
                "native_grid_count": r["dump"]["grid_count"],
                "native_BV_frames": [{k: v for k, v in f.items() if k != "aligned_atoms"} for f in r["dump"]["orientations"]],
                "native_Sterimol_basis_rows": transverse(donor32-nc),
                "source_center": ref["center"], "source_neighbors": ref["neighbors"], "source_BV_frames": frames,
                "source_Sterimol_rotation_matrix_to_x": R.tolist(),
                "source_near_definition": "sum of negative-z octants; far is positive-z octants",
                "native_zero_planes": "assigned positive; regional occupied fractions independently normalized",
                "source_zero_planes": "strict-open octants in contemporary Morfeus; source grid .001 has no zero planes",
                "atom_order": "all three planes for extrema; native totals average three; historical totals retain final original-order neighbor"}, separators=(",", ":")) + "\n")
            n += 1
    save(out / "complete.json", {"N": n, "native_frames": 3*n, "source_frames": 3*n,
        "data_sha256": sha(out / "all_frames.jsonl.gz"), "plan_sha256": sha(out / "plan.json")})


if __name__ == "__main__":
    main()
