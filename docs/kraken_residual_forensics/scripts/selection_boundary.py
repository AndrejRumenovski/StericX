"""Locate the floating-point grid decision behind ligand 390 selection."""
import copy
import json

import numpy as np

from analyze import read_jsonl
from freeze import OUT, save, sha
from minimal_reference import REF
from run_native import native, PROCESSES


def main():
    out = OUT / "selection_boundary_v2"
    out.mkdir(exist_ok=False)
    save(out / "plan.json", {"script_sha256": sha(__file__),
        "case": "KRAKEN:390:38849", "question": "Which occupancy decision creates the one-cell difference and changes argmin?",
        "scope": "Arithmetic diagnosis, no production change or target-based conformer selection."})
    req = next(r for r in read_jsonl(OUT / "campaigns/recovered_precision/requests.jsonl.gz") if r["id"] == "KRAKEN:390:38849")
    req = copy.deepcopy(req)
    xyz = np.array([a["position"] for a in req["atoms"]])
    v = np.sum(xyz[req["donor"]] - xyz[req["neighbors"]], axis=0)
    req["center"] = (xyz[req["donor"]] + 2.28*v/np.linalg.norm(v)).tolist()
    plan = json.loads((OUT / "conventions/plan.json").read_text())
    radii = next(r["value"] for r in plan["constant_audit"] if r["quantity"] == "BV_radii")
    for a in req["atoms"]:
        a["radius"] = radii[a["element"]]
    req["config"]["density"] = .001
    req["dump"] = True
    req["sterimol_only"] = False
    result = native(req)
    save(out / "native_request.json", req)
    save(out / "native_result.json", result)
    for p in PROCESSES:
        p.stdin.close()
        assert p.wait() == 0
    frame = result["dump"]["orientations"][-1]
    # Native dump supplies its actual transformed f32 atoms/radii, not an
    # independently re-created basis that could differ at the boundary.
    atoms = frame["aligned_atoms"]
    c32 = np.array([a["position"] for a in atoms], dtype=np.float32)
    r232 = np.array([a["radius_squared"] for a in atoms], dtype=np.float32)
    xyz = np.array([a["position"] for a in req["atoms"]])
    donor, neighbors = req["donor"], req["neighbors"]
    v = np.sum(xyz[donor] - xyz[neighbors], axis=0)
    center = xyz[donor] + 2.28 * v / np.linalg.norm(v)
    z = REF.unit(center - xyz[donor])
    w = xyz[neighbors[-1]] - center
    x = REF.unit(w - (w @ z) * z)
    basis = np.array([x, np.cross(z, x), z])
    include = np.array([a["element"] != "H" for a in req["atoms"]])
    c64 = ((xyz - center) @ basis.T)[include]
    r264 = (np.array([a["radius"] for a in req["atoms"]])[include] * 1.17)**2
    axis64 = np.linspace(-3.5, 3.5, 70)
    axis32 = -np.float32(3.5) + np.arange(70, dtype=np.float32) * (np.float32(7)/np.float32(69))
    def grid(axis):
        p = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), -1).reshape(-1, 3)
        return p, np.sum(p*p, axis=1) <= 3.5**2
    all64, inside64 = grid(axis64)
    all32, inside32 = grid(axis32)
    assert np.array_equal(inside64, inside32)
    p64, p32 = all64[inside64], all32[inside32]
    def occupancy(points, centers, radius2):
        best = np.full(len(points), np.inf, dtype=points.dtype)
        atom_idx = np.zeros(len(points), dtype=int)
        for i, (c, r2) in enumerate(zip(centers, radius2, strict=True)):
            d = points-c
            margin = (d[:, 0]*d[:, 0] + d[:, 1]*d[:, 1]) + d[:, 2]*d[:, 2] - r2
            better = margin < best
            best[better] = margin[better]
            atom_idx[better] = i
        return best <= 0, best, atom_idx
    m64, d64, a64 = occupancy(p64, c64, r264)
    m32, d32, a32 = occupancy(p32, c32, r232)
    cells = []
    for i in np.flatnonzero(m64 != m32):
        cells.append({"grid_index": int(i), "point64": p64[i].tolist(), "point32": p32[i].tolist(),
            "occupied64": bool(m64[i]), "occupied32": bool(m32[i]),
            "minimum_squared_distance_minus_radius_squared64": float(d64[i]),
            "minimum_squared_distance_minus_radius_squared32": float(d32[i]),
            "nearest_heavy_atom64": int(a64[i]), "nearest_heavy_atom32": int(a32[i]),
            "center64": c64[a64[i]].tolist(), "radius_squared64": float(r264[a64[i]]),
            "center32": c32[a32[i]].tolist(), "radius_squared32": float(r232[a32[i]])})
    sphere = 4*np.pi*3.5**3/3
    v32 = float(m32.mean()*sphere)
    v64 = float(m64.mean()*sphere)
    assert abs(v32 - frame["buried_volume"]) < 1e-5
    save(out / "complete.json", {"grid_points": len(p64), "grid_masks_identical": True,
        "occupied32": int(m32.sum()), "occupied64": int(m64.sum()),
        "native_dump_volume": frame["buried_volume"], "independent_count_volume32": v32,
        "independent_count_volume64": v64, "cell_volume": sphere/len(p64),
        "boundary_cells": cells,
        "classification": "Finite f32 grid/atom boundary decision, amplified by discontinuous argmin selection; no incorrect union-of-balls operation established.",
        "minimal_witness": "Each differing cell and its nearest sphere is sufficient to reproduce the occupancy divergence; coordinates/radius/margins retained above."})
    print((out / "complete.json").read_text())


if __name__ == "__main__":
    main()
