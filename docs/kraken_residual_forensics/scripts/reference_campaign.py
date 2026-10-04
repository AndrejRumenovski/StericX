"""Fresh Morfeus reference and independent analytic support on identical exports."""
from __future__ import annotations

import argparse
import concurrent.futures
import datetime
import gzip
import json
import math
import os
import time

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import numpy as np
from morfeus import BuriedVolume, Pyramidalization, Sterimol
from morfeus.utils import get_radii

from analyze import read_jsonl
from freeze import OUT, save, sha


def analytic_support(xyz, donor, center, radii):
    axis = xyz[donor] - center
    axis /= np.linalg.norm(axis)
    d = xyz - center
    axial = d @ axis
    radial = np.sqrt(np.maximum(0, np.sum(d * d, axis=1) - axial * axial))
    return {"sterimol_l": float(np.max(axial + radii) + .4),
            "sterimol_b5": float(np.max(radial + radii))}


def calculate(item):
    req, primary = item
    elems = [a["element"] for a in req["atoms"]]
    xyz = np.array([a["position"] for a in req["atoms"]])
    donor = req["donor"]
    neighbors = primary["primary_neighbors"]
    v = np.sum(xyz[donor] - xyz[neighbors], axis=0)
    center = xyz[donor] + 2.28 * v / np.linalg.norm(v)
    augmented = np.vstack((xyz, center))
    augmented_elements = elems + ["Pd"]
    metal = len(augmented)
    frames = []
    for neighbor in neighbors:
        bv = BuriedVolume(augmented_elements, augmented.copy(), metal, excluded_atoms=[metal],
            z_axis_atoms=[donor + 1], xz_plane_atoms=[neighbor + 1],
            density=.001, radius=3.5, radii_type="bondi", radii_scale=1.17, include_hs=False)
        bv.octant_analysis()
        q = list(bv.quadrants["buried_volume"].values())
        o = list(bv.octants["buried_volume"].values())
        frames.append({"plane": neighbor, "buried_volume": float(bv.buried_volume),
            "qvbur_min": float(min(q)), "qvbur_max": float(max(q)),
            "max_delta_qvbur": float(max(abs(q[j] - q[j-1]) for j in range(4))),
            "ovbur_min": float(min(o)), "ovbur_max": float(max(o)),
            "near_vbur": float(sum(o[4:])), "far_vbur": float(sum(o[:4])),
            "quadrants": list(map(float, q)), "octants": list(map(float, o))})
    volume = {k: float(min(f[k] for f in frames) if k in ("qvbur_min", "ovbur_min") else max(f[k] for f in frames))
              for k in ("qvbur_min", "qvbur_max", "max_delta_qvbur", "ovbur_min", "ovbur_max")}
    volume.update({k: frames[-1][k] for k in ("buried_volume", "near_vbur", "far_vbur")})
    volume["percent_buried_volume"] = volume["buried_volume"] * 100 / (4 * math.pi * 3.5**3 / 3)
    radii = np.array(get_radii(elems, radii_type="bondi"))
    radii[radii == 1.2] = 1.09
    ster = Sterimol(augmented_elements, augmented.copy(), metal, donor + 1,
                   radii=list(radii) + [1.63], n_rot_vectors=3600)
    pyr = Pyramidalization(xyz, donor + 1, elements=elems, neighbor_indices=[n + 1 for n in req["neighbors"]])
    return {"id": req["id"], "center": center.tolist(), "neighbors": neighbors,
        "buried_volume": volume, "frames": frames,
        "sterimol_coordination": {"l": float(ster.L_value), "b1": float(ster.B_1_value), "b5": float(ster.B_5_value)},
        "pyramidalization": {"pyr_p": float(pyr.P), "pyr_alpha": float(pyr.alpha)},
        "minimal_analytic": analytic_support(xyz, donor, center, radii)}


def protect(item):
    try:
        return calculate(item)
    except Exception as e:
        return {"id": item[0]["id"], "error": repr(e)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--name")
    parser.add_argument("--reuse-prefix", help="Interrupted gzip stream; preserve source, reuse only complete ordered JSON records")
    args = parser.parse_args()
    name = args.name or ("morfeus_primary" + (f"_probe{args.limit}" if args.limit else ""))
    out = OUT / "campaigns" / name
    out.mkdir(exist_ok=False)
    requests = list(read_jsonl(OUT / "frozen/requests.jsonl.gz"))
    refs = {r["id"]: r for r in read_jsonl(OUT / "frozen/prior_morfeus_analytic.jsonl.gz")}
    if args.limit:
        requests = requests[:args.limit]
    prefix = []
    recovery_error = None
    if args.reuse_prefix:
        try:
            for row in read_jsonl(args.reuse_prefix):
                assert row["id"] == requests[len(prefix)]["id"]
                prefix.append(row)
        except (EOFError, gzip.BadGzipFile) as error:
            recovery_error = repr(error)
    save(out / "plan.json", {"created_utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "script_sha256": sha(__file__), "N": len(requests), "workers": args.workers,
        "conventions_sha256": sha(OUT / "conventions/plan.json"),
        "reused_prefix": {"path": args.reuse_prefix, "sha256": sha(args.reuse_prefix), "complete_rows": len(prefix), "terminal_read_error": recovery_error} if args.reuse_prefix else None,
        "selection": "All frozen cohort conformers in input order (prefix only if probe)",
        "inputs": "Exactly the frozen SDF decimal coordinates promoted to float64; not original unrounded geometry",
        "independence": "Direct Morfeus calls; minimal L/B5 equations do not call native geometry code or use published targets."})
    start = time.monotonic()
    errors = []
    with gzip.open(out / "observations.jsonl.gz", "wt", compresslevel=3) as stream, concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        for row in prefix:
            stream.write(json.dumps(row, separators=(",", ":"), allow_nan=False) + "\n")
            if "error" in row:
                errors.append(row)
        for i, result in enumerate(pool.map(protect, ((r, refs[r["id"]]) for r in requests[len(prefix):]), chunksize=8), len(prefix) + 1):
            stream.write(json.dumps(result, separators=(",", ":"), allow_nan=False) + "\n")
            if "error" in result:
                errors.append(result)
            if i % 1000 == 0 or i == len(requests):
                print(f"{i}/{len(requests)}; {time.monotonic()-start:.1f}s; errors={len(errors)}", flush=True)
    save(out / "complete.json", {"N": len(requests), "errors": errors, "elapsed_seconds": time.monotonic()-start,
        "observations_sha256": sha(out / "observations.jsonl.gz"), "plan_sha256": sha(out / "plan.json")})


if __name__ == "__main__":
    main()
