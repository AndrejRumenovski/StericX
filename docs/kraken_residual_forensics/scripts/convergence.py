"""Predeclared density ladder on complete representative/adverse ensembles."""
import concurrent.futures
import copy
import datetime
import gzip
import json
import time

import numpy as np
import pandas as pd

from analyze import read_jsonl
from freeze import OUT, save, sha
from run_native import PROCESSES, native


def calculate(item):
    req, reference, radii, density = item
    req = copy.deepcopy(req)
    req.update(op="geometry_topology", center=reference["primary_center"],
               neighbors=reference["primary_neighbors"], dump=True)
    req["config"]["density"] = density
    for atom in req["atoms"]:
        atom["radius"] = radii[atom["element"]]
    r = native(req)
    if "error" in r.get("buried_volume", {}):
        return {"id": req["id"], "density": density, "error": r["buried_volume"]["error"]}
    frames = [{k: v for k, v in f.items() if k != "aligned_atoms"} for f in r["dump"]["orientations"]]
    return {"id": req["id"], "density": density, "grid_count": r["dump"]["grid_count"],
            "native_mean": r["buried_volume"], "frames": frames,
            "historical_last": {**r["buried_volume"], **{k: frames[-1][k] for k in ("buried_volume", "near_vbur", "far_vbur")}}}


def main():
    out = OUT / "convergence"
    out.mkdir(exist_ok=False)
    ligands = pd.read_csv(OUT / "baseline/ligand_provenance.csv")
    baseline = pd.read_csv(OUT / "baseline/rankings/max_delta_qvbur_min.csv")
    final = pd.read_csv(OUT / "campaigns/primary_grid/historical_last_plane/rankings/max_delta_qvbur_min.csv")
    ids = set(baseline.molecule_id.head(3)) | set(ligands.groupby("donor_environment").molecule_id.min()) | set(final.molecule_id.head(10))
    requests = [r for r in read_jsonl(OUT / "frozen/requests.jsonl.gz") if int(r["id"].split(":")[1]) in ids]
    refs = {r["id"]: r for r in read_jsonl(OUT / "frozen/prior_morfeus_analytic.jsonl.gz")}
    plan = json.loads((OUT / "conventions/plan.json").read_text())
    densities = plan["convergence_densities_A3"]
    radii = next(v["value"] for v in plan["constant_audit"] if v["quantity"] == "BV_radii")
    save(out / "plan.json", {"created_utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "source_plan_sha256": sha(OUT / "conventions/plan.json"), "script_sha256": sha(__file__),
        "selection": "all conformers of baseline headline top3, minimum-ID ligand in every donor environment, and final-profile headline top10",
        "ligand_ids": sorted(map(int, ids)), "request_ids": [r["id"] for r in requests],
        "densities": densities, "N_ligands": len(ids), "N_conformers": len(requests),
        "interpretation": "finest density estimate and successive-grid empirical envelope; not a rigorous convergence bound; no density selected by agreement with targets"})
    start = time.monotonic()
    with gzip.open(out / "observations.jsonl.gz", "wt", compresslevel=3) as dst, concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        tasks = ((r, refs[r["id"]], radii, d) for d in densities for r in requests)
        for i, r in enumerate(pool.map(calculate, tasks), 1):
            dst.write(json.dumps(r, separators=(",", ":")) + "\n")
            if i % 100 == 0:
                print(i, "of", len(requests) * len(densities), f"{time.monotonic()-start:.1f}s", flush=True)
    for p in PROCESSES:
        p.stdin.close()
        assert p.wait() == 0
        assert not p.stderr.read()
    save(out / "complete.json", {"N": len(requests) * len(densities), "elapsed_seconds": time.monotonic()-start,
        "observations_sha256": sha(out / "observations.jsonl.gz"), "plan_sha256": sha(out / "plan.json")})


if __name__ == "__main__":
    main()
