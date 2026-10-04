"""Independent exact support-envelope/union-of-balls reference on adverse cases.

The equations were sealed in the prior scientific audit before this experiment.
They do not call StericX. B1 enumerates analytic candidate angles instead of
sampling. The finite-grid volume reference has explicit strict-open octants.
"""
import concurrent.futures
import datetime
import functools
import gzip
import json
import os
import sys
import time

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import numpy as np
import pandas as pd

from analyze import read_jsonl
from freeze import OUT, ROOT, save, sha

sys.path.insert(0, str(ROOT / "docs/scientific_remediation/geometry"))
from sealed_reference import load_reference

REF = load_reference()


@functools.lru_cache(maxsize=1)
def grid():
    return REF.grid(3.5, .001)


def calculate(item):
    req, primary, radii = item
    xyz = np.array([a["position"] for a in req["atoms"]])
    center = np.array(primary["primary_center"])
    rs = np.array([radii[a["element"]] for a in req["atoms"]])
    ster_rs = rs.copy()
    ster_rs[ster_rs == 1.2] = 1.09
    ster = REF.sterimol_exact(xyz, ster_rs, center, req["donor"])
    ster["l"] += .4
    pyr = REF.pyramidalization_ref(req, xyz)["independent"]
    include = np.array([a["element"] != "H" for a in req["atoms"]])
    frames = []
    z = REF.unit(center - xyz[req["donor"]])
    for n in primary["primary_neighbors"]:
        plane = xyz[n] - center
        x = REF.unit(plane - (plane @ z) * z)
        y = np.cross(z, x)
        basis = np.array([x, y, z])
        aligned = (xyz - center) @ basis.T
        volumes = REF.volume_on_grid(grid(), aligned[include], rs[include] * 1.17, 3.5)
        frames.append({"plane": n, "basis": basis.tolist(), **volumes})
    return {"id": req["id"], "center": center.tolist(), "frames": frames,
            "sterimol_coordination": ster, "pyramidalization": pyr}


def protect(item):
    try:
        return calculate(item)
    except Exception as error:
        return {"id": item[0]["id"], "error": repr(error)}


def main():
    out = OUT / "minimal_reference"
    out.mkdir(exist_ok=False)
    selection = set()
    reasons = {}
    for stage in (OUT / "baseline", OUT / "campaigns/primary_grid/historical_last_plane"):
        for p in (stage / "rankings").glob("*.csv"):
            rows = pd.read_csv(p).head(100)
            for r in rows.to_dict("records"):
                for c in ("min_conformer_id", "max_conformer_id", "selected_conformer_id"):
                    if pd.notna(r[c]):
                        key = f"KRAKEN:{r['molecule_id']}:{int(r[c])}"
                        selection.add(key)
                        reasons.setdefault(key, []).append(f"{stage.relative_to(OUT)}:{p.stem}:{c}")
    requests = [r for r in read_jsonl(OUT / "frozen/requests.jsonl.gz") if r["id"] in selection]
    assert len(requests) == len(selection)
    primary = {r["id"]: r for r in read_jsonl(OUT / "frozen/prior_morfeus_analytic.jsonl.gz")}
    plan = json.loads((OUT / "conventions/plan.json").read_text())
    radii = next(v["value"] for v in plan["constant_audit"] if v["quantity"] == "BV_radii")
    save(out / "plan.json", {"created_utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "script_sha256": sha(__file__), "reference_source_sha256": sha(REF.__file__), "N": len(requests),
        "selection": "Every min/max/selected conformer for top100 rows of every one of 56 metrics, at baseline and source-convention profile",
        "reasons": reasons, "convention_plan_sha256": sha(OUT / "conventions/plan.json"),
        "independence": "Sealed pre-remediation mathematical reference; B1 exact support candidate enumeration; independent float64 grid and union-of-balls; no SUT calls"})
    start = time.monotonic()
    errors = []
    with gzip.open(out / "observations.jsonl.gz", "wt", compresslevel=3) as dst, concurrent.futures.ProcessPoolExecutor(max_workers=3) as pool:
        for i, r in enumerate(pool.map(protect, ((r, primary[r["id"]], radii) for r in requests), chunksize=8), 1):
            dst.write(json.dumps(r, separators=(",", ":"), allow_nan=False) + "\n")
            if "error" in r:
                errors.append(r)
            if i % 100 == 0:
                print(i, "of", len(requests), f"{time.monotonic()-start:.1f}s", "errors", len(errors), flush=True)
    save(out / "complete.json", {"N": len(requests), "elapsed_seconds": time.monotonic()-start, "errors": errors,
        "observations_sha256": sha(out / "observations.jsonl.gz"), "plan_sha256": sha(out / "plan.json")})


if __name__ == "__main__":
    main()
