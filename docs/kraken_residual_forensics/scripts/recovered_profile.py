"""Replay justified XYZ precision recovery; preserve every frozen ensemble ID."""
import concurrent.futures
import copy
import gzip
import json
import math
import time

import numpy as np
import pandas as pd

from analyze import FIELDS, read_jsonl, references, summarize, write_analysis
from freeze import OUT, save, sha
from reference_campaign import calculate as morfeus_calculate
from minimal_reference import calculate as minimal_calculate
from run_native import PROCESSES, native
from sampling_campaign import PROCS, observe as sample


def calculate(item):
    req, radii, neighbors = item
    xyz = np.array([a["position"] for a in req["atoms"]])
    v = np.sum(xyz[req["donor"]] - xyz[neighbors], axis=0)
    center = xyz[req["donor"]] + 2.28 * v / np.linalg.norm(v)
    given = copy.deepcopy(req)
    given.update(center=center.tolist(), op="geometry_topology", neighbors=neighbors, dump=True)
    given["config"]["density"] = .001
    for a in given["atoms"]:
        a["radius"] = radii[a["element"]]
    result = native(given)
    last = result["dump"]["orientations"][-1]
    for f in ("buried_volume", "near_vbur", "far_vbur"):
        result["buried_volume"][f] = last[f]
    result["buried_volume"]["percent_buried_volume"] = last["buried_volume"] * 100 / (4 * math.pi * 3.5**3 / 3)
    for a in given["atoms"]:
        if a["radius"] == 1.2:
            a["radius"] = 1.09
    given["n_rot_vectors"] = 3600
    result["sterimol_coordination"] = sample(given)
    reference = morfeus_calculate((req, {"primary_neighbors": neighbors}))
    minimal = minimal_calculate((req, {"primary_neighbors": neighbors, "primary_center": center.tolist()}, radii))
    return {"id": req["id"], "native": result, "morfeus": reference, "minimal": minimal}


def main():
    out = OUT / "campaigns/recovered_precision"
    out.mkdir(exist_ok=False)
    receipt = json.loads((OUT / "recovered_xyz/complete.json").read_text())
    assert not receipt["errors"] and not receipt["incompatible"], receipt
    changed = set(receipt["changed"])
    requests = list(read_jsonl(OUT / "frozen/requests.jsonl.gz"))
    responses = {r["id"]: r for r in read_jsonl(OUT / "recovered_xyz/responses.jsonl.gz") if r["id"] in changed}
    primary = {r["id"]: r for r in read_jsonl(OUT / "frozen/prior_morfeus_analytic.jsonl.gz")}
    plan = json.loads((OUT / "conventions/plan.json").read_text())
    radii = next(r["value"] for r in plan["constant_audit"] if r["quantity"] == "BV_radii")
    updates = []
    with gzip.open(out / "requests.jsonl.gz", "wt") as stream:
        for req in requests:
            if req["id"] in changed:
                new = copy.deepcopy(req)
                for a, p in zip(new["atoms"], responses[req["id"]]["coordinates"], strict=True):
                    a["position"] = p
                updates.append(new)
            else:
                new = req
            stream.write(json.dumps(new, separators=(",", ":")) + "\n")
    save(out / "plan.json", {"source_recovery_manifest_sha256": sha(OUT / "recovered_xyz/complete.json"),
        "script_sha256": sha(__file__), "N": len(requests), "changed_N": len(updates),
        "changed_IDs": sorted(changed), "input_sha256": sha(out / "requests.jsonl.gz"),
        "profile": "Source-convention BV including historical last plane; endpoint-inclusive native 3600-direction Sterimol; source precision used wherever API supplies it, irrespective of effect on target agreement.",
        "unchanged_reference_source_sha256": sha(OUT / "campaigns/morfeus_primary_resumed/observations.jsonl.gz")})
    start = time.monotonic()
    records = []
    with gzip.open(out / "changed_observations.jsonl.gz", "wt") as dst, concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        for r in pool.map(calculate, ((r, radii, primary[r["id"]]["primary_neighbors"]) for r in updates)):
            records.append(r)
            dst.write(json.dumps(r, separators=(",", ":")) + "\n")
    for p in PROCESSES + PROCS:
        p.stdin.close()
        assert p.wait() == 0 and not p.stderr.read()
    replacement = {r["id"]: r for r in records}
    with gzip.open(out / "morfeus_observations.jsonl.gz", "wt") as dst:
        for r in read_jsonl(OUT / "campaigns/morfeus_primary_resumed/observations.jsonl.gz"):
            dst.write(json.dumps(replacement[r["id"]]["morfeus"] if r["id"] in replacement else r, separators=(",", ":")) + "\n")
    seen = set()
    with gzip.open(out / "minimal_observations.jsonl.gz", "wt") as dst:
        for r in read_jsonl(OUT / "minimal_reference/observations.jsonl.gz"):
            seen.add(r["id"])
            dst.write(json.dumps(replacement[r["id"]]["minimal"] if r["id"] in replacement else r, separators=(",", ":")) + "\n")
        for key in sorted(set(replacement) - seen):
            dst.write(json.dumps(replacement[key]["minimal"], separators=(",", ":")) + "\n")
    base = pd.read_csv(OUT / "campaigns/primary_sampling_3600/analysis/conformers.csv.gz", float_precision="round_trip").set_index("id")
    for record in records:
        for f, (op, key, _, _) in FIELDS.items():
            base.loc[record["id"], f] = record["native"][op][key]
    base = base.reset_index()
    published, _ = references()
    lig = pd.read_csv(OUT / "baseline/ligand_provenance.csv")
    comparison, failures = summarize(base, published, lig)
    write_analysis(comparison, out / "analysis", plot=False)
    base.to_csv(out / "analysis/conformers.csv.gz", index=False)
    save(out / "analysis/failures.json", failures)
    save(out / "complete.json", {"N": len(requests), "changed_N": len(records), "elapsed_seconds": time.monotonic()-start,
        "native_conformers_sha256": sha(out / "analysis/conformers.csv.gz"),
        "morfeus_observations_sha256": sha(out / "morfeus_observations.jsonl.gz"),
        "plan_sha256": sha(out / "plan.json")})
    print(len(records), "precision-recovered conformers replayed; full cohort retained", flush=True)


if __name__ == "__main__":
    main()
