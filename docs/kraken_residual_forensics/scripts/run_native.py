"""Run fixed, independently sourced input profiles through frozen StericX APIs.

No native source is modified. Every ensemble member is retained. The historical
last-plane view is reported separately from the native permutation-symmetric mean.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import copy
import datetime
import gzip
import json
import os
import subprocess
import threading
import time

import numpy as np
import pandas as pd

from analyze import FIELDS, read_jsonl, references, summarize, write_analysis
from freeze import OUT, save, sha

LOCAL = threading.local()
PROCESSES = []


def native(req):
    if not hasattr(LOCAL, "proc"):
        LOCAL.proc = subprocess.Popen([str(OUT / "frozen/current_observer")], stdin=subprocess.PIPE,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
            env={**os.environ, "RAYON_NUM_THREADS": "1"})
        PROCESSES.append(LOCAL.proc)
    p = LOCAL.proc
    p.stdin.write(json.dumps(req, separators=(",", ":")) + "\n")
    p.stdin.flush()
    line = p.stdout.readline()
    if not line:
        raise RuntimeError("native observer terminated: " + p.stderr.read())
    r = json.loads(line)
    assert r["id"] == req["id"]
    return r


def run(item):
    source, primary, stage, radii = item
    req = copy.deepcopy(source)
    req["op"] = "geometry_topology"
    if stage != "baseline_replay":
        xyz = np.array([a["position"] for a in req["atoms"]])
        neighbors = primary["primary_neighbors"]
        assert set(neighbors) == set(req["neighbors"]), req["id"]
        v = np.zeros(3)
        for n in neighbors:
            v += xyz[req["donor"]] - xyz[n]
        center = xyz[req["donor"]] + 2.28 * v / np.linalg.norm(v)
        assert np.allclose(center, primary["primary_center"], atol=1e-12, rtol=0)
        req["center"] = center.tolist()
        req["neighbors"] = neighbors
    if stage in ("primary_radii", "primary_grid"):
        for a in req["atoms"]:
            a["radius"] = radii[a["element"]]
    if stage == "primary_grid":
        req["config"]["density"] = .001
        req["dump"] = True
    result = native(req)
    if stage in ("primary_radii", "primary_grid"):
        ster_req = copy.deepcopy(req)
        ster_req["sterimol_only"] = True
        ster_req["dump"] = False
        for a in ster_req["atoms"]:
            if a["radius"] == 1.2:
                a["radius"] = 1.09
        ster_result = native(ster_req)
        result["sterimol_coordination"] = ster_result["sterimol_coordination"]
    compact = {"id": result["id"], "center": result["center"],
               "neighbors": req["neighbors"], "config": req["config"],
               "buried_volume": result["buried_volume"], "sterimol_coordination": result["sterimol_coordination"],
               "pyramidalization": result["pyramidalization"]}
    if "dump" in result:
        compact["frames"] = {k: v for k, v in result["dump"].items() if k != "orientations"}
        compact["frames"]["orientations"] = [{k: v for k, v in frame.items() if k != "aligned_atoms"}
                                               for frame in result["dump"]["orientations"]]
    return compact


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=["baseline_replay", "raw_center", "primary_radii", "primary_grid"], required=True)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()
    plan = json.loads((OUT / "conventions/plan.json").read_text())
    radii = next(r["value"] for r in plan["constant_audit"] if r["quantity"] == "BV_radii")
    name = args.stage + (f"_probe{args.limit}" if args.limit else "")
    out = OUT / "campaigns" / name
    out.mkdir(parents=True, exist_ok=False)
    frozen = json.loads((OUT / "frozen/manifest.json").read_text())
    assert sha(OUT / "frozen/current_observer") == frozen["files"]["current_observer"]["sha256"]
    save(out / "plan.json", {"created_utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "stage": args.stage, "workers": args.workers, "limit": args.limit,
        "convention_plan_sha256": sha(OUT / "conventions/plan.json"),
        "frozen_manifest_sha256": sha(OUT / "frozen/manifest.json"), "script_sha256": sha(__file__),
        "native_source_changes": False, "sampling_limit": "Native B1 uses 360 directions; retained and measured, not silently labeled 3600."})
    refs = {r["id"]: r for r in read_jsonl(OUT / "frozen/prior_morfeus_analytic.jsonl.gz")}
    requests = list(read_jsonl(OUT / "frozen/requests.jsonl.gz"))
    if args.limit:
        requests = requests[:args.limit]
    started = time.monotonic()
    results = []
    with gzip.open(out / "observations.jsonl.gz", "wt", compresslevel=3) as dst, concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        for i, r in enumerate(pool.map(run, ((r, refs[r["id"]], args.stage, radii) for r in requests)), 1):
            dst.write(json.dumps(r, separators=(",", ":"), allow_nan=False) + "\n")
            results.append(r)
            if i % 1000 == 0 or i == len(requests):
                print(f"{name}: {i}/{len(requests)}; {time.monotonic()-started:.1f}s", flush=True)
    for p in PROCESSES:
        p.stdin.close()
        assert p.wait() == 0
        assert not p.stderr.read()
    save(out / "complete.json", {"N": len(results), "elapsed_seconds": time.monotonic()-started,
        "observations_sha256": sha(out / "observations.jsonl.gz"), "plan_sha256": sha(out / "plan.json")})
    if args.limit:
        return
    base = pd.read_csv(OUT / "baseline/conformers.csv.gz", float_precision="round_trip")
    lig = pd.read_csv(OUT / "baseline/ligand_provenance.csv", float_precision="round_trip")
    published, _ = references()
    views = ["native_mean", "historical_last_plane"] if args.stage == "primary_grid" else ["native_mean"]
    for view in views:
        rows = []
        for r in results:
            row = {"id": r["id"]}
            for f, (op, key, _, _) in FIELDS.items():
                row[f] = r[op].get(key)
            if view == "historical_last_plane":
                last = r["frames"]["orientations"][-1]
                for f in ("buried_volume", "near_vbur", "far_vbur"):
                    row[f] = last[f]
                row["percent_buried_volume"] = last["buried_volume"] * 100 / (4 * np.pi * 3.5 ** 3 / 3)
            rows.append(row)
        conformers = base.drop(columns=list(FIELDS)).merge(pd.DataFrame(rows), on="id", validate="one_to_one")
        comparison, errors = summarize(conformers, published, lig)
        analysis = out / view
        met = write_analysis(comparison, analysis, plot=False)
        conformers.to_csv(analysis / "conformers.csv.gz", index=False)
        save(analysis / "failures.json", errors)
        if args.stage == "baseline_replay":
            values = conformers[list(FIELDS)].to_numpy()
            original = base[list(FIELDS)].to_numpy()
            assert np.array_equal(values, original)
            save(out / "frozen_output_equality.json", {"rows": len(base), "fields": len(FIELDS), "bitwise_JSON_numeric_equality": True})
        print(view, json.dumps(next(m for m in met if m["metric"] == "max_delta_qvbur_min")), flush=True)


if __name__ == "__main__":
    main()
