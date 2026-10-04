"""Validate optional sampling against frozen defaults and analytic envelopes."""
import concurrent.futures
import copy
import gzip
import json
import math
import os
import subprocess
import threading
import time

import numpy as np
import pandas as pd

from analyze import FIELDS, read_jsonl, references, summarize, write_analysis
from freeze import OUT, ROOT, save, sha

LOCAL = threading.local()
PROCS = []
BINARY = OUT / "adapter/target/release/stericx-kraken-sampling-observer"


def observe(req):
    if not hasattr(LOCAL, "proc"):
        LOCAL.proc = subprocess.Popen([str(BINARY)], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        PROCS.append(LOCAL.proc)
    p = LOCAL.proc
    p.stdin.write(json.dumps(req, separators=(",", ":")) + "\n")
    p.stdin.flush()
    row = json.loads(p.stdout.readline())
    assert row["id"] == req["id"] and "error" not in row, row
    return row


def calculate(item):
    req, old, profile, reference, radii, exact = item
    baseline = copy.deepcopy(req)
    baseline["center"] = old["center"]
    default = observe(baseline)
    expected = old["sterimol_coordination"]
    assert [default[k] for k in ("l", "b1", "b5")] == [expected[k] for k in ("l", "b1", "b5")]
    assert default["bits"] == [int(expected["_bits"][k], 16) for k in ("l", "b1", "b5")]
    primary = copy.deepcopy(req)
    primary["center"] = reference["primary_center"]
    for a in primary["atoms"]:
        value = radii[a["element"]]
        a["radius"] = 1.09 if value == 1.2 else value
    coarse = observe(primary)
    assert [coarse[k] for k in ("l", "b1", "b5")] == [profile["sterimol_coordination"][k] for k in ("l", "b1", "b5")]
    primary["n_rot_vectors"] = 3600
    fine = observe(primary)
    assert fine["l"] == coarse["l"] and fine["b5"] == coarse["b5"]
    result = {"id": req["id"], "default": default, "profile_360": coarse, "profile_3600": fine,
              "default_matches_frozen_bits": True, "L_B5_unchanged": True}
    if exact is not None:
        xyz = np.array([a["position"] for a in primary["atoms"]])
        xyz32 = xyz.astype(np.float32).astype(float)
        center = np.array(primary["center"])
        center32 = center.astype(np.float32).astype(float)
        u = xyz[req["donor"]] - center
        u32 = xyz32[req["donor"]] - center32
        u /= np.linalg.norm(u)
        u32 /= np.linalg.norm(u32)
        rs = np.array([a["radius"] for a in primary["atoms"]])
        max_distance = float(np.linalg.norm(xyz-center, axis=1).max())
        input_bound = float(np.linalg.norm(xyz32-xyz, axis=1).max() + np.linalg.norm(center32-center)
            + max_distance * np.linalg.norm(u32-u) + np.max(abs(rs.astype(np.float32).astype(float)-rs)))
        arithmetic_bound = float(64 * np.finfo(np.float32).eps * (max_distance + rs.max() + 1))
        scan_bound = max_distance * math.pi / 3599
        delta = fine["b1"] - exact["sterimol_coordination"]["b1"]
        result["independent_B1"] = {"exact": exact["sterimol_coordination"]["b1"],
            "coarse_error": coarse["b1"]-exact["sterimol_coordination"]["b1"], "fine_error": delta,
            "input_rounding_bound": input_bound, "arithmetic_guard": arithmetic_bound,
            "analytic_scan_bound": scan_bound,
            "within_analytic_bound": bool(-(input_bound+arithmetic_bound) <= delta <= scan_bound+input_bound+arithmetic_bound)}
    return result


def main():
    out = OUT / "campaigns/primary_sampling_3600"
    out.mkdir(exist_ok=False)
    assert (OUT / "minimal_reference/complete.json").exists()
    save(out / "plan.json", {"candidate_plan_sha256": sha(OUT / "candidates/angular_sampling/plan_before_edit.json"),
        "script_sha256": sha(__file__), "binary_sha256": sha(BINARY),
        "production_source_sha256": sha(ROOT / "src/geometry/sterimol.rs"),
        "fixed_n_rot_vectors": 3600, "N_conformers": 31611,
        "acceptance_basis": "Default bit equality and unchanged L/B5 across whole corpus; fixed mathematical scan bound on all 5301 adverse analytic-reference cases; no R2 criterion"})
    old = {r["id"]: r for r in read_jsonl(OUT / "frozen/current_predictions.jsonl.gz")}
    profile = {r["id"]: r for r in read_jsonl(OUT / "campaigns/primary_grid/observations.jsonl.gz")}
    reference = {r["id"]: r for r in read_jsonl(OUT / "frozen/prior_morfeus_analytic.jsonl.gz")}
    exact = {r["id"]: r for r in read_jsonl(OUT / "minimal_reference/observations.jsonl.gz")}
    plan = json.loads((OUT / "conventions/plan.json").read_text())
    radii = next(r["value"] for r in plan["constant_audit"] if r["quantity"] == "BV_radii")
    requests = list(read_jsonl(OUT / "frozen/requests.jsonl.gz"))
    start = time.monotonic()
    rows = []
    with gzip.open(out / "observations.jsonl.gz", "wt", compresslevel=3) as dst, concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
        for i, r in enumerate(pool.map(calculate, ((r, old[r["id"]], profile[r["id"]], reference[r["id"]], radii, exact.get(r["id"])) for r in requests)), 1):
            rows.append(r)
            dst.write(json.dumps(r, separators=(",", ":")) + "\n")
            if i % 2000 == 0:
                print(i, "of", len(requests), f"{time.monotonic()-start:.1f}s", flush=True)
    for p in PROCS:
        p.stdin.close()
        assert p.wait() == 0 and not p.stderr.read()
    checks = [r["independent_B1"] for r in rows if "independent_B1" in r]
    save(out / "acceptance_checks.json", {"N": len(rows), "default_bit_mismatches": 0, "L_B5_changes": 0,
        "independent_B1_cases": len(checks), "bound_failures": [r for r in rows if "independent_B1" in r and not r["independent_B1"]["within_analytic_bound"]],
        "maximum_coarse_B1_error": max(abs(r["coarse_error"]) for r in checks),
        "maximum_fine_B1_error": max(abs(r["fine_error"]) for r in checks),
        "observations_sha256": sha(out / "observations.jsonl.gz")})
    base = pd.read_csv(OUT / "campaigns/primary_grid/historical_last_plane/conformers.csv.gz", float_precision="round_trip")
    updates = pd.DataFrame([{"id": r["id"], **{"sterimol_" + k: r["profile_3600"][k] for k in ("l", "b1", "b5")}} for r in rows])
    conformers = base.drop(columns=["sterimol_l", "sterimol_b1", "sterimol_b5"]).merge(updates, on="id", validate="one_to_one")
    published, _ = references()
    lig = pd.read_csv(OUT / "baseline/ligand_provenance.csv")
    comp, failures = summarize(conformers, published, lig)
    result = out / "analysis"
    met = write_analysis(comp, result, plot=False)
    conformers.to_csv(result / "conformers.csv.gz", index=False)
    save(result / "failures.json", failures)
    save(out / "complete.json", {"N": len(rows), "elapsed_seconds": time.monotonic()-start,
        "plan_sha256": sha(out / "plan.json"), "acceptance_checks_sha256": sha(out / "acceptance_checks.json")})


if __name__ == "__main__":
    main()
