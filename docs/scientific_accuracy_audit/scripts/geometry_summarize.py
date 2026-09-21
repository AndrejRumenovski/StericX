"""Full-precision per-field residuals and stratified audits; no outlier removal."""

import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1] / "geometry"
C = {
    x["id"]: x
    for x in json.loads((ROOT / "inputs.json").read_text())["cases"]
    + json.loads((ROOT / "supplement_inputs.json").read_text())["cases"]
}
S = {
    x["id"]: x
    for x in map(json.loads, (ROOT / "sut_outputs.jsonl").read_text().splitlines())
}
R = {
    x["id"]: x
    for x in map(
        json.loads, (ROOT / "reference_results.jsonl").read_text().splitlines()
    )
}


def finite(v):
    return isinstance(v, (int, float)) and math.isfinite(v)


def stats(pairs):
    y = np.array([p[0] for p in pairs])
    x = np.array([p[1] for p in pairs])
    e = y - x
    if len(e) == 0:
        return {"n": 0}
    ss = sum((x - x.mean()) ** 2)
    slope, intercept = (
        np.polyfit(x, y, 1) if len(x) > 1 and ss > 1e-20 else (None, None)
    )
    return {
        "n": len(e),
        "mae": float(abs(e).mean()),
        "rmse": float(np.sqrt((e * e).mean())),
        "max_absolute_error": float(abs(e).max()),
        "median_absolute_error": float(np.median(abs(e))),
        "mean_signed_error": float(e.mean()),
        "r2_identity": float(1 - sum(e * e) / ss) if ss > 1e-20 else None,
        "slope": float(slope) if slope is not None else None,
        "intercept": float(intercept) if intercept is not None else None,
    }


rows = []
failures = []
bins = []


def add(id, component, sut, ref, reference):
    for k, v in ref.items():
        if finite(v) and finite(sut.get(k)):
            d = sut[k] - v
            rows.append(
                {
                    "id": id,
                    "category": C[id]["category"],
                    "component": component,
                    "descriptor": k,
                    "reference": reference,
                    "sut": sut[k],
                    "reference_value": v,
                    "signed_error": d,
                    "absolute_error": abs(d),
                    "relative_error": abs(d / v) if abs(v) > 1e-15 else None,
                }
            )


for id, c in C.items():
    s = S[id]
    r = R[id]
    for key in ["bond", "coordination"]:
        sr = r.get("sterimol", {}).get(key, {})
        ss = s.get("sterimol_bond" if key == "bond" else "sterimol_dummy_raw", {})
        for ref in ["independent", "morfeus"]:
            add(id, "sterimol_" + key, ss, sr.get(ref, {}), ref)
    for ref in ["independent", "morfeus"]:
        add(
            id,
            "pyramidalization",
            s.get("pyramidalization", {}),
            r.get("pyramidalization", {}).get(ref, {}),
            ref,
        )
        add(
            id,
            "buried_volume",
            s.get("buried_volume", {}),
            r.get("buried_volume", {}).get(ref, {}),
            ref,
        )
    for k in ["buried_volume", "center", "pyramidalization"]:
        if isinstance(s.get(k), dict) and "error" in s[k]:
            failures.append(
                {
                    "id": id,
                    "category": c["category"],
                    "component": k,
                    "source": "SUT",
                    "error": s[k]["error"],
                }
            )
    for k in ["sterimol", "pyramidalization", "buried_volume"]:
        if "error" in r.get(k, {}):
            failures.append(
                {
                    "id": id,
                    "category": c["category"],
                    "component": k,
                    "source": "reference",
                    "error": r[k]["error"],
                }
            )
    obs = s.get("dump", {}).get("orientations", [])
    byplane = {o["plane"]: o for o in obs}
    for orientation in r.get("buried_volume", {}).get("orientations", []):
        o = byplane.get(orientation["plane"])
        if o is None:
            continue
        for ref in ["independent", "morfeus"]:
            for label in ["quadrants", "octants"]:
                for index, (a, b) in enumerate(
                    zip(o[label], orientation[ref][label], strict=True)
                ):
                    if finite(a) and finite(b):
                        bins.append(
                            {
                                "id": id,
                                "category": c["category"],
                                "plane": o["plane"],
                                "region": label,
                                "index": index,
                                "reference": ref,
                                "sut": a,
                                "reference_value": b,
                                "error": a - b,
                            }
                        )
for name, data in [
    ("descriptor_comparisons", rows),
    ("bin_comparisons", bins),
    ("failures", failures),
]:
    with (ROOT / (name + ".csv")).open("w") as fp:
        w = csv.DictWriter(fp, fieldnames=list(data[0]) if data else [])
        w.writeheader()
        w.writerows(data)
metrics = {}
for population, allowed in [
    ("all", set(C)),
    ("base", {k for k, v in C.items() if v["category"] == "base"}),
    (
        "molecular_smiles",
        {k for k, v in C.items() if v["category"] == "base" and "smiles" in v},
    ),
]:
    groups = defaultdict(list)
    for r in rows:
        if r["id"] in allowed:
            groups[(r["component"], r["descriptor"], r["reference"])].append(
                (r["sut"], r["reference_value"])
            )
    metrics[population] = {"/".join(k): stats(v) for k, v in groups.items()}
metrics["individual_bins"] = {}
for category in ["all", "base", "convergence"]:
    for region in ["quadrants", "octants"]:
        for ref in ["independent", "morfeus"]:
            metrics["individual_bins"][f"{category}/{region}/{ref}"] = stats(
                [
                    (r["sut"], r["reference_value"])
                    for r in bins
                    if (category == "all" or r["category"] == category)
                    and r["region"] == region
                    and r["reference"] == ref
                ]
            )
(ROOT / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
invariance = []
for category in ["rigid", "permutation", "precision"]:
    groups = defaultdict(list)
    for id, c in C.items():
        if c["category"] != category:
            continue
        for component in [
            "sterimol_bond",
            "sterimol_dummy_raw",
            "pyramidalization",
            "buried_volume",
        ]:
            s = S[id].get(component, {})
            base = S[c["parent"]].get(component, {})
            for d, a in base.items():
                if finite(a) and finite(s.get(d)):
                    groups[(c["parent"], component, d)].append(s[d] - a)
    for (parent, component, d), deltas in groups.items():
        x = np.array(deltas)
        invariance.append(
            {
                "campaign": category,
                "parent": parent,
                "component": component,
                "descriptor": d,
                "n": len(x),
                "max_absolute_deviation": float(abs(x).max()),
                "mean_signed_deviation": float(x.mean()),
                "mean_absolute_deviation": float(abs(x).mean()),
                "standard_deviation": float(x.std()),
            }
        )
with (ROOT / "invariance.csv").open("w") as fp:
    w = csv.DictWriter(fp, fieldnames=list(invariance[0]))
    w.writeheader()
    w.writerows(invariance)
# Connectivity compared to independently declared graphs, not a second radius heuristic.
connectivity = []
for id, c in C.items():
    if c["category"] not in ["base", "radii", "chemical_edge"]:
        continue
    observed = sorted(t["index"] for t in S[id]["bonded_neighbors"])
    expected = (
        None
        if c.get("sterimol_only")
        else c.get("expected_bonded_neighbors", c.get("neighbors"))
    )
    connectivity.append(
        {
            "id": id,
            "donor": c["donor"],
            "expected_neighbors": expected,
            "observed_neighbors": observed,
            "matches_graph": sorted(expected) == observed
            if expected is not None
            else None,
            "reason": c.get(
                "chemical_expectation",
                "explicit SMILES graph"
                if "smiles" in c
                else "synthetic geometric intent; chemical stability not asserted",
            ),
        }
    )
(ROOT / "connectivity.json").write_text(json.dumps(connectivity, indent=2) + "\n")
print("Rows", len(rows), "bins", len(bins), "failures", len(failures))
for k, v in metrics["molecular_smiles"].items():
    if k.endswith("independent"):
        print(k, v)
print("Worst all rows")
print(sorted(rows, key=lambda x: x["absolute_error"], reverse=True)[:8])
