"""Analytic coordinate-rounding bounds, not target-fitted sensitivity sweeps.

Conditional assumptions: nearest rounding to 4 decimals, same atom identity,
same complete ensemble, same donor/neighbors, primary raw-vector center/radii.
These cannot bound unknown historical ensemble membership or software changes.
"""
import datetime
import gzip
import json
import math

import numpy as np
import pandas as pd

from analyze import read_jsonl, references
from freeze import OUT, save, sha


def intervals(req, reference):
    xyz = np.array([a["position"] for a in req["atoms"]])
    donor = req["donor"]
    vectors = xyz[req["neighbors"]] - xyz[donor]
    lengths = np.linalg.norm(vectors, axis=1)
    unit = vectors / lengths[:, None]
    coordinate_error = math.sqrt(3) * .00005
    bond_error = 2 * coordinate_error
    du = 2 * bond_error / (lengths - bond_error)
    angles = []
    for k in range(3):
        i, j = [s for s in range(3) if s != k]
        cross = np.cross(unit[i], unit[j])
        norm = np.linalg.norm(cross)
        cross_error = du[i] + du[j]
        normal_error = min(2., 2 * cross_error / (norm - cross_error)) if norm > cross_error else 2.
        cosine = abs(float(unit[k] @ cross / norm))
        eps = normal_error + du[k]
        low = math.acos(min(1., cosine + eps))
        high = math.acos(max(0., cosine - eps))
        sign = float((unit[i] + unit[j]) @ unit[k])
        sign_error = du[i] + du[j] + 2 * du[k]
        if sign - sign_error > 0:
            angles.append((-high, -low))
        elif sign + sign_error < 0:
            angles.append((low, high))
        else:
            angles.append((-high, high))
    alpha_lower, alpha_upper = np.rad2deg(np.mean(angles, axis=0))
    det = abs(float(np.linalg.det(unit)))
    det_lower = max(0., det - du.sum())
    det_upper = min(1., det + du.sum())
    if alpha_lower >= 0:
        p_lower, p_upper = det_lower, det_upper
    elif alpha_upper < 0:
        p_lower, p_upper = 2 - det_upper, 2 - det_lower
    else:
        p_lower, p_upper = det_lower, 2 - det_lower
    raw_norm = float(np.linalg.norm(np.sum(-vectors, axis=0)))
    raw_error = 6 * coordinate_error
    d_axis = min(2., 2 * raw_error / (raw_norm - raw_error)) if raw_norm > raw_error else 2.
    max_distance = float(np.linalg.norm(xyz - xyz[donor], axis=1).max())
    support_bound = 2 * coordinate_error + max_distance * d_axis
    sampled_b1_bound = max_distance * math.pi / 3599
    ster = reference["sterimol_primary"]
    bounds = {
        "pyr_p": [p_lower, p_upper], "pyr_alpha": [float(alpha_lower), float(alpha_upper)],
        "sterimol_l": [ster["sterimol_l"] - support_bound, ster["sterimol_l"] + support_bound],
        "sterimol_b5": [ster["sterimol_b5"] - support_bound, ster["sterimol_b5"] + support_bound],
        "sterimol_b1": [ster["sterimol_b1"] - support_bound - sampled_b1_bound,
                        ster["sterimol_b1"] + support_bound + sampled_b1_bound],
    }
    # Guard far below the conservative analytic margins, not a fitted tolerance.
    bounds = {k: [v[0] - 1e-10, v[1] + 1e-10] for k, v in bounds.items()}
    return {"id": req["id"], "molecule_id": int(req["id"].split(":")[1]), "bounds": bounds,
        "unit_bond_error_bounds": du.tolist(), "raw_sum_norm_A": raw_norm,
        "raw_axis_error_bound": d_axis, "support_error_bound_A": support_bound,
        "b1_single_scan_error_bound_A": sampled_b1_bound, "alpha_sign_stable": bool(alpha_lower >= 0 or alpha_upper < 0)}


def main():
    out = OUT / "geometry_bounds"
    out.mkdir(exist_ok=False)
    save(out / "plan.json", {"created_utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "script_sha256": sha(__file__), "coordinate_halfwidth_A": .00005,
        "assumptions": "same atomic identity, same ensemble membership, same selected donor/neighbors, nearest SDF rounding, primary raw-vector center and documented radii; B1 scan phase free within 3600-point gap bound",
        "derivation": {
            "coordinate": "||delta r|| <= sqrt(3)*0.00005 = h; ||delta bond|| <= 2h; ||delta normalized bond|| <= 4h/(bond length-2h)",
            "P": "Telescoping determinant of three unit vectors gives |delta |det|| <= sum(unit-vector bounds). Use separate signed alpha interval to establish P=|det| or 2-|det| branch.",
            "alpha": "||delta cross(u_i,u_j)|| <= du_i+du_j; normalization bound 2epsilon/(cross norm-epsilon); bound absolute normal dot product by du_k+dnormal, then monotonic acos. Bisector sign bound du_i+du_j+2du_k; ambiguous sign uses both branches.",
            "support": "raw vector sum perturbation <= 6h; axis change <= 12h/(raw sum norm-6h). L/B5 envelope changes <= 2h + max donor-distance * axis change. Translation to Pd adds constant 2.28 to L and no transverse displacement.",
            "B1": "support is max-distance Lipschitz in azimuth; nearest sample angular gap <= pi/3599. A sampled value brackets continuous B1 within this one-sided bound; permit phase changes conservatively on both ends.",
            "ensemble": "min/max interval endpoints reduced monotonically; delta=[max lower-min upper,max upper-min lower]; unknown vburminconf membership uses union of conformer intervals.",
            "conditional_R2_upper_bound": "1 - sum(distance(published, admissible interval)^2)/SST; allows every output to move independently, hence an optimistic upper bound, not an achievable fitted result."},
        "arithmetic": "conservative analytic inequalities evaluated in float64 with 1e-10 outward guard; not a formal interval-arithmetic proof of the implementation"})
    refs = {r["id"]: r for r in read_jsonl(OUT / "frozen/prior_morfeus_analytic.jsonl.gz")}
    groups = {}
    with gzip.open(out / "conformer_intervals.jsonl.gz", "wt") as stream:
        for req in read_jsonl(OUT / "frozen/requests.jsonl.gz"):
            r = intervals(req, refs[req["id"]])
            stream.write(json.dumps(r, separators=(",", ":")) + "\n")
            groups.setdefault(r["molecule_id"], []).append(r)
    published, _ = references()
    names = {"pyr_p": "pyr_P", "pyr_alpha": "pyr_alpha", "sterimol_l": "sterimol_L", "sterimol_b1": "sterimol_B1", "sterimol_b5": "sterimol_B5"}
    rows = []
    for mid, g in groups.items():
        for name, prop in names.items():
            b = np.array([r["bounds"][name] for r in g])
            for red in ("min", "max", "delta", "vburminconf"):
                if red == "min":
                    low, high = b.min(axis=0)
                elif red == "max":
                    low, high = b.max(axis=0)
                elif red == "delta":
                    low, high = max(0., b[:, 0].max()-b[:, 1].min()), b[:, 1].max()-b[:, 0].min()
                else:
                    low, high = b[:, 0].min(), b[:, 1].max()
                target = published[mid][prop][red]
                rows.append({"molecule_id": mid, "metric": name + "_" + red, "kraken": target,
                    "lower": low, "upper": high, "unavoidable_abs_error_with_rounding_only": max(0., low-target, target-high),
                    "published_outside_interval": bool(target < low or target > high)})
    df = pd.DataFrame(rows)
    df.to_csv(out / "ligand_intervals.csv", index=False)
    result = []
    for name, g in df.groupby("metric"):
        sst = float(np.sum((g.kraken-g.kraken.mean())**2))
        floor = float(np.sum(g.unavoidable_abs_error_with_rounding_only**2))
        result.append({"metric": name, "N": len(g), "outside_interval": int(g.published_outside_interval.sum()),
            "SST": sst, "rounding_only_SSE_lower_bound": floor, "conditional_R2_upper_bound": 1-floor/sst,
            "RMSE_lower_bound": math.sqrt(floor/len(g))})
    pd.DataFrame(result).to_csv(out / "conditional_ceilings.csv", index=False)
    save(out / "complete.json", {"ligands": len(groups), "conformers": sum(map(len, groups.values())),
        "metrics": result, "plan_sha256": sha(out / "plan.json"), "intervals_sha256": sha(out / "ligand_intervals.csv")})
    print(pd.DataFrame(result).to_string(index=False))


if __name__ == "__main__":
    main()
