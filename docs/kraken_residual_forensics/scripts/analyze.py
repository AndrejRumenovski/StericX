"""Complete-cohort residual analysis; diagnostic regressions never alter outputs."""
from __future__ import annotations

import argparse
import gzip
import json
import math
import tarfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from freeze import OUT, sha, save

FIELDS = {
    "buried_volume": ("buried_volume", "buried_volume", "vbur_vbur", "angstrom^3"),
    "qvbur_min": ("buried_volume", "qvbur_min", "vbur_qvbur_min", "angstrom^3"),
    "qvbur_max": ("buried_volume", "qvbur_max", "vbur_qvbur_max", "angstrom^3"),
    "max_delta_qvbur": ("buried_volume", "max_delta_qvbur", "vbur_max_delta_qvbur", "angstrom^3"),
    "ovbur_min": ("buried_volume", "ovbur_min", "vbur_ovbur_min", "angstrom^3"),
    "ovbur_max": ("buried_volume", "ovbur_max", "vbur_ovbur_max", "angstrom^3"),
    "near_vbur": ("buried_volume", "near_vbur", "vbur_near_vbur", "angstrom^3"),
    "far_vbur": ("buried_volume", "far_vbur", "vbur_far_vbur", "angstrom^3"),
    "sterimol_l": ("sterimol_coordination", "l", "sterimol_L", "angstrom"),
    "sterimol_b1": ("sterimol_coordination", "b1", "sterimol_B1", "angstrom"),
    "sterimol_b5": ("sterimol_coordination", "b5", "sterimol_B5", "angstrom"),
    "pyr_p": ("pyramidalization", "pyr_p", "pyr_P", "dimensionless"),
    "pyr_alpha": ("pyramidalization", "pyr_alpha", "pyr_alpha", "degrees"),
    "percent_buried_volume": ("buried_volume", "percent_buried_volume", "vbur_vbur", "percentage points (derived reference)"),
}


def read_jsonl(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as stream:
        yield from map(json.loads, stream)


def metrics(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    assert np.isfinite(x).all() and np.isfinite(y).all()
    d = y - x
    sse = float(d @ d)
    sst = float(np.sum((x - x.mean()) ** 2))
    slope = float(np.dot(x - x.mean(), y - y.mean()) / sst) if sst else None
    return {"N": len(x), "R2": 1 - sse / sst if sst else None,
            "pearson_r": float(np.corrcoef(x, y)[0, 1]) if len(x) > 1 and np.ptp(x) and np.ptp(y) else None,
            "MAE": float(np.abs(d).mean()), "RMSE": float(np.sqrt(np.mean(d ** 2))),
            "median_absolute_error": float(np.median(np.abs(d))),
            "maximum_absolute_error": float(np.abs(d).max()),
            "slope": slope, "intercept": float(y.mean() - slope * x.mean()) if slope is not None else None,
            "mean_signed_error": float(d.mean()), "SSE": sse, "SST": sst,
            "target_0_99999_SSE_budget": sst * 1e-5}


def references():
    refs, molecules = {}, {}
    with tarfile.open(OUT / "frozen/inputs_and_reference.tar.gz") as archive:
        for member in archive:
            if member.name.startswith("primary/published_dft/"):
                mid = int(Path(member.name).stem)
                refs[mid] = {r["property"]: r for r in json.load(archive.extractfile(member))}
            elif member.name.startswith("primary/molecule/"):
                mid = int(Path(member.name).stem)
                molecules[mid] = json.load(archive.extractfile(member))
    return refs, molecules


def features(molecules):
    from rdkit import Chem
    pt = Chem.GetPeriodicTable()
    rows = []
    inventory = {r["id"]: r for r in read_jsonl(OUT / "frozen/inventory.jsonl.gz")}
    for req in read_jsonl(OUT / "frozen/requests.jsonl.gz"):
        inv = inventory[req["id"]]
        xyz = np.array([a["position"] for a in req["atoms"]])
        elements = [a["element"] for a in req["atoms"]]
        donor = xyz[req["donor"]]
        radii = [pt.GetRvdw(e) for e in elements]
        rows.append({"id": req["id"], "molecule_id": inv["molecule_id"],
                     "conformer_id": inv["conformer_id"], "atom_count": len(xyz),
                     "heavy_atom_count": sum(e != "H" for e in elements),
                     "P_H_count": inv["P_H_count"], "donor_environment": inv["topology_category"],
                     "molecular_weight": sum(pt.GetAtomicWeight(e) for e in elements),
                     "donor_max_distance_A": float(np.linalg.norm(xyz - donor, axis=1).max()),
                     "steric_bulk_radius_A": float(np.max(np.linalg.norm(xyz - donor, axis=1) + radii))})
    c = pd.DataFrame(rows)
    lig = c.groupby("molecule_id").agg(atom_count=("atom_count", "max"),
        heavy_atom_count=("heavy_atom_count", "max"), P_H_count=("P_H_count", "max"),
        donor_environment=("donor_environment", "first"), molecular_weight=("molecular_weight", "max"),
        donor_max_distance_A=("donor_max_distance_A", "max"), steric_bulk_radius_A=("steric_bulk_radius_A", "max"),
        conformer_count=("conformer_id", "count")).reset_index()
    classes = []
    for row in lig.to_dict("records"):
        mid = row["molecule_id"]
        mol = Chem.MolFromSmiles(molecules[mid]["smiles"])
        donors = [a for a in mol.GetAtoms() if a.GetSymbol() == "P"] if mol else []
        aryl = sum(n.GetIsAromatic() for a in donors for n in a.GetNeighbors())
        classes.append(f"{row['donor_environment']};aryl_neighbors={aryl}")
    lig["ligand_family"] = classes
    lig["provenance_class"] = "F"
    lig["provenance_detail"] = "Public DFT SDF export; historical full-precision geometry/ensemble identity unverified"
    lig["exact_available_input_for_algorithm_comparison"] = True
    lig["exact_historical_input_verified"] = False
    lig["historical_conformer_energy"] = None
    lig["historical_weights"] = None
    lig["historical_coordinate_RMSD_A"] = None
    lig["historical_provenance_status"] = "UNRESOLVABLE FROM AVAILABLE DATA"
    return c, lig


def summarize(conformers, refs, lig):
    rows, failures = [], []
    for mid, g in conformers.groupby("molecule_id", sort=True):
        count = int(lig.set_index("molecule_id").loc[mid, "conformer_count"])
        assert len(g) == count
        for descriptor, (_, _, prop, unit) in FIELDS.items():
            values = g[descriptor].to_numpy(float)
            for red in ("min", "max", "delta", "vburminconf"):
                published = refs[mid].get(prop, {}).get(red)
                if published is None or not np.isfinite(values).all() or (red == "vburminconf" and not np.isfinite(g.buried_volume).all()):
                    failures.append({"molecule_id": int(mid), "descriptor": descriptor, "reduction": red,
                                     "reason": "missing reference or incomplete calculated ensemble"})
                    continue
                imin, imax = int(np.argmin(values)), int(np.argmax(values))
                idx = {"min": imin, "max": imax, "vburminconf": int(np.argmin(g.buried_volume))}.get(red)
                value = float(values.max() - values.min()) if red == "delta" else float(values[idx])
                ref = float(published)
                if descriptor == "percent_buried_volume":
                    ref *= 100 / (4 * math.pi * 3.5 ** 3 / 3)
                rows.append({"molecule_id": int(mid), "metric": descriptor + "_" + red,
                    "descriptor": descriptor, "reduction": red, "units": unit,
                    "kraken": ref, "stericx": value, "residual": value - ref,
                    "absolute_residual": abs(value - ref),
                    "min_conformer_id": int(g.iloc[imin].conformer_id), "max_conformer_id": int(g.iloc[imax].conformer_id),
                    "selected_conformer_id": int(g.iloc[idx].conformer_id) if idx is not None else None})
    return pd.DataFrame(rows).merge(lig, on="molecule_id", validate="many_to_one"), failures


def write_analysis(comparison, out, plot=True):
    out.mkdir(exist_ok=False)
    (out / "rankings").mkdir()
    (out / "plots").mkdir()
    comparison.to_csv(out / "all_residuals.csv.gz", index=False)
    metric_rows, top_rows, strata, associations = [], [], [], []
    for name, g in comparison.groupby("metric", sort=True):
        m = metrics(g.kraken, g.stericx)
        metric_rows.append({"metric": name, "units": g.iloc[0].units, **m})
        ranked = g.sort_values(["absolute_residual", "molecule_id"], ascending=[False, True]).copy()
        ranked.insert(0, "rank", np.arange(1, len(g) + 1))
        ranked.to_csv(out / "rankings" / f"{name}.csv", index=False)
        for n in (10, 25, 50, 100):
            top = ranked.head(n)
            sse = float((top.residual ** 2).sum())
            top_rows.append({"metric": name, "top_n": n, "SSE": sse, "fraction_total_SSE": sse / m["SSE"] if m["SSE"] else 0,
                             "ligand_ids": top.molecule_id.tolist()})
        for column in ("ligand_family", "P_H_count", "donor_environment", "provenance_class"):
            for value, group in g.groupby(column):
                strata.append({"metric": name, "stratifier": column, "class": str(value), **metrics(group.kraken, group.stericx)})
        features = ("kraken", "molecular_weight", "atom_count", "heavy_atom_count", "P_H_count", "conformer_count", "steric_bulk_radius_A", "donor_max_distance_A")
        for f in features:
            for response in ("residual", "absolute_residual"):
                x, y = g[f], g[response]
                associations.append({"metric": name, "feature": f, "response": response,
                    "pearson_r": float(x.corr(y)) if x.nunique() > 1 and y.nunique() > 1 else None,
                    "spearman_r": float(x.corr(y, method="spearman")) if x.nunique() > 1 and y.nunique() > 1 else None})
        if plot:
            fig, axes = plt.subplots(3, 4, figsize=(17, 11))
            axes = axes.ravel()
            axes[0].scatter(g.kraken, g.stericx, s=5, alpha=.4)
            span = [min(g.kraken.min(), g.stericx.min()), max(g.kraken.max(), g.stericx.max())]
            axes[0].plot(span, span, "k--", lw=1)
            axes[0].set(xlabel="Kraken", ylabel="StericX", title=f"N={len(g)}; R²={m['R2']:.8f}")
            axes[1].scatter(g.kraken, g.residual, s=5, alpha=.4)
            axes[1].axhline(0, c="k", lw=.5)
            axes[1].set(xlabel="Kraken descriptor magnitude", ylabel="Residual (StericX − Kraken)")
            axes[2].hist(g.residual, bins=60)
            axes[2].set(xlabel="Residual", ylabel="Ligands")
            family = g.groupby("donor_environment").absolute_residual.mean().sort_values()
            axes[3].barh(family.index, family.values)
            axes[3].tick_params(axis="y", labelsize=6)
            axes[3].set(xlabel="Mean absolute error", title="Donor environment / family")
            for a, f in zip(axes[4:11], features[1:]):
                a.scatter(g[f], g.residual, s=5, alpha=.4)
                a.axhline(0, c="k", lw=.5)
                a.set(xlabel=f, ylabel="Residual")
            cats = sorted(g.ligand_family.unique())
            axes[11].scatter(g.ligand_family.map({c: i for i, c in enumerate(cats)}), g.absolute_residual, s=5, alpha=.4)
            axes[11].set(xlabel="Ligand family index (labels in CSV)", ylabel="Absolute residual")
            fig.suptitle(name + " — " + g.iloc[0].units)
            fig.tight_layout()
            fig.savefig(out / "plots" / f"{name}.png", dpi=110)
            plt.close(fig)
    pd.DataFrame(metric_rows).to_csv(out / "metrics.csv", index=False)
    pd.DataFrame(strata).to_csv(out / "stratified_metrics.csv", index=False)
    pd.DataFrame(associations).to_csv(out / "residual_associations.csv", index=False)
    save(out / "top_10_25_50_100.json", top_rows)
    save(out / "metrics.json", metric_rows)
    return metric_rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()
    frozen = OUT / "frozen"
    manifest = json.loads((frozen / "manifest.json").read_text())
    for name, data in manifest["files"].items():
        assert sha(frozen / name) == data["sha256"], name
    refs, molecules = references()
    cf, lig = features(molecules)
    rows = []
    for r in read_jsonl(frozen / "current_predictions.jsonl.gz"):
        row = {"id": r["id"]}
        for f, (op, key, _, _) in FIELDS.items():
            row[f] = r[op].get(key)
        rows.append(row)
    conformers = cf.merge(pd.DataFrame(rows), on="id", validate="one_to_one")
    comparison, failures = summarize(conformers, refs, lig)
    out = OUT / "baseline"
    result = write_analysis(comparison, out, plot=not args.no_plots)
    conformers.to_csv(out / "conformers.csv.gz", index=False)
    lig.to_csv(out / "ligand_provenance.csv", index=False)
    save(out / "failures.json", failures)
    save(out / "coverage.json", {"ligands": len(lig), "conformers": len(conformers), "metrics": len(result),
        "minimum_metric_N": min(r["N"] for r in result), "failures": len(failures),
        "historical_exact_geometry_subset_N": 0, "available_exact_export_algorithm_subset_N": len(lig),
        "provenance_class_note": "F means exact historical input is unverified, not that the export is known to be a different geometry. A/B/C/D/E require affirmative evidence."})
    save(out / "manifest.json", {"frozen_manifest_sha256": sha(frozen / "manifest.json"),
        "script_sha256": sha(Path(__file__)), "files": {str(p.relative_to(out)): sha(p) for p in sorted(out.rglob("*")) if p.is_file()},
        "R2_definition": "1-SSE/SST; Pearson correlation and diagnostic slope/intercept reported separately; no calibration."})
    print(json.dumps({"coverage": json.loads((out / "coverage.json").read_text()),
                      "headline": next(m for m in result if m["metric"] == "max_delta_qvbur_min")}))


if __name__ == "__main__":
    main()
