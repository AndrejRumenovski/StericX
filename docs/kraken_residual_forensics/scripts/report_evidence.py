"""Tabulate final evidence without changing any predictions or cohort members."""
import json

import numpy as np
import pandas as pd

from analyze import metrics
from freeze import OUT, save, sha


def csv(name):
    return pd.read_csv(OUT / name, float_precision="round_trip")


def main():
    out = OUT / "report_evidence"
    out.mkdir(exist_ok=False)
    dossiers = csv("results/all_residual_dossiers.csv.gz")
    final = csv("results/final_native/metrics.csv")
    rows = []
    for r in final.to_dict("records"):
        rows.append({"provenance_class": "F", **r})
        for label in "ABCDE":
            rows.append({"provenance_class": label, "metric": r["metric"], "N": 0,
                         "R2": None, "status": "No verified members; not scored"})
    pd.DataFrame(rows).to_csv(out / "provenance_class_metrics.csv", index=False)
    # Retain all 84 members. ID-block summaries diagnose provenance, not a
    # proposed filtered dataset, and no filtered R2 is computed.
    reference = csv("results/morfeus_conformers.csv.gz")
    ligand = reference[reference.molecule_id == 369]
    assert len(ligand) == 84
    records = []
    for name, group in [("all_84", ligand), ("lower_ID_block_44", ligand[ligand.conformer_id < 64409]),
                        ("higher_ID_block_40", ligand[ligand.conformer_id >= 64409])]:
        for field in ["buried_volume", "far_vbur", "near_vbur", "max_delta_qvbur"]:
            records.append({"block": name, "N": len(group), "descriptor": field,
                "minimum": group[field].min(), "maximum": group[field].max(),
                "member_IDs": list(map(int, group.conformer_id))})
    save(out / "ligand_369.json", {"records": records,
        "status": "Two present-day ID blocks support ensemble-membership mismatch; original retention records unavailable. No member removed."})
    # All requested top-k investigations use the same full-population denominator.
    ranked = []
    for stage, path in [("baseline", "baseline/all_residuals.csv.gz"), ("final", "results/final_native/all_residuals.csv.gz")]:
        data = csv(path)
        for metric, g in data.groupby("metric"):
            g = g.sort_values(["absolute_residual", "molecule_id"], ascending=[False, True])
            total = float(np.sum(g.residual**2))
            for k in [3, 10, 25, 50, 100]:
                h = g.head(k)
                ranked.append({"stage": stage, "metric": metric, "top_k": k, "full_N": len(g),
                    "top_k_SSE": float(np.sum(h.residual**2)), "full_SSE": total,
                    "fraction_of_full_SSE": float(np.sum(h.residual**2))/total,
                    "boundary_absolute_error": h.absolute_residual.iloc[-1]})
    pd.DataFrame(ranked).to_csv(out / "worst_case_concentration.csv", index=False)
    # Association with measured present-day geometric proxies, never fabricated
    # historical RMSD. Historical RMSD is null for the entire cohort.
    assoc = []
    for metric, g in dossiers.groupby("metric"):
        for feature in ["maximum_intracohort_RMSD_A", "maximum_native_coordinate_RMSD_A", "maximum_frame_axis_angle_deg"]:
            for response in ["residual", "absolute_residual"]:
                valid = g[[feature, response]].dropna()
                assoc.append({"metric": metric, "feature": feature, "response": response, "N": len(valid),
                    "pearson_r": valid[feature].corr(valid[response]),
                    "spearman_r": valid[feature].corr(valid[response], method="spearman"),
                    "interpretation": "Present-day proxy association; does not measure disagreement with historical geometry."})
    pd.DataFrame(assoc).to_csv(out / "geometric_proxy_associations.csv", index=False)
    stage_families = []
    for (metric, family), g in dossiers.groupby(["metric", "ligand_family"]):
        for stage in [c for c in g.columns if c.startswith("value_")]:
            stage_families.append({"metric": metric, "ligand_family": family,
                "stage": stage.removeprefix("value_"), **metrics(g.kraken, g[stage])})
    pd.DataFrame(stage_families).to_csv(out / "stage_family_metrics.csv", index=False)
    save(out / "complete.json", {"script_sha256": sha(__file__), "ligands": 1541,
        "metrics": len(final), "full_comparison_rows": len(dossiers),
        "exact_historical_N": 0, "no_predictions_changed": True})


if __name__ == "__main__":
    main()
