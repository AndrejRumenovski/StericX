"""Score both immutable reference versions explicitly; never switch silently."""
import math

import numpy as np
import pandas as pd

from analyze import FIELDS, metrics, write_analysis
from consolidate import STAGES, REASONS
from freeze import OUT, save, sha


def csv(path):
    return pd.read_csv(OUT / path, float_precision="round_trip")


def main():
    out = OUT / "historical_target"
    out.mkdir(exist_ok=False)
    save(out / "plan.json", {"script_sha256": sha(__file__),
        "reference": "frozen/historical_reference.csv",
        "reference_sha256": sha(OUT / "frozen/historical_reference.csv"),
        "reason": "Explicit historical-reference comparison after identifying substantive API drift for 821 and 1036. Previously frozen API-target results remain unchanged. No computed descriptor changes.",
        "selection": "All original 1541 IDs, every metric, every stage; no reference chosen by goodness of fit."})
    hist = csv("frozen/historical_reference.csv").set_index("Unnamed: 0")
    def replace_reference(d):
        values = []
        for row in d.itertuples():
            prop = FIELDS[row.descriptor][2]
            value = hist.loc[row.molecule_id, prop+"_"+row.reduction]
            if row.descriptor == "percent_buried_volume":
                value *= 100/(4*math.pi*3.5**3/3)
            values.append(value)
        d = d.copy()
        d["api_reference"] = d.kraken
        d["kraken"] = values
        d["reference_version_difference"] = d.api_reference-d.kraken
        d["residual"] = d.stericx-d.kraken
        d["absolute_residual"] = abs(d.residual)
        return d
    stages = []
    for stage, path in STAGES.items():
        d = replace_reference(csv(path+"/all_residuals.csv.gz"))
        d.to_csv(out / (stage+"_residuals.csv.gz"), index=False)
        for metric, g in d.groupby("metric"):
            stages.append({"stage": stage, "scientific_reason": REASONS[stage], "metric": metric,
                           **metrics(g.kraken, g.stericx)})
        if stage in ("baseline", "recovered_XYZ_precision"):
            write_analysis(d, out / ("baseline" if stage == "baseline" else "final"), plot=True)
    pd.DataFrame(stages).to_csv(out / "stage_metrics.csv", index=False)
    d = replace_reference(csv("results/all_residual_dossiers.csv.gz"))
    d["component_reference_version"] = d.api_reference-d.kraken
    d["component_published_reconstruction"] = d.morfeus-d.kraken
    bounded = d.lower.notna()
    d.loc[bounded, "unavoidable_abs_error_with_rounding_only"] = np.maximum(
        0, np.maximum(d.loc[bounded, "lower"]-d.loc[bounded, "kraken"],
                      d.loc[bounded, "kraken"]-d.loc[bounded, "upper"]))
    d.loc[bounded, "published_outside_interval"] = d.loc[bounded, "unavoidable_abs_error_with_rounding_only"] > 0
    d["remaining_evidence_class"] = np.where(abs(d.component_published_reconstruction) <= 1e-8,
        "published_rounding_precision_compatible", np.where(d.published_outside_interval.fillna(False),
        "beyond_coordinate_rounding_under_fixed_ensemble", np.where(d.descriptor.isin(["pyr_p", "pyr_alpha", "sterimol_l", "sterimol_b1", "sterimol_b5"]),
        "coordinate_rounding_or_historical_provenance_unresolved", "volume_historical_provenance_or_rounding_unresolved")))
    d["reference_precision_note"] = "1e-8 describes precision compatibility for numerically unchanged reference entries, NOT maximum dataset-version discrepancy. 821 and 1036 have substantive API/SI differences. No acceptance threshold."
    d.to_csv(out / "all_residual_dossiers.csv.gz", index=False)
    (out / "dossiers").mkdir()
    variance, partitions, bounds, morfeus = [], [], [], []
    for metric, g in d.groupby("metric"):
        ranked = g.sort_values(["absolute_residual", "molecule_id"], ascending=[False, True]).head(100)
        (out / "dossiers" / (metric+"_top100.json")).write_text(ranked.to_json(orient="records", indent=2, double_precision=15)+"\n")
        a = (g.stericx-g.morfeus).to_numpy()
        b = (g.morfeus-g.kraken).to_numpy()
        sst = float(np.sum((g.kraken-g.kraken.mean())**2))
        sse = float(np.sum(g.residual**2))
        variance.append({"metric": metric, "N": len(g), "SST": sst, "native_vs_published_SSE": sse,
            "same_input_algorithm_SSE": float(a@a), "published_reconstruction_SSE": float(b@b), "cross_term": float(2*(a@b)),
            "identity_error": sse-float(a@a+b@b+2*(a@b))})
        for category, h in g.groupby("remaining_evidence_class"):
            v = float(np.sum(h.residual**2))
            partitions.append({"metric": metric, "evidence_class": category, "N": len(h), "SSE": v, "SSE_fraction": v/sse, "lost_R2": v/sst})
        if g.lower.notna().all():
            v = float(np.sum(g.unavoidable_abs_error_with_rounding_only**2))
            bounds.append({"metric": metric, "N": len(g), "outside_interval": int(g.published_outside_interval.sum()),
                "SST": sst, "rounding_only_SSE_lower_bound": v, "conditional_R2_upper_bound": 1-v/sst, "RMSE_lower_bound": math.sqrt(v/len(g))})
        morfeus.append({"metric": metric, **metrics(g.kraken, g.morfeus)})
    pd.DataFrame(variance).to_csv(out / "remaining_residual_variance_components.csv", index=False)
    pd.DataFrame(partitions).to_csv(out / "remaining_evidence_class_SSE.csv", index=False)
    pd.DataFrame(bounds).to_csv(out / "conditional_ceilings.csv", index=False)
    pd.DataFrame(morfeus).to_csv(out / "morfeus_vs_published_metrics.csv", index=False)
    drift = d.loc[abs(d.reference_version_difference) > 1e-7,
        ["molecule_id", "metric", "kraken", "api_reference", "reference_version_difference", "stericx", "morfeus"]]
    drift.to_csv(out / "substantive_reference_version_differences.csv", index=False)
    save(out / "complete.json", {"N": 1541, "comparison_rows": len(d), "metrics": d.metric.nunique(),
        "substantive_reference_version_ligands": list(map(int, sorted(drift.molecule_id.unique()))),
        "substantive_reference_version_rows": len(drift), "no_predictions_changed": True,
        "erratum": "Earlier results/all_residual_dossiers.csv.gz reference_precision_note incorrectly called 1e-8 the maximum observed SI/CSV difference. It is only a numerical-precision compatibility label. This audit records the substantive 821/1036 version differences explicitly; all earlier numerical results are preserved."})


if __name__ == "__main__":
    main()
