"""Assemble descriptor-wise evidence, causal components and complete dossiers."""
import gzip
import json
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analyze import FIELDS, metrics, read_jsonl, references, summarize, write_analysis
from freeze import OUT, save, sha

STAGES = {
    "baseline": "baseline",
    "raw_center": "campaigns/raw_center/native_mean",
    "primary_radii": "campaigns/primary_radii/native_mean",
    "primary_grid_mean": "campaigns/primary_grid/native_mean",
    "historical_frame_reduction": "campaigns/primary_grid/historical_last_plane",
    "sampling_3600": "campaigns/primary_sampling_3600/analysis",
    "recovered_XYZ_precision": "campaigns/recovered_precision/analysis",
}
REASONS = {
    "baseline": "Frozen validated native defaults on the historical cohort",
    "raw_center": "Original source raw-bond sum replaces unit-bond sum in reproduction inputs",
    "primary_radii": "Original source Morfeus Bondi table and Sterimol Paton H substitution",
    "primary_grid_mean": "Original documented 0.001 density; native public frame mean retained",
    "historical_frame_reduction": "Historical source last-plane scalars; explicit reproduction view, not a general default",
    "sampling_3600": "Source-requested 3600 endpoint-inclusive directions through optional native API",
    "recovered_XYZ_precision": "Use maximum available original API precision for all 31611 IDs; 27 changed",
}


def csv(path):
    return pd.read_csv(path, float_precision="round_trip")


def minimal_values(r):
    frames = r["frames"]
    result = dict(r["pyramidalization"])
    result.update({"sterimol_"+k: v for k, v in r["sterimol_coordination"].items()})
    q = [v for f in frames for v in f["quadrants"]]
    o = [v for f in frames for v in f["octants"]]
    result.update(qvbur_min=min(q), qvbur_max=max(q), ovbur_min=min(o), ovbur_max=max(o),
        max_delta_qvbur=max(abs(f["quadrants"][i]-f["quadrants"][(i+3)%4]) for f in frames for i in range(4)))
    result.update({k: frames[-1][k] for k in ("buried_volume", "near_vbur", "far_vbur")})
    result["percent_buried_volume"] = result["buried_volume"] * 100 / (4*math.pi*3.5**3/3)
    return {"id": r["id"], **result}


def main():
    out = OUT / "results"
    out.mkdir(exist_ok=False)
    (out / "dossiers").mkdir()
    stages = []
    for stage, directory in STAGES.items():
        for row in json.loads((OUT / directory / "metrics.json").read_text()):
            stages.append({"stage": stage, "scientific_reason": REASONS[stage], **row})
    pd.DataFrame(stages).to_csv(out / "stage_metrics.csv", index=False)
    final = OUT / STAGES["recovered_XYZ_precision"]
    native = csv(final / "conformers.csv.gz")
    native_comparison = csv(final / "all_residuals.csv.gz")
    published, _ = references()
    lig = csv(OUT / "baseline/ligand_provenance.csv")
    reference_rows = []
    for r in read_jsonl(OUT / "campaigns/recovered_precision/morfeus_observations.jsonl.gz"):
        assert "error" not in r, r
        reference_rows.append({"id": r["id"], **{f: r[op][key] for f, (op, key, _, _) in FIELDS.items()}})
    reference = native.drop(columns=list(FIELDS)).merge(pd.DataFrame(reference_rows), on="id", validate="one_to_one")
    reference.to_csv(out / "morfeus_conformers.csv.gz", index=False)
    comparison, failures = summarize(reference, published, lig)
    assert not failures
    comparison.rename(columns={"stericx": "morfeus"}).to_csv(out / "morfeus_vs_published.csv.gz", index=False)
    reference_metrics = [{"metric": name, **metrics(g.kraken, g.stericx)} for name, g in comparison.groupby("metric")]
    pd.DataFrame(reference_metrics).to_csv(out / "morfeus_vs_published_metrics.csv", index=False)
    n = native.set_index("id")
    m = reference.set_index("id").loc[n.index]
    assert len(n) == len(m) == 31611 and native.molecule_id.nunique() == 1541
    direct = []
    for f in FIELDS:
        direct.append({"descriptor": f, "comparison": "native_vs_morfeus_same_export", **metrics(m[f], n[f])})
    small = pd.DataFrame([minimal_values(r) for r in read_jsonl(OUT / "minimal_reference_final/observations.jsonl.gz")]).set_index("id")
    for f in FIELDS:
        direct.append({"descriptor": f, "comparison": "native_vs_minimal_same_export", **metrics(small[f], n.loc[small.index, f])})
        direct.append({"descriptor": f, "comparison": "morfeus_vs_minimal_same_export", **metrics(small[f], m.loc[small.index, f])})
    pd.DataFrame(direct).to_csv(out / "exact_available_input_algorithm_metrics.csv", index=False)
    difference = n[list(FIELDS)] - m[list(FIELDS)]
    difference.to_csv(out / "native_minus_morfeus_per_conformer.csv.gz")
    small.to_csv(out / "minimal_conformers.csv.gz")
    save(out / "exact_historical_input_subset.json", {"N": 0, "R2": None,
        "reason": "No affirmative full historical geometry/ensemble provenance recovered. Low error is not proof of exact historical identity.",
        "separate_available_input_algorithm_comparison": {"ligands": 1541, "conformers": 31611,
            "reference": "Morfeus on identical exported inputs; minimal independent reference on all final top100 extremizers and adverse baseline cases",
            "metrics": "exact_available_input_algorithm_metrics.csv"}})
    # All observations remain; no score is recomputed on a trimmed population.
    write_analysis(native_comparison, out / "final_native", plot=True)
    nc = native_comparison.set_index(["molecule_id", "metric"])
    mc = comparison.set_index(["molecule_id", "metric"]).loc[nc.index]
    algorithms = []
    for metric in sorted(native_comparison.metric.unique()):
        a = nc.xs(metric, level="metric")
        b = mc.xs(metric, level="metric").loc[a.index]
        algorithms.append({"metric": metric, **metrics(b.stericx, a.stericx)})
    pd.DataFrame(algorithms).to_csv(out / "available_input_algorithm_ligand_metrics.csv", index=False)
    joined = native_comparison.copy()
    keys = ["molecule_id", "metric"]
    for stage, directory in STAGES.items():
        vals = csv(OUT / directory / "all_residuals.csv.gz")[keys+['stericx']].rename(columns={"stericx": "value_"+stage})
        joined = joined.merge(vals, on=keys, validate="one_to_one")
    refcols = mc.reset_index()[keys+['stericx', 'min_conformer_id', 'max_conformer_id', 'selected_conformer_id']].rename(columns={
        "stericx": "morfeus", "min_conformer_id": "morfeus_min_conformer_id", "max_conformer_id": "morfeus_max_conformer_id", "selected_conformer_id": "morfeus_selected_conformer_id"})
    joined = joined.merge(refcols, on=keys, validate="one_to_one")
    component_names = []
    names = list(STAGES)
    for a, b in zip(names[:-1], names[1:]):
        name = "component_" + b
        joined[name] = joined["value_"+a] - joined["value_"+b]
        component_names.append(name)
    joined["component_same_input_algorithm"] = joined.stericx - joined.morfeus
    joined["component_published_reconstruction"] = joined.morfeus - joined.kraken
    component_names += ["component_same_input_algorithm", "component_published_reconstruction"]
    telescoping = joined[component_names].sum(axis=1) - (joined.value_baseline-joined.kraken)
    assert np.max(abs(telescoping)) < 1e-11
    bounds = csv(OUT / "geometry_bounds/ligand_intervals.csv")
    joined = joined.merge(bounds[keys+["lower", "upper", "unavoidable_abs_error_with_rounding_only", "published_outside_interval"]], on=keys, how="left", validate="one_to_one")
    provenance = csv(OUT / "provenance_analysis/ligand_provenance.csv")
    joined = joined.merge(provenance[["molecule_id", "maximum_intracohort_RMSD_A", "maximum_native_coordinate_RMSD_A", "maximum_frame_axis_angle_deg", "connectivity_mismatch_count", "historical_exact_geometry_verified"]], on="molecule_id", validate="many_to_one")
    joined["historical_status"] = "UNRESOLVABLE FROM AVAILABLE DATA"
    joined["remaining_evidence_class"] = np.where(abs(joined.component_published_reconstruction) <= 1e-8,
        "published_precision_compatible", np.where(joined.published_outside_interval.fillna(False),
        "beyond_coordinate_rounding_under_fixed_ensemble", np.where(joined.descriptor.isin(["pyr_p", "pyr_alpha", "sterimol_l", "sterimol_b1", "sterimol_b5"]),
        "coordinate_rounding_or_historical_provenance_unresolved", "volume_historical_provenance_or_rounding_unresolved")))
    joined["reference_precision_note"] = "1e-8 compatibility label follows maximum observed SI/CSV precision difference; not a causal attribution or acceptance tolerance"
    joined.to_csv(out / "all_residual_dossiers.csv.gz", index=False)
    covariance = []
    variance = []
    sse_partition = []
    for metric, g in joined.groupby("metric"):
        sst = float(np.sum((g.kraken-g.kraken.mean())**2))
        a = g.component_same_input_algorithm.to_numpy()
        b = g.component_published_reconstruction.to_numpy()
        variance.append({"metric": metric, "N": len(g), "SST": sst, "native_vs_published_SSE": float(np.sum(g.residual**2)),
            "same_input_algorithm_SSE": float(a@a), "published_reconstruction_SSE": float(b@b), "cross_term": float(2*(a@b)),
            "identity_error": float(np.sum(g.residual**2)-(a@a+b@b+2*(a@b))),
            "warning": "Components covary; their SSEs cannot be added without cross-term. Unavailable geometry/conformer/reference causes cannot be uniquely partitioned."})
        for category, h in g.groupby("remaining_evidence_class"):
            sse = float(np.sum(h.residual**2))
            sse_partition.append({"metric": metric, "evidence_class": category, "N": len(h), "SSE": sse,
                "SSE_fraction": sse/float(np.sum(g.residual**2)) if np.sum(g.residual**2) else 0, "lost_R2": sse/sst})
        for i, aa in enumerate(component_names):
            for bb in component_names[i:]:
                covariance.append({"metric": metric, "component_a": aa, "component_b": bb,
                    "SSE_term": float(np.dot(g[aa], g[bb]) * (1 if aa == bb else 2))})
        ranked = g.sort_values(["absolute_residual", "molecule_id"], ascending=[False, True]).head(100)
        # Typed JSON nulls, never fabricated zeros for missing historical inputs.
        with (out / "dossiers" / f"{metric}_top100.json").open("w") as dst:
            dst.write(ranked.to_json(orient="records", indent=2, double_precision=15) + "\n")
    pd.DataFrame(variance).to_csv(out / "remaining_residual_variance_components.csv", index=False)
    pd.DataFrame(covariance).to_csv(out / "all_stage_SSE_cross_terms.csv", index=False)
    pd.DataFrame(sse_partition).to_csv(out / "remaining_evidence_class_SSE.csv", index=False)
    # Resolve vburminconf selection amplification without assigning historical IDs.
    selections = []
    for mid, gn in native.groupby("molecule_id"):
        gm = reference[reference.molecule_id == mid]
        ni = gn.loc[gn.buried_volume.idxmin()]
        mi = gm.loc[gm.buried_volume.idxmin()]
        same_ref = gm[gm.conformer_id == ni.conformer_id].iloc[0]
        for field in FIELDS:
            selections.append({"molecule_id": int(mid), "descriptor": field, "native_selected_id": int(ni.conformer_id),
                "morfeus_selected_id": int(mi.conformer_id), "selected_IDs_equal": int(ni.conformer_id) == int(mi.conformer_id),
                "native_value": ni[field], "reference_at_native_selection": same_ref[field], "reference_at_reference_selection": mi[field],
                "algorithm_at_fixed_conformer": ni[field]-same_ref[field],
                "selection_component": same_ref[field]-mi[field], "historical_selected_id": None})
    pd.DataFrame(selections).to_csv(out / "vburminconf_selection_forensics.csv.gz", index=False)
    save(out / "complete.json", {"ligands": 1541, "conformers": 31611, "descriptor_fields": len(FIELDS), "metrics": len(reference_metrics),
        "independent_minimal_conformers": len(small), "comparison_rows": len(joined), "failures": failures,
        "maximum_telescoping_error": float(np.max(abs(telescoping))),
        "script_sha256": sha(__file__), "native_source": str(final.relative_to(OUT)),
        "same_input_metrics_sha256": sha(out / "exact_available_input_algorithm_metrics.csv"),
        "residual_dossiers_sha256": sha(out / "all_residual_dossiers.csv.gz")})
    print(json.dumps(json.loads((out / "complete.json").read_text()), indent=2))


if __name__ == "__main__":
    main()
