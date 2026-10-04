"""Summarize the fixed grid ladder without fitting a convergence law to targets."""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analyze import FIELDS, read_jsonl, references
from freeze import OUT, save, sha


def main():
    out = OUT / "convergence/analysis"
    out.mkdir(exist_ok=False)
    fields = list(FIELDS)[:8]
    rows = []
    for r in read_jsonl(OUT / "convergence/observations.jsonl.gz"):
        assert "error" not in r
        for view in ("native_mean", "historical_last"):
            rows.append({"id": r["id"], "molecule_id": int(r["id"].split(":")[1]),
                "density": r["density"], "view": view, "grid_count": r["grid_count"],
                **{f: r[view][f] for f in fields}})
    cf = pd.DataFrame(rows)
    cf.to_csv(out / "conformers.csv.gz", index=False)
    refs, _ = references()
    reduced = []
    for (mid, density, view), group in cf.groupby(["molecule_id", "density", "view"]):
        for field in fields:
            prop = FIELDS[field][2]
            vals = group[field].to_numpy()
            for red in ("min", "max", "delta", "vburminconf"):
                value = {"min": vals.min(), "max": vals.max(), "delta": vals.max()-vals.min(),
                    "vburminconf": vals[int(np.argmin(group.buried_volume))]}[red]
                published = refs[int(mid)][prop][red]
                reduced.append({"molecule_id": mid, "density": density, "view": view, "metric": field+"_"+red,
                    "N_conformers": len(group), "calculated": value, "published": published,
                    "residual": value-published, "grid_count": int(group.grid_count.iloc[0])})
    df = pd.DataFrame(reduced)
    df.to_csv(out / "ligand_grid_ladder.csv.gz", index=False)
    rows = []
    for (mid, metric, view), g in df.groupby(["molecule_id", "metric", "view"]):
        g = g.sort_values("density", ascending=False)
        values = g.calculated.to_numpy()
        source = float(g.loc[g.density == .001, "calculated"].iloc[0])
        pub = float(g.published.iloc[0])
        rows.append({"molecule_id": mid, "metric": metric, "view": view,
            "finest_estimate": values[-1], "last_grid_change": abs(values[-1]-values[-2]),
            "last_two_steps_envelope": max(abs(values[-1]-values[-2]), abs(values[-2]-values[-3])),
            "source_grid_value": source, "published": pub,
            "source_grid_abs_residual": abs(source-pub), "finest_abs_residual": abs(values[-1]-pub),
            "published_closer_to_source_grid_than_finest": abs(source-pub) < abs(values[-1]-pub),
            "rigorous_bound": False})
    estimates = pd.DataFrame(rows)
    estimates.to_csv(out / "finest_estimates.csv", index=False)
    headline = df[(df.metric == "max_delta_qvbur_min") & (df.view == "historical_last")]
    ids = [1796, 1455, 1454, 1905, 163, 873, 1214, 542, 1487]
    fig, axes = plt.subplots(3, 3, figsize=(13, 10))
    for ax, mid in zip(axes.flat, ids):
        g = headline[headline.molecule_id == mid].sort_values("density", ascending=False)
        ax.semilogx(g.density, g.calculated, "o-", label="StericX primary profile")
        ax.axhline(g.published.iloc[0], c="black", ls="--", label="Published Kraken")
        ax.invert_xaxis()
        ax.set(title=f"Ligand {mid}", xlabel="density (Å³/nominal point), finer →", ylabel="min(max adjacent ΔQ), Å³")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "headline_grid_convergence.png", dpi=150)
    plt.close(fig)
    save(out / "complete.json", {"conformer_grid_records": len(cf)//2,
        "ligands": int(df.molecule_id.nunique()), "densities": sorted(map(float, df.density.unique()), reverse=True),
        "source_grid_closer_headline_count": int(estimates[(estimates.metric == 'max_delta_qvbur_min') & (estimates.view == 'historical_last')].published_closer_to_source_grid_than_finest.sum()),
        "limits": "Finest grid is an estimate, successive differences are empirical sensitivity not certified intervals. No target-fitted Richardson order, grid selection or production default change.",
        "script_sha256": sha(__file__), "data_sha256": sha(out / "ligand_grid_ladder.csv.gz")})


if __name__ == "__main__":
    main()
