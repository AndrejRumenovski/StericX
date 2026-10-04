"""Readable standalone family diagnostics; all ligands and classes remain."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from freeze import OUT, save, sha


def main():
    out = OUT / "family_plots_v2"
    out.mkdir(exist_ok=False)
    data = pd.read_csv(OUT / "results/all_residual_dossiers.csv.gz", float_precision="round_trip")
    for metric, g in data.groupby("metric"):
        families = sorted(g.ligand_family.unique())
        fig, axes = plt.subplots(1, 2, figsize=(15, max(12, .26*len(families))), sharey=True)
        labels = []
        for i, family in enumerate(families):
            h = g[g.ligand_family == family]
            labels.append(f"{family} (N={len(h)})")
            jitter = ((h.molecule_id.to_numpy()*17 % 101)/100 - .5)*.5
            axes[0].scatter(h.residual, i+jitter, s=8, alpha=.55)
            before = np.mean(abs(h.value_baseline-h.kraken))
            axes[1].barh(i-.17, before, height=.3, color="#adbac5", label="Baseline" if i == 0 else None)
            axes[1].barh(i+.17, h.absolute_residual.mean(), height=.3, color="#186ea3", label="Source profile" if i == 0 else None)
        axes[0].set_yticks(np.arange(len(families)), labels, fontsize=8)
        axes[0].axvline(0, color="black", linewidth=.5)
        axes[0].set_xlabel("Final residual: StericX − Kraken")
        axes[1].set_xlabel("Mean absolute error (all members)")
        axes[1].legend()
        axes[0].invert_yaxis()
        fig.suptitle(f"{metric} | N=1541 | {g.iloc[0].units}\nFamilies defined by donor-neighbor elements and aromatic-neighbor count")
        fig.tight_layout(rect=(0, 0, 1, .965))
        fig.savefig(out / (metric+".png"), dpi=110)
        plt.close(fig)
    save(out / "complete.json", {"script_sha256": sha(__file__), "metrics": data.metric.nunique(), "N_per_metric": 1541})


if __name__ == "__main__":
    main()
