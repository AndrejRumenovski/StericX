"""Untrimmed agreement and residual plots for every frozen Kraken metric."""
from pathlib import Path
import hashlib,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

A=Path(__file__).resolve().parents[1]
source=A/'kraken/analysis/comparisons.csv'
frame=pd.read_csv(source)
metrics=json.loads((A/'kraken/analysis/metrics.json').read_text())
out=A/'visuals/kraken_all_metrics';out.mkdir(parents=True,exist_ok=True)
index=['# All Kraken comparison residual plots','','Each plot contains every complete paired ligand for the named descriptor/reduction. Red points identify its twenty largest absolute residuals; none are excluded. The agreement line is the identity line, not a fitted correction. Full values and coverage remain in the Kraken analysis tables.','']
for name,g in frame.groupby('metric',sort=True):
    m=metrics[name];unit=str(g.units.iloc[0]);x=g.reference.to_numpy();y=g.sut.to_numpy();r=y-x
    worst=np.argsort(abs(r))[-min(20,len(r)):]
    fig,ax=plt.subplots(1,2,figsize=(10.5,4.1),constrained_layout=True)
    for a,z in zip(ax,[y,r]):
        a.scatter(x,z,s=8,alpha=.42,color='#2878a5',rasterized=True)
        a.scatter(x[worst],z[worst],s=13,alpha=.8,color='#c24b3a',rasterized=True)
        a.set_xlabel(f'Published reference ({unit})');a.grid(alpha=.13)
    lo=min(x.min(),y.min());hi=max(x.max(),y.max())
    ax[0].plot([lo,hi],[lo,hi],color='black',lw=.8)
    ax[0].set_ylabel(f'StericX ({unit})')
    ax[1].axhline(0,color='black',lw=.8);ax[1].set_ylabel(f'StericX − reference ({unit})')
    fig.suptitle(f"{name} · N={len(g)}\nMAE={m['MAE']:.6g}; RMSE={m['RMSE']:.6g}; max |error|={m['maximum_absolute_error']:.6g} {unit}",fontsize=11)
    fig.savefig(out/f'{name}.png',dpi=165);plt.close(fig)
    index.append(f'- [{name}]({name}.png) — N={len(g)}, units {unit}')
(out/'INDEX.md').write_text('\n'.join(index)+'\n')
(out/'manifest.json').write_text(json.dumps({'input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'plot_count':frame.metric.nunique(),'outlier_policy':'all rows retained; twenty largest residuals highlighted, never trimmed','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
print(f'Wrote {frame.metric.nunique()} plots')
