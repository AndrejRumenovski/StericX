"""Export standalone scientific figures from the admitted numerical tables."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(__file__).resolve().parent
data=json.loads((root/'final_measurements.json').read_text())
pareto=json.loads((root/'accuracy_pareto.json').read_text())
plt.rcParams.update({'font.size':10,'svg.fonttype':'none'})
fig,ax=plt.subplots(figsize=(7.4,4.6),layout='constrained')
threads=[1,2,4,6,12]
for side,label in [('baseline','Frozen baseline'),('candidate','Final engine')]:
    ax.plot(threads,[data['scaling'][str(t)][side]['files_per_second'] for t in threads],'o-',label=label)
ax.set(xlabel='Worker threads',ylabel='Files per second',title='10,000-file descriptor batch; median of eight trials',xticks=threads)
ax.grid(alpha=.25);ax.legend();fig.savefig(root/'throughput_scaling.svg');fig.savefig(root/'throughput_scaling.png',dpi=180);plt.close(fig)
cases=sorted({r['case'] for r in pareto['rows']})
fig,axes=plt.subplots(6,3,figsize=(14,20),layout='constrained')
for ax,case in zip(axes.flat,cases,strict=True):
    rows=[r for r in pareto['rows'] if r['case']==case]
    for prefix,color,marker in [('certified','tab:blue','o'),('reference','tab:orange','s'),('midpoint','tab:green','^')]:
        selected=[r for r in rows if r['method'].startswith(prefix)]
        for r in selected:
            x=r['ns']['median']/1e6
            if r['certified_error'] is not None:
                ax.scatter(x,max(r['certified_error'],1e-16),c=color,marker=marker,s=38)
            ax.scatter(x,max(r['observed_error_interval'][1],1e-16),edgecolors=color,facecolors='none',marker=marker,s=40)
        frontier=sorted([r for r in rows if r['certified_pareto']],key=lambda r:r['ns']['median'])
    ax.plot([r['ns']['median']/1e6 for r in frontier],[max(r['certified_error'],1e-16) for r in frontier],color='tab:blue',alpha=.45)
    ax.set(xscale='log',yscale='log',title=case.replace('/reference',''),xlabel='Public API latency (ms)',ylabel='Error upper bound (Å³)')
    ax.grid(alpha=.2)
    for threshold in [.1,.05,.01]:ax.axhline(threshold,color='gray',lw=.5,alpha=.3)
fig.suptitle('Per-input accuracy/cost frontier\nFilled: intrinsic certificate; open: observed error upper bound. Blue: adaptive chord; orange: 3D boxes; green: midpoint.\nMolecular observed error uses a separate same-algorithm 0.001 certificate; coaxial references are independent.',fontsize=13)
fig.savefig(root/'accuracy_cost_pareto.svg');fig.savefig(root/'accuracy_cost_pareto.png',dpi=140);plt.close(fig)
