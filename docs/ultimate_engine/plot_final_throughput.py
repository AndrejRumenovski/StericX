#!/usr/bin/env python3
"""Standalone throughput plot: full retained sample ranges, no inferred confidence bounds."""
import hashlib,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
p=HERE/'final_measurements.json';d=json.loads(p.read_text())
threads=[1,2,4,6,12]
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,2,figsize=(11,4.9))
for side,label,color in [('baseline','Frozen original','#64748b'),('candidate','Final engine','#12645b')]:
    r=[d['cli'][f'descriptors_10000/t{t}'][side] for t in threads]
    y=[1e13/x['median_ns'] for x in r]
    low=[v-1e13/x['maximum_ns'] for x,v in zip(r,y)]
    high=[1e13/x['minimum_ns']-v for x,v in zip(r,y)]
    axes[0].errorbar(threads,y,yerr=[low,high],fmt='o-',capsize=4,label=label,color=color)
    axes[1].plot(threads,[100*d['scaling'][str(t)][side]['parallel_efficiency'] for t in threads],
                 'o-',label=label,color=color)
axes[0].set(ylabel='Files per second',title='10,000-file complete CLI workflow',ylim=(0,31000))
axes[1].set(ylabel='Parallel efficiency (%)',title='Scaling relative to each one-thread run',ylim=(0,110))
for ax in axes:
    ax.set(xlabel='Worker threads',xticks=threads);ax.grid(True,alpha=.2);ax.legend(frameon=False,loc='lower right')
fig.suptitle('StericX controlled original–final comparison',fontsize=14)
fig.text(.5,.02,'Eight alternating trials per side; error bars span every observed trial. Exact scientific outputs preserved.\nRyzen 5 5600G, six physical cores / twelve logical CPUs. Startup, parsing and output included; file caches warm.',ha='center',fontsize=9)
fig.tight_layout(rect=(0,.13,1,.94))
fig.savefig(HERE/'final_throughput.svg');fig.savefig(HERE/'final_throughput.png',dpi=160)
(HERE/'final_throughput_receipt.json').write_text(json.dumps({'source':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
    'scope':'Reciprocal of median elapsed time; all-sample range, not a confidence interval; efficiency relative to one-thread throughput'},indent=2)+'\n')
