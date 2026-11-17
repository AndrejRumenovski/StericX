#!/usr/bin/env python3
"""Standalone scientific figures from retained exploratory samples, never invented timings."""
import argparse,hashlib,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--final',action='store_true');args=p.parse_args()
prefix='final_' if args.final else ''
paths=[HERE/(prefix+'n1_accuracy_measurements.json'),HERE/(prefix+'n2_accuracy_measurements.json')]
n1,n2=[json.loads(p.read_text()) for p in paths]
one=next(r for r in n1['rows'] if r['id']=='random_2_32')
two=[r for r in n2['rows'] if r['id'].startswith('three_overlapping_chain/')]
rows=[{'case':'32-disc Sterimol','method':r['method'],'time_us':r['median_ns']/1000,
       'error_upper':r['absolute_error_upper'],'certified_bound':None} for r in one['methods']]
rows += [{'case':'Three-ball volume','method':r['id'].split('/')[1],'time_us':r['median_ns']/1000,
          'error_upper':r['absolute_error_upper'],'certified_bound':r['result']['absolute_error_bound_angstrom3']} for r in two]
for row in rows:
    peers=[r for r in rows if r['case']==row['case']]
    row['empirical_pareto']=not any(r['time_us']<=row['time_us'] and r['error_upper']<=row['error_upper']
        and (r['time_us']<row['time_us'] or r['error_upper']<row['error_upper']) for r in peers)
receipt={'scope':('Final quiet-window medians; empirical upper residual to independent reference bracket, no universal optimality' if args.final else 'Empirical upper residual to independent reference bracket; exploratory medians, not final quiet-window timing or universal optimality'),
    'inputs':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'rows':rows}
(HERE/(prefix+'accuracy_cost_frontier.json')).write_text(json.dumps(receipt,indent=2)+'\n')
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,2,figsize=(12,5.3))
for ax,case in zip(axes,['32-disc Sterimol','Three-ball volume']):
    subset=[r for r in rows if r['case']==case]
    for row in subset:
        is_ref=row['method'].startswith('reference');is_event=row['method'].startswith('continuous')
        color='#b45309' if is_ref else '#12645b' if is_event else '#2563a6'
        ax.scatter(row['time_us'],max(row['error_upper'],1e-16),c=color,marker='D' if is_ref else 'o',s=45,zorder=3)
        if row['certified_bound'] is not None:
            ax.plot([row['time_us']]*2,[row['error_upper'],row['certified_bound']],c=color,alpha=.55)
            ax.scatter(row['time_us'],row['certified_bound'],facecolors='none',edgecolors=color,marker='s',s=45)
    frontier=sorted([r for r in subset if r['empirical_pareto']],key=lambda r:r['time_us'])
    ax.plot([r['time_us'] for r in frontier],[max(r['error_upper'],1e-16) for r in frontier],
        '--',c='#555555',linewidth=1,label='Empirical Pareto frontier')
    ax.set(xscale='log',yscale='log',xlabel='Native median runtime (µs)',title=case)
    ax.grid(True,which='major',alpha=.17)
axes[0].set_ylabel('B1 absolute residual upper estimate (Å)')
axes[1].set_ylabel('Volume absolute residual / bound (Å³)')
for method,label,offset in [('default_360_f32','Default 360',(-12,12)),('continuous_events_f64','Continuous events',(8,9)),('endpoint_scan_23041_f32','23,041 angles',(-87,-18))]:
    r=next(r for r in rows if r['case']=='32-disc Sterimol' and r['method']==method)
    axes[0].annotate(label,(r['time_us'],r['error_upper']),xytext=offset,textcoords='offset points',fontsize=9)
axes[1].plot([],[],'o',c='#2563a6',label='Midpoint slices: observed residual')
axes[1].plot([],[],'D',c='#b45309',label='Adaptive boxes: observed residual')
axes[1].plot([],[],'s',markerfacecolor='none',c='#b45309',label='Adaptive boxes: certified bound')
handles,labels=axes[1].get_legend_handles_labels()
fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.74,.09),ncol=2,fontsize=8)
fig.suptitle('Accuracy and computational cost on two frozen geometries',fontsize=14)
footer=('Controlled quiet-window samples. Same quantity and input within each panel.\nN1: independent high-precision bracket. N2: rational analytical enclosures. Residuals and certified bounds are distinct.' if args.final else 'Exploratory desktop samples. Same quantity and input within each panel.\nN1 reference is a high-precision bracket; N2 uses rational analytical enclosures. Final controlled timing is pending.')
fig.text(.5,.025,footer,ha='center',fontsize=9)
fig.tight_layout(rect=(0,.19,1,.95))
fig.savefig(HERE/(prefix+'accuracy_cost_frontier.svg'));fig.savefig(HERE/(prefix+'accuracy_cost_frontier.png'),dpi=160)
