"""Independent analytical profiles plus standalone and molecular interval checks."""
from pathlib import Path
import argparse,copy,json,subprocess
from validate_regions import ROOT,make_case,sha

def main():
 p=argparse.ArgumentParser();p.add_argument('--engine',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
 qs=[]
 for axis in range(3):
  for i,balls in enumerate([[(0.,1.)],[(-.75,1.),(.25,1.),(1.,.75)],[(2.,.5),(-2.,.5)],[(.5,.125)]]):
   radii=[.5,1.,1.5,2.5,3.5,1.5]
   refs=[make_case('ref',(0.,r),balls,axis,[0.,0.,0.]) for r in radii]
   qs.append(dict(id=f'analytic_{i}_axis{axis}',center=[0.,0.,0.],atoms=refs[0]['atoms'],radii=radii,tolerance=.02,budget=262144,controls=True,references=[r['references'] for r in refs]))
 molecules=json.loads((ROOT/'.stericx/occupancy_precision/molecular_chord_requests.json').read_text())
 for q in molecules:
  qs.append(dict(id=q['id'],center=q['region']['center'],atoms=q['atoms'],radii=[2.5,3.5,4.5,3.5],tolerance=.05,budget=1048576,controls=True))
 path=a.output/'requests.json';path.write_text(json.dumps(qs,indent=2)+'\n')
 with (a.output/'observed.jsonl').open('x') as out,(a.output/'process_metrics.txt').open('x') as err:subprocess.run(['/usr/bin/time','-v',str(a.engine.resolve()),str(path)],stdout=out,stderr=err,check=True)
 rs=list(map(json.loads,(a.output/'observed.jsonl').read_text().splitlines()));fail=[];unmet=[];points=0;analytic=0
 for q,row in zip(qs,rs,strict=True):
  assert q['id']==row['id'];p=row['profile'];work=0;queries=0
  for i,(e,c) in enumerate(zip(p['entries'],row['controls'],strict=True)):
   points+=1
   if e['result'] is None:fail.append([q['id'],i,'error',e]);continue
   r=e['result'];direct=copy.deepcopy(c['regions']);saved=copy.deepcopy(r)
   direct.pop('elapsed_seconds');saved.pop('elapsed_seconds')
   if direct!=saved:fail.append([q['id'],i,'standalone_identity'])
   if e['reused_from_index'] is None:work+=r['cells_evaluated'];queries+=r['geometric_evaluations']
   elif r['elapsed_seconds']!=0:fail.append([q['id'],i,'reuse_timing'])
   if not r['target_met']:unmet.append([q['id'],i,q['radii'][i]])
   total=c['total'];t=r['total']
   if max(total['lower_bound_angstrom3'],t['lower_bound_angstrom3'])>min(total['upper_bound_angstrom3'],t['upper_bound_angstrom3']):fail.append([q['id'],i,'total_disagreement'])
   if 'references' in q:
    for kind in ['total','quadrants','octants']:
     ms=[r[kind]] if kind=='total' else r[kind];refs=[q['references'][i][kind]] if kind=='total' else q['references'][i][kind]
     for m,ref in zip(ms,refs,strict=True):
      analytic+=1
      if not m['lower_bound_angstrom3']<=ref[0]<=ref[1]<=m['upper_bound_angstrom3']:fail.append([q['id'],i,kind,'containment'])
  if work!=p['cells_evaluated'] or queries!=p['geometric_evaluations']:fail.append([q['id'],'work_accounting'])
  unique=sorted({e['radius_angstrom']:e['result'] for e in p['entries'] if e['result']}.items())
  for (_,left),(_,right) in zip(unique,unique[1:]):
   for kind in ['total','quadrants','octants']:
    lm=[left[kind]] if kind=='total' else left[kind];rm=[right[kind]] if kind=='total' else right[kind]
    for l,r in zip(lm,rm,strict=True):
     if l['lower_bound_angstrom3']>r['upper_bound_angstrom3']:fail.append([q['id'],kind,'monotonicity_contradiction'])
 summary=dict(profiles=len(qs),requested_radii=points,independent_analytical_integrals=analytic,failures=fail,unmet=unmet,inputs_sha256=sha(path),observed_sha256=sha(a.output/'observed.jsonl'),engine_sha256=sha(a.engine),controller_sha256=sha(Path(__file__)),scope='Exact rational analytical profiles; molecular independent-total overlap is consistency, not exact molecular truth. Standalone identity includes every bound, work counter and status. Profile timing excludes controls.')
 (a.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
 if fail:raise SystemExit(1)
if __name__=='__main__':main()
