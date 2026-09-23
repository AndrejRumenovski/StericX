"""Retain descriptor-level investigations; diagnostics never define thresholds."""
import csv,json,sys
from pathlib import Path
R=Path(__file__).resolve().parent;H=R.parent
sys.path.insert(0,str(H));from compare import metrics

def read(p):return json.loads(p.read_text())
def write(p,d):p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')

def main():
 ev=R/'evaluation_2';gate=read(ev/'scientific_gate.json');assert gate['passed']
 records=read(ev/'volume_records.json');analytic={r['filename']:r for r in read(H/'sterimol_diagnosis.json')['rows']}
 residuals=list(csv.DictReader((ev/'scientific_comparisons.csv').open()))
 per_molecule=[]
 for record in records:
  filename=record['filename'];rows=[r for r in residuals if r['filename']==filename]
  assert len(rows)==42
  entry=dict(filename=filename,donor=record['donor'],neighbors=record['neighbors'],mapping_hydrogen_frame_evidence='Independent exact-decimal topology and interval frame: matching mapping and all point decisions in three planes. Source and controlled diagnostic lanes retained.',descriptors={})
  for field in sorted({r['descriptor'] for r in rows}):
   ds={r['lane']:{k:float(r[k]) if r[k] else None for k in ['stericx','morfeus','absolute_difference','relative_difference']} for r in rows if r['descriptor']==field}
   d=dict(residuals=ds)
   if field in gate['admitted_descriptors']:
    d.update(status='admitted finite estimator',cause='Exactly accounted for by independently reconstructed f32/f64 constants and arithmetic; all individual occupancy predicates certified identical.',reference_native=record['native_model'][field],reference_morfeus=record['morfeus_model'][field],ambiguous_or_differing_points=0)
   elif field.startswith('sterimol'):
    k=field.rsplit('_',1)[-1];a=analytic[filename]
    d.update(status='unresolved; excluded',analytic_diagnostic={key:a[key][k] for key in ['analytic','stericx','morfeus','signed_stericx_error','signed_morfeus_error']})
    d['cause']=dict(l='Common continuous support formula; remaining discrepancy involves input/axis/floating reduction. Diagnostic analytic comparison is not a certified forward-error bound.',b1='Different angular phase and finite grids estimate the same continuous minimum; matching unique angular counts alone does not align phase. The continuous support-envelope diagnostic quantifies each estimator error separately, without setting a tolerance.',b5='Native radial maximum is analytic; morfeus angular support is sampled. Diagnostic analytic maximum separates angular underestimation from native precision residual; full certified input/axis propagation remains absent.')[k]
   else:
    d.update(status='unresolved; excluded',cause='Same continuous trivalent definition, differing coordinate precision and independent algebra/trigonometric evaluation. Controlled represented-input comparison isolates part of the discrepancy; matching rounded outputs is not a certified error budget.')
   entry['descriptors'][field]=d
  per_molecule.append(entry)
 write(R/'descriptor_investigations.json',dict(all_56_examined=True,thresholds_used=False,rows=per_molecule,reference_implementation_correction=read(R/'reference_revision.json'),first_evaluation_failures=read(R/'evaluation_1/scientific_gate.json')['failures']))
 # Actual admitted scalar path, including percentage calculated after volume averaging.
 fields=gate['admitted_descriptors'];out=[];summary={}
 for field in fields:
  a=[r['morfeus_model'][field] for r in records];b=[r['native_model'][field] for r in records];summary[field]=metrics(a,b)
  for r,x,y in zip(records,a,b,strict=True):out.append(dict(filename=r['filename'],descriptor=field,stericx=y,morfeus=x,absolute_difference=abs(y-x),relative_difference=abs((y-x)/x) if x else None))
 write(R/'admitted_scalar_comparisons.json',dict(metrics=summary,rows=out,scope='Actual admitted nine-output expressions; exact rational models were checked bitwise against freshly executed tools.'))
 with (R/'admitted_scalar_comparisons.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
 print('Investigated all 56 × 14 descriptors; preserved all raw lanes and 504 admitted scalar residuals.')

if __name__=='__main__':main()
