"""Independent upstream Sterimol and statistical cross-checks on new experimental inputs."""
from pathlib import Path
import sys,json,hashlib,math,inspect
import numpy as np
from scipy.stats import linregress
from morfeus import Sterimol
ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/'.stericx/experimental_validation';OUT=BASE/'chemical_references_v2'
sys.path.append(str(ROOT/'.stericx/paton/reference/venv/lib/python3.12/site-packages'))
from dbstep import sterics,calculator

def main():
 OUT.mkdir(exist_ok=False);rows=[];failures=[];inputs=BASE/'chemistry_preparation/native.jsonl'
 (OUT/'morfeus_calculate_source.py').write_text(inspect.getsource(Sterimol.calculate));(OUT/'dbstep_sterimol_source.py').write_text(inspect.getsource(sterics.get_classic_sterimol))
 for ligand in map(json.loads,inputs.read_text().splitlines()):
  result=ligand['result']
  for c,ref in zip(result['descriptors']['conformers'],result['precision']):
   geo=c['geometry']['value'];p=ref['result'];n=p['continuous_sterimol_raw_l'];coordinates=np.array([p['center']]+geo['coordinates_angstrom']);radii=np.array([0.]+geo['radii_angstrom']);relative=coordinates-coordinates[0];axis=relative[p['donor']+1].copy();axis/=np.linalg.norm(axis);seed=np.eye(3)[int(np.argmin(abs(axis)))];x=seed-np.dot(seed,axis)*axis;x/=np.linalg.norm(x);y=np.cross(axis,x);aligned=np.column_stack([relative@x,relative@y,relative@axis]);rho=float(max(np.linalg.norm(relative-np.outer(relative@axis,axis),axis=1)));scale=float(max(np.linalg.norm(relative,axis=1)+radii));margin=128*np.finfo(float).eps*scale
   row=dict(identifier=c['identifier'],native=n,margin=margin)
   try:
    m=Sterimol(['Ni']+geo['elements'],coordinates.copy(),dummy_index=1,attached_index=p['donor']+2,radii=radii,n_rot_vectors=3600);length,b5,b1,_=sterics.get_classic_sterimol(aligned,radii,['Bq']+geo['elements']);row.update(morfeus=dict(l_raw=float(m.L_value_uncorrected),b1=float(m.B_1_value),b5=float(m.B_5_value)),dbstep=dict(l_absolute_axis=float(length),b1=float(b1),b5=float(b5)),b1_sampling_bounds=dict(morfeus=rho*math.pi/3599,dbstep=rho*math.pi/360),morfeus_b5_sampling_bound=rho*2*math.sin(math.pi/(2*3599))**2)
    assert abs(row['morfeus']['l_raw']-n['l'])<=margin
    for lib in ['morfeus','dbstep']:assert -margin<=row[lib]['b1']-n['b1']<=row['b1_sampling_bounds'][lib]+margin
    assert -margin<=n['b5']-row['morfeus']['b5']<=row['morfeus_b5_sampling_bound']+margin
    assert abs(row['dbstep']['b5']-n['b5'])<=margin
    row['passed']=True
   except Exception as e:row.update(passed=False,error=repr(e));failures.append(c['identifier'])
   rows.append(row)
  print(ligand['ligand_id'],len(rows),len(failures),flush=True)
 # Independent univariate fits use scipy's centered covariance, not the analysis's SVD.
 analysis=json.loads((BASE/'chemical_analysis_v1/statistical_results.json').read_text());import pandas as pd
 design=pd.read_csv(BASE/'chemical_analysis_v1/primary_design.csv').set_index('ligand_id');checks=[]
 for model in ['A_electronic','B_single','C_ensemble_min','D_published']:
  r=analysis['results']['primary_as_distributed'][model]
  for fold in r['folds']:
   train=design.loc[fold['train_ids']];test=design.loc[fold['test_ids']];fit=linregress(train[model],train.response);pred=fit.intercept+fit.slope*test[model].to_numpy();native=np.array([next(v['predicted'] for v in r['predictions'] if v['ligand_id']==i) for i in fold['test_ids']]);difference=float(max(abs(pred-native)));assert difference<1e-11;checks.append(dict(model=model,fold=fold['group'],max_prediction_difference=difference))
 result=dict(complete=True,rows=rows,failures=failures,statistical_crosschecks=checks,scope='Matched exported f32 geometry and explicit radii; upstream angular sampling bounds distinguished from continuous extrema. DBSTEP absolute-axis L is recorded but not equated with donor-referenced L.',input_sha256=hashlib.sha256(inputs.read_bytes()).hexdigest())
 (OUT/'summary.json').write_text(json.dumps(result,indent=2)+'\n');assert not failures
if __name__=='__main__':main()
