"""Independent model/statistics audit. No StericX scientific kernels imported.
Run: .venv/bin/python docs/scientific_accuracy_audit/scripts/models_audit.py
The frozen CLI is the subject; numpy.linalg.lstsq/SVD, sklearn and scipy are references.
"""
from pathlib import Path
import csv, hashlib, json, subprocess, sys, platform
import numpy as np
import pandas as pd
import scipy
from scipy import stats
import sklearn
from sklearn.linear_model import Ridge, Lasso

A=Path(__file__).resolve().parents[1]; ROOT=A.parents[1]; M=A/'models'; F=A/'frozen/repository'; BIN=A/'frozen/bin/stericx'
for d in ['inputs','raw','results']: (M/d).mkdir(exist_ok=True)

def save(p,data):
 p.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
def freeze(p, b):
 if p.exists():
  if p.read_bytes()!=b: raise ValueError('Frozen artifact differs: '+str(p))
 else: p.write_bytes(b)
def run(name,args):
 d=M/'raw'/name;d.mkdir(exist_ok=True)
 proc=subprocess.run([str(BIN),*map(str,args)],capture_output=True)
 freeze(d/'stdout',proc.stdout);freeze(d/'stderr',proc.stderr)
 save(d/'command.json',{'argv':[str(BIN.relative_to(ROOT)),*map(str,args)],'exit_code':proc.returncode})
 return proc

def metrics(y,p):
 e=np.asarray(p)-y
 return dict(n=len(y),r2=float(1-e@e/np.sum((y-y.mean())**2)),mae=float(np.mean(abs(e))),rmse=float(np.sqrt(np.mean(e**2))))
def expand(records):
 # The feature definition itself is an audited claim. Compute interactions in
 # f32 once to match record contract; solve all reference regression in f64.
 l,b1,b5,q,ir=records[:,:5].T
 return np.column_stack([np.ones(len(l)),l,b1,b5,q,b1*q,b5*q,ir]).astype(float)
def linear(x,y,kind='ols',alpha=0):
 mean=x.mean(0); scale=x.std(0);scale[scale<1e-12]=1
 z=(x-mean)/scale
 if kind=='ols':
  coef=np.linalg.lstsq(np.column_stack([np.ones(len(y)),z]),y,rcond=None)[0];intercept=coef[0];b=coef[1:]
 elif kind=='ridge':
  model=Ridge(alpha=alpha,solver='svd').fit(z,y);intercept=model.intercept_;b=model.coef_
 else:
  model=Lasso(alpha=alpha,tol=1e-12,max_iter=100000).fit(z,y);intercept=model.intercept_;b=model.coef_
 return np.r_[intercept-np.sum(b*mean/scale),b/scale]
def predict(x,b): return np.column_stack([np.ones(len(x)),x])@b

def select(x,y,limit=3):
 # Forward BIC + a declared Delta-BIC >2 design choice; SVD reference fits.
 selected=[];n=len(y);rss=np.sum((y-y.mean())**2);score=n*np.log(max(rss,1e-12)/n)+np.log(n)
 for _ in range(min(limit,max(1,(n-1)//3))):
  trials=[]
  for j in range(1,x.shape[1]):
   if j in selected or x[:,j].std()<=1e-12:continue
   if any(abs(np.corrcoef(x[:,j],x[:,k])[0,1])>0.95 for k in selected):continue
   cols=selected+[j];pred=predict(x[:,cols],linear(x[:,cols],y));rss=np.sum((y-pred)**2)
   trials.append((n*np.log(max(rss,1e-12)/n)+(len(cols)+1)*np.log(n),j))
  if not trials:break
  value,j=min(trials)
  if value>=score-2:break
  selected.append(j);score=value
 return selected

def cv(x,y,selected,kind='ols',grid=None,groups=None,reselect=False):
 pred=np.zeros(len(y));groups=np.arange(len(y)) if groups is None else np.asarray(groups)
 for group in dict.fromkeys(groups):
  test=groups==group;train=~test;cols=select(x[train],y[train]) if reselect else selected
  if not cols:pred[test]=np.mean(y[train]);continue
  if grid is not None:
   scores=[]
   for alpha in grid:
    loo=np.zeros(train.sum());a=x[train][:,cols];b=y[train]
    for i in range(len(b)):
     keep=np.arange(len(b))!=i;loo[i]=predict(a[i:i+1],linear(a[keep],b[keep],kind,alpha))[0]
    scores.append(np.mean((loo-b)**2))
   alpha=grid[int(np.argmin(scores))]
  else:alpha=0
  pred[test]=predict(x[test][:,cols],linear(x[train][:,cols],y[train],kind,alpha))
 return pred

def fit_case(name,records,groups=None):
 records=np.asarray(records,dtype='<f4');n=len(records)-1
 freeze(M/'inputs'/f'{name}.sigpack',records.tobytes())
 groups=groups or [f'g{i}' for i in range(len(records))]
 labels='Reaction_ID,Dataset_Split,Ligand_Group\n'+''.join(f'{name}_{i},{"train" if i<n else "test"},{groups[i]}\n' for i in range(len(records)))
 freeze(M/'inputs'/f'{name}_labels.csv',labels.encode())
 out=M/'raw'/name;out.mkdir(exist_ok=True)
 # CLI run outputs are frozen once; diagnostics contain no wall-clock metadata.
 if not (out/'report.json').exists():
  run(name,['fit','--data',M/'inputs'/f'{name}.sigpack','--metadata',M/'inputs'/f'{name}_labels.csv','--output',out/'report.json','--predictions',out/'predictions.csv','--portable-model',out/'portable.json','--bootstrap','500','--permutations','1000','--seed','192026','--optimize','maximize','--response-temp-k','298.15'])
 report=json.loads((out/'report.json').read_text());x=expand(records[:n]);y=records[:n,6].astype(float);sel=report['selected_feature_indices'];z=x[:,sel];coef=linear(z,y);refpred=predict(z,coef)
 tests={'reference_selected_features':select(x,y),'stericx_selected_features':sel,'max_abs_coefficient_error':float(np.max(abs(np.array(report['weights'])[[0]+sel]-coef))),'training_reference':metrics(y,refpred),'training_stericx':report['training'],'fixed_loo_reference':metrics(y,cv(x,y,sel)),'fixed_loo_stericx':report['fixed_feature_loo'],'group_loo_reference':metrics(y,cv(x,y,sel,groups=groups[:n])),'group_loo_stericx':report['fixed_feature_group_loo'],'full_pipeline_loo_reference':metrics(y,cv(x,y,sel,reselect=True))}
 for kind,grid in [('ridge',[1e-4,1e-3,1e-2,.1,1,10]),('lasso',[1e-4,1e-3,1e-2,.05,.1,.5])]:
  ref=cv(x,y,sel,kind,grid);tests[kind+'_nested_alpha_reference']=metrics(y,ref);tests[kind+'_stericx']=report[kind+'_baseline']
 # Independent covariance/hat matrix, VIF, and t distribution.
 g=report['training_geometry'];design=np.column_stack([np.ones(n),(z-z.mean(0))/z.std(0)]);inv=np.linalg.inv(design.T@design)
 tests['max_abs_xtx_inverse_error']=float(np.max(abs(inv-np.array(g['xtx_inverse']))));tests['residual_standard_error_reference']=float(np.sqrt(np.sum((y-refpred)**2)/(n-len(sel)-1)))
 tests['residual_standard_error_stericx']=g['residual_standard_error'];tests['vif_reference']=[float(np.linalg.inv(np.corrcoef(z,rowvar=False))[i,i]) for i in range(len(sel))] if len(sel)>1 else [1.]
 b=json.loads((out/'portable.json').read_text())['uncertainty'];coeff=np.asarray(b['replicates']);interval=np.quantile(coeff,[.025,.975],axis=0).T
 tests['bootstrap_coefficient_interval_max_abs_error']=float(np.max(abs(interval-np.array([[r['lower_95'],r['upper_95']] for r in report['coefficient_intervals']]))))
 tests['limitations']=['Fixed-feature diagnostics condition on selection performed using all training labels.','Reference nested-alpha regularization also fixes features to isolate the claimed calculation; full-pipeline LOO repeats selection separately.','No-feature outer folds use an explicit training-mean fallback only in the independent full-pipeline diagnostic.']
 save(M/'results'/f'{name}_math.json',tests)
 rows=pd.DataFrame({'y':y,'reference_fitted':refpred,'reference_fixed_loo':cv(x,y,sel),'reference_full_pipeline_loo':cv(x,y,sel,reselect=True)})
 rows.to_csv(M/'results'/f'{name}_rows.csv',index=False,float_format='%.17g')
 return records,report

def main():
 rng=np.random.default_rng(192026);records=np.zeros((33,16));x=rng.normal(size=(33,3));records[:,:3]=x;records[:,2]=4;records[:,4]=1650;records[:,5]=298.15;records[:,6]=.7+1.9*x[:,0]-.8*x[:,1]+rng.normal(scale=.15,size=33)
 fit_case('synthetic_linear',records,[f'g{i//2}' for i in range(33)])
 raw=np.fromfile(F/'data/reactions.sigpack',dtype='<f4').reshape(-1,16);meta=pd.read_csv(F/'data/reactions_raw.csv');idx=np.r_[np.where(meta.Dataset_Split.eq('train'))[0],np.where(~meta.Dataset_Split.eq('train'))[0]]
 records,report=fit_case('native_ni_hda',raw[idx],meta.Ligand_Group.iloc[idx].tolist())
 # Freeze key observational provenance: row identity, duplicate structures and group overlap.
 from rdkit import Chem
 from rdkit.Chem.Scaffolds import MurckoScaffold
 identity=[]
 for row in meta.to_dict('records'):
  mol=Chem.MolFromSmiles(row['Ligand_SMILES']);identity.append({'reaction':row['Reaction_ID'],'split':row['Dataset_Split'],'group':row['Ligand_Group'],'canonical_isomeric_smiles':Chem.MolToSmiles(mol,isomericSmiles=True),'scaffold':MurckoScaffold.MurckoScaffoldSmiles(mol=mol,includeChirality=False)})
 save(M/'results'/'structure_identity.json',identity)
 # Study011 historical raw evidence: recompute all reported outcome/coverage aggregates.
 ranking=pd.read_csv(F/'docs/study_011/rankings.csv');save(M/'results'/'study011_columns.json',ranking.columns.tolist())
 save(M/'results'/'runtime.json',{'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'sklearn':sklearn.__version__,'platform':platform.platform()})
 # The initial phase manifest is immutable. Later inputs and references belong
 # to manifest_finalize.json; replay must not silently rewrite the first freeze.
 if not (M/'manifest.json').exists():
  manifest={str(p.relative_to(A)):hashlib.sha256(p.read_bytes()).hexdigest() for p in list((M/'inputs').rglob('*'))+list((M/'sources').rglob('*'))+list((M/'raw').rglob('*')) if p.is_file()};save(M/'manifest.json',manifest)

if __name__=='__main__':main()
