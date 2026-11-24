"""Frozen-design experimental analysis. Rust supplies every molecular descriptor.
Python independently checks reductions, extracts tables and fits statistical models.
All outer-fold predictions (including poor models) are retained; no outcome-driven
feature choice, geometry selection, or replacement of the registered endpoint.
"""
from pathlib import Path
import json,hashlib,math,itertools,warnings
import numpy as np,pandas as pd
from scipy import stats
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/'.stericx/experimental_validation';OUT=BASE/'chemical_analysis_v2'
FEATURES=['sterimol_l','sterimol_b1','sterimol_b5','percent_buried_volume','max_delta_qvbur','qvbur_min','qvbur_max','ovbur_min','ovbur_max','pyr_p','pyr_alpha']
FAMILIES={401:'cage',498:'spiro',723:'benzoxaphosphole',724:'benzoxaphosphole',1057:'benzoxaphosphole',1058:'benzoxaphosphole',785:'BINOL',2062:'BINOL',2063:'BINOL',2064:'BINOL',2067:'BINOL'}
EE={401:4,498:8,723:88,724:91,785:68,1057:79,1058:68,2062:24,2063:18,2064:3,2067:32,906:8}
YIELD={401:69,498:35,723:63,724:95,785:23,1057:48,1058:16,2062:43,2063:70,2064:61,2067:12,906:99}
R=0.0019872042586408316;T=353.15

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return [json.loads(x) for x in p.read_text().splitlines()]
def values(c):
 assert all(c[k]['status']=='success' for k in ['sterimol','buried_volume','pyramidalization'])
 return {**{'sterimol_'+k:v for k,v in c['sterimol']['value'].items()},**c['buried_volume']['value']['summary'],**c['pyramidalization']['value']}
def summary(mat):return {kind:getattr(np,func)(mat,axis=0) for kind,func in [('minimum','min'),('maximum','max'),('mean','mean'),('standard_deviation','std')]}
def metrics(y,p):
 e=p-y;return dict(n=len(y),mae=float(np.mean(abs(e))),rmse=float(np.sqrt(np.mean(e*e))),r2=None if np.var(y)==0 else float(1-np.sum(e*e)/np.sum((y-y.mean())**2)),bias=float(e.mean()))
def splits(groups):return [(np.where(groups!=g)[0],np.where(groups==g)[0],str(g)) for g in sorted(set(groups))]
def ols(x,y,xt):
 design=np.column_stack([np.ones(len(y)),x]);test=np.column_stack([np.ones(len(xt)),xt]);coef,_,rank,_=np.linalg.lstsq(design,y,rcond=None);pred=test@coef;df=len(y)-rank
 sigma2=float(np.sum((y-design@coef)**2)/df) if df>0 else float('nan');cov=np.linalg.pinv(design.T@design);width=stats.t.ppf(.975,df)*np.sqrt(sigma2*(1+np.einsum('ij,jk,ik->i',test,cov,test)))
 return pred,width,dict(coefficients=coef.tolist(),rank=int(rank),residual_df=int(df),condition_number=float(np.linalg.cond(design)))
def evaluate(label,x,y,groups,ids,method='ols'):
 pred=np.full(len(y),np.nan);lo=pred.copy();hi=pred.copy();folds=[]
 for train,test,g in splits(groups):
  row=dict(group=g,train_ids=ids[train].tolist(),test_ids=ids[test].tolist())
  try:
   if method=='ols':
    p,width,detail=ols(x[train],y[train],x[test]);lo[test]=p-width;hi[test]=p+width;row.update(detail)
   else:
    inner=[]
    for alpha in [.01,.1,1,10,100,1000]:
     errors=[]
     for tr,te,_ in splits(groups[train]):
      model=make_pipeline(StandardScaler(),Ridge(alpha=alpha));model.fit(x[train][tr],y[train][tr]);errors.extend((model.predict(x[train][te])-y[train][te])**2)
     inner.append(dict(alpha=alpha,mse=float(np.mean(errors))))
    best=min(inner,key=lambda r:r['mse'])['alpha'];model=make_pipeline(StandardScaler(),Ridge(alpha=best));model.fit(x[train],y[train]);p=model.predict(x[test]);row.update(inner_scores=inner,chosen_alpha=best)
   pred[test]=p;row['metrics']=metrics(y[test],p)
  except Exception as e:row['error']=repr(e)
  folds.append(row)
 result=dict(model=label,method=method,n_features=x.shape[1],parameter_count=x.shape[1]+1,folds=folds,predictions=[dict(ligand_id=int(i),family=str(g),observed=float(yv),predicted=None if not np.isfinite(pv) else float(pv),conditional_95_prediction_interval=None if not np.isfinite(lv) else [float(lv),float(hv)]) for i,g,yv,pv,lv,hv in zip(ids,groups,y,pred,lo,hi)],failures=[r for r in folds if 'error'in r])
 if np.all(np.isfinite(pred)):
  result['metrics']=metrics(y,pred);result['fold_rmse_sd']=float(np.std([r['metrics']['rmse'] for r in folds],ddof=1));result['fold_rmse_range']=[min(r['metrics']['rmse'] for r in folds),max(r['metrics']['rmse'] for r in folds)]
  if method=='ols':result['conditional_interval_coverage']=float(np.mean((y>=lo)&(y<=hi)));result['conditional_interval_mean_width']=float(np.mean(hi-lo))
 return result

def paired_uncertainty(y,b,c,groups):
 rng=np.random.default_rng(20261124);clusters=sorted(set(groups));draws=[];independent=[]
 for _ in range(10000):
  chosen=rng.choice(clusters,len(clusters),replace=True);ix=np.concatenate([np.where(groups==g)[0] for g in chosen]);draws.append(float(np.sqrt(np.mean((c[ix]-y[ix])**2))-np.sqrt(np.mean((b[ix]-y[ix])**2))))
  ix=rng.integers(0,len(y),len(y));independent.append(float(np.sqrt(np.mean((c[ix]-y[ix])**2))-np.sqrt(np.mean((b[ix]-y[ix])**2))))
 losses=np.array([np.sum((c[groups==g]-y[groups==g])**2-(b[groups==g]-y[groups==g])**2) for g in clusters]);observed=float(losses.sum());permuted=[float(np.dot(signs,losses)) for signs in itertools.product([-1,1],repeat=len(losses))]
 return dict(delta_rmse_c_minus_b=float(np.sqrt(np.mean((c-y)**2))-np.sqrt(np.mean((b-y)**2))),family_bootstrap_95_interval=np.quantile(draws,[.025,.975]).tolist(),ligand_bootstrap_95_interval=np.quantile(independent,[.025,.975]).tolist(),exact_signflip_one_sided_p=float(np.mean(np.array(permuted)<=observed+1e-14)),family_loss_differences=dict(zip(clusters,losses.tolist())),signflip_values=permuted,bootstrap_draws=10000,seed=20261124,interpretation='Conditional resampling of fixed cross-fitted predictions, not a full refit bootstrap. Four clusters, dependent CV training sets and tiny chemical sample limit uncertainty calibration; ligand bootstrap ignores family dependence.')
def main():
 OUT.mkdir(exist_ok=False);published=load(BASE/'chemistry_preparation/native.jsonl');generated=load(BASE/'generated_v2/native.jsonl');source=pd.read_csv(BASE/'sources/nihda.csv').set_index('Unnamed: 0');stereo=json.loads((BASE/'stereochemistry_v1/audit.json').read_text());pstatus={f"KRAKEN:{r['ligand_id']}:{r['conformer_id']}":r['source_p_comparison'] for r in stereo['records']}
 rows=[];ligands=[];reduction_failures=[];precision=[];ensembles={}
 for r in published+generated:
  assert 'error' not in r,r
  cs=r['result']['descriptors']['conformers'];native=r['result']['sensitivity']['statistics'];vals=[values(c) for c in cs];mat=np.array([[v[f] for f in FEATURES] for v in vals]);mid=int(r['ligand_id']);lane=r['lane'];ens=dict(id=r['id'],ligand_id=mid,lane=lane,n_conformers=len(cs),first_identifier=cs[0]['identifier'],features={})
  for j,f in enumerate(FEATURES):
   got=native['descriptors'][f];col=mat[:,j]
   for key,expected in [('minimum',min(col)),('maximum',max(col)),('mean',math.fsum(col)/len(col)),('standard_deviation',float(np.std(col)))]:
    if not math.isclose(got[key],expected,rel_tol=2e-13,abs_tol=2e-13):reduction_failures.append(dict(id=r['id'],feature=f,key=key,native=got[key],reference=expected))
   ens['features'][f]={k:got[k] for k in ['minimum','maximum','mean','standard_deviation']};ens['features'][f].update(first=float(col[0]),range=float(np.ptp(col)))
  for pop in r['result']['population_analyses']:
   temp=pop['temperature_k'];energies=np.array([c['energy']['value'] for c in cs]);w=np.exp(-(energies-energies.min())/(R*temp));w/=sum(w)
   for j,f in enumerate(FEATURES):
    got=pop['statistics']['descriptors'][f]['weighted_mean'];expect=float(w@mat[:,j]);assert math.isclose(got,expect,rel_tol=2e-12,abs_tol=2e-12),(f,got,expect);ens['features'][f]['weighted_'+str(temp)]=got
   ens.setdefault('energy_analysis',[]).append(dict(temperature_k=temp,energies_kcal_mol=energies.tolist(),relative_energies_kcal_mol=(energies-energies.min()).tolist(),weights=w.tolist(),effective_states=float(1/sum(w*w)),max_weight=float(max(w))))
  for c,v in zip(cs,vals):rows.append(dict(lane=lane,ligand_id=mid,identifier=c['identifier'],**v,energy_kcal_mol=None if c['energy'] is None else c['energy']['value'],source_p_comparison=pstatus.get(c['identifier'],'generated_separate_lane')))
  ligands.append(ens);ensembles[r['id']]=ens
  if lane=='published_dft_library':
   for c,v,p in zip(cs,vals,r['result']['precision']):
    assert 'error' not in p
    ref=p['result'];co=ref['continuous_sterimol_raw_l'];intervals=[];grid=[];total=[]
    byplane={x['plane_atom_index']:x for x in c['buried_volume']['value']['orientations']}
    for orient in ref['certified_regions']:
     cr=orient['volumes'];assert cr['target_met'];q=cr['quadrants'];total.append(cr['total']['volume_angstrom3']);plane=orient['plane_atom_index'];grid.extend(abs(a['volume_angstrom3']-b) for a,b in zip(q,byplane[plane]['quadrants']))
     for k in range(4):
      a,b=q[k],q[(k+1)%4];lo=a['lower_bound_angstrom3']-b['upper_bound_angstrom3'];hi=a['upper_bound_angstrom3']-b['lower_bound_angstrom3'];intervals.append((0 if lo<=0<=hi else min(abs(lo),abs(hi)),max(abs(lo),abs(hi))))
    lower=max(a for a,b in intervals);upper=max(b for a,b in intervals);estimate=(lower+upper)/2
    precision.append(dict(ligand_id=mid,identifier=c['identifier'],l_grid_minus_continuous=v['sterimol_l']-(co['l']+.4),b1_grid_minus_continuous=v['sterimol_b1']-co['b1'],b5_grid_minus_continuous=v['sterimol_b5']-co['b5'],delta_grid=v['max_delta_qvbur'],delta_continuous=estimate,delta_lower=lower,delta_upper=upper,delta_grid_minus_continuous=v['max_delta_qvbur']-estimate,max_quadrant_grid_difference=max(grid),mean_continuous_total_volume=float(np.mean(total))))
 assert not reduction_failures,reduction_failures
 pd.DataFrame(rows).to_csv(OUT/'conformer_descriptors.csv',index=False);pd.DataFrame(precision).to_csv(OUT/'numerical_comparison.csv',index=False);(OUT/'ensemble_statistics.json').write_text(json.dumps(ligands,indent=2)+'\n')
 obs=[]
 for mid in [723,1057,1058,370,724,401,785,906,2063,2064,2067,2062,498]:
  ee=EE.get(mid);obs.append(dict(ligand_id=mid,ee_percent=ee,yield_percent=YIELD.get(mid),apparent_abs_ddg_kcal_mol=None if ee is None else R*T*math.log((100+ee)/(100-ee)),author_csv_ddg=None if mid not in source.index else (None if pd.isna(source.loc[mid,'ddG_abs']) else float(source.loc[mid,'ddG_abs'])),eligible_primary=mid in FAMILIES,exclusion='unquantified ee/trace yield' if mid==370 else 'missing verified structure and author predictor row' if mid==906 else None,family=FAMILIES.get(mid),source='Cadge et al. 2025 SI Table S3, 353.15 K, GC-FID ee; yield is not modeled'))
 pd.DataFrame(obs).to_csv(OUT/'experimental_observations.csv',index=False)
 ids=np.array(sorted(FAMILIES));groups=np.array([FAMILIES[i] for i in ids]);y=np.array([R*T*math.log((100+EE[i])/(100-EE[i])) for i in ids]);raw={}
 for kind in ['first','minimum','maximum','mean']:
  raw[kind]=np.array([[ensembles[f'published_{i}']['features'][f][kind] for f in FEATURES] for i in ids])
 j=FEATURES.index('max_delta_qvbur');features={'A0_train_mean':np.zeros((len(ids),0)),'A_electronic':source.loc[ids,['nbo_P_boltz']].to_numpy(),'B_single':raw['first'][:,j:j+1],'C_ensemble_min':raw['minimum'][:,j:j+1],'D_published':source.loc[ids,['vbur_max_delta_qvbur_min']].to_numpy()};features['E_combined']=np.column_stack([features['D_published'],features['C_ensemble_min']]);results={}
 for label,target,keep,splitgroups in [('primary_as_distributed',y,np.ones(len(ids),bool),groups),('author_ddg',source.loc[ids,'ddG_abs'].to_numpy(),np.ones(len(ids),bool),groups),('exclude_2064',y,ids!=2064,groups),('ligand_wise',y,np.ones(len(ids),bool),ids.astype(str)),('exclude_stereo_401_1058',y,~np.isin(ids,[401,1058]),groups)]:
  result={name:evaluate(name,x[keep],target[keep],splitgroups[keep],ids[keep]) for name,x in features.items()}
  b=np.array([r['predicted'] for r in result['B_single']['predictions']]);c=np.array([r['predicted'] for r in result['C_ensemble_min']['predictions']]);result['paired_uncertainty']=paired_uncertainty(target[keep],b,c,splitgroups[keep]);results[label]=result
 # Filtered published lane: exclude 401; 1058's first matching source stereoconformer.
 filtered={k:v.copy() for k,v in features.items()};mid=1058;cs=next(r for r in published if r['ligand_id']==mid)['result']['descriptors']['conformers'];vv=[values(c)['max_delta_qvbur'] for c in cs if pstatus[c['identifier']]=='match'];ix=list(ids).index(mid);filtered['B_single'][ix,0]=vv[0];filtered['C_ensemble_min'][ix,0]=min(vv);filtered['E_combined'][ix,1]=min(vv);keep=ids!=401
 results['filtered_published_stereo']={name:evaluate(name,x[keep],y[keep],groups[keep],ids[keep]) for name,x in filtered.items()}
 # Secondary multivariate panels use inner family splits, scaling wholly inside folds.
 panels={'B_single_11':raw['first'],'C_mean_11':raw['mean'],'C_minmaxmean_33':np.column_stack([raw[k] for k in ['minimum','maximum','mean']])}
 results['nested_ridge']={name:evaluate(name,x,y,groups,ids,method='ridge') for name,x in panels.items()}
 # Generated lanes are distinct preparation/approximate-weighting sensitivities.
 results['generated_separate']={}
 for seed in [20261124,20261125]:
  chosen=[ensembles[f'generated_{i}_{seed}'] for i in ids];fs={kind:np.array([[r['features']['max_delta_qvbur'][kind]] for r in chosen]) for kind in ['first','minimum','mean','weighted_298.15','weighted_353.15','weighted_400.0']}
  results['generated_separate'][str(seed)]={name:evaluate(name,x,y,groups,ids) for name,x in fs.items()}
 # No inference is drawn from choosing the best of these correlations.
 corr=[]
 for name,x in {**features,**{f'single_{f}':raw['first'][:,i:i+1] for i,f in enumerate(FEATURES)},**{f'mean_{f}':raw['mean'][:,i:i+1] for i,f in enumerate(FEATURES)}}.items():
  if x.shape[1]==1:
   with warnings.catch_warnings():
    warnings.simplefilter('ignore');pear=stats.pearsonr(x[:,0],y);spear=stats.spearmanr(x[:,0],y)
   corr.append(dict(feature=name,pearson_r=None if not np.isfinite(pear.statistic) else float(pear.statistic),pearson_p_unadjusted=None if not np.isfinite(pear.pvalue) else float(pear.pvalue),spearman_rho=None if not np.isfinite(spear.statistic) else float(spear.statistic),spearman_p_unadjusted=None if not np.isfinite(spear.pvalue) else float(spear.pvalue)))
 pd.DataFrame(corr).to_csv(OUT/'descriptive_correlations.csv',index=False)
 design=pd.DataFrame(dict(ligand_id=ids,family=groups,response=y,**{k:v[:,0] for k,v in features.items() if v.shape[1]==1}));design.to_csv(OUT/'primary_design.csv',index=False)
 result=dict(design='STATISTICAL_PROTOCOL.md and pre-fit stereo amendment; primary as-distributed lane cannot be interpreted as uniformly verified experimental stereoisomers',results=results,random_forests_gradient_boosting='Omitted as pre-permitted: four scaffold clusters, outer training sets as small as six; multiple flexible learners would give poorly calibrated complexity comparisons. Nested ridge addresses multivariate regularization.',reduction_checks=dict(ligands=len(ligands),conformers=len(rows),features=len(FEATURES),failures=reduction_failures),inputs={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),BASE/'chemistry_preparation/native.jsonl',BASE/'generated_v2/native.jsonl',BASE/'stereochemistry_v1/audit.json',BASE/'sources/nihda.csv',ROOT/'docs/experimental_validation/STATISTICAL_PROTOCOL.md']})
 (OUT/'statistical_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({k:v['metrics'] for k,v in results['primary_as_distributed'].items() if k!='paired_uncertainty'},indent=2));print(json.dumps(results['primary_as_distributed']['paired_uncertainty'],indent=2))
if __name__=='__main__':main()
