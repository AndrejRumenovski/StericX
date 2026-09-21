"""Preserved statistical falsification witnesses; no production changes."""
from models_audit import *

def main():
 records=np.fromfile(M/'inputs/native_ni_hda.sigpack',dtype='<f4').reshape(-1,16);x=expand(records[:-1]);y=records[:-1,6].astype(float);sel=select(x,y);folds=[]
 for i in range(len(y)):
  keep=np.arange(len(y))!=i;chosen=select(x[keep],y[keep]);p=predict(x[i:i+1,chosen],linear(x[keep][:,chosen],y[keep]))[0] if chosen else y[keep].mean();folds.append({'held_out_index':i,'selected_columns':chosen,'fallback_train_mean':not bool(chosen),'actual':y[i],'prediction':float(p)})
 save(M/'results/selection_folds.json',folds)
 rng=np.random.default_rng(192026);observed=metrics(y,predict(x[:,sel],linear(x[:,sel],y)))['r2'];fixed=[];pipeline=[]
 for _ in range(2000):
  yp=rng.permutation(y);f=metrics(yp,predict(x[:,sel],linear(x[:,sel],yp)))['r2'];s=select(x,yp);p=metrics(yp,predict(x[:,s],linear(x[:,s],yp)))['r2'] if s else 0.;fixed.append(f);pipeline.append(p)
 save(M/'results/permutation_selection.json',{'observed_r2':observed,'samples':2000,'seed':192026,'fixed_features_p':(1+sum(f>=observed for f in fixed))/2001,'repeat_selection_p':(1+sum(f>=observed for f in pipeline))/2001,'fixed_null_r2':fixed,'pipeline_null_r2':pipeline,'note':'Independent NumPy permutations; MonteCarlo p-values use add-one correction. Fixed and pipeline share each permutation; not a reproduction of StericX RNG.'})
 # A constructed bootstrap ensemble disproves a distribution-free guarantee
 # that an interval-arithmetic box of marginal95%CIs is conservative.
 d=json.loads((M/'raw/synthetic_linear/portable.json').read_text());indices=d['uncertainty']['column_indices'];base=np.asarray(d['weights'],dtype=np.float32).astype(float)[indices];reps=np.tile(base,(1000,1))
 for start,column,sign in [(0,0,1),(24,0,-1),(48,1,1),(72,1,-1)]:reps[start:start+24,column]+=sign*100
 d['uncertainty']['replicates']=reps.tolist();d['uncertainty']['replicate_count']=1000;d['uncertainty']['requested_samples']=1000
 q=np.quantile(reps,[.025,.975],axis=0).T
 for row,bounds in zip(d['coefficient_intervals'],q):row['lower_95']=float(bounds[0]);row['upper_95']=float(bounds[1])
 inp=M/'inputs/marginal_ci_counterexample_model.json';freeze(inp,(json.dumps(d,indent=2)+'\n').encode());lib=M/'inputs/marginal_ci_counterexample.csv';freeze(lib,b'Reaction_ID,sterimol_l,sterimol_b1\nwitness,1,0\n');out=M/'raw/marginal_ci_counterexample'
 if not(out/'stdout').exists():run('marginal_ci_counterexample',['screen',inp,lib,'--format','json'])
 pred=reps[:,0]+reps[:,1];band=[sum(q[:2,0]),sum(q[:2,1])]
 save(M/'results/marginal_ci_counterexample.json',{'input':str(inp.relative_to(A)),'independent_marginal_band':band,'independent_joint_propagation_quantiles':np.quantile(pred,[.025,.975]).tolist(),'empirical_coverage_of_marginal_band':float(np.mean((pred>=band[0])&(pred<=band[1]))),'interpretation':'A mathematical counterexample for the unqualified conservative label. Artificial serialized bootstrap population, not an observed chemistry dataset; StericX mean-response bootstrap interval correctly propagates the joint replicates.'})

if __name__=='__main__':main()
