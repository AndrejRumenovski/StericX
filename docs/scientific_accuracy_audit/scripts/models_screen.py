"""Frozen CLI screening versus transparent independent arithmetic and SciPy."""
from models_audit import *
from scipy.spatial.distance import cdist

def main():
 model=M/'raw/synthetic_linear/portable.json'; d=json.loads(model.read_text());g=d['training_geometry'];w=np.asarray(d['weights'],dtype=np.float32).astype(float);sel=d['selected_feature_indices'];b=d['uncertainty'];beta=np.asarray(b['replicates']);means=np.asarray(g['means']);scales=np.asarray(g['scales']);train=np.asarray(g['standardized_training_points'])
 candidates=[('negative',-3.,0.),('high',2.,-1.),('tie_a',1.,.5),('tie_b',1.,.5),('middle',0.,0.),('extrapolated',8.,5.)]
 candidates += [('training_'+str(i),float(means[0]+z[0]*scales[0]),float(means[1]+z[1]*scales[1])) for i,z in enumerate(train[:3])]
 library=M/'inputs/screen_v2.csv';freeze(library,('Reaction_ID,sterimol_l,sterimol_b1,sterimol_b5\n'+''.join(f'{name},{l:.17g},{b1:.17g},4\n' for name,l,b1 in candidates)+'missing,,1,4\n').encode())
 for name,flags in [('default',[]),('ascending',['--ascending']),('descending',['--descending']),('diverse',['--diverse','--diversity-weight','.65']),('diverse_zero',['--diverse','--diversity-weight','0']),('domain_only',['--in-domain-only'])]:
  out=M/'raw'/('screen_v2_'+name)
  if not (out/'stdout').exists():run('screen_v2_'+name,['screen',model,library,'--format','json',*flags])
 reports={n:json.loads((M/'raw'/('screen_v2_'+n)/'stdout').read_text()) for n in ['default','ascending','descending','diverse','diverse_zero','domain_only']}
 # Locate independently with public output schema.
 report=reports['default'];save(M/'results'/'screen_schema.json',list(report))
 hits=report.get('hits',report.get('candidates',[]));rows=[];refs={}
 for name,l,b1 in candidates:
  row=np.array([1,l,b1,4,0,0,0,0],dtype=np.float32).astype(float);pred=float(row@w);z=(row[sel]-means)/scales;design=np.r_[1,z]
  h=float(design@np.linalg.inv(np.column_stack([np.ones(len(train)),train]).T@np.column_stack([np.ones(len(train)),train]))@design)
  t=stats.t.ppf(.975,g['observations']-g['parameters']);half=t*g['residual_standard_error']*np.sqrt(1+h);dist=float(cdist(z[None],train).min());mahal=float(np.sqrt(z@np.linalg.inv(np.cov(train,rowvar=False))@z));interval=np.quantile(beta@row[b['column_indices']],[.025,.975]);outside=any(row[c]<r['minimum'] or row[c]>r['maximum'] for c,r in zip(sel,d['applicability_domain']));verdict='extrapolation' if outside else ('sparse_interpolation' if dist>g['neighbor_calibration']['maximum'] else 'interpolation')
  refs[name]={'predicted':pred,'leverage':h,'pi_low':pred-half,'pi_high':pred+half,'nn':dist,'mahalanobis':mahal,'bootstrap_low':float(interval[0]),'bootstrap_high':float(interval[1]),'verdict':verdict}
 for hit in hits:
  name=hit.get('ligand_id',hit.get('ligand'));ref=refs[name];rows.append({'ligand':name,'stericx':hit,'independent':ref})
 save(M/'results'/'screen_rows.json',rows)
 save(M/'results'/'screen_ranking_reference.json',{'descending':sorted(refs,key=lambda k:(-refs[k]['predicted'],k)),'ascending':sorted(refs,key=lambda k:(refs[k]['predicted'],k)),'independent':refs})
 # Study011 is recalculated from individual frozen observations, not metrics.json.
 r=pd.read_csv(F/'docs/study_011/rankings.csv');err=r.predicted_ddg_kcal_mol-r.observed_ddg_kcal_mol;covered=(r.observed_ddg_kcal_mol>=r.prediction_interval_low)&(r.observed_ddg_kcal_mol<=r.prediction_interval_high)
 by={}
 for name,group in r.groupby('domain_verdict'):
  e=group.predicted_ddg_kcal_mol-group.observed_ddg_kcal_mol;c=(group.observed_ddg_kcal_mol>=group.prediction_interval_low)&(group.observed_ddg_kcal_mol<=group.prediction_interval_high);by[name]={'n':len(group),'mae':float(abs(e).mean()),'coverage':float(c.mean())}
 split=[]
 for key,group in r.groupby('split_id'):
  order=group.sort_values('predicted_rank');actual=group.sort_values('observed_ddg_kcal_mol',ascending=False);best=actual.iloc[0].reaction_id
  split.append({'split':key,'top1':int(order.iloc[0].reaction_id==best),'top2':len(set(order.head(2).reaction_id)&set(actual.head(2).reaction_id))/2,'best_rank':int(order.loc[order.reaction_id.eq(best),'predicted_rank'].iloc[0]),'rho':float(stats.spearmanr(group.predicted_ddg_kcal_mol,group.observed_ddg_kcal_mol).statistic)})
 save(M/'results'/'study011_recalculated.json',{'rows':len(r),'covered':int(covered.sum()),'coverage':float(covered.mean()),'mae':float(abs(err).mean()),'rmse':float(np.sqrt(np.mean(err**2))),'by_domain':by,'successful_splits':len(split),'failed_splits':57-len(split),'top1_all57':sum(s['top1'] for s in split)/57,'top2_all57':sum(s['top2'] for s in split)/57,'max_screen_fit_difference':float(abs(r.screen_minus_fit_prediction_kcal_mol).max()),'mean_rho_successful':float(np.mean([s['rho'] for s in split])),'split_metrics':split,'limitation':'Overlapping panels repeatedly predict the same11ligands; coverage is descriptive, not141independent calibrationobservations.'})

if __name__=='__main__':main()
