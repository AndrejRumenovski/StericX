"""Additional model variants, screening objectives and training partition checks."""
from models_audit import *

def main():
 d=json.loads((M/'raw/synthetic_linear/portable.json').read_text());lib=M/'inputs/screen_v2.csv';rankings={}
 for direction in ['minimize','maximize_magnitude','unspecified']:
  clone=json.loads(json.dumps(d));clone['inference']['response']['optimization']=direction;p=M/'inputs'/f'objective_{direction}.json';freeze(p,(json.dumps(clone,indent=2)+'\n').encode());name='objective_'+direction
  if not(M/'raw'/name/'stdout').exists():run(name,['screen',p,lib,'--format','json'])
  try:response=json.loads((M/'raw'/name/'stdout').read_text());rankings[direction]=[(x['ligand'],x['predicted_ddg_kcal_mol']) for x in response['hits']]
  except (json.JSONDecodeError,KeyError):rankings[direction]={'error':(M/'raw'/name/'stderr').read_text()}
 exclusions=M/'inputs/tested.csv';freeze(exclusions,b'Reaction_ID\nhigh\nnot_present\nhigh\n');name='exclude_tested'
 if not(M/'raw'/name/'stdout').exists():run(name,['screen',M/'raw/synthetic_linear/portable.json',lib,'--exclude-tested',exclusions,'--format','json'])
 save(M/'results/screen_objectives.json',rankings)
 refs=json.loads((M/'results/screen_ranking_reference.json').read_text());base=json.loads((M/'raw/screen_v2_default/stdout').read_text());pred={h['ligand']:h['predicted_ddg_kcal_mol'] for h in base['hits']};rows=json.loads((M/'results/screen_rows.json').read_text());errors={k:0. for k in ['prediction','leverage','pi','nn','mahalanobis','bootstrap']}
 for r in rows:
  a,b=r['stericx'],r['independent'];errors['prediction']=max(errors['prediction'],abs(a['predicted_ddg_kcal_mol']-b['predicted']));errors['leverage']=max(errors['leverage'],abs(a['leverage']-b['leverage']));errors['pi']=max(errors['pi'],abs(a['prediction_interval_low']-b['pi_low']),abs(a['prediction_interval_high']-b['pi_high']));errors['nn']=max(errors['nn'],abs(a['nearest_training_distance']-b['nn']));errors['mahalanobis']=max(errors['mahalanobis'],abs(a['mahalanobis_distance']-b['mahalanobis']));errors['bootstrap']=max(errors['bootstrap'],abs(a['uncertainty']['lower']-b['bootstrap_low']),abs(a['uncertainty']['upper']-b['bootstrap_high']))
 checks={}
 for kind in ['ascending','descending','diverse','diverse_zero','domain_only']:
  out=json.loads((M/f'raw/screen_v2_{kind}/stdout').read_text());checks[kind]={'ranking':[h['ligand'] for h in out['hits']],'max_prediction_change':max(abs(h['predicted_ddg_kcal_mol']-pred[h['ligand']]) for h in out['hits'])}
 save(M/'results/screen_check_summary.json',{'max_absolute_math_errors':errors,'normal_ranking':[h['ligand'] for h in base['hits']],'independent_descending':refs['descending'],'independent_ascending':refs['ascending'],'variants':checks,'missing_exclusions':base['excluded']})
 # Holdout perturbation tests isolation of fitting/scaling/selection from frozen rows.
 arr=np.fromfile(M/'inputs/synthetic_linear.sigpack',dtype='<f4').reshape(-1,16).copy();arr[-1,:7]=[1e6,-1e6,1e6,100,1e6,100,999];p=M/'inputs/holdout_perturbed.sigpack';freeze(p,arr.tobytes());name='holdout_perturbed';out=M/'raw'/name;out.mkdir(exist_ok=True)
 if not(out/'report.json').exists():run(name,['fit','--data',p,'--metadata',M/'inputs/synthetic_linear_labels.csv','--output',out/'report.json','--predictions',out/'predictions.csv','--bootstrap','500','--permutations','1000','--seed','192026'])
 original=json.loads((M/'raw/synthetic_linear/report.json').read_text());changed=json.loads((out/'report.json').read_text());save(M/'results/holdout_partition.json',{'training_report_identical':original==changed,'changed_fields':[k for k in original if original[k]!=changed.get(k)],'perturbation':'Onlyheld-out row descriptors,electronics,temperature,target changed byextremes. No training row or label changes.'})
 # Recompute historical model variants from their frozen perligand descriptors.
 source=pd.read_csv(M/'sources/ni_hda_reaction_data.csv',index_col=0);ids=[401,498,724,785,1057,1058,2062,2063,2064,2067];results={}
 for study,file,modelname in [('002','official_kraken_comparison.csv','native_model.json'),('003','official_kraken_lmo_comparison.csv','quantum_model.json')]:
  table=pd.read_csv(F/f'docs/study_{study}'/file).set_index('Source_ID');x=table.loc[ids,'max_delta_qvbur_min'].to_numpy();y=source.loc[ids,'ddG_abs'].to_numpy();coef=linear(x[:,None],y);pred=predict(x[:,None],coef);loo=np.array([predict(x[i:i+1,None],linear(np.delete(x,i)[:,None],np.delete(y,i)))[0] for i in range(len(y))]);stored=json.loads((F/f'docs/study_{study}'/modelname).read_text());hold=float(predict(table.loc[[723],['max_delta_qvbur_min']].to_numpy(),coef)[0]);results[study]={'training':metrics(y,pred),'fixed_loo':metrics(y,loo),'coefficient_reference':coef.tolist(),'coefficient_stored':[stored['intercept'],stored['slope']],'holdout723':hold,'holdout_absolute_error':abs(hold-source.loc[723,'ddG_abs']),'limitation':'Reference refits use frozen rounded perligand descriptors. Doesnotrevalidate theirgeometrygeneration; treat all723replays afterStudy001 ashistorical.'};pd.DataFrame({'source_id':ids,'x':x,'actual':y,'fitted':pred,'loo':loo}).to_csv(M/f'results/study_{study}_model_rows.csv',index=False,float_format='%.17g')
 save(M/'results/reaction_variants.json',results)

if __name__=='__main__':main()
