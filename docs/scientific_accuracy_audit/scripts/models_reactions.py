"""Primary-data reaction-model audit, independent fits and threshold arithmetic."""
from models_audit import *
import re
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import matthews_corrcoef, f1_score

def ni_hda():
 d=pd.read_csv(M/'sources/ni_hda_reaction_data.csv',index_col=0);x=pd.read_csv(M/'sources/ni_hda_kraken_features.csv',index_col=0);ids=[401,498,724,785,1057,1058,2062,2063,2064,2067];term='vbur_max_delta_qvbur_min';a=x.loc[ids,term].to_numpy();y=d.loc[ids,'ddG_abs'].to_numpy();coef=linear(a[:,None],y);p=predict(a[:,None],coef);loo=np.array([predict(a[i:i+1,None],linear(np.delete(a,i)[:,None],np.delete(y,i)))[0] for i in range(len(y))]);stored=json.loads((F/'docs/study_001/published_model.json').read_text());hold=float(predict(x.loc[[723],[term]].to_numpy(),coef)[0]);rows=pd.DataFrame({'source_id':ids,'x':a,'actual':y,'fitted':p,'loo':loo});rows.to_csv(M/'results/published_ni_hda_rows.csv',index=False,float_format='%.17g');save(M/'results/published_ni_hda.json',{'training':metrics(y,p),'fixed_feature_loo':metrics(y,loo),'coefficients_reference':coef.tolist(),'coefficients_stericx_study':[stored['intercept'],stored['slope']],'holdout_723_reference':hold,'holdout_723_actual':float(d.loc[723,'ddG_abs']),'absolute_holdout_error':abs(hold-d.loc[723,'ddG_abs']),'training_ids':ids,'feature':term,'primary_notebook_feature_search_count':len(x.columns)-1,'limitation':'The original authors searched the same10labels across descriptors. Treat historical fixedfeatureLOO as conditional diagnostic; this is not independent validation of that original search.'})

def read_tables():
 text=(M/'sources/crosscoupling_preprint_si.txt').read_text();headers=list(re.finditer(r'Table S(\d+)\. (?:Compiled|Tabulated)[^\n]*',text));blocks={int(h.group(1)):text[h.end():headers[i+1].start() if i+1<len(headers) else len(text)] for i,h in enumerate(headers)}
 rows=[];rejected=[]
 for table,reactions,tail in [(1,['I'],5),(2,['II','III','IV','V','RS1'],9),(9,['VII'],3),(10,['VIII'],3),(11,['IX'],3),(12,['X'],3),(15,['XI'],3),(16,['XII'],3)]:
  seen=set()
  for line in blocks[table].splitlines():
   fields=line.split()
   # Some ligand names wrap entirely above/below the numeric row. The numeric
   # row then contains only ID + values; requiring a name silently loses data.
   if len(fields)<tail+1 or not fields[0].isdigit():continue
   try:nums=list(map(float,fields[-tail:]))
   except ValueError:continue
   ident=int(fields[0]);vmin=nums[3] if tail in [5,9] else nums[1];ys=nums[4:] if tail in [5,9] else nums[2:]
   if not(10<=vmin<=80 and all(0<=v<=100 for v in ys)):
    rejected.append({'table':table,'line':line,'reason':'numeric_nondata_or_invalid_value'});continue
   if ident in seen:
    rejected.append({'table':table,'line':line,'reason':'duplicateID_in_extractedblock'});continue
   seen.add(ident)
   for reaction,y in zip(reactions,ys):rows.append({'reaction':reaction,'id':ident,'name':' '.join(fields[1:-tail]),'vbur_min_published_rounded':vmin,'yield':y,'table':table,'line':line})
 freeze(M/'inputs/crosscoupling_preprint_rows_v2.json',(json.dumps(rows,indent=2)+'\n').encode());save(M/'results/crosscoupling_parse_exceptions.json',rejected);return pd.DataFrame(rows)

def stump(x,y,balanced):
 # Enumerate all distinct adjacent cutpoints and minimize weighted Gini risk.
 # No sklearn or StericX tree is called for this reference calculation.
 weight=np.where(y, len(y)/(2*y.sum()) if balanced else 20.,len(y)/(2*(len(y)-y.sum())) if balanced else 1.)
 unique=np.unique(x);trials=[]
 for threshold in (unique[1:]+unique[:-1])/2:
  left=x<=threshold;loss=0;labels=[]
  for mask in [left,~left]:
   masses=np.bincount(y[mask].astype(int),weights=weight[mask],minlength=2);labels.append(int(np.argmax(masses)));total=masses.sum();loss+=total*(1-np.sum((masses/total)**2))
  pred=np.where(left,labels[0],labels[1]);trials.append((loss,float(threshold),pred,labels))
 loss,threshold,pred,labels=min(trials,key=lambda z:(z[0],z[1]))
 return {'threshold':threshold,'direction':'Right' if labels==[0,1] else 'Left' if labels==[1,0] else 'constant','accuracy':float(np.mean(pred==y)),'mcc':float(matthews_corrcoef(y,pred)),'f1':float(f1_score(y,pred)),'predictions':pred.tolist(),'gini_loss':float(loss)}

def crosscoupling():
 rows=read_tables();primary=pd.read_csv(M/'sources/ni_hda_kraken_features.csv',index_col=0);cuts={'I':10,'II':10,'III':5,'IV':5,'V':20,'RS1':10,'VII':10,'VIII':10,'IX':10,'X':5,'XI':30,'XII':10};summary=[];predrows=[]
 native_path=A/'kraken/analysis/percent_buried_volume_min_by_ligand.csv';freeze(M/'inputs/current_native_percent_buried_volume_min.csv',native_path.read_bytes());native=pd.read_csv(native_path,index_col='molecule_id')
 for name,r in rows.groupby('reaction',sort=False):
  # SI descriptor values are rounded to0.1. Retain that primary input; separately
  # compare the public fullprecisionKraken descriptor, without conflatingversions.
  for descriptor in ['si_rounded','public_kraken_fullprecision','current_stericx_fullprecision']:
   if descriptor=='si_rounded':x=r.vbur_min_published_rounded.to_numpy();keep=np.ones(len(r),bool)
   elif descriptor=='public_kraken_fullprecision':
    keep=r.id.isin(primary.index).to_numpy();x=100*primary.loc[r.id[keep],'vbur_vbur_min'].to_numpy()/(4*np.pi*3.5**3/3)
   else:
    keep=r.id.isin(native.index).to_numpy();x=native.loc[r.id[keep],'sut'].to_numpy()
   rr=r[keep];balanced=name in ['VII','VIII','IX','X','XI','XII'];cut=cuts[name]
   for convention in ['official_ge','stericx_study_gt']:
    y=(rr['yield'].to_numpy()>=cut if convention=='official_ge' else rr['yield'].to_numpy()>cut).astype(int);ref=stump(x,y,balanced);sk=DecisionTreeClassifier(max_depth=1,class_weight='balanced' if balanced else {0:1,1:20},random_state=0).fit(x[:,None],y);s=sk.predict(x[:,None]);summary.append({'reaction':name,'descriptor':descriptor,'label_convention':convention,'cutoff':cut,'n':len(y),'n_active':int(y.sum()),'boundary_count':int((rr['yield']==cut).sum()),'reference':{k:v for k,v in ref.items() if k!='predictions'},'sklearn_threshold':float(sk.tree_.threshold[0]),'sklearn_accuracy':float(np.mean(s==y)),'sklearn_mcc':float(matthews_corrcoef(y,s)),'reference_sklearn_prediction_disagreements':int(np.sum(s!=ref['predictions'])),'omitted_missing_reference_ids':r.id[~keep].tolist()})
    for i,(_,row) in enumerate(rr.iterrows()):predrows.append({'reaction':name,'id':int(row.id),'descriptor':descriptor,'label_convention':convention,'yield':row['yield'],'feature':float(x[i]),'actual_label':int(y[i]),'reference_prediction':int(ref['predictions'][i]),'sklearn_prediction':int(s[i])})
 save(M/'results/crosscoupling_reference_models.json',summary);pd.DataFrame(predrows).to_csv(M/'results/crosscoupling_reference_rows.csv',index=False,float_format='%.17g')

if __name__=='__main__':ni_hda();crosscoupling()
