"""Observe actual outer-fold feature selection failures, without hiding them."""
from models_audit import *

def main():
 original=pd.read_csv(M/'inputs/native_ni_hda_labels.csv');refs=json.loads((M/'results/selection_folds.json').read_text());rows=[]
 for i in range(10):
  name=f'native_outer_{i}';labels=original.copy();labels.loc[i,'Dataset_Split']='test';path=M/'inputs'/f'{name}_labels.csv';freeze(path,labels.to_csv(index=False).encode());out=M/'raw'/name;out.mkdir(exist_ok=True)
  if not(out/'command.json').exists():run(name,['fit','--data',M/'inputs/native_ni_hda.sigpack','--metadata',path,'--output',out/'report.json','--predictions',out/'predictions.csv','--bootstrap','0','--permutations','0','--seed','192026'])
  cmd=json.loads((out/'command.json').read_text());row={'fold':i,'exit_code':cmd['exit_code'],'stderr':(out/'stderr').read_text(),'reference':refs[i]}
  if(out/'report.json').exists():
   r=json.loads((out/'report.json').read_text());p=pd.read_csv(out/'predictions.csv');row.update(selected_columns=r['selected_feature_indices'],frozen_predictions=p.to_dict('records'))
  rows.append(row)
 save(M/'results/outer_fold_observations.json',rows)

if __name__=='__main__':main()
