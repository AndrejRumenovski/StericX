"""Disclosed post-primary identity sensitivity, retaining the registered analysis."""
from pathlib import Path
import json,numpy as np,pandas as pd
from analyze_chemistry import evaluate,paired_uncertainty,values
ROOT=Path(__file__).resolve().parents[2];B=ROOT/'.stericx/experimental_validation';d=pd.read_csv(B/'chemical_analysis_v2/primary_design.csv').set_index('ligand_id');ids=d.index.to_numpy();groups=d.family.to_numpy();y=d.response.to_numpy();features={name:(np.zeros((len(ids),0)) if name=='A0_train_mean' else d[[name]].to_numpy(copy=True)) for name in ['A0_train_mean','A_electronic','B_single','C_ensemble_min','D_published']};features['E_combined']=np.column_stack([features['D_published'],features['C_ensemble_min']]);result={}
keep=~np.isin(ids,[401,1058,2063]);result['exclude_401_1058_2063']={n:evaluate(n,x[keep],y[keep],groups[keep],ids[keep]) for n,x in features.items()}
pstatus={f"KRAKEN:{r['ligand_id']}:{r['conformer_id']}":r['source_p_comparison'] for r in json.loads((B/'stereochemistry_v1/audit.json').read_text())['records']};axial=next(r for r in json.loads((B/'conformer_sensitivity_v2/summary.json').read_text())['published_relative_axial_stereo'] if r['ligand_id']==2063);allowed={r['conformer_id'] for r in axial['records'] if np.sign(np.sin(r['dihedrals_rad'][0]))==np.sign(np.sin(axial['records'][0]['dihedrals_rad'][0]))}
for r in map(json.loads,(B/'chemistry_preparation/native.jsonl').read_text().splitlines()):
 mid=r['ligand_id']
 if mid not in [1058,2063]:continue
 cs=r['result']['descriptors']['conformers'];cs=[c for c in cs if pstatus[c['identifier']]=='match'] if mid==1058 else [c for c in cs if int(c['identifier'].split(':')[-1]) in allowed];v=[values(c)['max_delta_qvbur'] for c in cs];ix=list(ids).index(mid);features['B_single'][ix,0]=v[0];features['C_ensemble_min'][ix,0]=min(v);features['E_combined'][ix,1]=min(v)
keep=ids!=401;result['filtered_P_and_relative_axial']={n:evaluate(n,x[keep],y[keep],groups[keep],ids[keep]) for n,x in features.items()}
for lane,rs in result.items():
 b=np.array([r['predicted'] for r in rs['B_single']['predictions']]);c=np.array([r['predicted'] for r in rs['C_ensemble_min']['predictions']]);sub=[r['ligand_id'] for r in rs['B_single']['predictions']];rs['paired_uncertainty']=paired_uncertainty(d.loc[sub,'response'].to_numpy(),b,c,d.loc[sub,'family'].to_numpy())
(B/'chemical_analysis_v2/post_stereo_sensitivity.json').write_text(json.dumps(dict(results=result,note='2063 includes one opposite relative BINOL-axis sector. Its primary scalar features and absolute ee are invariant to a whole-molecule reflection; excluding it is conservative and does not establish that the remaining absolute atropisomer labels are experimentally verified. This post-primary sensitivity does not replace the registered analysis.'),indent=2)+'\n')
print({k:{m:v['metrics']['rmse'] for m,v in rs.items() if 'metrics'in v} for k,rs in result.items()})
