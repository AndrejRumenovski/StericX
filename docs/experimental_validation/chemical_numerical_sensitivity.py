"""Integrate certified regional precision into the same fixed experimental splits."""
import json,numpy as np,pandas as pd
from pathlib import Path
from analyze_chemistry import evaluate,paired_uncertainty
ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/'.stericx/experimental_validation';out=BASE/'chemical_analysis_v2';d=pd.read_csv(out/'primary_design.csv').set_index('ligand_id');p=pd.read_csv(out/'numerical_comparison.csv');ids=d.index.to_numpy();g=d.family.to_numpy();y=d.response.to_numpy();models={};features=[]
for mid in ids:
 rs=p[p.ligand_id==mid];first=rs.iloc[0];features.append(dict(ligand_id=int(mid),single=first.delta_continuous,single_lower=first.delta_lower,single_upper=first.delta_upper,ensemble_min=float(rs.delta_continuous.min()),ensemble_lower=float(rs.delta_lower.min()),ensemble_upper=float(rs.delta_upper.min())))
for name,col in [('B_certified_midpoint','single'),('C_certified_midpoint','ensemble_min'),('B_all_lower','single_lower'),('B_all_upper','single_upper'),('C_all_lower','ensemble_lower'),('C_all_upper','ensemble_upper')]:
 x=np.array([[r[col]] for r in features]);models[name]=evaluate(name,x,y,g,ids)
b=np.array([r['predicted'] for r in models['B_certified_midpoint']['predictions']]);c=np.array([r['predicted'] for r in models['C_certified_midpoint']['predictions']]);result=dict(models=models,features=features,paired_uncertainty=paired_uncertainty(y,b,c,g),scope='Same as-distributed molecules and family splits; only finite-grid versus explicit continuous regional descriptor changes. Endpoint fits at all-lower/all-upper inputs are sensitivity scenarios, NOT certified extrema of the fitted model over the descriptor box. The maximum/absolute-difference/minimum operations on per-conformer intervals are monotone interval operations; no volume certainty is transferred to molecular population or model uncertainty.')
(out/'certified_model_sensitivity.json').write_text(json.dumps(result,indent=2)+'\n');print({name:r['metrics']['rmse'] for name,r in models.items()})
