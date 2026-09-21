#!/usr/bin/env python3
"""Conservative analytic rounding bound for P=|det(unit bond directions)|."""
import itertools,json,math,datetime
import numpy as np
from pathlib import Path
from kraken_analyze import ROOT,save,sha
out=ROOT/'focused_diagnosis';ids={v['id'] for v in json.loads((out/'pyramidalization_coordinate_rounding.json').read_text())};rows=[]
for r in map(json.loads,(ROOT/'prepared_all/requests.jsonl').open()):
 if r['id'] not in ids:continue
 c=np.asarray([v['position'] for v in r['atoms']]);v=c[r['neighbors']]-c[r['donor']];norm=np.linalg.norm(v,axis=1);u=v/norm[:,None];eps=2*math.sqrt(3)*.00005;delta=2*eps/(norm-eps)
 checks=[]
 for i,j,k in [(0,1,2),(0,2,1),(1,2,0)]:
  value=float((u[i]+u[j])@u[k]);bound=float(delta[i]+delta[j]+2*delta[k]);cross=float(np.linalg.norm(np.cross(u[i],u[j])));checks.append({'pair':[i,j],'apex':k,'bisector_sign_numerator':value,'perturbation_bound':bound,'positive_alpha_guaranteed':value+bound<0,'pair_cross_norm':cross,'cross_norm_lower_bound':cross-float(delta[i]+delta[j])})
 valid=all(v['positive_alpha_guaranteed'] and v['cross_norm_lower_bound']>0 for v in checks)
 rows.append({'id':r['id'],'coordinate_rounding_half_width_A':.00005,'bond_vector_perturbation_bound_A':eps,'bond_lengths_A':norm.tolist(),'unit_vector_perturbation_bounds':delta.tolist(),'unperturbed_absolute_triple_product':float(abs(np.linalg.det(u))),'P_absolute_error_upper_bound':float(sum(delta)) if valid else None,'branch_checks':checks,'valid_positive_alpha_branch_bound':valid})
save(out/'pyramidalization_rounding_bound.json',{'created_utc':datetime.datetime.now(datetime.UTC).isoformat(),'script_sha256':sha(Path(__file__)),'input_sha256':sha(ROOT/'prepared_all/requests.jsonl'),'proof':'Each bond perturbation norm≤epsilon=2sqrt(3)delta. Unit vector change≤2epsilon/(bondlength−epsilon). Telescoping determinant with unit vectors bounds |Δabsdet|≤sum(unit-vector bounds). Signs of each bisector numerator remain negative by inner-product perturbation bound, and paircrossnorm stayspositive, so all signed alphas staypositive and P=absdet branch is fixed. This is a conservative mathematical bound for arbitrary continuous perturbations in the rounding box, not only corners. Assumes ordinary nearest rounding to4decimals and same atom/ensemble identity.','rows':rows})
print('N',len(rows),'branch guaranteed',sum(v['valid_positive_alpha_branch_bound'] for v in rows),'maxbound',max(v['P_absolute_error_upper_bound'] or 0 for v in rows))
for v in rows:
 if v['id'].split(':')[1]=='1290':print(v['id'],v['P_absolute_error_upper_bound'])
