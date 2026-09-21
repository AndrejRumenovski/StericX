"""Minimize two real Kraken alignment failures without changing SUT algorithms."""
import hashlib,json,subprocess
from pathlib import Path
import numpy as np
from geometry_reference import sterimol_exact
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'geometry/alignment';OUT.mkdir(exist_ok=True)
source=ROOT/'kraken/morfeus_outliers_attempt1/selected_inputs.jsonl'
wanted={'KRAKEN:963:48394':'l','KRAKEN:1075:49864':'b5'};requests=[];notes=[]
for line in source.read_text().splitlines():
 pair=json.loads(line);c=pair['request'];s=pair['sut'];id=c['id']
 if id not in wanted:continue
 xyz=np.array([a['position'] for a in s['atoms']]);rs=np.array([a['vdw_radius'] for a in s['atoms']]);center=np.array(s['center']);u=xyz[c['donor']]-center;u/=np.linalg.norm(u);rel=xyz-center
 axial=rel@u;transverse=rel-np.outer(axial,u)
 scores=axial+rs if wanted[id]=='l' else np.linalg.norm(transverse,axis=1)+rs
 active=int(np.argmax(scores));keep=[c['donor'],active]
 full=c|{'id':id+'__controlled','center':s['center'],'sterimol_only':True,'atoms':[{'element':a['element'],'position':a['position'],'radius':a['vdw_radius']} for a in s['atoms']]}
 minimum=full|{'id':id+'__minimal','donor':0,'reference':1,'neighbors':[1,1,1],'atoms':[full['atoms'][i] for i in keep]}
 requests.extend([full,minimum]);notes.append({'id':id,'descriptor':wanted[id],'active_atom_original':active,'retained_atom_indices':keep,'unit_direction':u.tolist(),'dot_with_positive_z':float(u[2]),'angle_from_positive_z_degrees':float(np.rad2deg(np.arccos(u[2])))})
req=OUT/'requests.jsonl';blob=''.join(json.dumps(r)+'\n' for r in requests)
if req.exists():assert req.read_text()==blob
else:req.write_text(blob)
output=OUT/'sut_outputs.jsonl';binary=ROOT/'frozen/bin/stericx-audit-observer'
if not output.exists():
 with req.open('rb') as src,output.open('wb') as dst,(OUT/'stderr.txt').open('wb') as err:subprocess.run([str(binary)],stdin=src,stdout=dst,stderr=err,check=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
(OUT/'freeze.json').write_text(json.dumps({'input_source_sha256':sha(source),'requests_sha256':sha(req),'sut_outputs_sha256':sha(output),'observer_sha256':sha(binary)},indent=2)+'\n')
sut={r['id']:r for r in map(json.loads,output.read_text().splitlines())};refs=[]
for c in requests:
 s=sut[c['id']];xyz=np.array([a['position'] for a in s['atoms']]);rs=np.array([a['vdw_radius'] for a in s['atoms']]);ref=sterimol_exact(xyz,rs,np.array(s['center']),c['donor'])
 refs.append({'id':c['id'],'independent':ref,'sut':s['sterimol_dummy_raw'],'errors':{k:s['sterimol_dummy_raw'][k]-ref[k] for k in ref}})
(OUT/'reference_results.json').write_text(json.dumps({'case_notes':notes,'comparison':refs,'operation':'glam0.30.10 Quat::from_rotation_arc snaps within abs(dot)>1-2*f32EPSILON; source documents about0.001 near-singularity accuracy. This approximate axis feeds axial/radial projections.'},indent=2)+'\n')
print(json.dumps(refs,indent=2))
