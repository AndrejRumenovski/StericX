"""Frozen real conformers, radii/center profiles and independent box cross-checks."""
import argparse,json,re,struct,subprocess,sys
from pathlib import Path
from campaign import ROOT,EVIDENCE,sha,save
sys.path.insert(0,str(ROOT/'scripts'))
from dbstep_resolution_reference import certified_union_volume
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--engine',type=Path,required=True);p.add_argument('--reference-engine',type=Path,required=True);p.add_argument('--output',required=True);a=p.parse_args()
out=EVIDENCE/a.output;out.mkdir(exist_ok=False)
fixtures_path=EVIDENCE/'baseline_observations/geometry/stdout.jsonl'
fixtures=list(map(json.loads,fixtures_path.read_text().splitlines()))
radii={}
# Adversarial fixtures can override radii. Only the first twelve chemical
# fixtures establish the configured element table below.
for f in fixtures[:12]:
    for atom in f['atoms']:
        e=atom['element'].upper();r=atom['vdw_radius'];assert radii.setdefault(e,r)==r
f32=lambda x:struct.unpack('f',struct.pack('f',float(x)))[0]
radii_source=EVIDENCE/'baseline/source/src/geometry/xyz.rs'
table=radii_source.read_text().split('fn van_der_waals_radius(')[1].split('pub fn covalent_radius(')[0]
for element,value in re.findall(r'"([A-Za-z]+)"\s*=>\s*([0-9.]+)',table):
    configured=f32(value);assert radii.setdefault(element.upper(),configured)==configured
catalog=ROOT/'docs/benchmarks/morfeus_current/input_manifest.json'
manifest=json.loads(catalog.read_text())
probe=ROOT/'.stericx/ultimate_engine/competitive_gate/probe.stdout'
frames=list(map(json.loads,probe.read_text().splitlines()))
requests=[];sources={str(catalog):sha(catalog),str(probe):sha(probe),str(fixtures_path):sha(fixtures_path),str(radii_source):sha(radii_source)}
for i,(row,frame) in enumerate(zip(manifest['corpus'],frames,strict=True)):
    path=catalog.parent/row['filename'];assert sha(path)==row['sha256'];sources[str(path)]=sha(path)
    lines=path.read_text().splitlines();count=int(lines[0]);all_atoms=[]
    for line in lines[2:2+count]:
        element,*xyz=line.split();element=element.upper()
        assert element in radii,(path,element)
        all_atoms.append({'center':[f32(x) for x in xyz[:3]],'radius':radii[element]*1.1699999570846558,'element':element})
    atoms=[{k:v for k,v in atom.items() if k!='element'} for atom in all_atoms if atom['element']!='H']
    for radius in [2.5,3.5,5.0]:
        requests.append({'id':f'{row["ligand"]}/{row["conformer"]}/R{radius}','region':{'center':frame['center'],'radius':radius},'atoms':atoms,'tolerance':.05,'budget':1048575})
    if i%4==0:
        requests.append({'id':f'conformer_{i}/explicit_shifted_center','region':{'center':[x+d for x,d in zip(frame['center'],[.125,-.25,.5])],'radius':3.5},'atoms':atoms,'tolerance':.05,'budget':1048575})
        requests.append({'id':f'conformer_{i}/include_hydrogens','region':{'center':frame['center'],'radius':3.5},'atoms':[{k:v for k,v in atom.items() if k!='element'} for atom in all_atoms],'tolerance':.05,'budget':1048575})
for f in fixtures[:12]:
    for radius in [2.5,5.0]:
        requests.append({'id':f['id']+f'/unscaled_R{radius}','region':{'center':f['center'],'radius':radius},'atoms':[{'center':atom['position'],'radius':atom['vdw_radius']} for atom in f['atoms'] if atom['element']!='H'],'tolerance':.05,'budget':1048575})
save(out/'requests.json',requests)
with (out/'observed.jsonl').open('x') as log,(out/'stderr').open('x') as err:
    subprocess.run([str(a.engine.resolve()),str(out/'requests.json')],stdout=log,stderr=err,check=True)
observed=list(map(json.loads,(out/'observed.jsonl').read_text().splitlines()))
reference_requests=[{'id':q['id'],'region':q['region'],'atoms':q['atoms'],'method':{'method':'reference','absolute_tolerance':.05,'max_boxes':65535},'samples':0} for q in requests]
save(out/'reference_requests.json',reference_requests)
with (out/'reference.jsonl').open('x') as log,(out/'reference.stderr').open('x') as err:
    subprocess.run([str(a.reference_engine.resolve()),str(out/'reference_requests.json')],stdout=log,stderr=err,check=True)
reference=list(map(json.loads,(out/'reference.jsonl').read_text().splitlines()));failures=[];rows=[]
for q,r,b in zip(requests,observed,reference,strict=True):
    old=b['result'];common=[max(r['lower'],old['lower_bound_angstrom3']),min(r['upper'],old['upper_bound_angstrom3'])]
    if common[0]>common[1]:failures.append([q['id'],'independent_method_disjoint'])
    if r['target_met']!=(r['absolute_error_bound']<=q['tolerance']):failures.append([q['id'],'target_flag'])
    if r['cells']>q['budget']:failures.append([q['id'],'budget'])
    rows.append({'id':q['id'],'new':r,'old':old,'common':common})
# A separate exact-predicate/rational-volume Python octree on eight preselected cases.
independent=[]
for i in [0,7,19,35,57,91,137,195]:
    q=requests[i];r=observed[i]
    ref=certified_union_volume([s['center'] for s in q['atoms']],[s['radius'] for s in q['atoms']],q['region']['radius'],center=q['region']['center'],max_depth=6)
    overlap=max(r['lower'],ref['lower'])<=min(r['upper'],ref['upper'])
    if not overlap:failures.append([q['id'],'rational_octree_disjoint'])
    independent.append({'id':q['id'],'reference':ref,'overlap':overlap})
    save(out/'independent.json',independent)
receipt={'complete':True,'cases':len(requests),'failures':failures,'target_met':sum(r['target_met'] for r in observed),'unmet':[r['id'] for r in observed if not r['target_met']],
    'scope':'56 frozen real conformers, three radii, explicit shifted centers, hydrogens and unscaled radii. Separate 3D-box and rational octree overlap is corroboration, not exact molecular truth. Timing under validation load is not performance evidence.',
    'configured_radii_f32':radii,'sources':sources,'engine_sha256':sha(a.engine),'reference_engine_sha256':sha(a.reference_engine),'controller_sha256':sha(Path(__file__)),'rows':rows,'independent':independent}
save(out/'summary.json',receipt);print(json.dumps({'cases':len(requests),'target_met':receipt['target_met'],'failures':failures,'unmet':receipt['unmet']}))
if failures:raise SystemExit(1)
