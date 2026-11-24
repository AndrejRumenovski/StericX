"""Geometry/energy/search sensitivity, without selecting structures from experimental outcomes."""
from pathlib import Path
import json,collections,numpy as np,sys
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from generate_conformers import mdl_molecule
from audit_stereochemistry import p_geometry
from analyze_chemistry import values
ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/'.stericx/experimental_validation';OUT=BASE/'conformer_sensitivity_v2'
def main():
 OUT.mkdir(exist_ok=False);stats=json.loads((BASE/'chemical_analysis_v2/ensemble_statistics.json').read_text());byid={r['id']:r for r in stats};native=[json.loads(x) for x in (BASE/'chemistry_preparation/native.jsonl').read_text().splitlines()];rng=np.random.default_rng(20261124);rows=[];stereo=[];generated_native={r['id']:r for r in map(json.loads,(BASE/'generated_v2/native.jsonl').read_text().splitlines())}
 for r in native:
  mid=r['ligand_id'];cs=r['result']['descriptors']['conformers'];v=np.array([c['buried_volume']['value']['summary']['max_delta_qvbur'] for c in cs]);coverage=[]
  for n in sorted(set([1,min(2,len(v)),min(4,len(v)),min(8,len(v)),len(v)])):
   estimates=np.array([min(rng.choice(v,n,replace=False)) for _ in range(1000)]);coverage.append(dict(sampled_conformers=n,minimum_descriptor_5_50_95_quantiles=np.quantile(estimates,[.05,.5,.95]).tolist(),probability_sample_min_equals_full_min=float(np.mean(estimates==min(v)))))
  gen=[byid[f'generated_{mid}_{seed}'] for seed in [20261124,20261125]];embedded=byid.get(f'embedded_{mid}_20261124');features={}
  for f in gen[0]['features']:
   a,b=[g['features'][f] for g in gen];pairs={}
   if embedded is not None:
    pre={c['identifier'].removesuffix(':unoptimized'):values(c)[f] for c in generated_native[f'embedded_{mid}_20261124']['result']['descriptors']['conformers']}
    pairs={c['identifier']:values(c)[f]-pre[c['identifier']] for c in generated_native[f'generated_{mid}_20261124']['result']['descriptors']['conformers'] if c['identifier'] in pre}
   temp=[a['weighted_'+str(t)] for t in [298.15,353.15,400.]];features[f]=dict(published_range=byid[r['id']]['features'][f]['range'],generated_seed_mean_difference=abs(a['mean']-b['mean']),generated_seed_min_difference=abs(a['minimum']-b['minimum']),temperature_weighted_mean_range=max(temp)-min(temp),optimization_mean_difference=None if not pairs else float(np.mean(list(pairs.values()))),optimization_matched_pairs=len(pairs),published_minus_generated_mean=byid[r['id']]['features'][f]['mean']-a['mean'])
  rows.append(dict(ligand_id=mid,published_conformers=len(cs),generated_counts=[g['n_conformers'] for g in gen],coverage=coverage,features=features,generated_energy_analyses=[g['energy_analysis'] for g in gen]))
  paths=[ROOT/c['source'] for c in cs];ref=mdl_molecule(paths[0]);axes=[]
  for b in ref.GetBonds():
   i,j=b.GetBeginAtomIdx(),b.GetEndAtomIdx()
   if b.GetBondType()==Chem.BondType.SINGLE and b.IsInRing() and ref.GetAtomWithIdx(i).GetIsAromatic() and ref.GetAtomWithIdx(j).GetIsAromatic():axes.append([min(n.GetIdx() for n in ref.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx()!=j),i,j,min(n.GetIdx() for n in ref.GetAtomWithIdx(j).GetNeighbors() if n.GetIdx()!=i)])
  rec=[]
  for p in paths:
   m=mdl_molecule(p);rec.append(dict(conformer_id=int(p.stem),dihedrals_rad=[rdMolTransforms.GetDihedralRad(m.GetConformer(),*ids) for ids in axes]))
  stereo.append(dict(ligand_id=mid,axes=axes,records=rec,sign_pattern_counts={str(k):v for k,v in collections.Counter(tuple(int(np.sign(np.sin(x))) for x in r['dihedrals_rad']) for r in rec).items()}))
 a=json.loads((BASE/'generated_v2/audit.json').read_text());mismatches=[]
 for r in a['records']:
  if r['rejection_reasons']:continue
  m=mdl_molecule(ROOT/r['optimized_path']);_,label,_=p_geometry(m);expected=r['encoded_source_p']
  if expected in ['R','S'] and label!=expected:mismatches.append(dict(ligand_id=r['ligand_id'],seed=r['seed'],id=r['embedded_id'],expected=expected,actual=label))
 result=dict(rows=rows,published_relative_axial_stereo=stereo,generated_encoded_p_mismatches=mismatches,generated_embedding_count=len(a['records']),generated_retained_count=sum(not r['rejection_reasons'] for r in a['records']),generated_rejections=dict(collections.Counter(x for r in a['records'] for x in r['rejection_reasons'])),notes='Coverage experiment samples recorded library conformers without replacement and does not fit outcomes. Relative axial sign is not an absolute experimental assignment. Optimization differences compare matched retained pairs; missing pairs are not replaced.')
 (OUT/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['rows','published_relative_axial_stereo']},indent=2));print('axial',[(r['ligand_id'],r['sign_pattern_counts']) for r in stereo]);assert not mismatches
if __name__=='__main__':main()
