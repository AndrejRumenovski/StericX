#!/usr/bin/env python3
"""Freeze topology-based requests without consulting SUT outcomes or residuals."""
from __future__ import annotations
import collections,datetime,hashlib,json
from pathlib import Path
from rdkit import Chem
ROOT=Path(__file__).resolve().parents[1]
CONFIG={'sphere_radius':3.5,'density':0.01,'center_distance':2.28,'radii_scale':1.17,'include_hydrogens':False}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def emit(stream,value):stream.write(json.dumps(value,separators=(',',':'),allow_nan=False)+'\n')
def main():
    out=ROOT/'prepared_all';out.mkdir(exist_ok=False)
    ligands=json.loads((ROOT/'primary/ligands_manifest.json').read_text())['ligands']
    manifests={r['conformer_id']:r for r in json.loads((ROOT/'primary/conformers_manifest.json').read_text())['conformers']}
    counts=collections.Counter();topology_counts=collections.Counter();input_hashes={}
    with (out/'requests.jsonl').open('w') as requests,(out/'inventory.jsonl').open('w') as inventory:
        for ligand in ligands:
            mid=ligand['molecule_id']
            if not ligand.get('dft_available'):
                emit(inventory,{'molecule_id':mid,'status':'no_DFT_geometry_available','acquisition_record':ligand});counts['ligands_no_dft']+=1;continue
            counts['ligands_with_dft']+=1
            for cid in ligand['conformer_ids']:
                row={'id':f'KRAKEN:{mid}:{cid}','molecule_id':mid,'conformer_id':cid,'smiles':ligand['smiles']}
                path=ROOT/'primary/conformers'/f'{cid}.body'
                row['sdf_path']=str(path.relative_to(ROOT));row['sdf_sha256']=sha(path)
                assert row['sdf_sha256']==manifests[cid]['sha256']
                input_hashes[row['sdf_path']]=row['sdf_sha256']
                counts['conformers_present']+=1
                try:
                    mol=Chem.MolFromMolBlock(path.read_text(),removeHs=False,sanitize=False,strictParsing=True)
                    if mol is None or mol.GetNumConformers()!=1:raise ValueError('SDF did not yield exactly one coordinate set')
                    copied=Chem.Mol(mol)
                    row['sanitization_flag']=str(Chem.SanitizeMol(copied,catchErrors=True))
                    phosphorus=[a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol()=='P']
                    topological=[i for i in phosphorus if mol.GetAtomWithIdx(i).GetDegree()<=3]
                    row['phosphorus_atoms']=[{'index':i,'degree':mol.GetAtomWithIdx(i).GetDegree(),'charge':mol.GetAtomWithIdx(i).GetFormalCharge(),'neighbors':list(n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors())} for i in phosphorus]
                    if len(phosphorus)==1:donor=phosphorus[0]
                    elif len(topological)==1:donor=topological[0]
                    else:raise ValueError(f'independent_topology_ambiguous_donor:{phosphorus}')
                    neighbors=sorted(n.GetIdx() for n in mol.GetAtomWithIdx(donor).GetNeighbors())
                    heavy=[i for i in neighbors if mol.GetAtomWithIdx(i).GetSymbol()!='H']
                    if not neighbors:raise ValueError('independent_topology_no_neighbor')
                    reference=(heavy or neighbors)[0]
                    coordinate=mol.GetConformer()
                    atoms=[{'element':a.GetSymbol(),'position':list(coordinate.GetAtomPosition(a.GetIdx()))} for a in mol.GetAtoms()]
                    labels=[mol.GetAtomWithIdx(i).GetSymbol() for i in neighbors]
                    category=f"P_degree_{len(neighbors)}_"+'-'.join(sorted(labels))
                    row.update(status='submitted',donor=donor,reference=reference,neighbors=neighbors,neighbor_elements=labels,topology_category=category,atom_count=len(atoms),P_H_count=labels.count('H'),chemical_three_coordinate=len(neighbors)==3)
                    emit(requests,{'op':'geometry','id':row['id'],'atoms':atoms,'donor':donor,'reference':reference,'neighbors':neighbors,'config':CONFIG})
                    counts['requests']+=1;topology_counts[category]+=1
                except Exception as error:
                    row.update(status='independent_preparation_failure',error=repr(error));counts['preparation_failures']+=1
                emit(inventory,row)
    manifest={'created_utc':datetime.datetime.now(datetime.UTC).isoformat(),'script_sha256':sha(Path(__file__)),'config':CONFIG,'counts':dict(counts),'topology_counts':dict(topology_counts),'eligibility':'SDF bond graph only: if one P, attempt it regardless of degree; if several P, require a unique P with <=3 neighbors. Reference is first bonded non-H in SDF atom order, or first bonded H when no heavy neighbor exists, so PH3 is attempted and any SUT restriction is recorded. No StericX cutoff or result is consulted. All failed preparations remain in inventory; SUT operation failures remain in raw outputs.','scientific_population':'Primary comparisons require all source conformers and operation outputs for each ligand. Three-coordinate chemistry is recorded, never used to delete an unfavorable result.','raw_geometry_sha256':input_hashes,'prepared_files':{name:sha(out/name) for name in ['requests.jsonl','inventory.jsonl']},'input_manifests':{name:sha(ROOT/'primary'/name) for name in ['ligands_manifest.json','conformers_manifest.json','acquisition_summary.json']}}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'counts':dict(counts),'topology_counts':dict(topology_counts)},indent=2))
if __name__=='__main__':main()
