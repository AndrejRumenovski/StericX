"""Audit original molecular identities/order before requesting Rust descriptors."""
from pathlib import Path
import hashlib,json
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
ROOT=Path(__file__).resolve().parents[2]
SRC=ROOT/'.stericx/experimental_validation/structures'
OUT=ROOT/'.stericx/experimental_validation/chemistry_preparation'

def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def canonical(m,stereo=True):return Chem.MolToSmiles(Chem.RemoveHs(m),isomericSmiles=stereo)
def signature(m):
 return dict(elements=[a.GetSymbol() for a in m.GetAtoms()],bonds=sorted((min(b.GetBeginAtomIdx(),b.GetEndAtomIdx()),max(b.GetBeginAtomIdx(),b.GetEndAtomIdx()),str(b.GetBondType())) for b in m.GetBonds()))

def main():
 OUT.mkdir(parents=True,exist_ok=False);records=[];requests=[]
 for mid in [401,498,723,724,785,906,1057,1058,2062,2063,2064,2067]:
  mp=SRC/str(mid)/'molecule.json'
  try:metadata=json.loads(mp.read_text())
  except Exception as e:records.append(dict(id=mid,error=f'metadata unavailable: {e}'));continue
  if not metadata.get('conformers_id'):records.append(dict(id=mid,error='no conformer IDs'));continue
  source=Chem.MolFromSmiles(metadata['smiles']);assert source is not None
  origin_signature=None;files=[];ligand=[]
  for cid in sorted(metadata['conformers_id']):
   p=mp.parent/f'{cid}.sdf';r=dict(ligand_id=mid,conformer_id=cid,path=str(p.relative_to(ROOT)),sha256=sha(p),geometry_provenance='original MolSSI Kraken DFT-library V2000 export, retrieved 2026-11-24; not newly optimized; exported coordinates rounded to 4 decimals',energy=None)
   try:
    m=Chem.MolFromMolBlock(p.read_text(),removeHs=False,sanitize=True);assert m is not None
    sig=signature(m)
    if origin_signature is None:origin_signature=sig
    Chem.AssignStereochemistryFrom3D(m,confId=0,replaceExistingTags=True)
    mirror=Chem.Mol(m)
    for atom in mirror.GetAtoms():atom.InvertChirality()
    source_smiles=canonical(source);observed=canonical(m)
    stereo='exact_encoded_match' if observed==source_smiles else ('encoded_enantiomer' if canonical(mirror)==source_smiles else 'unresolved_or_partially_encoded')
    donors=[a.GetIdx() for a in m.GetAtoms() if a.GetSymbol()=='P'];assert len(donors)==1
    donor=donors[0];neighbors=sorted(n.GetIdx() for n in m.GetAtomWithIdx(donor).GetNeighbors());assert len(neighbors)==3
    r.update(formula=rdMolDescriptors.CalcMolFormula(m),formal_charge=Chem.GetFormalCharge(m),fragments=len(Chem.GetMolFrags(m)),atom_count=m.GetNumAtoms(),explicit_hydrogens=sum(a.GetAtomicNum()==1 for a in m.GetAtoms()),donor_index=donor,neighbors=neighbors,source_smiles=source_smiles,export_smiles=observed,stereo_comparison=stereo,connectivity_matches_source=canonical(source,False)==canonical(m,False),ordered_graph_matches_first=sig==origin_signature,chiral_centers_from_3d=Chem.FindMolChiralCenters(m,includeUnassigned=True),source_chiral_centers=Chem.FindMolChiralCenters(source,includeUnassigned=True),atom_mapping='identity within retained original export ordering',metal_atoms=[a.GetSymbol() for a in m.GetAtoms() if a.GetAtomicNum() in [26,27,28,29,44,45,46,47,76,77,78,79]])
    assert r['formal_charge']==0 and r['fragments']==1 and r['connectivity_matches_source']
    assert r['ordered_graph_matches_first']
    files.append(dict(path=str(p.relative_to(ROOT)),identifier=f'KRAKEN:{mid}:{cid}',energy=None))
   except Exception as e:r['error']=repr(e)
   ligand.append(r)
  records.extend(ligand)
  # Complete ensembles only; incomplete ligands remain in the audit and do not silently shrink.
  if len(files)==len(metadata['conformers_id']):requests.append(dict(id=f'published_{mid}',ligand_id=mid,lane='published_dft_library',files=files,temperature_k=353.15,certified=True))
 summary=dict(records=records,eligible_ligands=[q['ligand_id'] for q in requests],errors=[r for r in records if 'error' in r],stereo_policy='Encoded stereo and 3D assignment compared explicitly. Axial/P stereochemistry not encoded or not perceived remains unresolved; primary scalar steric features and absolute ee are reflection-invariant, but a mismatched diastereomer would not be equivalent.')
 (OUT/'audit.json').write_text(json.dumps(summary,indent=2)+'\n');(OUT/'requests.json').write_text(json.dumps(requests,indent=2)+'\n')
 print(json.dumps(dict(eligible=summary['eligible_ligands'],records=len(records),errors=summary['errors'],stereo={s:sum(r.get('stereo_comparison')==s for r in records) for s in sorted({r['stereo_comparison'] for r in records if 'stereo_comparison'in r})}),indent=2))
if __name__=='__main__':main()
