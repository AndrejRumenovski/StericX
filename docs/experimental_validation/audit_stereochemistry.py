"""Audit P handedness with retained atom mapping and independently embedded source stereo.
The virtual fourth point is a diagnostic construction, never a molecular atom used
for descriptors. Its removal changes RDKit's implicit-neighbor convention; labels
are assigned only after removal, and checked against stereo-constrained embeddings.
"""
from pathlib import Path
import json,hashlib,collections,numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D
ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/'.stericx/experimental_validation';OUT=BASE/'stereochemistry_v1'
def p_geometry(m):
 e=Chem.RWMol(m);ip=next(a.GetIdx() for a in e.GetAtoms() if a.GetAtomicNum()==15);pos=e.GetConformer().GetPositions();ns=sorted(a.GetIdx() for a in e.GetAtomWithIdx(ip).GetNeighbors());vecs=pos[ns]-pos[ip]
 i=e.AddAtom(Chem.Atom(0));e.AddBond(ip,i,Chem.BondType.SINGLE);e.GetAtomWithIdx(ip).SetFormalCharge(1);e.GetConformer().SetAtomPosition(i,Point3D(*(pos[ip]-(vecs/np.linalg.norm(vecs,axis=1)[:,None]).sum(axis=0))));e.UpdatePropertyCache(strict=False);Chem.AssignStereochemistryFrom3D(e,replaceExistingTags=True);e.RemoveAtom(i);e.GetAtomWithIdx(ip).SetFormalCharge(0)
 for atom in e.GetAtoms():
  if atom.HasProp('_CIPCode'):atom.ClearProp('_CIPCode')
 Chem.AssignStereochemistry(e,force=True,cleanIt=False)
 return e,dict(Chem.FindMolChiralCenters(e,includeUnassigned=True)).get(ip),float(np.linalg.det(vecs))
def main():
 OUT.mkdir(exist_ok=False);rows=[];refs=[]
 for mid in [401,498,723,724,785,1057,1058,2062,2063,2064,2067]:
  d=BASE/'structures'/str(mid);source=Chem.MolFromSmiles(json.loads((d/'molecule.json').read_text())['smiles']);paths=sorted(d.glob('*.sdf'),key=lambda p:int(p.stem));ref=Chem.MolFromMolFile(str(paths[0]),removeHs=False);match=ref.GetSubstructMatch(source);inv={v:i for i,v in enumerate(match)};ip=next(a.GetIdx() for a in ref.GetAtoms() if a.GetAtomicNum()==15);ns=sorted(n.GetIdx() for n in ref.GetAtomWithIdx(ip).GetNeighbors());expected=dict(Chem.FindMolChiralCenters(source,includeUnassigned=True)).get(inv[ip]);check=None
  if expected in ['R','S']:
   embed=Chem.AddHs(source);opts=AllChem.ETKDGv3();opts.randomSeed=20261124;opts.numThreads=1;ids=list(AllChem.EmbedMultipleConfs(embed,10,opts));dets=[]
   for cid in ids:
    xyz=embed.GetConformer(cid).GetPositions();dets.append(float(np.linalg.det(xyz[[inv[n] for n in ns]]-xyz[inv[ip]])))
   assert len(ids)==10 and len(set(np.sign(dets)))==1
   check=float(np.sign(dets[0]));refs.append(dict(id=mid,expected=expected,source_to_export_atom_mapping=match,source_embedding_determinants=dets))
   w=Chem.SDWriter(str(OUT/f'source_stereo_reference_{mid}.sdf'))
   for cid in ids:w.write(embed,confId=cid)
   w.close()
  for path in paths:
   m=Chem.MolFromMolFile(str(path),removeHs=False);diagnostic,assigned,det=p_geometry(m)
   decision='unspecified_or_nonstereogenic_in_source' if check is None else ('match' if np.sign(det)==check else 'mismatch')
   if check is not None:assert (assigned==expected)==(decision=='match')
   rows.append(dict(ligand_id=mid,conformer_id=int(path.stem),source_p_label=expected,geometric_p_label=assigned,ordered_neighbor_determinant_angstrom3=det,source_p_comparison=decision))
 out=dict(method='Ordered neighbor determinant calibrated against 10 stereo-constrained source-SMILES ETKDGv3 embeddings for each encoded P center; agrees with virtual-fourth-point diagnostic after removal and native implicit-neighbor parity conversion.',records=rows,reference_embeddings=refs,counts={str(mid):dict(collections.Counter(r['source_p_comparison'] for r in rows if r['ligand_id']==mid)) for mid in sorted(set(r['ligand_id'] for r in rows))},limitations='Axial/spiro absolute configurations are not fully represented by the supplied SMILES; these remain linked by source identity and relative geometry, not claimed experimentally verified absolute stereo.')
 (OUT/'audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out['counts'],indent=2))
if __name__=='__main__':main()
