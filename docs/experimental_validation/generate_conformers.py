"""Separate reproducible ETKDGv3/MMFF94 preparation experiment; no published energies inferred."""
from pathlib import Path
import json,hashlib,numpy as np,sys
from rdkit import Chem,rdBase
from rdkit.Chem import AllChem,rdMolAlign,rdMolTransforms
from audit_stereochemistry import p_geometry
ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/'.stericx/experimental_validation';OUT=BASE/'generated_v2'
SEEDS=[20261124,20261125]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def determinant(m,center,neighbors,cid=0):
 p=m.GetConformer(cid).GetPositions();return float(np.linalg.det(p[neighbors[:3]]-p[center]))
def mdl_molecule(path):
 m=Chem.MolFromMolFile(str(path),removeHs=False,sanitize=False)
 Chem.SanitizeMol(m,sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
 Chem.SetAromaticity(m,Chem.AromaticityModel.AROMATICITY_MDL);Chem.SetConjugation(m);Chem.SetHybridization(m)
 return m
def main():
 OUT.mkdir(exist_ok=False);records=[];requests=[]
 for mid in [401,498,723,724,785,1057,1058,2062,2063,2064,2067]:
  src=BASE/'structures'/str(mid);meta=json.loads((src/'molecule.json').read_text());path=sorted(src.glob('*.sdf'),key=lambda p:int(p.stem))[0];ref=Chem.MolFromMolFile(str(path),removeHs=False);source=Chem.AddHs(Chem.MolFromSmiles(meta['smiles']));mapping=source.GetSubstructMatch(ref,useChirality=False);assert len(mapping)==ref.GetNumAtoms()
  source_labels={list(mapping).index(i):label for i,label in Chem.FindMolChiralCenters(source,includeUnassigned=True) if label in ['R','S']};ref=mdl_molecule(path);source=Chem.Mol(ref);Chem.AssignStereochemistryFrom3D(source,replaceExistingTags=True)
  for idx,label in source_labels.items():
   for tag in [Chem.ChiralType.CHI_TETRAHEDRAL_CW,Chem.ChiralType.CHI_TETRAHEDRAL_CCW]:
    source.GetAtomWithIdx(idx).SetChiralTag(tag);Chem.AssignStereochemistry(source,force=True,cleanIt=True)
    if dict(Chem.FindMolChiralCenters(source,includeUnassigned=True)).get(idx)==label:break
   else:raise ValueError('cannot transfer source-encoded center')
  ip=next(a.GetIdx() for a in source.GetAtoms() if a.GetAtomicNum()==15);pn=sorted(a.GetIdx() for a in source.GetAtomWithIdx(ip).GetNeighbors());source_p=dict(Chem.FindMolChiralCenters(source,includeUnassigned=True)).get(ip)
  diag,ref_label,ref_det=p_geometry(ref)
  # Map the diagnostic tag through the unchanged export neighbor order; source molecules
  # have a different bond insertion order, so use CIP-target iteration rather than raw tags.
  if source_p=='?' and ref_label in ['R','S']:
   for tag in [Chem.ChiralType.CHI_TETRAHEDRAL_CW,Chem.ChiralType.CHI_TETRAHEDRAL_CCW]:
    source.GetAtomWithIdx(ip).SetChiralTag(tag);Chem.AssignStereochemistry(source,force=True,cleanIt=True)
    if dict(Chem.FindMolChiralCenters(source,includeUnassigned=True)).get(ip)==ref_label:break
   else:raise ValueError('cannot transfer unspecified P geometry')
  axial=[]
  for b in ref.GetBonds():
   i,j=b.GetBeginAtomIdx(),b.GetEndAtomIdx()
   if b.GetBondType()==Chem.BondType.SINGLE and b.IsInRing() and ref.GetAtomWithIdx(i).GetIsAromatic() and ref.GetAtomWithIdx(j).GetIsAromatic():
    a=min(n.GetIdx() for n in ref.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx()!=j);d=min(n.GetIdx() for n in ref.GetAtomWithIdx(j).GetNeighbors() if n.GetIdx()!=i);ids=[a,i,j,d];angle=rdMolTransforms.GetDihedralRad(ref.GetConformer(),*ids);axial.append((ids,float(np.sign(np.sin(angle)))))
  # Unencoded spiro centers: preserve handedness of ordered heavy-neighbor vectors.
  spiro=[]
  for a in ref.GetAtoms():
   ns=sorted(n.GetIdx() for n in a.GetNeighbors() if n.GetAtomicNum()>1)
   if a.GetAtomicNum()==6 and len(ns)==4 and a.IsInRing():spiro.append((a.GetIdx(),ns,float(np.sign(determinant(ref,a.GetIdx(),ns)))))
  for seed in SEEDS:
   work=Chem.Mol(source);params=AllChem.ETKDGv3();params.randomSeed=seed;params.numThreads=1;params.enforceChirality=True;params.pruneRmsThresh=0.25;params.onlyHeavyAtomsForRMS=True;params.useSymmetryForPruning=True;params.maxIterations=2000
   ids=list(AllChem.EmbedMultipleConfs(work,numConfs=20,params=params));folder=OUT/str(mid)/str(seed);folder.mkdir(parents=True);accepted=[];pre=[];retained=[];props=AllChem.MMFFGetMoleculeProperties(work,mmffVariant='MMFF94');assert props is not None
   def relative_ok(cid):
    failures=[]
    for ix,sign in axial:
     if float(np.sign(np.sin(rdMolTransforms.GetDihedralRad(work.GetConformer(cid),*ix))))!=sign:failures.append('opposite_source_relative_axial_sector')
    for center,ns,sign in spiro:
     if float(np.sign(determinant(work,center,ns,cid)))!=sign:failures.append('opposite_source_relative_spiro_handedness')
    return failures
   for cid in ids:
    prepath=folder/f'{cid}_embedded.sdf';prepath.write_text(Chem.MolToMolBlock(work,confId=cid)+'\n$$$$\n');before_fail=relative_ok(cid)
    ff=AllChem.MMFFGetMoleculeForceField(work,props,confId=cid);initial=ff.CalcEnergy();status=ff.Minimize(maxIts=2000,forceTol=1e-4,energyTol=1e-6);optimized=ff.CalcEnergy();postpath=folder/f'{cid}_optimized.sdf';postpath.write_text(Chem.MolToMolBlock(work,confId=cid)+'\n$$$$\n');fail=relative_ok(cid)
    if status!=0:fail.append('optimizer_not_converged')
    rms=None
    if retained:
     heavy=Chem.RemoveHs(work);rms=min(rdMolAlign.GetBestRMS(Chem.Mol(heavy),heavy,prbId=cid,refId=prev) for prev in retained)
     if rms<0.25:fail.append('postoptimization_duplicate_heavy_rms_below_0.25_A')
    # Energy evaluated again on the actual four-decimal SDF used by the Rust engine.
    saved=mdl_molecule(postpath);sp=AllChem.MMFFGetMoleculeProperties(saved,mmffVariant='MMFF94');energy=AllChem.MMFFGetMoleculeForceField(saved,sp).CalcEnergy()
    row=dict(ligand_id=mid,seed=seed,embedded_id=cid,embedded_path=str(prepath.relative_to(ROOT)),optimized_path=str(postpath.relative_to(ROOT)),embedded_sha256=sha(prepath),optimized_sha256=sha(postpath),embedding_relative_failures=before_fail,rejection_reasons=fail,optimizer_status=status,initial_energy_kcal_mol=initial,optimized_energy_kcal_mol=optimized,serialized_geometry_energy_kcal_mol=energy,minimum_previous_heavy_rms=rms,source_to_original_export_mapping=list(mapping),encoded_source_p=source_p,relative_axial_constraints=axial,relative_spiro_constraints=spiro)
    records.append(row)
    if not fail:
     retained.append(cid);ident=f'generated:{mid}:{seed}:{cid}';accepted.append(dict(path=str(postpath.relative_to(ROOT)),identifier=ident,energy=dict(value=energy,unit='kcal_mol',kind='force_field',model=f'RDKit {rdBase.rdkitVersion} MMFF94 gas-phase single-point on optimized 4-decimal SDF; no entropy/solvent',provenance=str(postpath.relative_to(ROOT)),temperature_k=None,degeneracy=1)))
     if not before_fail:pre.append(dict(path=str(prepath.relative_to(ROOT)),identifier=ident+':unoptimized',energy=None))
   if accepted:requests.append(dict(id=f'generated_{mid}_{seed}',ligand_id=mid,lane=f'generated_mmff94_seed_{seed}',files=accepted,certified=False))
   if seed==SEEDS[0] and pre:requests.append(dict(id=f'embedded_{mid}_{seed}',ligand_id=mid,lane=f'generated_unoptimized_seed_{seed}',files=pre,certified=False))
   print(mid,seed,'embedded',len(ids),'retained',len(accepted),flush=True)
   (OUT/'audit.json').write_text(json.dumps(dict(preparation='Original SDF bond orders; MDL aromaticity and hybridization before ETKDG, explicit source R/S labels, indexed ring biaryl and spiro handedness filters; MMFF uses its own internal aromaticity',rdkit_version=rdBase.rdkitVersion,seeds=SEEDS,requested_per_seed=20,records=records,requests=requests,complete=False),indent=2)+'\n')
 audit=dict(preparation='Original SDF bond orders; MDL aromaticity and hybridization before ETKDG, explicit source R/S labels, indexed ring biaryl and spiro handedness filters; MMFF uses its own internal aromaticity',rdkit_version=rdBase.rdkitVersion,seeds=SEEDS,requested_per_seed=20,records=records,requests=requests,complete=True,limitations='Small gas-phase search; relative P/axial/spiro choices follow one published geometry where source SMILES lacks absolute labels. 401 follows its encoded source P label, deliberately different from published-library coordinates. RMS pruning and unequal search coverage do not define measured state degeneracies or solution populations.')
 (OUT/'audit.json').write_text(json.dumps(audit,indent=2)+'\n');(OUT/'requests.json').write_text(json.dumps(requests,indent=2)+'\n')
if __name__=='__main__':main()
