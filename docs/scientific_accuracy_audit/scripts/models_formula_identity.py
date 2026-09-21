"""Frozen Study 007 formula-helper counterexample; not a historical remapping.

RDKit supplies distinct explicit chemical graphs and an independent molecular
formula. The observed helpers are compiled directly from the frozen study AST.
No production file or scientific kernel is modified.
"""
from models_audit import *
import ast
from collections import Counter
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdMolDescriptors


def main():
    directory=M/'inputs/formula_identity';directory.mkdir(exist_ok=True)
    cases=[]
    for name,smiles in [('n_propylphosphine','CCCP'),('methylethylphosphine','CCPC')]:
        mol=Chem.AddHs(Chem.MolFromSmiles(smiles))
        assert AllChem.EmbedMolecule(mol,randomSeed=20260921)==0
        mol.SetProp('_Name',name)
        xyz=directory/f'{name}.xyz';sdf=directory/f'{name}.sdf'
        freeze(xyz,Chem.MolToXYZBlock(mol).encode())
        freeze(sdf,(Chem.MolToMolBlock(mol)+'\n$$$$\n').encode())
        donor=next(a for a in mol.GetAtoms() if a.GetAtomicNum()==15)
        cases.append({'name':name,'input_smiles':smiles,'canonical_smiles':Chem.MolToSmiles(Chem.RemoveHs(mol)),
                      'rdkit_formula':rdMolDescriptors.CalcMolFormula(mol),
                      'donor_bonded_elements':sorted(a.GetSymbol() for a in donor.GetNeighbors()),
                      'atoms':[{'index':a.GetIdx(),'element':a.GetSymbol()} for a in mol.GetAtoms()],
                      'bonds':[{'begin':b.GetBeginAtomIdx(),'end':b.GetEndAtomIdx(),'order':b.GetBondTypeAsDouble()} for b in mol.GetBonds()],
                      'xyz':str(xyz.relative_to(A)),'sdf':str(sdf.relative_to(A))})
    graph_path=directory/'explicit_graphs.json';freeze(graph_path,(json.dumps(cases,indent=2)+'\n').encode())
    source=F/'studies/study_007_crosscoupling.py'
    phase={'scope':'Input graphs and coordinates frozen before invoking the formula helpers; observation is helper-level, not a rerun of historical mappings.',
           'rdkit_version':rdBase.rdkitVersion,'embedding_seed':20260921,
           'frozen_study_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
           'inputs':{str(p.relative_to(A)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(directory.iterdir()) if p.is_file()}}
    freeze(M/'inputs/formula_identity_manifest_before_execution.json',(json.dumps(phase,indent=2)+'\n').encode())
    tree=ast.parse(source.read_text());names={'_xyz_formula','_sdf_formula'}
    definitions=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names]
    assert {n.name for n in definitions}==names
    namespace={'Counter':Counter,'Path':Path}
    exec(compile(ast.Module(body=definitions,type_ignores=[]),str(source),'exec'),namespace)
    comparisons=[]
    for left in cases:
        for right in cases:
            a=namespace['_xyz_formula']((A/left['xyz']).read_text())
            b=namespace['_sdf_formula'](A/right['sdf'])
            comparisons.append({'xyz_graph':left['name'],'sdf_graph':right['name'],
                                'xyz_formula_helper_result':a,'sdf_formula_helper_result':b,
                                'formula_guard_rejects':a!=b,
                                'independent_canonical_graphs_equal':left['canonical_smiles']==right['canonical_smiles']})
    result={'cases':cases,'comparisons':comparisons,
            'finding':'Both distinct explicit graphs have C3H9P. The formula helpers accept the cross-isomer pairing; donor heavy-neighbor counts differ (one versus two). Formula comparison cannot guarantee identity or reject constitutional isomers.',
            'limitation':'This synthetic helper-level witness does not show that any of the eighteen historical mappings was wrong. Other lookup/presence conditions are not chemical-identity tests and were not rerun here.'}
    payload=(json.dumps(result,indent=2)+'\n').encode()
    freeze(M/'raw/formula_identity_helper_observation.json',payload)
    save(M/'results/formula_identity_counterexample.json',result)


if __name__=='__main__':main()
