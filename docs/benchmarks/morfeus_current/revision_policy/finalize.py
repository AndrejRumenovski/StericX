"""Post-campaign documentation/reproduction workflow; seals are made after exit."""
import json,os,shutil,subprocess,sys,tarfile,hashlib
from pathlib import Path
R=Path(__file__).resolve().parent;H=R.parent;ROOT=H.parents[2]
def read(p):return json.loads(p.read_text())
def write(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(name,cmd):
 with (R/(name+'.stdout')).open('wb') as out,(R/(name+'.stderr')).open('wb') as err:
  p=subprocess.run(cmd,cwd=ROOT,stdout=out,stderr=err)
 write(R/(name+'.command.json'),dict(argv=cmd,returncode=p.returncode))
 p.check_returncode();print(name+' completed',flush=True)
def main():
 assert read(R/'admitted_campaign/completion.json')['passed']
 for name in ['postflight','investigate']:run(name,[sys.executable,str(R/(name+'.py'))])
 scratch=H.parent/'.morfeus_policy_reproduction_smoke_20260922';assert not scratch.exists()
 run('reproduction-smoke',[sys.executable,str(H/'reproduce.py'),str(scratch),'--prepare-only'])
 gate=read(scratch/'revision_policy/evaluation_2/scientific_gate.json');assert gate['passed']
 # Make regenerable symlink aliases portable; source/file identities do not change.
 links=[]
 for directory in [H,scratch]:
  for base,dirs,names in os.walk(directory,followlinks=False):
   dirs[:]=[n for n in dirs if n not in ['target','venv','__pycache__','repeated_10000']]
   for name in dirs+names:
    p=Path(base)/name
    if p.is_symlink():
     target=p.resolve();old=os.readlink(p);new=os.path.relpath(target,p.parent)
     p.unlink();p.symlink_to(new,target_is_directory=target.is_dir())
     if directory==H:links.append(dict(path=str(p.relative_to(H)),old=old,new=new))
 archive=R/'reproduction-smoke.tar.gz';assert not archive.exists()
 with tarfile.open(archive,'w:gz') as t:
  for base,dirs,names in os.walk(scratch,followlinks=False):
   dirs[:]=[n for n in dirs if n not in ['target','venv','__pycache__','repeated_10000']]
   for name in names:
    p=Path(base)/name;t.add(p,arcname=str(p.relative_to(scratch)),recursive=False)
   for name in dirs:
    p=Path(base)/name
    if p.is_symlink():t.add(p,arcname=str(p.relative_to(scratch)),recursive=False)
 write(R/'reproduction_smoke.json',dict(passed=True,scope='Clean sibling rebuild, pinned fresh Python environment, independent reference tests and all56 prospective scientific checks; no performance samples from smoke run',policy_sha256=gate['policy_sha256'],archive_sha256=sha(archive),archive='reproduction-smoke.tar.gz'))
 assert read(scratch/'freeze.json')['reproduced_from']==str(H)
 shutil.rmtree(scratch)
 write(R/'portable_symlinks.json',dict(note='Only alias paths changed after computation; pointed-to bytes unchanged.',symlinks=links))
 # Replace historical gate labels at entry points, retaining original bytes in prior archive.
 gate=read(R/'evaluation_2/scientific_gate.json')
 write(H/'scientific_gate.json',gate|dict(descriptors=gate['admitted_descriptors'],numerical_tolerance=None,scope='Prospectively frozen finite-lattice volume gate; see revision_policy/evaluation_2'))
 p=H/'pyramidalization/scientific_gate.json';old=read(p);old['exploratory_check_passed']=old.pop('passed');old['passed']=False;old['publication_status']='excluded_under_frozen_policy';old['reason']='Observed agreement after f32 rounding is not a prospectively derived certified error bound.';write(p,old)
 p=H/'pyramidalization/benchmark_results.json';old=read(p);old['publication_status']='withdrawn_exploratory_only';old['scientific_gate']=read(H/'pyramidalization/scientific_gate.json');write(p,old)
 p=H/'scientific_summary.json';old=read(p);old['scientific_gate']='Original exploratory residual report. Current admission is revision_policy/evaluation_2/scientific_gate.json under frozen EQUIVALENCE_POLICY.md.';write(p,old)
 for name in ['exact_volume_gate.py','pyramid_benchmark.py','finish_measurements.py']:
  p=H/name;text=p.read_text();p.write_text('if __name__ == "__main__":\n    raise SystemExit("Historical exploratory runner retired. Use reproduce.py and the frozen EQUIVALENCE_POLICY.md.")\n'+text)
 (H/'report.py').write_text('"""Regenerate only the prospectively admitted report; all gates are checked."""\nfrom pathlib import Path\nimport runpy\nif __name__ == "__main__":\n    runpy.run_path(str(Path(__file__).resolve().parent/"revision_policy/publish.py"),run_name="__main__")\n')
 shutil.copy2(R/'verify_revision.py',H/'verify.py')
 run('publication',[sys.executable,str(R/'publish.py')])
 run('whitespace-check',['git','diff','--check'])
 write(R/'finalization.json',dict(passed=True,needs_final_evidence_seal=True))
 print('Postflight, clean reproduction and publication passed. Seal only after these logs close.',flush=True)
if __name__=='__main__':main()
