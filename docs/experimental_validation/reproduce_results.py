"""Build fresh native adapters and reproduce archived numerical/chemical/model results."""
from pathlib import Path
import argparse,subprocess,os,json,hashlib,contextlib,shutil,sys,time
from profile_methods import strip
ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/'.stericx/experimental_validation'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();out=a.output.resolve();assert out.is_relative_to(BASE);out.mkdir(parents=True,exist_ok=False);receipt=dict(complete=False,checks=[],runs=[],source={str(p.relative_to(ROOT)):sha(p) for p in [*sorted((ROOT/'src').rglob('*.rs')),*sorted((ROOT/'docs/experimental_validation/observer/src').rglob('*.rs')),ROOT/'Cargo.lock',Path(__file__)]})
 def save():(out/'summary.json').write_text(json.dumps(receipt,indent=2)+'\n')
 def run(label,cmd):
  t=time.monotonic();cmd=list(map(str,cmd))
  with (out/(label+'.stdout')).open('wb') as stdout,(out/(label+'.stderr')).open('wb') as stderr:r=subprocess.run(cmd,cwd=ROOT,env=dict(os.environ,CARGO_BUILD_JOBS='2',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',RAYON_NUM_THREADS='2'),stdout=stdout,stderr=stderr)
  receipt['runs'].append(dict(label=label,command=cmd,returncode=r.returncode,seconds=time.monotonic()-t));save();print(label,r.returncode,flush=True);assert r.returncode==0;return out/(label+'.stdout')
 run('verify_restore',[ROOT/'.venv/bin/python',ROOT/'docs/experimental_validation/archive_chemical_data.py','restore'])
 target=Path('/tmp/stericx-experimental-reproduction');run('build',['cargo','build','--locked','--offline','--release','--manifest-path',ROOT/'docs/experimental_validation/observer/Cargo.toml','--target-dir',target]);binary=target/'release'
 cases=[('geometric','geometric',BASE/'numerical_frozen/analytical_requests.json',BASE/'numerical_frozen/geometric_observed.jsonl'),('regions','regions',BASE/'numerical_frozen/regional_requests.json',BASE/'numerical_frozen/regional_observed.jsonl'),('profiles','profiles',BASE/'numerical_frozen/profile_requests.json',BASE/'numerical_frozen/profile_observed.jsonl'),('published','chemistry',BASE/'chemistry_preparation/requests.json',BASE/'chemistry_preparation/native.jsonl'),('generated','chemistry',BASE/'generated_v2/requests.json',BASE/'generated_v2/native.jsonl')]
 for name,engine,request,expected in cases:
  observed=run(name,[binary/engine,request]);left=list(map(json.loads,expected.read_text().splitlines()));right=list(map(json.loads,observed.read_text().splitlines()));assert strip(left)==strip(right),name;receipt['checks'].append(dict(lane=name,rows=len(left),exact_without_time_or_allocations=True,request_sha256=sha(request),expected_sha256=sha(expected),observed_sha256=sha(observed),binary_sha256=sha(binary/engine)));save()
 # Independently rerun Python reductions and every locked fold using fresh Rust output.
 staging=out/'analysis_inputs';staging.mkdir()
 for directory,filename,source in [('chemistry_preparation','native.jsonl',out/'published.stdout'),('generated_v2','native.jsonl',out/'generated.stdout'),('stereochemistry_v1','audit.json',BASE/'stereochemistry_v1/audit.json'),('sources','nihda.csv',BASE/'sources/nihda.csv')]:
  d=staging/directory;d.mkdir(exist_ok=True);shutil.copyfile(source,d/filename)
 import analyze_chemistry as analysis
 analysis.BASE=staging;analysis.OUT=out/'analysis'
 with (out/'analysis.log').open('w') as log,contextlib.redirect_stdout(log):analysis.main()
 expected=json.loads((BASE/'chemical_analysis_v2/statistical_results.json').read_text());observed=json.loads((out/'analysis/statistical_results.json').read_text());expected.pop('inputs');observed.pop('inputs');assert expected==observed
 receipt['checks'].append(dict(lane='all_statistical_results',exact=True))
 for path,h in receipt['source'].items():assert sha(ROOT/path)==h,path
 receipt['complete']=True;receipt['all_checks_passed']=True;save();print('All fresh scientific replays match archived results.',flush=True)
if __name__=='__main__':main()
