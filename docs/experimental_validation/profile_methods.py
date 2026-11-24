"""Exploratory source-bound CPU/heap/instruction and parallel diagnostics.
Ordinary workstation load: these timings are descriptive, not controlled speed claims.
"""
from pathlib import Path
import json,os,subprocess,time,hashlib,statistics
ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/'.stericx/experimental_validation';OUT=BASE/'performance_v1';BIN=Path('/tmp/stericx-experimental-profile/release');FROZEN=Path('/tmp/stericx-experimental-baseline-target/release/certified')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def strip(v):
 if isinstance(v,dict):return {k:strip(x) for k,x in v.items() if k not in ['elapsed_seconds','elapsed_ns','allocations','geometry_preparation_seconds','preparation_seconds']}
 if isinstance(v,list):return [strip(x) for x in v]
 return v

def main():
 OUT.mkdir(exist_ok=False);qs=json.loads((ROOT/'.stericx/occupancy_precision/analytical_chord_requests.json').read_text());unmet=json.loads((ROOT/'docs/occupancy_precision/c1_validation.json').read_text())['analytical']['unmet'];original=[q for q in qs if q['id'] in unmet];expanded=[dict(q,budget=4*q['budget'],id=q['id']+'_budget4x') for q in original];save(OUT/'original.json',original);save(OUT/'expanded.json',expanded)
 # Twelve frozen molecular geometries, each with identical per-request tolerances.
 profiles=json.loads((BASE/'profiles_independent_v1/requests.json').read_text());molecules=[q for q in profiles if q['id'].startswith('molecular')];
 if not molecules:molecules=profiles[12:]
 regional=[dict(id=q['id'],region=dict(center=q['center'],radius=3.5),atoms=q['atoms'],tolerance=.05,budget=1048576) for q in molecules];assert len(regional)==12
 save(OUT/'regions.json',regional);save(OUT/'one_region.json',regional[:1]);oneprofile=dict(molecules[0],controls=False);save(OUT/'one_profile.json',[oneprofile]);save(OUT/'one_chain.json',[next(q for q in original if 'overlapping_chain' in q['id'])])
 receipt=dict(scope=__doc__,runs=[],complete=False,bindings={str(p):sha(p) for p in [Path(__file__),FROZEN,*[BIN/n for n in ['certified','geometric','regions','profiles','regional_batch']]]},host=dict(cpuinfo=Path('/proc/cpuinfo').read_text(),loadavg=Path('/proc/loadavg').read_text(),allowed_affinity=sorted(os.sched_getaffinity(0))),original_targets_unchanged=True)
 env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',RAYON_NUM_THREADS='1')
 def run(name,cmd,extra=None):
  folder=OUT/name;folder.mkdir();cmd=list(map(str,cmd));t=time.monotonic()
  with (folder/'stdout').open('wb') as out,(folder/'stderr').open('wb') as err:r=subprocess.run(['/usr/bin/time','-v','-o',str(folder/'process_metrics.txt'),*cmd],cwd=ROOT,env=dict(env,**(extra or {})),stdout=out,stderr=err)
  row=dict(name=name,command=cmd,extra_environment=extra,returncode=r.returncode,controller_wall_seconds=time.monotonic()-t,stdout_sha256=sha(folder/'stdout'),metrics=(folder/'process_metrics.txt').read_text());receipt['runs'].append(row);save(OUT/'summary.json',receipt);print(name,r.returncode,flush=True);return folder
 for repeat in range(3):
  for kind,engine in ([('frozen',FROZEN),('geometric',BIN/'geometric')] if repeat%2==0 else [('geometric',BIN/'geometric'),('frozen',FROZEN)]):run(f'original_{kind}_{repeat}',[engine,OUT/'original.json'])
 for kind,engine in [('frozen',FROZEN),('geometric',BIN/'geometric')]:run('expanded_'+kind,[engine,OUT/'expanded.json'])
 run('allocations_old',[BIN/'certified',OUT/'original.json']);run('allocations_regions',[BIN/'regions',OUT/'regions.json']);run('profile_with_duplicate',[BIN/'profiles',OUT/'one_profile.json'])
 parallel=[];expected=None
 for threads in [1,2,4,6]:
  times=[]
  for repeat in range(3):
   folder=run(f'parallel_{threads}_{repeat}',[BIN/'regional_batch',OUT/'regions.json',threads]);result=json.loads((folder/'stdout').read_text());normalized=strip(result['rows']);expected=expected or normalized;assert expected==normalized;times.append(result['elapsed_seconds'])
  parallel.append(dict(threads=threads,elapsed_seconds=times,median_seconds=statistics.median(times)))
 for row in parallel:row['exploratory_speedup']=parallel[0]['median_seconds']/row['median_seconds'];row['exploratory_efficiency']=row['exploratory_speedup']/row['threads']
 receipt['parallel']=parallel;save(OUT/'summary.json',receipt)
 vg=ROOT/'.stericx/profiling/tools/valgrind';vgenv={'VALGRIND_LIB':str(vg/'usr/libexec/valgrind')};vgb=vg/'usr/bin/valgrind.bin';receipt['bindings'][str(vgb)]=sha(vgb)
 for label,binary,requests in [('old_chain',FROZEN,'one_chain.json'),('new_chain',BIN/'geometric','one_chain.json'),('regions',BIN/'regions','one_region.json'),('profile',BIN/'profiles','one_profile.json')]:
  run('callgrind_'+label,[vgb,'--tool=callgrind','--callgrind-out-file='+str(OUT/('callgrind_'+label+'.out')),binary,OUT/requests],vgenv)
 run('memcheck_regions',[vgb,'--tool=memcheck','--error-exitcode=99','--leak-check=full','--errors-for-leak-kinds=definite',BIN/'regions',OUT/'one_region.json'],vgenv)
 run('perf_hardware_attempt',['perf','stat','-e','instructions,cycles',BIN/'geometric',OUT/'one_chain.json'])
 # Containment remains checked independently for every expanded request as well.
 for label,requests in [('original',original),('expanded',expanded)]:
  for method in ['frozen','geometric']:
   folder=OUT/(f'original_{method}_0' if label=='original' else f'expanded_{method}');rows=list(map(json.loads,(folder/'stdout').read_text().splitlines()));assert len(rows)==len(requests)
   for q,r in zip(requests,rows):
    lo,hi=q['independent_reference'];assert r['lower']<=lo<=hi<=r['upper'];assert r['cells']<=q['budget']
   receipt.setdefault('target_results',{})[label+'_'+method]=dict(count=len(rows),met=sum(r['target_met'] for r in rows),maximum_error=max(r['absolute_error_bound'] for r in rows))
 for p,h in receipt['bindings'].items():assert sha(Path(p))==h
 receipt['complete']=True;receipt['scientific_checks_passed']=True;receipt['known_hardware_counter_unavailability']='perf_event_paranoid=4; failure retained, no permissions changed.';save(OUT/'summary.json',receipt)
if __name__=='__main__':main()
