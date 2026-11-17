#!/usr/bin/env python3
"""Final A/B controller: gated quiet start, all outputs/trials retained, no outlier removal."""
import argparse,hashlib,json,os,statistics,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import profile_stericx as harness
from host_activity import measure, postflight
from benchmark_artifacts import compact_metrics, retain_fingerprint

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline',type=Path,required=True);p.add_argument('--candidate',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--preflight-only',action='store_true');a=p.parse_args()
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    journal={'scope':'Portable native CLI, equal outputs, fixed affinity, warmed file cache and fresh processes',
        'pairs':8,'warmups_per_side':1,'attempts':[],'complete':False,'preflight_only':a.preflight_only}
    def save():(out/'summary.json').write_text(json.dumps(journal,indent=2)+'\n')
    if not a.preflight_only:
        journal['host_preflight']=measure();save()
        if not journal['host_preflight']['quiet']:
            raise SystemExit('Host not quiet: no timing attempted; request quiet window after independent work.')
    manifest_path=ROOT/'.stericx/profiling/scientifically_validated_optimization/corrected_workloads_v2/workloads/manifest.json'
    manifest=json.loads(manifest_path.read_text());harness.verify_manifest(manifest)
    binaries={'baseline':a.baseline.resolve()/'native-target/release/stericx',
              'candidate':a.candidate.resolve()/'native-target/release/stericx'}
    bindings=[manifest_path,Path(__file__),Path(__file__).with_name('benchmark_artifacts.py'),Path(__file__).with_name('host_activity.py'),*binaries.values()]
    journal['bindings']={str(p):harness.digest(p) for p in bindings}
    launcher=harness.build_launcher(out/'launcher')
    env=dict(os.environ,LC_ALL='C',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
    env.pop('STERICX_PROFILE_PATH',None);expected={}
    (out/'fingerprints').mkdir()
    for threads,affinity in [(1,[2]),(2,[2,3]),(4,list(range(4))),(6,list(range(6))),(12,list(range(12)))]:
        for name in ['single_small','single_large','conformers_56','ensemble_sdf','descriptors_10000','parse_ensembles_1000','db_build_ensembles']:
            work=manifest['workloads'][name]
            for pair in ([-1] if a.preflight_only else range(-1,8)):
                sides=['baseline','candidate'] if pair<0 or pair%2==0 else ['candidate','baseline']
                for side in sides:
                    folder=out/f't{threads}'/name/f'{pair}_{side}'
                    argv=[str(binaries[side]),*[x.replace('{run_dir}',str(folder)) for x in work['argv']]]
                    attempt={'threads':threads,'affinity':affinity,'workload':name,'pair':pair,'warmup':pair<0,'side':side,'status':'started'}
                    journal['attempts'].append(attempt);save()
                    metrics=harness.measure(argv,dict(env,RAYON_NUM_THREADS=str(threads)),set(affinity),folder,launcher)
                    attempt['metrics']=compact_metrics(metrics,folder);save()
                    assert metrics.get('exit_code',metrics.get('returncode',0))==0,attempt
                    scientific=harness.check_result(work,folder)
                    assert scientific==expected.setdefault(name,scientific),(name,side,threads)
                    fingerprint=retain_fingerprint(scientific,out/'fingerprints'/(name+'.json'))
                    attempt.update(status='passed',scientific=fingerprint);save()
            print(threads,name,flush=True)
    harness.verify_manifest(manifest)
    for path,digest in journal['bindings'].items():assert harness.digest(Path(path))==digest
    if a.preflight_only:
        journal.update(complete=True,interpretation='Functional preflight only. Incidental metrics are not controlled performance evidence.');save();return
    journal['host_postflight']=postflight(out);journal['complete']=True
    journal['quiet_endpoints']=journal['host_postflight']['quiet'] and journal['host_preflight']['quiet']
    journal['results']={}
    for t in [1,2,4,6,12]:
        for name in expected:
            rows=[r for r in journal['attempts'] if r['threads']==t and r['workload']==name and not r['warmup']]
            vals={s:[r['metrics']['wall_ns'] for r in rows if r['side']==s] for s in binaries}
            stats={s:{'values_ns':v,'median_ns':statistics.median(v),'MAD_ns':statistics.median(abs(x-statistics.median(v)) for x in v),'minimum_ns':min(v),'maximum_ns':max(v)} for s,v in vals.items()}
            stats['speedup']=statistics.median(vals['baseline'])/statistics.median(vals['candidate'])
            journal['results'][f'{name}/t{t}']=stats
    save()

if __name__=='__main__':main()
