"""Prospective common-grid admission of a fixed twelve-conformer radius profile."""
import gzip,importlib.util,json,os,subprocess,sys
from fractions import Fraction as Q
from pathlib import Path
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import numpy as np
from campaign import ROOT,EVIDENCE,sha,save
import radius_worker as worker
out=EVIDENCE/'competitive_radius_gate';out.mkdir(exist_ok=False)
reference_source=ROOT/'docs/benchmarks/morfeus_current/revision_policy/reference.py'
catalog=ROOT/'docs/benchmarks/morfeus_current/input_manifest.json';corpus=json.loads(catalog.read_text())['corpus']
selection=list(range(0,48,4));paths=[catalog.parent/corpus[i]['filename'] for i in selection]
native=EVIDENCE/'competitive_radius_build/radius_workflow'
receipt={'complete':False,'passed':False,'selection':selection,'radii':worker.RADII,'density':'0.01*(R/3.5)^3; grid side count 32 at every radius','scope':'Matched finite-grid radius profile; distinct from a constant-density profile. DBSTEP occupancy kernel receives common points, not its default grid. Twelve prospectively selected conformers; do not generalize to all radii.','bindings':{str(p):sha(p) for p in [Path(__file__),Path(worker.__file__),reference_source,catalog,native,*paths]},'rows':[],'failures':[]}
save(out/'summary.json',receipt)
refs={}
for radius in worker.RADII:
    text=reference_source.read_text();needle='Q(7*(2*int(t)-31),62)';assert text.count(needle)==1
    text=text.replace(needle,f"Q('{radius}')*Q(2*int(t)-31,31)")
    text=text.replace("r3=Q('3.5')**3",f"r3=Q('{radius}')**3")
    path=out/f'reference_{radius}.py';path.write_text(text);spec=importlib.util.spec_from_file_location(f'radius_reference_{radius}',path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);refs[radius]=module;receipt['bindings'][str(path)]=sha(path)
command=[str(native),'0',*map(str,paths)];r=subprocess.run(command,capture_output=True,check=True);(out/'native.json').write_bytes(r.stdout);(out/'native.stderr').write_bytes(r.stderr);observed=json.loads(r.stdout)['rows'];assert len(observed)==len(paths)
for index,(path,public) in enumerate(zip(paths,observed,strict=True)):
    data=subprocess.run([str(native),'probe',str(path)],capture_output=True,check=True);(out/f'probe_{index}.jsonl.gz').write_bytes(gzip.compress(data.stdout,mtime=0));(out/f'probe_{index}.stderr').write_bytes(data.stderr)
    probes=list(map(json.loads,data.stdout.splitlines()));assert len(probes)==3
    external={tool:worker.calculate((tool,str(path)),details=True) for tool in ['morfeus','dbstep']};issues=[]
    for ri,(radius,probe) in enumerate(zip(worker.RADII,probes,strict=True)):
        ref=refs[radius];proof=ref.certify(path);n=probe['probe'];nf,nv,nvolume=ref.reduce_counts(proof['frames'],24);mf,mv,mvolume=ref.reduce_counts(proof['frames'],53)
        def check(value,kind):
            if not value:issues.append({'radius':radius,'issue':kind})
        check(public['profile'][ri]['values']==nv,'native_reduction');check(n['sphere_volume']==nvolume,'native_sphere_volume')
        for tool,(row,frames) in external.items():check(row['profile'][ri]['values']==mv,tool+'_reduction')
        indices=np.array([p['index'] for p in n['points']]);ids=indices[:,0]*1024+indices[:,1]*32+indices[:,2];order=np.argsort(ids);check(np.array_equal(ids[order],ref.IDS),'native_grid')
        check(np.array_equal(np.array([p['region'] for p in n['points']])[order],ref.REGIONS),'native_regions')
        for fi,(reference,native_frame) in enumerate(zip(proof['frames'],n['frames'],strict=True)):
            check(not np.any(reference['unknown']),f'unknown_{fi}');check(np.array_equal(np.array(native_frame['occupied'])[order],reference['occupied']),f'native_mask_{fi}')
            for tool,(_,frames) in external.items():
                frame=frames[ri][fi]
                def ids_of(points):
                    ijk=np.rint((points+radius)*31/(2*radius)).astype(int);return ijk[:,0]*1024+ijk[:,1]*32+ijk[:,2]
                check(np.array_equal(np.sort(ids_of(frame['points'])),ref.IDS),f'{tool}_grid_{fi}')
                check(np.array_equal(np.isin(ref.IDS,ids_of(frame['occupied'])),reference['occupied']),f'{tool}_mask_{fi}')
    receipt['rows'].append({'path':str(path),'native':public,'morfeus':external['morfeus'][0],'dbstep':external['dbstep'][0],'issues':issues})
    receipt['failures'] += [{'path':str(path),**issue} for issue in issues];save(out/'summary.json',receipt);print(index+1,len(issues),flush=True)
receipt['complete']=True;receipt['passed']=not receipt['failures'];receipt['point_decisions']=len(paths)*3*3*len(next(iter(refs.values())).IDS)
for path,digest in receipt['bindings'].items():assert sha(Path(path))==digest
save(out/'summary.json',receipt)
if not receipt['passed']:raise SystemExit(1)
