"""Fresh untimed all-corpus evaluation of the frozen policy; never changes it."""
import gzip,hashlib,json,os,shutil,subprocess,sys
from pathlib import Path
from datetime import datetime,timezone
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
import numpy as np
import reference as independent
from fractions import Fraction as Q
R=Path(__file__).resolve().parent;H=R.parent;ROOT=H.parents[2]
sys.path.insert(0,str(H))
import compare
import morfeus_batch as mb
from morfeus import BuriedVolume


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def policy():
 frozen=json.loads((R/'policy_freeze.json').read_text());assert sha(H/frozen['policy'])==frozen['sha256'];return frozen


def reference_tool(path):
 elements,xyz,donor,ns,center,radii=mb.geometry(str(path))
 frames=[]
 for plane in ns:
  z=center-xyz[donor];z/=np.linalg.norm(z);v=xyz[plane]-center;x=v-np.dot(v,z)*z;x/=np.linalg.norm(x);y=np.cross(z,x)
  aligned=(xyz-center)@np.array([x,y,z]).T
  b=BuriedVolume(['H',*elements],np.vstack([np.zeros(3),aligned]),1,radii=np.r_[0.,radii],radius=3.5,density=.01,include_hs=False);b.octant_analysis()
  def ids(points):
   idx=np.rint((points+3.5)*31/7).astype(int)
   assert np.all((idx>=0)&(idx<=31))
   assert np.all(np.abs((points+3.5)*31/7-idx)<.5) # index identification, not science tolerance
   return idx[:,0]*1024+idx[:,1]*32+idx[:,2]
  grid=ids(b._sphere.points);buried=ids(b._buried_points)
  octs=[float(b.octants['buried_volume'][i]) for i in [0,1,2,3,7,6,5,4]]
  frames.append(dict(plane=plane,grid=grid,buried=buried,buried_volume=float(b.buried_volume),quadrants=[octs[i]+octs[i+4] for i in range(4)],octants=octs,near_vbur=sum(octs[4:]),far_vbur=sum(octs[:4])))
 return dict(donor=donor,neighbors=ns,frames=frames,values=mb.calculate(str(path)))


def main():
 os.sched_setaffinity(0,{2})
 frozen=policy();out=R/'evaluation_2';assert not out.exists();out.mkdir();(out/'raw').mkdir()
 # Seal numerical rules, reference/probe implementation and software before evaluating.
 freeze=json.loads((H/'freeze.json').read_text())
 for name,digest in freeze['source_sha256'].items():assert sha(ROOT/name)==digest
 for name in ['inputs','build','adapter']:(out/name).symlink_to((H/name).resolve(),target_is_directory=True)
 for name in ['input_manifest.json','freeze.json']:shutil.copy2(H/name,out/name)
 binary=H/'adapter/target/release/policy-grid-observer'
 shutil.copy2(binary,R/'policy-grid-observer')
 files=[R/'reference.py',R/'evaluate.py',R/'selftest.py',R/'selftest.json',R/'grid_append.rs',R/'policy-grid-observer',H/'adapter/policy_copied_volume.rs',H/'adapter/src/bin/policy-grid-observer.rs',H/'compare.py',H/'morfeus_batch.py',H/'build/stericx']
 assert (H/'adapter/policy_copied_volume.rs').read_bytes()==(ROOT/'src/geometry/buried_volume.rs').read_bytes()+(R/'grid_append.rs').read_bytes()
 receipt=dict(started_utc=datetime.now(timezone.utc).isoformat(),policy=frozen,files={str(p.relative_to(H)):sha(p) for p in files},classification_complete=False)
 write(out/'execution_freeze.json',receipt)
 compare.HERE=out;compare.main() # Fresh full14-descriptor residuals in all3 lanes.
 obs=[json.loads(x) for x in (out/'raw/stericx-observer.stdout').read_text().splitlines()]
 manifest=json.loads((H/'input_manifest.json').read_text())
 paths=[H/row['filename'] for row in manifest['corpus']]
 cmd=[str(R/'policy-grid-observer'),*map(str,paths)]
 p=subprocess.run(cmd,capture_output=True);write(out/'raw/probe-command.json',dict(argv=cmd,returncode=p.returncode));(out/'raw/probe.stdout.gz').write_bytes(gzip.compress(p.stdout));(out/'raw/probe.stderr').write_bytes(p.stderr);p.check_returncode()
 native=[json.loads(x) for x in p.stdout.splitlines()];assert len(native)==len(paths)==56
 records=[];failures=[]
 for number,(path,n,o) in enumerate(zip(paths,native,obs,strict=True)):
  ref=independent.certify(path);m=reference_tool(path);issues=[];frames=[];arrays={}
  def issue(kind,**kw):issues.append(dict(kind=kind,**kw))
  for key in ['donor','neighbors']:
   if n[key]!=ref[key] or m[key]!=ref[key]:issue('mapping',field=key,native=n[key],morfeus=m[key],reference=ref[key])
  idx=np.array([q['index'] for q in n['points']]);ni=idx[:,0]*1024+idx[:,1]*32+idx[:,2];order=np.argsort(ni)
  native_grid_ok=len(ni)==len(set(ni.tolist()))==len(independent.IDS) and np.array_equal(ni[order],independent.IDS)
  if not native_grid_ok:issue('native_grid_membership')
  if native_grid_ok and not np.array_equal(np.array([p['region'] for p in n['points']])[order],independent.REGIONS):issue('native_regions')
  nfmodel,nvalues,nvolume=independent.reduce_counts(ref['frames'],24)
  mfmodel,mvalues,mvolume=independent.reduce_counts(ref['frames'],53)
  if Q(n['sphere_volume'])!=Q(nvolume):issue('native_sphere_constant',actual=n['sphere_volume'],expected=nvolume)
  for k in nvalues:
   if Q(o[k])!=Q(nvalues[k]):issue('native_scalar_rounding',field=k,actual=o[k],expected=nvalues[k])
   if Q(m['values'][k])!=Q(mvalues[k]):issue('morfeus_scalar_rounding',field=k,actual=m['values'][k],expected=mvalues[k])
  for j,(rf,nf,mf,nmodel,mmodel) in enumerate(zip(ref['frames'],n['frames'],m['frames'],nfmodel,mfmodel,strict=True)):
   if rf['plane']!=nf['plane'] or rf['plane']!=mf['plane']:issue('plane_mapping',frame=j)
   mi=np.sort(mf['grid']);mgok=len(mi)==len(set(mi.tolist()))==len(independent.IDS) and np.array_equal(mi,independent.IDS)
   if not mgok:issue('morfeus_grid_membership',frame=j)
   nm=np.array(nf['occupied'])[order] if native_grid_ok else np.zeros(len(independent.IDS),dtype=bool)
   mm=np.isin(independent.IDS,mf['buried'])
   unknown=rf['unknown'];nd=(nm!=rf['occupied'])&~unknown;md=(mm!=rf['occupied'])&~unknown
   if np.any(unknown):issue('uncertified_reference_predicate',frame=j,point_ids=independent.IDS[unknown].tolist())
   for tool,mask in [('stericx',nd),('morfeus',md)]:
    if np.any(mask):issue('occupancy_disagreement',tool=tool,frame=j,plane=rf['plane'],point_ids=independent.IDS[mask].tolist(),clearance_intervals=list(zip(rf['clearance_lo'][mask].tolist(),rf['clearance_hi'][mask].tolist(),strict=True)))
   for tool,actual,model in [('stericx',nf,nmodel),('morfeus',mf,mmodel)]:
    for key,value in model.items():
     av=actual[key];av=av if isinstance(av,list) else [av];ev=value if isinstance(value,list) else [value]
     if any(Q(a)!=Q(b) for a,b in zip(av,ev,strict=True)):issue('plane_rounding',tool=tool,frame=j,field=key,actual=actual[key],expected=value)
   for key,value in [('reference',rf['occupied']),('uncertified',unknown),('native',nm),('morfeus',mm),('clearance_lo',rf['clearance_lo']),('clearance_hi',rf['clearance_hi'])]:arrays[f'frame{j}_{key}']=value
   frames.append(dict(plane=rf['plane'],uncertified=int(unknown.sum()),native_differences=int(nd.sum()),morfeus_differences=int(md.sum()),reference_counts=rf['occupied_counts'],reference_populations=rf['grid_populations'],native_model=nmodel,morfeus_model=mmodel))
  np.savez_compressed(out/f'reference-{number:02d}.npz',lattice_ids=independent.IDS,**arrays)
  record=dict(filename=manifest['corpus'][number]['filename'],donor=ref['donor'],neighbors=ref['neighbors'],frames=frames,native_model=nvalues,morfeus_model=mvalues,issues=issues)
  records.append(record);failures.extend(dict(filename=record['filename'],**q) for q in issues)
  write(out/'volume_records.json',records)
  print(json.dumps(dict(completed=number+1,issues=len(issues),filename=record['filename'])),flush=True)
 policy()
 for p in files:assert sha(p)==receipt['files'][str(p.relative_to(H))]
 result=dict(completed_utc=datetime.now(timezone.utc).isoformat(),policy_sha256=frozen['sha256'],structures=56,frames=168,points_per_frame=len(independent.IDS),certified_point_decisions=56*3*len(independent.IDS),passed=not failures,failures=failures,
             admitted_descriptors=list(nvalues) if not failures else [],excluded=['sterimol_l','sterimol_b1','sterimol_b5','pyr_p','pyr_alpha','combined','continuum_volume_accuracy'],prior_timings_admissible=False)
 write(out/'scientific_gate.json',result)
 print(json.dumps({k:v for k,v in result.items() if k!='failures'},indent=2),flush=True)

if __name__=='__main__':main()
