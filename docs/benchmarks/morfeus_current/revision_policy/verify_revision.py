"""Standard-library artifact verification for the prospectively admitted revision.

This checks retained certificates and measurements; it does not recompute geometry.
"""
import argparse,ast,gzip,hashlib,json,os,struct,tarfile,zipfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
if HERE.name=='revision_policy':HERE=HERE.parent
R=HERE/'revision_policy';C=R/'admitted_campaign'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
EXCLUDED={'target','venv','__pycache__','repeated_10000'}
EXCLUDED_FILES={'evidence_manifest.json','evidence_manifest.json.sha256','verification.json'}

def inventory():
 files=[];links=[]
 for base,dirs,names in os.walk(HERE,followlinks=False):
  dirs[:]=[d for d in dirs if d not in EXCLUDED]
  for name in sorted(dirs+names):
   p=Path(base)/name;rel=p.relative_to(HERE).as_posix()
   if p.is_symlink():links.append(dict(path=rel,target=os.readlink(p)));continue
   if not p.is_file() or rel in EXCLUDED_FILES:continue
   files.append(dict(path=rel,bytes=p.stat().st_size,sha256=sha(p)))
 return dict(files=sorted(files,key=lambda r:r['path']),symlinks=sorted(links,key=lambda r:r['path']))


def npy(z,name):
 data=z.read(name+'.npy');assert data[:6]==b'\x93NUMPY'
 if data[6]==1:length=struct.unpack('<H',data[8:10])[0];offset=10
 else:length=struct.unpack('<I',data[8:12])[0];offset=12
 header=ast.literal_eval(data[offset:offset+length].decode());return header,data[offset+length:]


def verify():
 manifest=read(HERE/'evidence_manifest.json');assert sha(HERE/'evidence_manifest.json')==(HERE/'evidence_manifest.json.sha256').read_text().strip()
 actual=inventory();assert manifest['files']==actual['files'],'Missing, modified or extra published file'
 assert manifest['symlinks']==actual['symlinks']
 revision=read(R/'revision.json');assert sha(R/'prior_publication.tar.gz')==revision['prior_archive_sha256']
 with tarfile.open(R/'prior_publication.tar.gz') as t:
  data=t.extractfile('evidence_manifest.json').read();assert hashlib.sha256(data).hexdigest()==revision['prior_manifest_sha256']
  prior=json.loads(data);expected={x['path']:x for x in prior['files']}
  for member in t:
   if member.name in expected:
    data=t.extractfile(member).read();row=expected.pop(member.name)
    assert len(data)==row['bytes'] and hashlib.sha256(data).hexdigest()==row['sha256']
  assert not expected
 policy=read(R/'policy_freeze.json');assert sha(HERE/'EQUIVALENCE_POLICY.md')==policy['sha256']
 gate=read(R/'evaluation_2/scientific_gate.json');assert gate['passed'] and not gate['failures']
 assert gate['policy_sha256']==policy['sha256']
 assert gate['structures']==56 and gate['frames']==168 and gate['certified_point_decisions']==2588544
 raw=read(C/'raw_timings.json');assert raw['status']=='complete' and not raw['exclusions']
 assert policy['frozen_utc']<read(R/'evaluation_1/execution_freeze.json')['started_utc']<read(R/'evaluation_2/execution_freeze.json')['started_utc']<gate['completed_utc']<raw['started_utc']
 assert read(R/'reference_revision.json')['policy_unchanged']
 for name,digest in read(R/'evaluation_2/execution_freeze.json')['files'].items():assert sha(HERE/name)==digest,name
 assert sha(HERE/'build/stericx')==read(HERE/'freeze.json')['executable_sha256']
 corpus=read(HERE/'input_manifest.json');assert len(corpus['corpus'])==56
 for row in corpus['corpus']:assert sha(HERE/row['filename'])==row['sha256']
 for i,row in enumerate(corpus['repeated_10000']):assert row['corpus_index']==i%56 and row['sha256']==corpus['corpus'][i%56]['sha256']
 assert len(corpus['repeated_10000'])==10000
 records=read(R/'evaluation_2/volume_records.json');assert len(records)==56 and all(not r['issues'] for r in records)
 for i,row in enumerate(records):
  with zipfile.ZipFile(R/'evaluation_2'/f'reference-{i:02d}.npz') as z:
   for j in range(3):
    hdr,ref=npy(z,f'frame{j}_reference');assert hdr['descr']=='|b1' and hdr['shape']==(15408,)
    assert npy(z,f'frame{j}_native')[1]==ref==npy(z,f'frame{j}_morfeus')[1]
    assert not any(npy(z,f'frame{j}_uncertified')[1])
  assert all(f['uncertified']==f['native_differences']==f['morfeus_differences']==0 for f in row['frames'])
 pre=read(C/'preflight.json');assert pre['passed']
 for name,digest in pre['identities'].items():
  p=C/name
  if name=='adapter/target/release/volume-batch':p=HERE/'build/volume-batch'
  assert sha(p)==digest,name
 for tool,key in [('stericx','native_model'),('morfeus','morfeus_model')]:assert pre['expected'][tool]==[r[key] for r in records]
 expected={(n,w,t,p) for n in [56,1000,10000] for w in [1,2,4,6] for t in ['stericx','morfeus'] for p in [-1,0,1,2,3]}
 seen=set()
 for s in raw['samples']:
  key=(s['structures'],s['workers'],s['tool'],s['pair']);assert key not in seen;seen.add(key)
  assert s['status']=='passed' and s['metrics']['returncode']==0
  stem=C/s['raw_prefix'];data=gzip.decompress(stem.with_suffix('.stdout.gz').read_bytes())
  assert hashlib.sha256(data).hexdigest()==s['stdout_sha256'] and sha(stem.with_suffix('.stderr'))==s['stderr_sha256']
  assert read(stem.with_suffix('.metrics.json'))==s['metrics']
  rows=json.loads(data);assert len(rows)==s['structures']
  for i,r in enumerate(rows):
   assert {k:r[k] for k in gate['admitted_descriptors']}==pre['expected'][s['tool']][i%56]
   fmt='<f' if s['tool']=='stericx' else '<d'
   for field in gate['admitted_descriptors']:
    assert struct.pack(fmt,r[field])==struct.pack(fmt,pre['expected'][s['tool']][i%56][field])
 assert seen==expected
 result=read(HERE/'benchmark_results.json');assert result['publication_status']=='admitted_after_frozen_policy' and not result['prior_timings_used']
 assert all(result['acceptance_gates'].values())
 for r in result['summaries']:
  assert r['paired_speedup']['minimum']>1
  assert r['tools']['stericx']['wall_ns']['maximum']<r['tools']['morfeus']['wall_ns']['minimum']
 memory=read(C/'memory_results.json');assert memory['complete'] and len(memory['results'])==8
 check=dict(passed=True,published_files=len(manifest['files']),policy_sha256=policy['sha256'],scientific_structures=56,certified_point_decisions=2588544,new_timing_launches=120,new_measured_launches=96,new_warmups=24,memory_runs=8,prior_archive_verified=True,scope='Artifact integrity, certified masks/model evidence, new timing completeness/output agreement and policy-before-evaluation-before-timing ordering; not a new geometry calculation or timing rerun')
 write(HERE/'verification.json',check);print(json.dumps(check,indent=2))


def main():
 parser=argparse.ArgumentParser();parser.add_argument('--seal',action='store_true');args=parser.parse_args()
 if args.seal:
  old=HERE/'evidence_manifest.json'
  assert sha(old)==read(R/'revision.json')['prior_manifest_sha256'],'Only replace the explicitly archived predecessor'
  current=inventory();current.update(revision='prospective_policy',predecessor_sha256=sha(old),policy_sha256=read(R/'policy_freeze.json')['sha256'])
  write(old,current);(HERE/'evidence_manifest.json.sha256').write_text(sha(old)+'\n')
 verify()

if __name__=='__main__':main()
