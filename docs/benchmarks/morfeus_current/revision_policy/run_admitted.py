"""Fresh paired campaign, callable only after the frozen interval policy passes.

Reuse the archived measurement runner; all original timings remain untouched.
"""
import hashlib,json,os,shutil,sys
from pathlib import Path
R=Path(__file__).resolve().parent;H=R.parent;ROOT=H.parents[2]
sys.path.insert(0,str(H))
import benchmark as b
C=R/'admitted_campaign'


def verify():
 policy=b.read(R/'policy_freeze.json');assert b.sha(H/'EQUIVALENCE_POLICY.md')==policy['sha256']
 evaluation=b.read(R/'evaluation_2/scientific_gate.json')
 assert evaluation['passed'] and not evaluation['failures'] and evaluation['policy_sha256']==policy['sha256']
 frozen=b.read(H/'freeze.json')
 for name,digest in frozen['source_sha256'].items():assert b.sha(ROOT/name)==digest
 assert b.sha(H/'build/stericx')==frozen['executable_sha256']
 execution=b.read(R/'evaluation_2/execution_freeze.json')
 for name,digest in execution['files'].items():assert b.sha(H/name)==digest,name
 gate=b.read(C/'scientific_gate.json')
 assert gate['evaluation_sha256']==b.sha(R/'evaluation_2/scientific_gate.json')
 return frozen,gate


def main():
 # Refuse to set up, let alone time, an unadmitted descriptor.
 ev=b.read(R/'evaluation_2/scientific_gate.json');assert ev['passed'] and not ev['failures']
 assert not C.exists();C.mkdir();(C/'raw').mkdir();(C/'build').mkdir()
 for name in ['inputs','adapter','environment']:(C/name).symlink_to((H/name).resolve(),target_is_directory=True)
 for name in ['freeze.json','input_manifest.json','benchmark.py','morfeus_batch.py']:shutil.copy2(H/name,C/name)
 shutil.copy2(H/'build/stericx',C/'build/stericx')
 for name in ['stericx-observer.stdout','morfeus-reference.json']:shutil.copy2(R/'evaluation_2/raw'/name,C/'raw'/name)
 gate=ev|dict(descriptors=ev['admitted_descriptors'],numerical_tolerance=None,evaluation_sha256=b.sha(R/'evaluation_2/scientific_gate.json'),policy_freeze=b.read(R/'policy_freeze.json'))
 b.write(C/'scientific_gate.json',gate)
 b.HERE=C;b.ROOT=ROOT;b.verify_freeze=verify
 verify();b.prepare()
 # Check actual timed-driver preflight against independent rational predictions,
 # in addition to checking every worker configuration against the fresh comparison.
 rows=b.read(R/'evaluation_2/volume_records.json');pre=b.read(C/'preflight.json')
 for tool,key in [('stericx','native_model'),('morfeus','morfeus_model')]:
  assert pre['expected'][tool]==[r[key] for r in rows]
 pre['identities'].update({os.path.relpath(R/'run_admitted.py',C):b.sha(R/'run_admitted.py'),os.path.relpath(H/'EQUIVALENCE_POLICY.md',C):b.sha(H/'EQUIVALENCE_POLICY.md'),os.path.relpath(R/'evaluation_2/scientific_gate.json',C):b.sha(R/'evaluation_2/scientific_gate.json')})
 b.write(C/'preflight.json',pre)
 b.run();b.summarize()
 # Separate RSS collection occurs only after uninstrumented timed pairs finish.
 import memory
 memory.HERE=C;memory.verify_freeze=verify
 memory.main()
 verify()
 b.write(C/'completion.json',dict(completed_utc=b.now(),passed=True,policy_sha256=gate['policy_sha256'],no_prior_timing_reused=True))


if __name__=='__main__':main()
