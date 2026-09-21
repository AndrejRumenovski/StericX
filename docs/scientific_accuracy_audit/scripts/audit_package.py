"""Seal/verify audit artifacts without changing any production code or captures.

The final inventory detects later changes. Separate checks below anchor it to
the original pre-audit snapshot and pre-comparison captures, so a new inventory
is not mistaken for proof that historical artifacts never changed.
"""
from pathlib import Path
from datetime import datetime, UTC
import argparse, hashlib, json, subprocess

A=Path(__file__).resolve().parents[1]
REPO=A.parents[1]
FINAL=A/'manifest_final.json'
EXCLUDE={'manifest_final.json','reproducibility/package_verification.json'}

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()

def check_original():
    failures=[];checks=[]
    def check(p,expected,label):
        ok=p.is_file() and sha(p)==expected
        checks.append({'path':str(p.relative_to(A)) if p.is_relative_to(A) else str(p),'check':label,'ok':ok})
        if not ok:failures.append(checks[-1])
    initial=json.loads((A/'manifest_initial.json').read_text())
    for row in [initial['binary'],*initial['repository_files'].values(),*initial['reference_tools']['morfeus-ml']['python_source_files'].values()]:check(A/row['frozen_path'],row['sha256'],'initial freeze')
    for rel in ['observer/manifest.json','kinetics/inputs_manifest.json','kinetics/sut_manifest.json']:
        for row in json.loads((A/rel).read_text())['files']:check(A/row['path'],row['sha256'],rel)
    for rel in ['models/manifest.json']:
        for name,digest in json.loads((A/rel).read_text()).items():
            if isinstance(digest,str) and len(digest)==64:
                if name=='models/sources/downloads.json' and sha(A/name)!=digest:
                    # An acquisition journal grew by two records after its first
                    # snapshot. Verify recovered ORIGINAL bytes, then establish
                    # that all original records remain unchanged in the journal.
                    original=A/'models/sources/downloads_initial_manifest_snapshot.json'
                    check(original,digest,'documented append-only acquisition journal amendment')
                    old=json.loads(original.read_text());current=json.loads((A/name).read_text())
                    if not isinstance(old,list) or current[:len(old)]!=old:
                        failures.append({'check':'acquisition journal unchanged original prefix','ok':False})
                    if not (A/'models/provenance_amendment_downloads.json').is_file():
                        failures.append({'check':'acquisition journal amendment must be documented','ok':False})
                else:check(A/name,digest,rel)
    g=json.loads((A/'geometry/sut_freeze.json').read_text())
    check(A/'geometry/sut_outputs.jsonl',g['output_sha256'],'geometry pre-comparison capture')
    check(A/'geometry/requests.jsonl',g['request_sha256'],'geometry frozen requests')
    for name,digest in g['input_hashes'].items():check(A/'geometry'/name,digest,'geometry frozen input')
    for folder in ['focused','alignment','numerical']:
        base=A/'geometry'/folder;p=base/'freeze.json'
        if not p.exists():continue
        d=json.loads(p.read_text())
        for key,name in [('input_sha256','inputs.json'),('request_sha256','requests.jsonl'),('output_sha256','sut_outputs.jsonl'),('requests_sha256','requests.jsonl'),('sut_outputs_sha256','sut_outputs.jsonl')]:
            if key in d:check(base/name,d[key],str(p.relative_to(A)))
        for name,digest in d.items():
            if name.endswith(('.json','.jsonl')) and isinstance(digest,str):check(base/name,digest,str(p.relative_to(A)))
    rotation=A/'geometry/alignment/rotation'
    d=json.loads((rotation/'sut_freeze.json').read_text())
    check(rotation/'requests.jsonl',d['requests_sha256'],'rotation frozen requests')
    check(rotation/'sut_outputs.jsonl',d['sut_sha256'],'rotation frozen observations')
    check(A/'frozen/bin/stericx-audit-observer',d['observer_sha256'],'rotation observer identity')
    for base in [A/'kraken/full_analytic_reference',A/'kraken/morfeus_outliers',A/'kraken/delta_outliers']:
        d=json.loads((base/'manifest_frozen.json').read_text())
        expected=d.get('raw_sha256',d.get('reference_raw_sha256'))
        if expected is None:raise ValueError(f'Unrecognized reference manifest {base}: {list(d)}')
        check(base/'reference_raw.jsonl',expected,'frozen independent reference')
        for key,name in [('plan_sha256','plan_frozen.json'),('selection_manifest_sha256','selection_frozen.json'),('selection_sha256','selection_frozen.json')]:
            if key in d:check(base/name,d[key],'independent reference plan/selection')
    edge=json.loads((A/'models/inputs/evaluate_edges_manifest_before_execution.json').read_text())
    for name,digest in edge['inputs'].items():check(A/name,digest,'evaluate edge inputs frozen before execution')
    check(A/'frozen/bin/stericx',edge['binary_sha256'],'evaluate edge executable identity')
    identity=json.loads((A/'models/inputs/formula_identity_manifest_before_execution.json').read_text())
    for name,digest in identity['inputs'].items():check(A/name,digest,'formula-identity inputs frozen before helper execution')
    check(A/'frozen/repository/studies/study_007_crosscoupling.py',identity['frozen_study_sha256'],'formula-identity frozen study implementation')
    for rel in ['models/manifest_finalize.json','kraken/manifest_final.json']:
        consolidated=json.loads((A/rel).read_text())
        for name,digest in consolidated['files'].items():check(A/name,digest,rel+' post-analysis consolidation (not a pre-comparison freeze)')
    for metadata in (A/'kraken/primary').rglob('*.metadata.json'):
        row=json.loads(metadata.read_text())
        if row.get('sha256'):
            body=metadata.with_name(metadata.name.replace('.metadata.json','.body'))
            check(body,row['sha256'],'primary API response acquisition hash')
    d=json.loads((A/'kraken/sut/manifest_frozen.json').read_text())
    for name,row in d['raw_output_files'].items():check(A/'kraken/sut'/name,row['sha256'],'Kraken pre-comparison capture')
    check(A/'kraken/prepared_all/requests.jsonl',d['request_sha256'],'Kraken frozen requests')
    # Key names are retained from the original capture producer.
    for key,name in [('stdout_sha256','stdout.jsonl'),('output_sha256','stdout.jsonl'),('stderr_sha256','stderr.txt')]:
        if key in d:check(A/'kraken/sut'/name,d[key],'Kraken pre-comparison capture')
    tracked=subprocess.run(['git','diff','--name-only','HEAD','--'],cwd=REPO,text=True,capture_output=True,check=True).stdout.splitlines()
    if tracked:failures.append({'check':'production tree unchanged','changed_tracked_files':tracked})
    head=subprocess.run(['git','rev-parse','HEAD'],cwd=REPO,text=True,capture_output=True,check=True).stdout.strip()
    if head!=initial['git_commit']:failures.append({'check':'current commit remains frozen commit','current':head,'frozen':initial['git_commit']})
    return {'checked':len(checks),'failures':failures,'tracked_changes':tracked,'initial_commit':initial['git_commit']}

def files():
    return sorted(p for p in A.rglob('*') if p.is_file() and p.relative_to(A).as_posix() not in EXCLUDE and '__pycache__' not in p.parts and p.suffix!='.pyc')

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('mode',choices=['seal','verify','check-original']);args=parser.parse_args()
    original=check_original()
    if original['failures']:
        print(json.dumps(original,indent=2));raise SystemExit('Frozen evidence or production tree differs; investigate before sealing')
    if args.mode=='check-original':print(json.dumps(original,indent=2));return
    if args.mode=='seal':
        required=['SCIENTIFIC_ACCURACY_AUDIT.md','CLAIMS.md','COMPLETION_CHECKLIST.md','REPRODUCE.md','geometry/REPORT.md','models/REPORT.md','kraken/REPORT.md','kinetics/KINETICS_CONFORMERS.md']
        missing=[s for s in required if not (A/s).is_file()]
        if missing:raise SystemExit('Required artifacts missing: '+', '.join(missing))
        inventory=[{'path':p.relative_to(A).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)} for p in files()]
        FINAL.write_text(json.dumps({'schema':1,'sealed_utc':datetime.now(UTC).isoformat(),'initial_commit':original['initial_commit'],'policy':'Final immutable-artifact inventory; original phase hashes checked independently. Python bytecode caches and this self-referential inventory/verification log are excluded. New replay logs may be added later without changing frozen evidence.','original_checks':original,'files':inventory},indent=2)+'\n')
    final=json.loads(FINAL.read_text());failures=[]
    for row in final['files']:
        p=A/row['path']
        if not p.is_file() or sha(p)!=row['sha256']:failures.append(row['path'])
    result={'verified_utc':datetime.now(UTC).isoformat(),'manifest_sha256':sha(FINAL),'files_checked':len(final['files']),'bytes_checked':sum(row['bytes'] for row in final['files']),'original_checks':original,'final_hash_mismatches':failures,'passed':not failures}
    dest=A/'reproducibility/package_verification.json';dest.parent.mkdir(exist_ok=True);dest.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    if failures:raise SystemExit(1)

if __name__=='__main__':main()
