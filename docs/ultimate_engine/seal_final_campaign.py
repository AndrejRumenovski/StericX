#!/usr/bin/env python3
"""Seal final reports and raw evidence; verify historical scientific hashes unchanged."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
BASE=ROOT/'.stericx/ultimate_engine'
OUT=HERE/'final_campaign_evidence.json'
PRODUCTION=['src','tests','examples','Cargo.toml','Cargo.lock','pyproject.toml','python']

def digest(p):
    with p.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def production_check():
    diff=subprocess.check_output(['git','diff','f9c952e2','--',*PRODUCTION],cwd=ROOT)
    assert not diff,'Production source differs from final validated engine'
    return True
def quiet(d):
    if isinstance(d,dict):
        if 'quiet' in d:assert d['quiet'] is True
        for v in d.values():quiet(v)
    elif isinstance(d,list):
        for v in d:quiet(v)
def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--verify',action='store_true');a=parser.parse_args()
    production_check()
    if a.verify:
        value=read(OUT)
        for r in value['records']:
            p=ROOT/r['path'];assert p.stat().st_size==r['bytes'] and digest(p)==r['sha256'],str(p)
        print(json.dumps({'verified':True,'records':len(value['records']),'bytes':sum(r['bytes'] for r in value['records']),
                          'seal_sha256':digest(OUT),'production_diff_empty':True}));return
    assert not OUT.exists(),'Never overwrite a previous seal'
    original=read(HERE/'final_scientific_evidence.json');assert len(original['records'])==65
    for r in original['records']:
        p=ROOT/r['path'];assert digest(p)==r['sha256'],r['path']
    assert original['legacy_reference_preservation_pass']
    paths={ROOT/r['path'] for r in original['records']}
    for folder in [BASE/'final_window_6',BASE/'sterimol_competitive_preflight',BASE/'sterimol_competitive_timing']:
        paths.update(p for p in folder.rglob('*') if p.is_file() and not p.is_symlink())
    # Include all campaign documentation/source, not interpreter caches or build products.
    paths.update(p for p in HERE.rglob('*') if p.is_file() and not p.is_symlink()
                 and not any(part in ('__pycache__','target') for part in p.relative_to(HERE).parts) and p!=OUT)
    for name in ['final_window_6_guard.json','final_window_6_service.json','final_window_6.log',
                 'sterimol_competitive_timing_guard.json','sterimol_competitive_timing_service.json','sterimol_competitive_timing.log']:
        paths.add(BASE/name)
    build=BASE/'sterimol_competitive_build'
    for name in ['receipt.json','started.json','build.log','cases.json','sterimol_workflow','observer/Cargo.toml',
                 'observer/Cargo.lock','observer/src/bin/sterimol_workflow.rs']:paths.add(build/name)
    for number in range(1,6):
        for suffix in ['_guard.json','_service.json','.log']:
            p=BASE/f'final_window_{number}{suffix}'
            if p.is_file():paths.add(p)
        folder=BASE/f'final_window_{number}'
        if folder.exists():paths.update(folder.rglob('summary.json'));paths.update(folder.rglob('receipt.json'))
    archive=read(BASE/'final_window_6/controller_archive/receipt.json')
    replacements={r['original_path']:ROOT/r['archive_path'] for r in archive['files']}
    for r in archive['files']:assert digest(ROOT/r['archive_path'])==r['sha256']
    service=read(BASE/'final_window_6_service.json')
    for p,sha in service['controller_hashes'].items():assert digest(replacements[p])==sha
    receipt=read(BASE/'final_window_6/receipt.json');assert receipt['complete'];quiet(receipt)
    assert len(receipt['runs'])==8 and all(r['status']=='completed' and r['returncode']==0 for r in receipt['runs'])
    # Bind timing sources and binaries to their pre-run receipts, using archived
    # controller bytes when an explicitly post-run report helper was updated.
    def bindings(d):
        if isinstance(d,dict):
            for key,val in d.items():
                if isinstance(key,str) and key.startswith(str(ROOT)) and isinstance(val,str) and re.fullmatch('[0-9a-f]{64}',val):
                    p=Path(key)
                    if digest(p)!=val:
                        p=replacements.get(key,p)
                    assert digest(p)==val,key
                    paths.add(p)
                bindings(val)
        elif isinstance(d,list):
            for v in d:bindings(v)
    names=['cli_ab/summary.json','kernel_ab/summary.json','kernel_profile/journal.json','cli_profile/journal.json',
           'python_ab/summary.json','python_native_bridge/summary.json','accuracy_cost/summary.json','competitive/summary.json']
    for name in names:
        d=read(BASE/'final_window_6'/name);assert d['complete'];quiet(d);bindings(d)
    sg=read(BASE/'sterimol_competitive_preflight/summary.json');st=read(BASE/'sterimol_competitive_timing/summary.json')
    assert sg['complete'] and st['complete'];quiet(st);bindings(sg);bindings(st)
    assert sg['admitted']=={'sampled361':True,'continuous':True,'morfeus':False,'dbstep':False}
    for name in ['final_window_6_guard.json','sterimol_competitive_timing_guard.json']:
        g=read(BASE/name);assert g['resumed'] and g['returncode']==0 and g['identity_after']['state'] not in ('T','t')
    active=[]
    for p in Path('/proc').iterdir():
        if p.name.isdigit():
            try:args=(p/'cmdline').read_bytes().split(b'\0')
            except OSError:continue
            if any(arg.endswith((b'/quiet_process_guard.py',b'/run_performance_window.py',b'/sterimol_competitive.py')) for arg in args):active.append(p.name)
    assert not active,active
    # Validate local links in every requested current report and final tables/audit.
    reports=['BASELINE.md','SCIENTIFIC_AUDIT.md','PROFILING_REPORT.md','NUMERICAL_ACCURACY.md',
             'VALIDATION_REPORT.md','COMPETITIVE_BENCHMARKS.md','FINAL_REPORT.md','FINAL_TABLES.md','COMPLETION_AUDIT.md']
    for name in reports:
        for target in re.findall(r'\]\(([^)]+)\)',(HERE/name).read_text()):
            if target.startswith(('http:','https:','#')):continue
            p=HERE/target.split('#')[0]
            assert p==OUT or p.exists(),(name,target)
    records=[]
    for p in sorted(paths):
        assert p.is_file() and not p.is_symlink(),str(p)
        records.append({'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':digest(p)})
    value={'scope':'Final completed campaign on available hardware; unavailable counters/platforms and failed external Sterimol targets explicitly retained.',
           'sealed_at_unix':time.time(),'head_before_report_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
           'production_commit':'f9c952e2','production_diff_empty':True,'original_scientific_records_reverified':65,
           'complete_final_stages':8,'all_required_host_gates_passed':True,'ui_restore_receipts_passed':True,
           'active_measurement_processes':active,'sterimol_admission':sg['admitted'],'records':records,
           'limitations':['Hardware counters/real bandwidth unavailable','No physical ARM/macOS/Windows performance',
                          'No controlled disk-cache-cold timing','No equal-accuracy speed ratio for failed external Sterimol methods',
                          'No universal numerical accuracy or fastest-engine claim']}
    with OUT.open('x') as f:json.dump(value,f,indent=2);f.write('\n')
    print(json.dumps({'sealed':True,'records':len(records),'bytes':sum(r['bytes'] for r in records),'production_diff_empty':True}))

if __name__=='__main__':main()
