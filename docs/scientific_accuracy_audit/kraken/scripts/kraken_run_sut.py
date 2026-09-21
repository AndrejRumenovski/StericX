#!/usr/bin/env python3
"""Run frozen public-API observer; freeze raw bytes before any comparisons."""
from __future__ import annotations
import datetime,hashlib,json,os,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2,sort_keys=True)+'\n')
def main():
    binary=ROOT.parent/'frozen/bin/stericx-audit-observer'
    assert sha(binary)=='9a8e778c661cc26f83e1252827a4452e8eb88b74f743fb5a853a24ee2101d957'
    prepared=ROOT/'prepared_all';manifest=json.loads((prepared/'manifest.json').read_text())
    for name,digest in manifest['prepared_files'].items():assert sha(prepared/name)==digest
    out=ROOT/'sut';out.mkdir(exist_ok=False)
    env=dict(os.environ,LC_ALL='C',TZ='UTC',RAYON_NUM_THREADS='1',OMP_NUM_THREADS='1')
    record={'started_utc':datetime.datetime.now(datetime.UTC).isoformat(),'argv':[str(binary)],'binary_sha256':sha(binary),'request_sha256':sha(prepared/'requests.jsonl'),'preparation_manifest_sha256':sha(prepared/'manifest.json'),'observer_manifest_sha256':sha(ROOT.parent/'observer/manifest.json'),'script_sha256':sha(Path(__file__)),'environment':{k:env[k] for k in ['LC_ALL','TZ','RAYON_NUM_THREADS','OMP_NUM_THREADS']},'scope':'Fresh full-precision public API observations on independently prepared primary DFT geometries; not an independent reference implementation.'}
    save(out/'manifest_started.json',record)
    with (prepared/'requests.jsonl').open('rb') as inp,(out/'stdout.jsonl').open('wb') as stdout,(out/'stderr.txt').open('wb') as stderr:
        result=subprocess.run([str(binary)],stdin=inp,stdout=stdout,stderr=stderr,env=env)
    record.update(returncode=result.returncode,finished_utc=datetime.datetime.now(datetime.UTC).isoformat(),raw_output_files={name:{'sha256':sha(out/name),'bytes':(out/name).stat().st_size} for name in ['stdout.jsonl','stderr.txt']},frozen_before_reference_comparison=True)
    save(out/'manifest_frozen.json',record)
    print(json.dumps(record,indent=2))
    if result.returncode:raise RuntimeError('Observer failed; retain frozen partial raw outputs')
if __name__=='__main__':main()
