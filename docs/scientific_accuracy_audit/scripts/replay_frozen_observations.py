"""Re-execute frozen SUT JSONL requests and compare complete output byte streams.

This checks reproducibility only, not scientific correctness. Original captures
are never overwritten. Failed replays are retained for investigation.
"""
from pathlib import Path
from datetime import datetime, UTC
import argparse, hashlib, json, os, subprocess, tempfile, time

A=Path(__file__).resolve().parents[1]
LANES={
    'geometry':(A/'geometry/requests.jsonl',A/'geometry/sut_outputs.jsonl'),
    'geometry_focused':(A/'geometry/focused/requests.jsonl',A/'geometry/focused/sut_outputs.jsonl'),
    'geometry_alignment':(A/'geometry/alignment/requests.jsonl',A/'geometry/alignment/sut_outputs.jsonl'),
    'geometry_alignment_rotation':(A/'geometry/alignment/rotation/requests.jsonl',A/'geometry/alignment/rotation/sut_outputs.jsonl'),
    'geometry_numerical':(A/'geometry/numerical/requests.jsonl',A/'geometry/numerical/sut_outputs.jsonl'),
    'kraken':(A/'kraken/prepared_all/requests.jsonl',A/'kraken/sut/stdout.jsonl'),
}
def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lanes',nargs='+',choices=list(LANES),default=list(LANES))
    args=parser.parse_args()
    binary=A/'frozen/bin/stericx-audit-observer'
    root=A/'reproducibility';root.mkdir(exist_ok=True)
    now=datetime.now(UTC);output=root/('replay_'+now.strftime('%Y%m%dT%H%M%SZ')+'.json')
    env={**os.environ,'LC_ALL':'C','TZ':'UTC','RAYON_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
    records=[]
    for name in args.lanes:
        requests,expected=LANES[name]
        if not requests.exists() or not expected.exists():raise SystemExit(f'Missing lane input/capture: {name}: {requests} / {expected}')
        print(f'Replaying {name}',flush=True)
        started=time.monotonic()
        with tempfile.NamedTemporaryFile(dir=root,prefix=f'{name}_',suffix='.jsonl',delete=False) as actual, tempfile.NamedTemporaryFile(dir=root,prefix=f'{name}_',suffix='.stderr',delete=False) as errors, requests.open('rb') as src:
            result=subprocess.run(['timeout','600',str(binary)],stdin=src,stdout=actual,stderr=errors,env=env)
            actual_path=Path(actual.name);error_path=Path(errors.name)
        row={'lane':name,'started_utc':now.isoformat(),'requests':str(requests.relative_to(A)),'request_sha256':digest(requests),'expected':str(expected.relative_to(A)),'expected_sha256':digest(expected),'actual_sha256':digest(actual_path),'binary_sha256':digest(binary),'returncode':result.returncode,'stderr':error_path.read_text(),'seconds':time.monotonic()-started}
        row['exact_bytes']=row['returncode']==0 and row['actual_sha256']==row['expected_sha256']
        if row['exact_bytes']:actual_path.unlink();error_path.unlink()
        else:row['preserved_actual']=str(actual_path.relative_to(A));row['preserved_stderr']=str(error_path.relative_to(A))
        records.append(row)
        output.write_text(json.dumps({'purpose':'Full byte reproducibility, not independent scientific validation','lanes':records},indent=2)+'\n')
        print(name,'MATCH' if row['exact_bytes'] else 'MISMATCH',f"{row['seconds']:.3f}s",flush=True)
    if not all(row['exact_bytes'] for row in records):raise SystemExit(1)
    print(output)

if __name__=='__main__':main()
