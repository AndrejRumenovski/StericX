"""Compare source-bound Rust observations with independent frozen references."""
from pathlib import Path
import argparse, hashlib, json
ROOT=Path(__file__).resolve().parents[2]

def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def main():
 p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args()
 requests=ROOT/'.stericx/occupancy_precision/analytical_chord_requests.json'
 frozen=ROOT/'.stericx/occupancy_precision/c1_public_analytical.jsonl'
 qs=json.loads(requests.read_text());rs=[json.loads(s) for s in (a.run/'observed.jsonl').read_text().splitlines()]
 failed_ids=json.loads((ROOT/'docs/occupancy_precision/c1_validation.json').read_text())['analytical']['unmet']
 old={r['id']:r for r in map(json.loads,frozen.read_text().splitlines())}
 assert len(qs)==len(rs)==167
 failures=[];cases=[]
 for q,r in zip(qs,rs,strict=True):
  assert q['id']==r['id']
  lo,hi=q['independent_reference']
  valid=r['lower']<=lo<=hi<=r['upper'] and r['lower']<=r['volume']<=r['upper']
  valid &= r['target_met']==(r['absolute_error_bound']<=q['tolerance']) and r['cells']<=q['budget']
  if not valid:failures.append(q['id'])
  if q['id'] in failed_ids:
   cases.append(dict(request=q,baseline=old[q['id']],candidate=r,valid=valid))
 summary=dict(cases=len(qs),invalid=failures,original_failed_count=len(failed_ids),
   resolved=sum(c['candidate']['target_met'] for c in cases),
   remaining=[c['request']['id'] for c in cases if not c['candidate']['target_met']],
   original_cases=cases,inputs_sha256=sha(requests),baseline_sha256=sha(frozen),
   observed_sha256=sha(a.run/'observed.jsonl'),
   source={str(p.relative_to(ROOT)):sha(p) for p in sorted((ROOT/'src').rglob('*.rs'))})
 (a.run/'validation.json').write_text(json.dumps(summary,indent=2)+'\n')
 print(json.dumps({k:v for k,v in summary.items() if k not in ['original_cases','source']},indent=2))
 if failures:raise SystemExit(1)
if __name__=='__main__':main()
