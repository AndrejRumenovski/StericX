"""Preserve every previously independently validated certificate after optimization."""
import argparse,json,subprocess
from pathlib import Path
from campaign import EVIDENCE,sha,save
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--engine',type=Path,required=True);p.add_argument('--output',required=True);a=p.parse_args()
out=EVIDENCE/a.output;out.mkdir(exist_ok=False)
lanes=[('analytical','analytical_chord_requests.json','c1_public_analytical.jsonl'),
       ('extended','c1_public_extended/requests.json','c1_public_extended/observed.jsonl'),
       ('coaxial','c1_coaxial/requests.json','c1_coaxial/observed.jsonl'),
       ('molecular','c1_molecular_broad_v3/requests.json','c1_molecular_broad_v3/observed.jsonl'),
       ('expanded_budget','c1_extended_budget/requests.json','c1_extended_budget/observed.jsonl')]
receipt={'engine_sha256':sha(a.engine),'controller_sha256':sha(Path(__file__)),'scope':'Exact JSON scientific-field comparison, including signed zero representations; only diagnostic elapsed_ns and allocator diagnostics are excluded. Reference containment results transfer through exact equality, not assumed similarity.','lanes':[],'differences':[],'complete':False}
def science(r):return json.dumps({k:v for k,v in r.items() if k not in ['elapsed_ns','allocations']},sort_keys=True,separators=(',',':'))
for name,q,old in lanes:
    q=EVIDENCE/q;old=EVIDENCE/old;observed=out/(name+'.jsonl')
    with observed.open('x') as log,(out/(name+'.stderr')).open('x') as err:subprocess.run([str(a.engine.resolve()),str(q)],stdout=log,stderr=err,check=True)
    rows=0
    for original,now in zip(old.read_text().splitlines(),observed.read_text().splitlines(),strict=True):
        before,after=json.loads(original),json.loads(now);rows+=1
        if science(before)!=science(after):receipt['differences'].append({'lane':name,'id':before.get('id'),'before':before,'after':after})
    receipt['lanes'].append({'name':name,'rows':rows,'input':str(q),'input_sha256':sha(q),'reference_output':str(old),'reference_sha256':sha(old),'output_sha256':sha(observed)});save(out/'summary.json',receipt);print(name,rows,len(receipt['differences']),flush=True)
receipt['complete']=True;receipt['rows']=sum(x['rows'] for x in receipt['lanes']);save(out/'summary.json',receipt)
if receipt['differences']:raise SystemExit(1)
