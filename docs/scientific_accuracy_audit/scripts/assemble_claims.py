"""Assemble individual, eight-field claim records without merging verdicts."""
from pathlib import Path
from collections import Counter
import json,re,runpy

A=Path(__file__).resolve().parents[1]
VALID={'VERIFIED','VERIFIED WITH NUMERICAL LIMIT','SUPPORTED','UNCERTAIN','INCORRECT','TERMINOLOGY ISSUE','OUT OF SCOPE'}
claims=[]
def add(row,component):
    keys={
      'source':['source','primary_or_reference_source'],
      'convention':['expected_convention','exact_convention_expected','convention'],
      'implementation':['implementation','implementation_location'],
      'method':['independent_method','independent_validation_method','method'],
    }
    out={k:row[k] for k in ['id','claim','result','status','confidence','limitations']}
    for key,aliases in keys.items():
        out[key]=next(row[a] for a in aliases if a in row)
    if 'basis' in row:out['basis']=row['basis']
    if 'claim_origin' in row:out['basis']=row['claim_origin']
    out['component']=component
    if out['status'] not in VALID:raise ValueError(out)
    if not all(str(out[k]).strip() for k in ['claim','source','convention','implementation','method','result','confidence','limitations']):raise ValueError(out)
    claims.append(out)

for row in json.loads((A/'geometry/claims.json').read_text()):add(row,'Geometry')
namespace=runpy.run_path(str(A/'scripts/kinetics_claims.py'))
for ident,claim,source,convention,location,method,status,result,confidence,limitations in namespace['claims']:
    add(dict(id=ident,claim=claim,source=source,convention=convention,implementation=location,method=method,status=status,result=result,confidence=confidence,limitations=limitations,basis=namespace['basis'].get(ident,'')),'Conformers and kinetics')
for row in json.loads((A/'models/results/claims.json').read_text()):add(row,'Models, statistics, uncertainty and screening')

text=(A/'claims_kraken.md').read_text()
for part in re.split(r'(?m)^## ',text)[1:]:
    heading,body=part.split('\n',1)
    match=re.match(r'(K\d+)\s*[—–-]\s*(.*)',heading)
    if not match:continue
    fields={k.strip():v.strip() for k,v in re.findall(r'(?m)^(?:-|\d+\.) \*\*([^*]+):\*\*\s*(.*)$',body)}
    def field(*names):
        return next(fields[n] for n in names if n in fields)
    result_key=next(k for k in fields if k.startswith('Result'))
    result=result_key+' — '+fields[result_key]
    status=next((s for s in sorted(VALID,key=len,reverse=True) if re.search(r'\b'+re.escape(s)+r'\b',result)),None)
    if status is None:raise ValueError((match.group(1),result))
    add(dict(id=match.group(1),claim=field('Claim','StericX claim'),source=field('Primary source','Primary/reference source'),convention=field('Exact convention','Exact convention expected','Expected convention'),implementation=field('Implementation location').strip('`'),method=field('Independent validation method','Independent validation'),status=status,result=fields[result_key],confidence=field('Confidence'),limitations=field('Limitations')),'Kraken reproduction')

ids=[r['id'] for r in claims]
if len(ids)!=len(set(ids)):raise ValueError('Duplicate claim IDs')
counts=Counter(r['status'] for r in claims)
out=['# Scientific claim inventory','',f'{len(claims)} claims and explicitly identified audit propositions, each with a separate scientific classification. Frozen system: commit `6393aafe0d983e504baf8abc1e18dc2a0f0d40e7`.','',
'This inventory examines the README, scientific/validation studies, code comments, model documentation and tests as claims to investigate. They are not scientific reference truth. Some entries test a stronger interpretation or a numerical-domain property rather than quote an explicit promise; those distinctions and existing documentation caveats are retained. A high-confidence negative finding is not evidence that every normal workload fails.','',
'Classification uses the requested seven categories. **VERIFIED** requires independent support for the stated scope; **VERIFIED WITH NUMERICAL LIMIT** adds a measured precision/discretization limit; **SUPPORTED** has favorable but incomplete evidence; **UNCERTAIN** has missing or conflicting evidence; **INCORRECT** has a constructive contradiction; **TERMINOLOGY ISSUE** concerns scientific/statistical wording; **OUT OF SCOPE** cannot be established computationally from available evidence.','',
'The counts below count propositions, not independent experiments or a percentage of scientific correctness. Related claims can share evidence without being merged into one PASS.','',
'| Classification | Claims |','|---|---:|']
for status in ['VERIFIED','VERIFIED WITH NUMERICAL LIMIT','SUPPORTED','UNCERTAIN','INCORRECT','TERMINOLOGY ISSUE','OUT OF SCOPE']:out.append(f'| {status} | {counts[status]} |')
out+=['','Evidence reports: [overview](SCIENTIFIC_ACCURACY_AUDIT.md), [geometry](geometry/REPORT.md), [kinetics/conformers](kinetics/KINETICS_CONFORMERS.md), [Kraken](kraken/REPORT.md), [models](models/REPORT.md). Machine-readable records: [claims.json](claims.json).','']
previous=None
for c in claims:
    if c['component']!=previous:out += [f"## {c['component']}",''];previous=c['component']
    out += [f"### {c['id']} — {c['claim']}",'',f"1. **StericX claim / audited proposition:** {c['claim']} "+c.get('basis',''),f"2. **Primary/reference source:** {c['source']}",f"3. **Exact convention expected:** {c['convention']}",f"4. **Implementation location:** `{c['implementation']}`",f"5. **Independent validation method:** {c['method']}",f"6. **Result:** **{c['status']}** — {c['result']}",f"7. **Confidence:** {c['confidence']}",f"8. **Limitations:** {c['limitations']}",'']
(A/'CLAIMS.md').write_text('\n'.join(out))
(A/'claims.json').write_text(json.dumps({'claims':claims,'status_counts':dict(counts)},indent=2,ensure_ascii=False)+'\n')
print(json.dumps({'claims':len(claims),'status_counts':dict(counts)},indent=2))
