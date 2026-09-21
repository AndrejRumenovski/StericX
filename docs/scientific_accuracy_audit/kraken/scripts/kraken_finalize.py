#!/usr/bin/env python3
"""Freeze final Kraken audit derivatives and verify preserved raw-reference hashes."""
import datetime,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(1024*1024):h.update(b)
 return h.hexdigest()
for phase,key in [('morfeus_outliers','reference_raw_sha256'),('full_analytic_reference','raw_sha256'),('delta_outliers','raw_sha256')]:
 m=json.loads((ROOT/phase/'manifest_frozen.json').read_text());assert sha(ROOT/phase/'reference_raw.jsonl')==m[key]
dossiers=json.loads((ROOT/'delta_interpretation/all_top20_outlier_dossiers.json').read_text());delta=[v for k,rows in dossiers.items() if k.endswith('_delta') for v in rows];assert len(dossiers)==56 and sum(map(len,dossiers.values()))==1120 and len(delta)==280;assert all(len(v['input_ids'])==2 and v['independent_complete_ensemble'] for v in delta)
paths=[ROOT/'REPORT.md',ROOT/'REPLAY_NOTES.md',ROOT.parent/'claims_kraken.md',ROOT/'report_manifest.json']
for folder in ['scripts','interpretation','focused_diagnosis','energy_source_search','delta_interpretation']:
 paths.extend(p for p in (ROOT/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts)
paths.extend(ROOT/folder/name for folder,name in [('analysis','manifest.json'),('sut','manifest_frozen.json'),('prepared_all','manifest.json'),('morfeus_outliers','manifest_frozen.json'),('full_analytic_reference','manifest_frozen.json'),('delta_outliers','manifest_frozen.json'),('primary_si','manifest.json')])
result={'created_utc':datetime.datetime.now(datetime.UTC).isoformat(),'scope':'Final report,15claim records,reproduction scripts,focused analyses,and complete1120outlier dossier withboth280delta extrema. Primary corpus and raw SUT/reference hashes remain in preceding phase manifests.','status':'Scientific audit completed to available evidence; historical source energies/selection records remain explicitly UNCERTAIN. No production code changed.','verification':{'raw_reference_hashes_valid':True,'metric_groups':56,'all_outlier_entries':1120,'delta_entries':280,'all_delta_entries_have_two_inputs_and_complete_reference':True},'files':{str(p.relative_to(ROOT.parent)):sha(p) for p in sorted(set(paths))}}
(ROOT/'manifest_final.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print('Frozen final Kraken manifest',len(result['files']),'artifacts')
