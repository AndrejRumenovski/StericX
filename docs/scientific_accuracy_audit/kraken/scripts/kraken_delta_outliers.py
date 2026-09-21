#!/usr/bin/env python3
"""Complete independent BV ensembles for every top20 delta-outlier ligand."""
import concurrent.futures,datetime,json
from pathlib import Path
from kraken_morfeus_outliers import ROOT,run,save,sha,emit
OUT=ROOT/'delta_outliers'
def main():
 OUT.mkdir(exist_ok=False)
 top=json.loads((ROOT/'analysis/top20_by_descriptor.json').read_text())
 mids={r['molecule_id'] for name,rows in top.items() if name.endswith('_delta') and not name.startswith(('sterimol','pyr')) for r in rows}
 requests={r['id']:r for r in map(json.loads,(ROOT/'prepared_all/requests.jsonl').open()) if int(r['id'].split(':')[1]) in mids}
 existing={r['id']:r for r in map(json.loads,(ROOT/'morfeus_outliers/reference_raw.jsonl').open()) if r['id'] in requests}
 todo=[]
 with (OUT/'selected_inputs.jsonl').open('w') as dst:
  for s in map(json.loads,(ROOT/'sut/stdout.jsonl').open()):
   if s['id'] not in requests:continue
   pair={'request':requests[s['id']],'sut':s};emit(dst,pair)
   if s['id'] not in existing:todo.append(pair)
 save(OUT/'selection_frozen.json',{'created_utc':datetime.datetime.now(datetime.UTC).isoformat(),'ligands':sorted(mids),'N':len(requests),'new_calculations':len(todo),'reused_previously_frozen':len(existing),'reuse_manifest_sha256':sha(ROOT/'morfeus_outliers/manifest_frozen.json'),'selected_input_sha256':sha(OUT/'selected_inputs.jsonl'),'script_sha256':sha(Path(__file__)),'helper_sha256':sha(Path(__file__).with_name('kraken_morfeus_outliers.py')),'source_SUT_sha256':sha(ROOT/'sut/stdout.jsonl'),'selection':'Every available conformer of every ligand among top20 absolute delta residuals for all9reported buried-volume fields. Reuse exact immutable independent record if already computed; calculate all missing records. All four untuned variants as original diagnostic campaign; no filtering or tuning.'})
 print('Frozen',len(mids),'ligands',len(requests),'conformers',len(todo),'new',flush=True)
 with (OUT/'reference_raw.jsonl').open('w') as dst:
  for key in sorted(existing):emit(dst,existing[key])
  with concurrent.futures.ProcessPoolExecutor(max_workers=3) as pool:
   for i,row in enumerate(pool.map(run,todo,chunksize=1),1):
    emit(dst,row)
    if i%25==0:print('Completed',i,'/',len(todo),flush=True)
 save(OUT/'manifest_frozen.json',{'created_utc':datetime.datetime.now(datetime.UTC).isoformat(),'record_count':len(requests),'raw_sha256':sha(OUT/'reference_raw.jsonl'),'selection_sha256':sha(OUT/'selection_frozen.json'),'frozen_before_interpretation':True})
if __name__=='__main__':main()
