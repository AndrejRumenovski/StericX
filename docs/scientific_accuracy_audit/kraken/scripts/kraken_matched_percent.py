#!/usr/bin/env python3
"""Reproduce the original diagnostic sample's percent-volume comparison."""
import math,json
from kraken_analyze import ROOT,metrics,save,sha
refs={v['id']:v for v in map(json.loads,(ROOT/'morfeus_outliers/reference_raw.jsonl').open())};rows=[]
for pair in map(json.loads,(ROOT/'morfeus_outliers/selected_inputs.jsonl').open()):
 s=pair['sut'];r=refs[s['id']].get('matched_sut_001',{});bv=s['buried_volume']
 if 'error' in bv or 'error' in r:continue
 ref=r['orientations'][0]['buried_volume']*100/(4*math.pi*3.5**3/3);val=bv['percent_buried_volume'];rows.append({'id':s['id'],'reference_percent':ref,'sut_percent':val,'difference_percentage_points':val-ref})
m=metrics([v['reference_percent'] for v in rows],[v['sut_percent'] for v in rows]);m['maximum_case']=max(rows,key=lambda r:abs(r['difference_percentage_points']))
save(ROOT/'interpretation/matched_default_percent.json',{'metrics':m,'rows':rows,'source_manifest_sha256':sha(ROOT/'morfeus_outliers/manifest_frozen.json'),'equation':'100×first_orientation_occupied_volume/(4π3.5³/3), versus original SUT percent output.'})
