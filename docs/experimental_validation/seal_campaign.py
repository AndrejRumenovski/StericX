"""Audit completion gates and seal scientific sources, evidence and deliverables."""
from pathlib import Path
import json,hashlib,subprocess,datetime,re
ROOT=Path(__file__).resolve().parents[2];DOC=Path(__file__).resolve().parent;BASE=ROOT/'.stericx/experimental_validation'
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def check_hash(p,h):assert sha(p)==h,str(p)
def main():
 required=['BASELINE.md','DATASET_REVIEW.md','DESCRIPTOR_REPRODUCTION.md','CONFORMER_STUDY.md','STATISTICAL_PROTOCOL.md','STATISTICAL_RESULTS.md','REGIONAL_MATHEMATICS.md','CERTIFIED_RADIUS_PROFILES.md','NUMERICAL_FAILURE_ANALYSIS.md','PERFORMANCE.md','VALIDATION_REPORT.md','REPRODUCIBILITY.md','FINAL_REPORT.md','WORK_LOG.md']
 for name in required:assert (DOC/name).stat().st_size>100,name
 baseline=read(DOC/'baseline_manifest.json');changed=[]
 for name,r in baseline['records'].items():
  if sha(ROOT/name)!=r['sha256']:changed.append(name)
 allowed={'src/geometry/buried_volume.rs','src/geometry/continuous_volume.rs','src/geometry/continuous_volume/certified.rs','src/geometry/mod.rs'}
 assert set(changed)<=allowed,changed
 for name,r in baseline['retained_evidence'].items():check_hash(ROOT/name,r['sha256'])
 check_hash(BASE/'baseline/production_source.tar',baseline['source_archive_sha256']);check_hash(BASE/'baseline/tracked_index.z',baseline['tracked_index_sha256'])
 engineering=read(BASE/'regression/engineering_v2/summary.json');assert engineering['complete'] and engineering['all_commands_succeeded']
 for row in engineering['runs']:check_hash(BASE/'regression/engineering_v2'/(row['name']+'.log'),row['log_sha256'])
 exact=read(BASE/'regression/exact_comparison/summary.json');assert exact['exact'] and exact['total_rows']==128021 and len(exact['lanes'])==13
 numerical=read(DOC/'geometric_subnormal_validation.json');assert numerical['resolved']==18 and not numerical['invalid'] and not numerical['remaining']
 regional=read(DOC/'regional_validation.json');assert regional['cases']==96 and not regional['failures'] and len(regional['unmet'])==4
 profile=read(DOC/'profile_validation.json');assert not profile['failures'] and not profile['unmet'] and profile['requested_radii']==120
 reproduced=read(BASE/'reproduction_v1/summary.json');assert reproduced['complete'] and reproduced['all_checks_passed']
 for name,h in reproduced['source'].items():check_hash(ROOT/name,h)
 refs=read(BASE/'chemical_references_v2/summary.json');assert len(refs['rows'])==135 and not refs['failures'] and len(refs['statistical_crosschecks'])==16
 for lane in ['references_inferred','references_topology']:
  r=read(BASE/'regression'/lane/'result.json');assert r['commands_completed'] and all(x['returncode']==0 for x in r['runs'])
 volume=read(BASE/'regression/volume_reference/complete.json');assert volume['complete'] and volume['complete_source_stream_validated']
 allbins=read(BASE/'regression/all_bins/complete.json');assert (allbins['rows'],allbins['frames'],allbins['bins'])==(31721,95163,1141956)
 assert read(BASE/'regression/all_bins/metrics.json')==read(ROOT/'.stericx/occupancy_precision/final_science/all_bins/metrics.json')
 performance=read(BASE/'performance_v1/summary.json');assert performance['complete'] and performance['target_results']['expanded_geometric']['met']==18
 vg=read(BASE/'valgrind_v2/summary.json');assert vg['complete'] and all(r['returncode']==0 and not r['preload_failed'] for r in vg['rows']);assert read(BASE/'valgrind_v2/output_comparison.json')['exact_normalized_certificates']
 assert read(BASE/'profile_allocations_v1/summary.json')['exact_certificates']
 for name,meta in read(DOC/'data/manifest.json').items():check_hash(DOC/'data'/name,meta['sha256'])
 results=read(DOC/'results/statistical_results.json');primary=results['results']['primary_as_distributed'];assert primary['B_single']['metrics']['n']==11 and not primary['B_single']['failures'] and not primary['C_ensemble_min']['failures'];assert primary['paired_uncertainty']['family_bootstrap_95_interval'][1]>0
 assert read(DOC/'results/correction_table_comparison.json')['same_whitespace_normalized']
 assert not read(DOC/'results/conformer_sensitivity.json')['generated_encoded_p_mismatches']
 # Scientific fields of fresh reference metrics must equal admitted historical fields.
 comparison=[]
 for lane in ['references_inferred','references_topology']:
  for name in ['kraken/analysis/coverage.json','kraken/interpretation/same_center_analytic_metrics.json','kraken/interpretation/full_primary_conformers_metrics.json','kraken/interpretation/matched_default_percent.json']:
   current=BASE/'regression'/lane/name;previous=ROOT/'.stericx/occupancy_precision/final_science'/lane/name;a,b=read(current),read(previous)
   if name.endswith('matched_default_percent.json'):a.pop('source_manifest_sha256');b.pop('source_manifest_sha256')
   assert a==b,(lane,name);comparison.append(dict(lane=lane,path=name,scientific_fields_exact=True,current_sha256=sha(current),previous_sha256=sha(previous)))
 files=set(p for p in DOC.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name not in ['campaign_manifest.json','campaign_manifest.sha256'])
 files.update((ROOT/'src').rglob('*.rs'));files.update(ROOT/'tests'/name for name in ['certified_geometric.rs','certified_regions.rs','certified_radius_profiles.rs','fixtures/geometric_sum_oracle.json'])
 for directory in ['regression/engineering_v2','regression/exact_comparison','regression/references_inferred','regression/references_topology','regression/volume_reference','regression/all_bins','reproduction_v1','valgrind_v2','profile_allocations_v1','geometric_subnormal_v1','geometric_graph_v1','regional_independent_v1','profiles_independent_v1','chemical_references_v2','performance_v1','sources']:
  d=BASE/directory
  # Hash provenance/completion summaries and direct evidence, avoid copied reference trees.
  for p in d.iterdir():
   if p.is_file():files.add(p)
 files.add(BASE/'regression/upstream_rdkit_api_v2.log')
 manifest=dict(complete=True,created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),baseline_commit=baseline['production_commit'],accepted_parent_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),baseline_existing_files_changed=changed,baseline_retained_evidence_unchanged=len(baseline['retained_evidence']),historical_outputs_exact=exact['total_rows'],scientific_reference_replay=comparison,claims=dict(experimental_hypothesis='insufficient evidence; favorable minimum-descriptor point estimate only',historical_targets_resolved=18,regional_unmet=regional['unmet'],hardware_counters='unavailable; no controlled wall-time claims'),records={str(p.relative_to(ROOT)):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(files)})
 (DOC/'campaign_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');(DOC/'campaign_manifest.sha256').write_text(sha(DOC/'campaign_manifest.json')+'\n');print(json.dumps(dict(complete=True,files=len(files),baseline_changed=changed,historical_rows=128021,original_targets_resolved=18),indent=2))
if __name__=='__main__':main()
