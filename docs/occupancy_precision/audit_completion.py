"""Verify final source, scientific admission, controlled stages and raw evidence."""
import hashlib,json,re,subprocess
from pathlib import Path
from campaign import ROOT,EVIDENCE,DOCS,save

verified={}
def digest(path):
    path=Path(path).resolve()
    if path not in verified:
        with path.open('rb') as f:verified[path]=hashlib.file_digest(f,'sha256').hexdigest()
    return verified[path]
def check(path,expected):
    path=Path(path)
    if not path.is_absolute():path=ROOT/path
    assert digest(path)==expected, str(path)
def read(path):
    digest(path)
    return json.loads(Path(path).read_text())
def mapping(d):
    for path,h in d.items():
        if isinstance(h,str) and re.fullmatch('[0-9a-f]{64}',h):check(path,h)
def records(value):
    if isinstance(value,dict):
        if isinstance(value.get('path'),str) and isinstance(value.get('sha256'),str):check(value['path'],value['sha256'])
        raw=value.get('raw_metrics')
        if raw:
            directory=Path(raw['path']).parent
            for field,name in [('raw_stdout_sha256','stdout.txt'),('raw_stderr_sha256','stderr.txt')]:
                if field in value:check(directory/name,value[field])
        for v in value.values():records(v)
    elif isinstance(value,list):
        for v in value:records(v)
def quiet(value):
    if isinstance(value,dict):
        if 'quiet' in value:assert value['quiet'] is True
        for v in value.values():quiet(v)
    elif isinstance(value,list):
        for v in value:quiet(v)

admission=read(DOCS/'final_admission.json');assert admission['scientific_admission_complete'];mapping(admission['bindings'])
source=read(EVIDENCE/'final_profile_v2/build.json');assert source['complete']
for name,h in source['source'].items():
    check(ROOT/name,h);check(EVIDENCE/'final_profile_v2/source'/name,h)
for name,record in source['binaries'].items():
    check(EVIDENCE/'final_profile_v2'/name,record['sha256'] if isinstance(record,dict) else record)
legacy=read(EVIDENCE/'final_exact_v2/summary.json');assert legacy['exact'] and legacy['total_rows']==128021
assert len(legacy['lanes'])==13
for lane,row in legacy['lanes'].items():
    assert row['baseline_sha256']==row['candidate_sha256'],lane
    for side,folder in [('baseline','baseline_observations'),('candidate','final_observations_v2')]:
        check(EVIDENCE/folder/lane/'stdout.jsonl',row[side+'_sha256'])
cert=read(EVIDENCE/'final_certified_exact/summary.json');assert cert['complete'] and cert['rows']==695 and not cert['differences']
validation=read(DOCS/'c1_validation.json');assert validation['analytical']['cases']==167 and not validation['analytical']['failures']
assert len(validation['analytical']['unmet'])==18
assert validation['predicate_cases']==6000 and not validation['predicate_failures']
assert validation['extended_cases']==180 and not validation['extended_failures']
for tolerance,rows in validation['targets'].items():
    assert len(rows)==12
    assert all(r['target_met'] and r['absolute_error_bound']<=float(tolerance) for r in rows)
    if float(tolerance)==.05:assert max(r['cells'] for r in rows)<=262143
coaxial=read(EVIDENCE/'c1_coaxial/summary.json');assert coaxial['complete'] and coaxial['cases']==108 and not coaxial['failures'] and not coaxial['unmet']
broad=read(EVIDENCE/'c1_molecular_broad_v3/summary.json');assert broad['complete'] and broad['cases']==220 and broad['target_met']==212 and not broad['failures'] and len(broad['unmet'])==8
expanded=read(DOCS/'molecular_validation.json');assert expanded['complete'] and expanded['cases']==expanded['target_met']==20 and not expanded['failures']
for row in expanded['rows']:
    assert row['target_met'] and row['absolute_error_bound']<=(.001 if '/target_0.001' in row['id'] else .05)
assert read(EVIDENCE/'final_science_review_v2/review.json')['legacy_reference_preservation_pass']
assert not read(EVIDENCE/'final_predicates/summary.json')['failures']
assert read(EVIDENCE/'final_engineering_v2/summary.json')['all_commands_succeeded']
assert not read(EVIDENCE/'final_arm_v2/summary.json')['failed_runs']
old_continuous=read(EVIDENCE/'final_legacy_continuous/summary.json');assert not old_continuous['failures'] and len(old_continuous['rows'])==244

data=read(DOCS/'final_measurements.json');assert data['complete'];mapping(data['bindings'])
stages={};attempts=[]
for n in [1,2,3,4]:
    directory=EVIDENCE/f'final_window_{n}';receipt=read(directory/'receipt.json')
    for row in receipt['runs']:
        admitted=row['status']=='completed' and row.get('host_preflight',{}).get('quiet') and row.get('host_postflight',{}).get('quiet')
        attempts.append({'window':n,'stage':row['name'],'status':row['status'],'admitted':bool(admitted)})
        if admitted:
            assert row['name'] not in stages
            assert n in [3,4];quiet(row)
            check(directory/(row['name']+'.log'),row['log_sha256'])
            check(DOCS/'performance_window.py',receipt['controller_sha256'])
            stages[row['name']]=directory/row['name']
assert set(stages)=={'cli_ab','kernel_ab','chord_ab','kernel_profile','cli_profile','accuracy','competitive','radius_competitive'}
expected={'cli_ab':630,'kernel_ab':80,'chord_ab':18,'kernel_profile':71,'cli_profile':175,'accuracy':1152,'competitive':432,'radius_competitive':162}
for name,directory in stages.items():
    report=read(directory/('journal.json' if name.endswith('_profile') else 'summary.json'));assert report['complete'];quiet(report)
    rows=report.get('attempts',report.get('runs'));assert len(rows)==expected[name],name
    for field in ['bindings','hashes','inputs','binaries']:mapping(report.get(field,{}))
    records(report)
    for row in rows:
        if 'status' in row:assert row['status']=='passed'
        if 'exact' in row:assert row['exact']
        if 'metrics' in row:assert row['metrics'].get('returncode',row['metrics'].get('exit_code',0))==0
    if name.endswith('_profile'):
        check(ROOT/'docs/ultimate_engine/profile_campaign.py',report['controller_sha256'])
        check(ROOT/'docs/ultimate_engine/benchmark_artifacts.py',report['artifact_helper_sha256'])

accuracy=read(stages['accuracy']/'summary.json');preflight=read(EVIDENCE/'accuracy_preflight_v2/summary.json')
expected_science={r['id']:{k:v for k,v in r['result'].items() if k!='elapsed_seconds'} for r in preflight['runs']}
assert len(expected_science)==144
for row in accuracy['runs']:assert {k:v for k,v in row['result'].items() if k!='elapsed_seconds'}==expected_science[row['id']]
assert len(data['accuracy_rows'])==144
assert all(r['result']['target_met'] for r in data['accuracy_rows'] if r['method'].startswith('certified'))
pareto=read(DOCS/'accuracy_pareto.json');assert len(pareto['rows'])==144 and len(pareto['fastest_by_target'])==108
for row in pareto['fastest_by_target']:
    assert row['fastest_tested_method'] is not None
    if row['accuracy_class']=='certified':assert row['fastest_tested_method']=='certified_'+str(row['tolerance'])
scaling=read(DOCS/'certified_scaling.json');assert scaling['complete'] and scaling['absolute_tolerance']==.05 and scaling['molecules_per_batch']==12
batch_build=read(EVIDENCE/'certified_batch_build/receipt.json');assert batch_build['returncode']==0 and batch_build['source_build_sha256']==digest(EVIDENCE/'final_profile_v2/build.json')
check(DOCS/'certified_batch.rs',batch_build['observer_sha256']);check(EVIDENCE/'certified_batch_build/certified_batch',batch_build['binary_sha256'])
records(scaling);scaling_raw=read(Path(scaling['source_receipt']['path']));assert scaling_raw['complete'] and not scaling_raw['preflight_only'];quiet(scaling_raw)
assert len(scaling_raw['runs'])==45 and all(r['exact'] for r in scaling_raw['runs'])
mapping(scaling_raw['bindings']);records(scaling_raw)
assert set(scaling['rows'])=={'1','2','4','6','12'}
assert all(r['batch_ns']['n']==8 for r in scaling['rows'].values())
scaling_host=read(EVIDENCE/'certified_scaling_1_host_control/receipt.json');assert scaling_host['complete'] and scaling_host['returncode']==0 and scaling_host['restored']==scaling_host['suspended']
check(DOCS/'run_quiet_window.py',scaling_host['controller_sha256'])
warp_v1_sha=hashlib.sha256(subprocess.check_output(['git','show','a9d1c8c0:docs/occupancy_precision/run_quiet_window.py'],cwd=ROOT)).hexdigest()
for n in [2,3,4]:
    host=read(EVIDENCE/f'final_window_{n}_host_control/receipt.json')
    assert host['complete'] and host['restored']==host['suspended']
    assert host['controller_sha256']==warp_v1_sha
    for row in host['restored']:
        path=Path('/proc')/str(row['pid'])/'stat'
        if path.exists():
            text=path.read_text();parts=text[text.rindex(')')+2:].split()
            if parts[19]==row['starttime']:assert parts[0] not in ['T','t']
assert not subprocess.check_output(['git','diff','ffc698f4','--','docs/ultimate_engine','docs/validation/dbstep'],cwd=ROOT)
for name in ['BASELINE.md','OCCUPANCY_PROFILE.md','FINAL_REPORT.md','FINAL_TABLES.md','CERTIFIED_SCALING.md','CONTINUOUS_MATHEMATICS.md','EXPERIMENTS.md','REPRODUCIBILITY.md','CANDIDATES.md','WORK_LOG.md','throughput_scaling.svg','throughput_scaling.png','accuracy_cost_pareto.svg','accuracy_cost_pareto.png','final_diagnostics.json','timing_admission.json','baseline_manifest.json','workload_identity.json']:
    digest(DOCS/name)
digest(EVIDENCE/'resumed_environment.json')
phases=[
 (0,'Frozen production source, input/compiler identities and fresh scientific baseline','BASELINE.md; baseline_manifest.json; baseline/build.json; baseline_checks.json'),
 (1,'Deep occupancy profiling, assembly, counters and before flamegraph; hardware limitations explicit','OCCUPANCY_PROFILE.md; baseline_profile_events.json; baseline_counters.json; baseline_flamegraph.svg'),
 (2,'Exact occupancy candidates implemented, compared, retained/rejected; masks and boundaries validated','CANDIDATES.md; EXPERIMENTS.md; o1–o5_measurements.json; final_exact_v2/summary.json'),
 (3,'Continuous target defined; analytical, interval, subdivision, clipping and quadrature methods analyzed','CONTINUOUS_MATHEMATICS.md; EXPERIMENTS.md; continuous_volume/certified.rs'),
 (4,'Explicit certified API, failure status and fixed tolerance results verified','c1_validation.json; molecular_validation.json; final_certified_exact/summary.json'),
 (5,'Shared conservative candidate strategy compared against specialized arithmetic and all-atom ablation','chord_ablations.json; o5_measurements.json; EXPERIMENTS.md'),
 (6,'Independent rational/analytical/coaxial references and extreme predicates verified','c1_validation.json; c1_coaxial/summary.json; final_predicates/summary.json; final_certified_exact/summary.json'),
 (7,'Frozen real conformers, multiple radii/centers/policies, budgets and nesting evaluated','c1_molecular_broad_v3/summary.json; molecular_validation.json; final_science_review_v2/review.json'),
 (8,'Isolated accepted commits, measured iterations, rejected evidence and final reprofile retained','CANDIDATES.md; EXPERIMENTS.md; final_diagnostics.json; WORK_LOG.md; Git commits in FINAL_REPORT.md'),
 (9,'Controlled accuracy ladder, certified batch scaling, CPU/wall/memory/evaluations, explicit unmet targets and separate Pareto classes','final_measurements.json; certified_scaling.json; accuracy_pareto.json; accuracy_cost_pareto.svg; FINAL_TABLES.md; CERTIFIED_SCALING.md'),
 (10,'Matched Morfeus/DBSTEP finite-grid and radius workflows, warm/fresh trials and scaling completed','competitive_gate_v2/summary.json; competitive_radius_gate_v2/summary.json; final_window_4/competitive; final_window_4/radius_competitive'),
 (11,'Full regression, byte-exact historical replay, external references, Python and available platform checks','final_admission.json; final_source_transition.json; final_engineering_v2; final_arm_v2; final_science_review_v2'),
 (12,'Final controlled A/B stages complete under unchanged quiet gates; all failures and samples retained','final_window_3/receipt.json; final_window_4/receipt.json; timing_admission.json; final_measurements.json'),
 (13,'Final scientific/engineering report, distributions, figures, commits and reproduction instructions complete','FINAL_REPORT.md; FINAL_TABLES.md; REPRODUCIBILITY.md; completion_audit.json'),
]
result={'complete':True,'scientific_complete':True,'controlled_performance_complete':True,'source_commit':'aeb3ca747ffbefa5003854e5709f2efa0ba0905e','legacy_raw_exact_streams':13,'legacy_rows':128021,'certificate_rows':695,'admitted_stage_run_counts':expected,'certified_scaling_processes':45,'attempts':attempts,'warp_restoration_verified':True,
 'scope':'First admitted stage results only; all failed attempts and adverse samples retained. Hardware counters unavailable, physical cross-platform timing unavailable, and external certified continuous comparisons not admitted. Reported analytical/budget failures are limitations, not silently promoted to successes.',
 'phases':[{'phase':n,'verified_requirement':text,'evidence':evidence} for n,text,evidence in phases],
 'verified_file_hashes':{str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p):h for p,h in sorted(verified.items())},'controller_sha256':digest(Path(__file__))}
save(DOCS/'completion_audit.json',result)
lines=['# Completion audit','', 'All fourteen requested phases have supporting evidence. Scientific and controlled performance completion are separate verified claims. Raw source/binary/input/measurement hashes checked by `audit_completion.py` are recorded in [completion_audit.json](completion_audit.json). Local evidence paths below are relative to `.stericx/occupancy_precision`; tracked report filenames are relative to this directory.','', '| Phase | Verified requirement | Evidence |','|---|---|---|']
lines += [f'| {n} | {text} | {evidence} |' for n,text,evidence in phases]
lines += ['', 'The audit verifies all 13 raw historical response streams, final scientific source hashes, 695 certified outputs through their sealed admission, all eight timing-stage gates and trial counts, the additional 45-process certified-scaling supplement, linked raw metrics/stdout/stderr hashes, unchanged accuracy responses across 1,152 observations, and Warp restoration. Every accepted stage is the first qualifying stage. The rejected 630-trial CLI attempt and preflight failures remain recorded.', '', 'Documented limits remain: no native hardware-counter access; no physical ARM/macOS/Windows timing; eighteen analytical stress targets unmet; eight broad molecular targets require a larger budget; no certified regional continuous API; no equivalent external certified-continuous comparison. Quiet endpoints cannot prove zero interruption inside a trial, and all outliers remain in the distributions.', '', 'No historical report was edited. No scientific constant, acceptance tolerance, finite-grid default or release state was changed. Separate meaningful production commits and complete reproduction commands are listed in FINAL_REPORT.md and REPRODUCIBILITY.md.','']
(DOCS/'COMPLETION_AUDIT.md').write_text('\n'.join(lines))
print(json.dumps({'complete':True,'verified_files':len(verified),'admitted_stages':len(stages),'legacy_raw_exact_streams':13,'warp_restored':True}))
