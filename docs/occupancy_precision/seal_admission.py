"""Record final scientific admission independently of controlled timing."""
import json, re
from pathlib import Path
from campaign import EVIDENCE, DOCS, ROOT, sha, save

paths=['final_source_transition.json','final_exact_v2/summary.json',
       'final_science_review_v2/review.json','final_certified_exact/summary.json',
       'final_predicates/summary.json','final_engineering_v2/summary.json',
       'final_arm_v2/summary.json','final_cli_preflight_v2/summary.json',
       'competitive_gate_v2/summary.json','competitive_radius_gate_v2/summary.json',
       'competitive_preflight_v2/summary.json','competitive_warm_preflight_v2/summary.json',
       'competitive_radius_preflight_v2/summary.json','competitive_memory_v2/summary.json',
       'accuracy_preflight_v2/summary.json','accuracy_memory_v2/summary.json',
       'final_valgrind_v2/journal.json','final_legacy_continuous/summary.json']
values={p:json.loads((EVIDENCE/p).read_text()) for p in paths}
assert values['final_exact_v2/summary.json']['exact']
assert values['final_source_transition.json']['finite_grid_source_unchanged']
assert values['final_science_review_v2/review.json']['legacy_reference_preservation_pass']
assert values['final_certified_exact/summary.json']['rows']==695 and not values['final_certified_exact/summary.json']['differences']
assert values['final_engineering_v2/summary.json']['all_commands_succeeded']
assert not values['final_arm_v2/summary.json']['failed_runs']
for p in paths:
    if 'gate_v2/' in p:assert values[p]['passed']
    if 'complete' in values[p]:assert values[p]['complete'],p
msrv=EVIDENCE/'final_engineering_v2/msrv.log'
count=sum(map(int,re.findall(r'test result: ok\. (\d+) passed;',msrv.read_text())))
assert count==372
save(DOCS/'final_admission.json',{
    'scientific_admission_complete':True,'controlled_performance_complete':False,
    'legacy_exact_rows':128021,'certificate_exact_rows':695,
    'rust_1_85_tests':count,'qemu_aarch64_tests':values['final_arm_v2/summary.json']['passed_tests'],
    'competitive_point_decisions':values['competitive_gate_v2/summary.json']['certified_point_decisions'],
    'radius_point_decisions':values['competitive_radius_gate_v2/summary.json']['point_decisions'],
    'scope':'Scientific admission only. Final timing requires separate successful host gates and complete controlled measurements. Known legacy reference disagreements remain preserved. Final-v1 external replay transfers through the explicit source transition and fresh final-v2 exact replay.',
    'bindings':{str((EVIDENCE/p).relative_to(ROOT)):sha(EVIDENCE/p) for p in paths},
    'controller_sha256':sha(Path(__file__)),
})
