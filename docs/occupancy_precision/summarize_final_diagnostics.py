"""Summarize simulator and allocator diagnostics without inferring native speed."""
import collections, json, re, statistics, sys
from pathlib import Path
from campaign import ROOT, EVIDENCE, DOCS, sha, save
sys.path.insert(0, str(ROOT / 'docs/ultimate_engine'))
from callgrind_flamegraph import parse

def profile(path):
    costs, total = parse(path)
    exclusive = collections.Counter()
    for context, n in costs.items():
        exclusive[context.split("'")[0]] += n
    return {'total_instructions': total, 'exclusive_top': exclusive.most_common(30)}

occpath = EVIDENCE / 'final_call_context_v2/callgrind.out'
contpath = EVIDENCE / 'final_continuous_profile/callgrind.out'
cachepath = EVIDENCE / 'final_valgrind_v2/cachegrind/cachegrind.out'
lines = cachepath.read_text().splitlines()
events = next(s.split()[1:] for s in lines if s.startswith('events:'))
counts = next(map(int, s.split()[1:]) for s in lines if s.startswith('summary:'))
assembly = EVIDENCE / 'final_call_context_v2/assembly.txt'
match = re.search(r'^([0-9a-f]+) <steric_x::geometry::buried_volume::occupied_prepared_avx2>:\n(.*?)(?=\n[0-9a-f]+ <)', assembly.read_text(), re.M | re.S)
assert match
kernel = match.group(0)
(EVIDENCE / 'final_call_context_v2/occupancy.asm').write_text(kernel)
assert not re.search(r'\bv(?:f|fn)m(?:add|sub)', kernel)
oldpath, newpath = [EVIDENCE / p / 'summary.json' for p in ['accuracy_memory', 'accuracy_memory_v2']]
old, new = [json.loads(p.read_text()) for p in [oldpath, newpath]]
assert old['complete'] and new['complete'] and len(old['runs']) == len(new['runs']) == 144
before = {r['id']: r for r in old['runs']}
rows = []
for row in new['runs']:
    prior = before[row['id']]
    assert {k:v for k,v in row['result'].items() if k != 'elapsed_seconds'} == {k:v for k,v in prior['result'].items() if k != 'elapsed_seconds'}
    if row['configuration'].startswith('certified'):
        a,b = prior['allocations'][0], row['allocations'][0]
        rows.append({'id':row['id'], 'configuration':row['configuration'],
                     'before':a, 'after':b, 'call_reduction_fraction':1-b['calls']/a['calls']})
groups={}
for configuration in ['certified_0.1','certified_0.05','certified_0.01']:
    molecular=[r for r in rows if r['configuration']==configuration and '/reference/' in r['id']]
    groups[configuration]={'cases':len(molecular),
        'median_call_reduction_fraction':statistics.median(r['call_reduction_fraction'] for r in molecular),
        'before_calls_range':[min(r['before']['calls'] for r in molecular),max(r['before']['calls'] for r in molecular)],
        'after_calls_range':[min(r['after']['calls'] for r in molecular),max(r['after']['calls'] for r in molecular)],
        'before_peak_live_max':max(r['before']['process_peak_live_bytes'] for r in molecular),
        'after_peak_live_max':max(r['after']['process_peak_live_bytes'] for r in molecular)}
result={'scope':'Modeled instructions/cache events and exact allocator diagnostics. Memory passes overlapped validation; their timings are not comparative performance evidence. No hardware cycle, IPC or bandwidth inference.',
        'occupancy':profile(occpath), 'continuous':profile(contpath),
        'cachegrind_simulated':dict(zip(events,counts)), 'occupancy_fma_instructions':0,
        'allocation_groups':groups,'allocation_rows':rows,
        'bindings':{str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),occpath,contpath,cachepath,assembly,oldpath,newpath]}}
save(DOCS / 'final_diagnostics.json',result)
print(json.dumps({k:result[k] for k in ['allocation_groups','occupancy_fma_instructions']},indent=2))
