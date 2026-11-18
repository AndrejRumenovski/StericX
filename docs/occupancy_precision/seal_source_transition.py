"""Bind expensive independent legacy replay to the final source by exact evidence."""
import json
from pathlib import Path
from campaign import ROOT, EVIDENCE, sha, save

old = json.loads((EVIDENCE / 'final_profile/build.json').read_text())
new = json.loads((EVIDENCE / 'final_profile_v2/build.json').read_text())
changes = {p: [old['source'].get(p), new['source'].get(p)]
           for p in old['source'].keys() | new['source'].keys()
           if old['source'].get(p) != new['source'].get(p)}
production = {p: v for p, v in changes.items() if p.startswith(('src/', 'python/rust/'))}
assert set(production) == {'src/geometry/continuous_volume/certified.rs'}, production
for p, digest in new['source'].items():
    assert sha(EVIDENCE / 'final_profile_v2/source' / p) == digest, p
    assert sha(ROOT / p) == digest, p
exact = EVIDENCE / 'final_exact_v2/summary.json'
comparison = json.loads(exact.read_text())
assert comparison['exact'] and comparison['total_rows'] == 128021
# The independently replayed certified binary was built from an identical
# scientific source snapshot. Verify files, rather than trusting its label.
cert_source = EVIDENCE / 'final_certified/source'
for p in new['source']:
    if p.startswith(('src/', 'python/rust/')):
        assert sha(cert_source / p) == new['source'][p], p
paths = ['final_profile/build.json', 'final_profile_v2/build.json',
         'final_scientific_build/manifest.json', 'final_scientific_build_v2/manifest.json',
         'final_exact_comparison/summary.json', 'final_exact_v2/summary.json',
         'final_science_review_v2/review.json', 'final_certified/certified_build.json',
         'final_certified_exact/summary.json', 'final_predicates/summary.json']
save(EVIDENCE / 'final_source_transition.json', {
    'finite_grid_source_unchanged': True,
    'changed_files': changes,
    'production_changes': production,
    'scope': 'Full independent legacy references were freshly evaluated on final-v1. '
             'Final-v2 changes only the separate certified API, reproduces all 128021 '
             'legacy responses, and has fresh engineering checks. Its identical '
             'certified source reproduces all 695 independently validated certificates. '
             'No claim that every external reference program was run twice.',
    'bindings': {p: sha(EVIDENCE / p) for p in paths},
})
