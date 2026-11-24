"""Compare profiler certificates, omitting named time/allocation diagnostics only."""
from pathlib import Path
import json
from profile_methods import strip
ROOT=Path(__file__).resolve().parents[2];B=ROOT/'.stericx/experimental_validation';p=B/'performance_v1';v=B/'valgrind_v2';checks=[]
for label,plain in [('old_chain',p/'original_frozen_0/stdout'),('new_chain',p/'original_geometric_0/stdout'),('regions',p/'allocations_regions/stdout'),('profile',p/'profile_with_duplicate/stdout')]:
 old=list(map(json.loads,plain.read_text().splitlines()));native=json.loads((v/f'callgrind_{label}.stdout').read_text());ref=next(r for r in old if r['id']==native['id']);assert strip(native)==strip(ref),label;checks.append(label)
for label in ['regions','profile']:
 a=json.loads((v/f'callgrind_{label}.stdout').read_text());c=json.loads((v/f'memcheck_{label}.stdout').read_text());assert strip(a)==strip(c)
(v/'output_comparison.json').write_text(json.dumps(dict(exact_normalized_certificates=True,cases=checks,omitted_fields=['elapsed_seconds','elapsed_ns','allocations','geometry_preparation_seconds','preparation_seconds'],first_check='The initial comparison reported only /profile/preparation_seconds, an unlisted timing field. Explicit normalization added; no numerical fields relaxed.'),indent=2)+'\n')
print('All profiler numerical outputs match exactly.')
