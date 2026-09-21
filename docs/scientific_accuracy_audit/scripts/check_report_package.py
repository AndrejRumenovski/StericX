"""Check report coverage and artifact references; this is not scientific validation."""
from pathlib import Path
from datetime import datetime, UTC
import ast
import json
import re
from urllib.parse import unquote

A = Path(__file__).resolve().parents[1]
failures = []
checks = []


def check(name, ok, detail=None):
    row = {'check': name, 'passed': bool(ok)}
    if detail is not None:
        row['detail'] = detail
    checks.append(row)
    if not ok:
        failures.append(row)


report = (A / 'SCIENTIFIC_ACCURACY_AUDIT.md').read_text()
check('Opening component table', report.startswith('| Component |'))
sections = re.findall(r'^## (\d+)\. ', report, flags=re.M)
check('All fourteen requested report sections', sections == list(map(str, range(1, 15))))
claims = json.loads((A / 'claims.json').read_text())['claims']
check('Unique claim IDs', len({r['id'] for r in claims}) == len(claims), len(claims))
required = ['claim', 'source', 'convention', 'implementation', 'method', 'result', 'confidence', 'limitations']
check('Every claim has all eight evidence fields', all(all(str(r.get(k, '')).strip() for k in required) for r in claims))
allowed = {'VERIFIED', 'VERIFIED WITH NUMERICAL LIMIT', 'SUPPORTED', 'UNCERTAIN', 'INCORRECT', 'TERMINOLOGY ISSUE', 'OUT OF SCOPE'}
check('Only requested claim classifications', all(r['status'] in allowed for r in claims))

metrics = json.loads((A / 'kraken/analysis/metrics.json').read_text())
fields = {'N', 'MAE', 'RMSE', 'maximum_absolute_error', 'median_absolute_error', 'R2_one_to_one', 'slope', 'intercept'}
check('All 56 Kraken metric records contain requested statistics', len(metrics) == 56 and all(fields <= set(v) for v in metrics.values()))
check('All 56 Kraken residual plots exist', all((A / 'visuals/kraken_all_metrics' / (k + '.png')).is_file() for k in metrics))
dossiers = json.loads((A / 'kraken/delta_interpretation/all_top20_outlier_dossiers.json').read_text())
check('Twenty outliers for every Kraken metric', set(dossiers) == set(metrics) and all(len(v) == 20 for v in dossiers.values()))
delta = [r for k, rows in dossiers.items() if k.endswith('_delta') for r in rows]
check('All 280 conformer-range dossiers map both extrema', len(delta) == 280 and all(len(r.get('native_extremizers', [])) == 2 and all(v.get('input_sha256') for v in r['native_extremizers']) for r in delta))
bv_delta = [r for r in delta if not r['descriptor'].startswith(('sterimol', 'pyr'))]
check('All 180 buried-volume range dossiers include complete independent ensembles', len(bv_delta) == 180 and all(len(r.get('independent_complete_ensemble', [])) == 4 and all('error' not in v and v['N_conformers'] == r['conformers'] for v in r['independent_complete_ensemble']) for r in bv_delta))
other = [r for k, rows in dossiers.items() if not k.endswith('_delta') for r in rows]
check('All 840 remaining outliers identify their native selected input', len(other) == 840 and all(r.get('input_id') and r.get('input_path') for r in other))

lanes = {}
for p in sorted((A / 'reproducibility').glob('replay_*.json')):
    for row in json.loads(p.read_text())['lanes']:
        lanes[row['lane']] = row
expected = {'geometry', 'geometry_focused', 'geometry_alignment', 'geometry_alignment_rotation', 'geometry_numerical', 'kraken'}
check('All six frozen native streams replayed exactly', expected <= set(lanes) and all(lanes[k]['exact_bytes'] and lanes[k]['returncode'] == 0 for k in expected), sorted(lanes))
model = json.loads((A / 'models/results/completion_checks.json').read_text())
check('Complete primary cross-coupling extraction and independent fitting', model['crosscoupling_rows'] == 746 and model['crosscoupling_fits'] == 72 and model['reference_sklearn_total_prediction_disagreement'] == 0 and not model['missing_crosscoupling_descriptor_ids'])

docs = ['SCIENTIFIC_ACCURACY_AUDIT.md', 'CLAIMS.md', 'COMPLETION_CHECKLIST.md', 'REPRODUCE.md', 'geometry/REPORT.md', 'geometry/METHODS.md', 'models/REPORT.md', 'models/DISCREPANCIES.md', 'kinetics/KINETICS_CONFORMERS.md', 'kraken/REPORT.md', 'visuals/kraken_all_metrics/INDEX.md']
links = 0
missing = []
for name in docs:
    p = A / name
    for target in re.findall(r'\[[^\]\n]*\]\(([^)\n]+)\)', p.read_text()):
        target = target.strip('<>').split('#', 1)[0]
        if not target or re.match(r'^[a-zA-Z]+:', target):
            continue
        links += 1
        if not (p.parent / unquote(target)).exists():
            missing.append({'document': name, 'target': target})
check('Authoritative local Markdown links resolve', not missing, {'links_checked': links, 'missing': missing})
scripts = list((A / 'scripts').glob('*.py')) + list((A / 'kraken/scripts').glob('*.py'))
invalid = []
for p in scripts:
    try:
        ast.parse(p.read_text(), filename=str(p))
    except SyntaxError as exc:
        invalid.append({'script': str(p.relative_to(A)), 'error': str(exc)})
check('Audit Python scripts parse', not invalid, {'scripts': len(scripts), 'invalid': invalid})
result = {'checked_utc': datetime.now(UTC).isoformat(), 'purpose': __doc__, 'checks': checks, 'passed': not failures}
(A / 'reproducibility/report_package_checks.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
if failures:
    raise SystemExit(1)
