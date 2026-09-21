"""Replay final Kraken interpretation in scratch space using immutable references."""
from pathlib import Path
from datetime import datetime, UTC
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile

A = Path(__file__).resolve().parents[1]
SOURCE = A / 'kraken'
scratch = Path(tempfile.mkdtemp(prefix='stericx-kraken-interpretation-'))
root = scratch / 'kraken'
root.mkdir()
shutil.copytree(SOURCE / 'scripts', root / 'scripts', ignore=shutil.ignore_patterns('__pycache__'))
for name in ['analysis', 'primary', 'prepared_all', 'sut', 'morfeus_outliers', 'full_analytic_reference', 'delta_outliers']:
    (root / name).symlink_to(SOURCE / name, target_is_directory=True)
records = []
for name in ['kraken_interpret.py', 'kraken_matched_percent.py', 'kraken_finish_delta.py', 'kraken_369_dossier.py']:
    result = subprocess.run([sys.executable, str(root / 'scripts' / name)], text=True, capture_output=True, timeout=600)
    records.append({'script': name, 'returncode': result.returncode, 'stdout': result.stdout, 'stderr': result.stderr})
    if result.returncode:
        break


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def normalized(x):
    if isinstance(x, dict):
        return {k: normalized(v) for k, v in x.items() if k != 'created_utc'}
    if isinstance(x, list):
        return [normalized(v) for v in x]
    return x


comparisons = []
for folder in ['interpretation', 'delta_interpretation']:
    for p in sorted((root / folder).glob('*')):
        if not p.is_file() or p.name in {'manifest.json', 'analytic_raw_frozen.json'}:
            continue
        relative = p.relative_to(root)
        expected = SOURCE / relative
        exact = expected.is_file() and sha(p) == sha(expected)
        same = exact
        if not exact and expected.is_file() and p.suffix == '.json':
            same = normalized(json.loads(p.read_text())) == normalized(json.loads(expected.read_text()))
        comparisons.append({'path': relative.as_posix(), 'exact_bytes': exact, 'equal_except_created_utc': same, 'expected_sha256': sha(expected) if expected.is_file() else None, 'replay_sha256': sha(p)})
passed = len(records) == 4 and all(r['returncode'] == 0 for r in records) and len(comparisons) >= 10 and all(r['equal_except_created_utc'] for r in comparisons)
record = {'created_utc': datetime.now(UTC).isoformat(), 'purpose': __doc__, 'scope': 'Reinterprets frozen independent numerical references; does not rerun Morfeus kernels or claim a new scientific reference.', 'scratch': str(scratch), 'commands': records, 'comparisons': comparisons, 'passed': passed}
target = A / 'reproducibility' / ('kraken_interpretation_' + datetime.now(UTC).strftime('%Y%m%dT%H%M%SZ') + '.json')
target.write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps({'passed': passed, 'files_compared': len(comparisons), 'record': str(target)}, indent=2))
if passed:
    shutil.rmtree(scratch)
else:
    raise SystemExit('Mismatch preserved in ' + str(scratch))
