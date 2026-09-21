"""Add a dated verification without overwriting the sealed scientific audit.

Run from the repository root with .venv/bin/python3. This checks current
identity/environment and independently recomputes report statistics from the
frozen row table. It does not claim to rerun independent descriptor kernels.
"""
from pathlib import Path
from datetime import datetime, UTC
from collections import defaultdict
import csv
import hashlib
import importlib.metadata
import importlib.util
import json
import math
import statistics
import subprocess
import sys

A = Path(__file__).resolve().parents[1]
REPO = A.parents[1]


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def command(*args):
    return subprocess.check_output(args, cwd=REPO, text=True).strip()


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def main():
    initial = json.loads((A / 'manifest_initial.json').read_text())
    sealed = A / 'manifest_final.json'
    stamp = datetime.now(UTC).strftime('%Y%m%dT%H%M%SZ')
    out = A / 'rechecks' / stamp
    out.mkdir(exist_ok=False)
    packages = {d.metadata['Name'].lower(): d.version for d in importlib.metadata.distributions()}
    reference_root = Path(importlib.util.find_spec('morfeus').origin).parent
    reference = initial['reference_tools']['morfeus-ml']['python_source_files']
    identity = {
        'checked_utc': datetime.now(UTC).isoformat(),
        'commit': command('git', 'rev-parse', 'HEAD'),
        'tracked_changes': command('git', 'diff', '--name-only', 'HEAD', '--'),
        'executable_sha256': sha(REPO / 'target/release/stericx'),
        'rustc_verbose': command('rustc', '-vV'),
        'cargo': command('cargo', '-V'),
        'python': sys.version,
        'packages': packages,
        'reference_source_hashes': {name: sha(reference_root / name) for name in reference},
        'sealed_manifest_sha256': sha(sealed),
    }
    identity['matches_original'] = (
        identity['commit'] == initial['git_commit']
        and not identity['tracked_changes']
        and identity['executable_sha256'] == initial['binary']['sha256']
        and identity['rustc_verbose'] == initial['rustc']['stdout'].strip()
        and identity['cargo'] == initial['cargo']['stdout'].strip()
        and identity['python'] == initial['python']['version']
        and packages == initial['python']['packages']
        and all(identity['reference_source_hashes'][k] == v['sha256'] for k, v in reference.items())
    )
    write(out / 'current_identity.json', identity)
    receipt = json.loads((A / 'reproducibility/package_verification.json').read_text())
    assert receipt['passed'] and receipt['manifest_sha256'] == sha(sealed)
    write(out / 'sealed_package_verification.json', receipt)

    source = A / 'kraken/analysis/comparisons.csv'
    grouped = defaultdict(list)
    with source.open() as f:
        for row in csv.DictReader(f):
            grouped[row['metric']].append((float(row['reference']), float(row['sut'])))
    published = json.loads((A / 'kraken/analysis/metrics.json').read_text())
    results = {}
    for key, rows in grouped.items():
        n = len(rows)
        xs, ys = zip(*rows)
        mx, my = math.fsum(xs) / n, math.fsum(ys) / n
        errors = [y - x for x, y in rows]
        absolute = list(map(abs, errors))
        sse = math.fsum(e * e for e in errors)
        sst = math.fsum((x - mx) ** 2 for x in xs)
        slope = math.fsum((x - mx) * (y - my) for x, y in rows) / sst if sst else None
        actual = {'N': n, 'MAE': math.fsum(absolute) / n,
                  'RMSE': math.sqrt(sse / n), 'maximum_absolute_error': max(absolute),
                  'median_absolute_error': statistics.median(absolute),
                  'R2_one_to_one': 1 - sse / sst if sst else None,
                  'slope': slope, 'intercept': my - slope * mx if slope is not None else None}
        differences = {k: abs(v - published[key][k]) if v is not None and published[key][k] is not None else None for k, v in actual.items()}
        agrees = all(v is None and published[key][k] is None or v is not None and published[key][k] is not None and math.isclose(v, published[key][k], rel_tol=1e-10, abs_tol=1e-12) for k, v in actual.items())
        results[key] = {'recomputed': actual, 'absolute_statistic_differences': differences, 'agrees_with_report': agrees}
    statistics_check = {'scope': 'Independent standard-library recomputation of all eight reported statistics for each frozen Kraken comparison; not a new descriptor reference.', 'input_sha256': sha(source), 'metric_count': len(results), 'row_count': sum(map(len, grouped.values())), 'all_agree': set(results) == set(published) and all(v['agrees_with_report'] for v in results.values()), 'results': results}
    write(out / 'statistics_recheck.json', statistics_check)
    summary = {'purpose': __doc__, 'identity_matches': identity['matches_original'], 'sealed_artifacts_pass': receipt['passed'], 'statistics_pass': statistics_check['all_agree'], 'scientific_production_code_changed': False, 'original_sealed_manifest_sha256': sha(sealed)}
    write(out / 'review_checks.json', summary)
    write(out / 'manifest.json', {'created_utc': datetime.now(UTC).isoformat(), 'original_sealed_manifest_sha256': sha(sealed), 'script_sha256': sha(Path(__file__)), 'files': {p.name: sha(p) for p in sorted(out.iterdir()) if p.is_file()}})
    print(json.dumps({'directory': str(out), **summary}, indent=2))
    if not all([summary['identity_matches'], summary['sealed_artifacts_pass'], summary['statistics_pass']]):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
