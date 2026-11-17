"""Source-bound occupancy/precision campaign utilities; never overwrite runs."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import time

ROOT = Path(__file__).resolve().parents[2]
EVIDENCE = ROOT / '.stericx/occupancy_precision'
DOCS = ROOT / 'docs/occupancy_precision'
TEST_FIXTURES = [
    'docs/scientific_remediation/geometry/focused_v1/requests.jsonl',
    'docs/scientific_remediation/geometry/planar_precision/requests.jsonl',
    'docs/scientific_remediation/geometry/topology_witnesses.jsonl',
    'docs/scientific_remediation/geometry/finite_mean/finite_mean_1.sdf',
    'docs/scientific_remediation/geometry/finite_mean/finite_mean_4.sdf',
    'docs/study_001/stericx_portable_model.json',
    'docs/study_001/stericx_model.json',
    'docs/study_001/stericx_frozen_predictions.csv',
    'docs/study_001/stericx_evaluation.json',
    'docs/study_007/crosscoupling_metrics.json',
    'docs/study_009/pd_crosscoupling_metrics.json',
]


def sha(path):
    return hashlib.file_digest(Path(path).open('rb'), 'sha256').hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + '\n')


def capture(argv):
    p = subprocess.run(argv, cwd=ROOT, text=True, capture_output=True)
    return {'argv': argv, 'returncode': p.returncode, 'stdout': p.stdout, 'stderr': p.stderr}


def freeze():
    out = EVIDENCE / 'baseline'
    out.mkdir(parents=True, exist_ok=False)
    paths = subprocess.check_output(['git', 'ls-files', '-z'], cwd=ROOT).decode().split('\0')
    selected = [p for p in paths if p and (
        p.startswith(('src/', 'tests/', 'data/', 'examples/', 'python/', 'scripts/',
                      'docs/ultimate_engine/', 'docs/validation/dbstep/'))
        or p in ('Cargo.toml', 'Cargo.lock', 'pyproject.toml'))]
    records = {p: {'sha256': sha(ROOT / p), 'bytes': (ROOT / p).stat().st_size}
               for p in selected}
    source = out / 'source'
    for p in selected:
        if p.startswith(('src/', 'tests/', 'examples/', 'python/', 'data/')) or p in ('Cargo.toml', 'Cargo.lock'):
            target = source / p
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ROOT / p, target)
    with tarfile.open(out / 'source.tar.gz', 'w:gz') as archive:
        archive.add(source, arcname='source')
    retained = [ROOT / '.stericx/ultimate_engine/baseline_observations/geometry/stdout.jsonl',
                ROOT / '.stericx/ultimate_engine/n2_molecular_validation/requests.json']
    # Freeze all historical response streams without changing their original locations.
    retained += sorted((ROOT / '.stericx/ultimate_engine/baseline_observations').glob('*/stdout.jsonl'))
    retained += sorted((ROOT / '.stericx/ultimate_engine').glob('n2*/requests.json'))
    workload = ROOT / '.stericx/profiling/scientifically_validated_optimization/corrected_workloads_v2/workloads/manifest.json'
    retained.append(workload)
    manifest = {
        'created_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        'commit': capture(['git', 'rev-parse', 'HEAD'])['stdout'].strip(),
        'initial_status': capture(['git', 'status', '--short']),
        'records': records,
        'retained_inputs': {str(p.relative_to(ROOT)): sha(p) for p in sorted(set(retained)) if p.is_file()},
        'source_archive_sha256': sha(out / 'source.tar.gz'),
        'host': {key: capture(argv) for key, argv in {
            'rustc': ['rustc', '-Vv'], 'cargo': ['cargo', '-V'], 'cpu': ['lscpu'],
            'kernel': ['uname', '-a'], 'perf': ['perf', 'stat', '-e',
                'cycles,instructions,branches,branch-misses,cache-references,cache-misses', '--', 'true'],
        }.items()},
        'environment': {key: os.environ.get(key) for key in
            ['RUSTFLAGS', 'CARGO_BUILD_TARGET', 'RAYON_NUM_THREADS', 'CARGO_PROFILE_RELEASE_LTO']},
        'compiler_contract': 'Portable release; thin LTO; one codegen unit; no fast math/FMA contraction or native target override. Observer debug=1. CLI RUSTFLAGS=-C debuginfo=1.',
        'protocol': {'threads': [1, 2, 4, 6, 12], 'warmups': 1,
                     'final_alternating_pairs': 8, 'absolute_tolerances_angstrom3': [.1, .05, .01],
                     'legacy_molecular_box_budget': 262143,
                     'final_host_gate': 'Same existing <=0.15 busy cores and <=0.05 I/O wait cores endpoint gate; no other apps stopped.'},
    }
    save(out / 'manifest.json', manifest)
    save(DOCS / 'baseline_manifest.json', manifest)
    print(out, flush=True)


def build(name, quick=False):
    out = EVIDENCE / name
    if name != 'baseline':
        out.mkdir(parents=True, exist_ok=False)
        for directory in ['src', 'tests', 'examples', 'python', 'data']:
            shutil.copytree(ROOT / directory, out / 'source' / directory,
                            ignore=shutil.ignore_patterns('target', '__pycache__', '.venv', '*.egg-info'))
        for file in ['Cargo.toml', 'Cargo.lock']:
            shutil.copy2(ROOT / file, out / 'source' / file)
        for file in TEST_FIXTURES:
            destination = out / 'source' / file
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ROOT / file, destination)
    source = out / 'source'
    observer = out / 'observer'
    shutil.copytree(ROOT / 'docs/ultimate_engine/observer', observer)
    manifest = observer / 'Cargo.toml'
    manifest.write_text(manifest.read_text().replace('path = "../../.."', f'path = "{source}"'))
    shutil.copy2(ROOT / '.stericx/ultimate_engine/baseline/observer/Cargo.lock', observer / 'Cargo.lock')
    target = Path('/tmp/stericx-occupancy-precision-target') / ('candidates' if quick else name)
    env = dict(os.environ, CARGO_BUILD_JOBS='2', CARGO_TARGET_DIR=str(target))
    env.pop('RUSTFLAGS', None)
    receipt = {'source': {str(p.relative_to(source)): sha(p) for p in sorted(source.rglob('*.rs'))},
               'observer': {str(p.relative_to(observer)): sha(p) for p in sorted(observer.rglob('*')) if p.is_file()},
               'commands': [], 'binaries': {}, 'complete': False}
    (out / 'bin').mkdir()

    def run(label, argv, extra_env=None):
        start = time.monotonic()
        with (out / f'{label}.log').open('x') as log:
            p = subprocess.run(list(map(str, argv)), env=dict(env, **(extra_env or {})),
                               cwd=source, stdout=log, stderr=subprocess.STDOUT)
        receipt['commands'].append({'label': label, 'argv': list(map(str, argv)),
                                   'returncode': p.returncode, 'seconds': time.monotonic()-start,
                                   'env': {'CARGO_TARGET_DIR': str(target), **(extra_env or {})}})
        save(out / 'build.json', receipt)
        p.check_returncode()
        print(name, label, 'complete', flush=True)

    def keep(src, dest):
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest)
        receipt['binaries'][str(dest.relative_to(out))] = sha(dest)

    bins = ['stericx-ultimate-observer'] if quick else ['stericx-ultimate-observer', 'continuous_volume', 'continuous', 'accuracy', 'volume_workflow']
    run('observer', ['cargo', 'build', '--locked', '--offline', '--release', '--manifest-path', observer / 'Cargo.toml', *(['--bin', bins[0]] if quick else ['--bins'])])
    for binary in bins:
        keep(target / 'release' / binary, out / 'bin' / ('observer' if binary == 'stericx-ultimate-observer' else binary))
    if quick:
        receipt['complete'] = True
        receipt['scope'] = 'Exploratory native observer only; CLI/allocation observers are built for final acceptance separately.'
        save(out / 'build.json', receipt)
        return
    run('observer_alloc', ['cargo', 'build', '--locked', '--offline', '--release', '--manifest-path', observer / 'Cargo.toml', '--bin', 'stericx-ultimate-observer', '--features', 'profiling'])
    keep(target / 'release/stericx-ultimate-observer', out / 'bin/observer_alloc')
    for label, flags in [('native', []), ('diagnostic', ['--features', 'profiling'])]:
        run(label, ['cargo', 'build', '--locked', '--offline', '--release', '--bin', 'stericx', *flags], {'RUSTFLAGS': '-C debuginfo=1'})
        keep(target / 'release/stericx', out / f'{label}-target/release/stericx')
    receipt['complete'] = True
    save(out / 'build.json', receipt)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['freeze', 'build'])
    parser.add_argument('--name', default='baseline')
    parser.add_argument('--quick', action='store_true')
    args = parser.parse_args()
    freeze() if args.action == 'freeze' else build(args.name, args.quick)
