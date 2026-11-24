"""Freeze scientific source, tracked identities and external evidence; refuse overwrite."""
from pathlib import Path
import hashlib
import importlib.metadata
import json
import os
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / '.stericx/experimental_validation/baseline'
DOC = ROOT / 'docs/experimental_validation'

def run(args):
    try:
        p = subprocess.run(args, cwd=ROOT, capture_output=True, text=True)
        return dict(argv=args, returncode=p.returncode, stdout=p.stdout, stderr=p.stderr)
    except OSError as e:
        return dict(argv=args, returncode=None, stdout='', stderr=str(e), available=False)

def sha(p):
    with p.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def main():
    OUT.mkdir(parents=True, exist_ok=False)
    head = run(['git', 'rev-parse', 'HEAD'])['stdout'].strip()
    tracked = subprocess.check_output(['git', 'ls-files', '-s', '-z'], cwd=ROOT)
    (OUT / 'tracked_index.z').write_bytes(tracked)
    paths = [x.split(b'\t', 1)[1].decode() for x in tracked.split(b'\0') if x]
    selected = [p for p in paths if p.startswith(('src/', 'tests/', 'data/', 'python/stericx/', 'python/src/', 'examples/')) or p in ('Cargo.toml', 'Cargo.lock', 'pyproject.toml', 'uv.lock')]
    selected += [p for p in paths if p.startswith(('docs/occupancy_precision/', 'docs/ultimate_engine/', 'docs/validation/dbstep/')) and p.endswith(('.md', '.json', '.py', '.rs')) and '/source/' not in p and '/baseline/' not in p]
    records = {}
    for p in sorted(set(selected)):
        path = ROOT / p
        if path.is_file():
            records[p] = dict(sha256=sha(path), bytes=path.stat().st_size)
    archived = OUT / 'production_source.tar'
    with archived.open('xb') as f:
        subprocess.run(['git', 'archive', head, 'src', 'tests', 'examples', 'Cargo.toml', 'Cargo.lock', 'python'], cwd=ROOT, stdout=f, check=True)
    retained = list((ROOT / '.stericx/ultimate_engine/baseline_observations').glob('*/stdout.jsonl'))
    retained += list((ROOT / '.stericx/occupancy_precision').glob('*analytical*'))
    retained += list((ROOT / '.stericx/ultimate_engine').glob('n2_validation*/requests.json'))
    manifest = dict(created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), production_commit=head,
        initial_status=run(['git', 'status', '--short']), tracked_index_sha256=sha(OUT/'tracked_index.z'),
        source_archive_sha256=sha(archived), records=records,
        retained_evidence={str(p.relative_to(ROOT)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(set(retained)) if p.is_file()},
        toolchain={key:run(args) for key,args in {
            'rustc':['rustc','-Vv'], 'cargo':['cargo','-V'], 'toolchains':['rustup','toolchain','list'],
            'targets':['rustup','target','list','--installed'], 'kernel':['uname','-a'], 'cpu':['lscpu'],
            'perf':['perf','stat','-e','cycles,instructions','--','true'], 'valgrind':['valgrind','--version'],
        }.items()}, python=sys.version, packages={d.metadata['Name']:d.version for d in importlib.metadata.distributions()},
        env={k:os.environ.get(k) for k in ['RUSTFLAGS','CARGO_BUILD_TARGET','RAYON_NUM_THREADS','CARGO_PROFILE_RELEASE_LTO']})
    text=json.dumps(manifest,indent=2,sort_keys=True)+'\n'
    (OUT/'manifest.json').write_text(text)
    (DOC/'baseline_manifest.json').write_text(text)
    print(json.dumps(dict(commit=head, tracked=len(paths), sha256_records=len(records), retained=len(manifest['retained_evidence']))))

if __name__ == '__main__':
    main()
