"""Freeze the system under investigation before independent audit comparisons."""
from __future__ import annotations
import datetime
import hashlib
import importlib.metadata
import json
import platform
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
AUDIT = REPO / 'docs/scientific_accuracy_audit'

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def command(args):
    r = subprocess.run(args, cwd=REPO, text=True, capture_output=True)
    return {'argv': args, 'returncode': r.returncode, 'stdout': r.stdout, 'stderr': r.stderr}

def main():
    target = AUDIT / 'manifest_initial.json'
    if target.exists():
        raise RuntimeError('Refusing to overwrite frozen initial audit manifest')
    tracked = subprocess.check_output(['git', 'ls-files', '-z'], cwd=REPO).decode().split('\0')
    files = {}
    for name in tracked:
        if not name or name.startswith(('docs/profiling/', 'docs/media/', 'docs/scientific_accuracy_audit/')):
            continue
        p = REPO / name
        wanted = name in ('Cargo.toml', 'Cargo.lock', 'pyproject.toml', 'uv.lock', 'README.md') or (
            name.startswith(('src/', 'data/', 'studies/', 'scripts/', 'tests/', 'docs/'))
            and p.suffix.lower() in ('.rs', '.py', '.c', '.csv', '.tsv', '.json', '.jsonl', '.md', '.xyz', '.sdf', '.sigpack', '.toml', '.lock'))
        if not wanted or not p.is_file():
            continue
        dst = AUDIT / 'frozen/repository' / name
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, dst)
        files[name] = {'bytes': p.stat().st_size, 'sha256': sha(p), 'frozen_path': str(dst.relative_to(AUDIT))}
    binary = REPO / 'target/release/stericx'
    frozen_binary = AUDIT / 'frozen/bin/stericx'
    frozen_binary.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(binary, frozen_binary)
    packages = {d.metadata['Name']: d.version for d in importlib.metadata.distributions() if d.metadata['Name']}
    import morfeus
    reference_root = Path(morfeus.__file__).parent
    reference_files = {}
    for source in reference_root.rglob('*.py'):
        rel = source.relative_to(reference_root)
        dst = AUDIT / 'frozen/reference_sources/morfeus_0.8.0' / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, dst)
        reference_files[str(rel)] = {'sha256': sha(source), 'frozen_path': str(dst.relative_to(AUDIT))}
    manifest = {
        'schema_version': 1,
        'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'git_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip(),
        'git_status_before_freeze': command(['git','status','--short']),
        'binary': {'original': str(binary), 'frozen_path': str(frozen_binary.relative_to(AUDIT)), 'bytes': binary.stat().st_size, 'sha256': sha(binary), 'version': command([str(frozen_binary),'--version'])},
        'python': {'version':sys.version,'executable':sys.executable,'platform':platform.platform(),'packages':dict(sorted(packages.items()))},
        'rustc': command(['rustc','-vV']), 'cargo': command(['cargo','-V']), 'uv':command(['uv','--version']),
        'repository_files':files,
        'reference_tools': {'morfeus-ml': {'version':importlib.metadata.version('morfeus-ml'),'python_source_files':reference_files}},
        'policy': 'This initial snapshot is immutable. Generated/downloaded audit inputs must be hashed in phase manifests before comparisons; raw StericX outputs must be frozen before independent-result interpretation. Production scientific code is not modified. Existing outputs/tests are observations, not scientific truth.'
    }
    target.write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
    print('Frozen',len(files),'repository inputs/source files,',len(reference_files),'reference source files and executable',manifest['binary']['sha256'])
if __name__ == '__main__':
    main()
