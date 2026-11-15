#!/usr/bin/env python3
"""Freeze current engine, then build source-bound native/allocation observers."""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('name')
args = parser.parse_args()
out = ROOT / '.stericx/ultimate_engine' / args.name
out.mkdir(exist_ok=False)
source = out / 'source'
source.mkdir()
for name in ('src', 'Cargo.toml', 'Cargo.lock'):
    p = ROOT / name
    if p.is_dir():
        shutil.copytree(p, source/name)
    else:
        shutil.copyfile(p, source/name)
observer = out/'observer'
shutil.copytree(ROOT/'docs/ultimate_engine/observer', observer)
p = observer/'Cargo.toml'
p.write_text(p.read_text().replace('path = "../../.."', f'path = "{source}"'))
shutil.copyfile(ROOT/'.stericx/ultimate_engine/baseline/observer/Cargo.lock', observer/'Cargo.lock')
manifest = {'source': {str(p.relative_to(source)): hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in sorted(source.rglob('*')) if p.is_file()}, 'builds': []}
(out/'bin').mkdir()
env = dict(os.environ, CARGO_BUILD_JOBS='2', CARGO_TARGET_DIR=str(out/'target'))
env.pop('RUSTFLAGS', None)
for name, extra in [('observer', []), ('observer_alloc', ['--features', 'profiling'])]:
    command = ['cargo', 'build', '--locked', '--offline', '--release',
               '--manifest-path', str(observer/'Cargo.toml'), '--bin', 'stericx-ultimate-observer', *extra]
    start = time.monotonic()
    with (out/f'{name}_build.log').open('w') as log:
        process = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT)
    row = {'command': command, 'returncode': process.returncode, 'seconds': time.monotonic()-start}
    manifest['builds'].append(row)
    (out/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    process.check_returncode()
    binary = out/'bin'/name
    shutil.copyfile(out/'target/release/stericx-ultimate-observer', binary)
    binary.chmod(0o755)
    row['binary_sha256'] = hashlib.sha256(binary.read_bytes()).hexdigest()
manifest['complete'] = True
(out/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
print(out)
