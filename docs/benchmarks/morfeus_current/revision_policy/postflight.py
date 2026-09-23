"""Post-timing freeze/input verification; never supplies performance measurements."""
import hashlib,importlib.metadata,json,sys
from pathlib import Path
from datetime import datetime,timezone
import morfeus
R=Path(__file__).resolve().parent;H=R.parent;ROOT=H.parents[2]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
frozen=read(H/'freeze.json')
for name,digest in frozen['source_sha256'].items():assert sha(ROOT/name)==digest
assert sha(H/'build/stericx')==frozen['executable_sha256']
package=Path(morfeus.__file__).parent
for name,digest in read(H/'environment/reference-sources.stdout').items():assert sha(package/name)==digest
versions={}
for line in (H/'environment/requirements.txt').read_text().splitlines():
 name,version=line.split('==');assert importlib.metadata.version(name)==version;versions[name]=version
assert sys.version_info[:3]==(3,12,13)
manifest=read(H/'input_manifest.json')
for row in manifest['corpus']:assert sha(H/row['filename'])==row['sha256']
for row in manifest['repeated_10000']:assert sha(H/'inputs/repeated_10000'/row['filename'])==row['sha256']
(R/'postflight.json').write_text(json.dumps(dict(passed=True,checked_utc=datetime.now(timezone.utc).isoformat(),original_inputs=56,repeated_inputs=10000,stericx_source_files=len(frozen['source_sha256']),morfeus_source_files=len(read(H/'environment/reference-sources.stdout')),versions=versions,production_unchanged=True),indent=2)+'\n')
print('Postflight passed: production, reference package, environment and all 10056 input files unchanged.')
