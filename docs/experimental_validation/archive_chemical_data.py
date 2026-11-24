"""Create/restore deterministic research bundles without rewriting existing inputs."""
from pathlib import Path
import argparse,gzip,json,hashlib
ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/'.stericx/experimental_validation';DEST=Path(__file__).with_name('data')
def sha(data):return hashlib.sha256(data).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('mode',choices=['archive','restore']);a=p.parse_args()
 if a.mode=='archive':
  DEST.mkdir(exist_ok=True)
  groups={'published_inputs':list((BASE/'structures').rglob('*.sdf'))+list((BASE/'structures').rglob('molecule.json'))+[BASE/'structures/receipts.json',BASE/'sources/nihda.csv',BASE/'chemistry_preparation/requests.json',BASE/'chemistry_preparation/audit.json',BASE/'stereochemistry_v1/audit.json'],'generated_inputs':list((BASE/'generated_v1').rglob('*.sdf'))+[BASE/'generated_v1/audit.json',BASE/'generated_v1/requests.json',BASE/'generated_v1/post_stereo_audit.json'],'native_results':[BASE/'chemistry_preparation/native.jsonl',BASE/'generated_v1/native.jsonl',BASE/'chemical_references_v2/summary.json']}
  groups['corrected_generated_inputs']=list((BASE/'generated_v2').rglob('*.sdf'))+[BASE/'generated_v2/audit.json',BASE/'generated_v2/requests.json',BASE/'generated_v2/native.jsonl']
  groups['analysis_and_provenance']=list((BASE/'chemical_analysis_v1').glob('*.json'))+list((BASE/'chemical_analysis_v2').glob('*.json'))+list((BASE/'sources').glob('*.json'))+list((BASE/'structures').rglob('*.json'))+[BASE/'conformer_sensitivity_v2/summary.json',BASE/'chemical_analysis_v1/protocol_at_analysis.md',BASE/'generated_v1/controller.py',BASE/'sources/table_s3_original.txt',BASE/'sources/table_s3_corrected.txt']
  groups['numerical_inputs']=list((BASE/'numerical_frozen').glob('*'))
  manifest={}
  for name,paths in groups.items():
   rows=[]
   for path in sorted(set(paths)):
    if not path.exists():raise FileNotFoundError(path)
    data=path.read_bytes();rows.append(dict(path=str(path.relative_to(ROOT)),sha256=sha(data),text=data.decode('utf8')))
   dest=DEST/(name+'.jsonl.gz');data=gzip.compress(('\n'.join(json.dumps(r,separators=(',',':')) for r in rows)+'\n').encode(),mtime=0)
   if dest.exists():assert dest.read_bytes()==data, 'Refuse to replace existing bundle: '+str(dest)
   else:dest.write_bytes(data)
   manifest[dest.name]=dict(files=len(rows),sha256=sha(dest.read_bytes()),uncompressed_files_bytes=sum(len(r['text'].encode()) for r in rows))
  (DEST/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 else:
  for name,meta in json.loads((DEST/'manifest.json').read_text()).items():
   path=DEST/name;assert sha(path.read_bytes())==meta['sha256']
   for row in map(json.loads,gzip.decompress(path.read_bytes()).splitlines()):
    dest=(ROOT/row['path']).resolve();assert dest.is_relative_to(BASE.resolve());data=row['text'].encode();assert sha(data)==row['sha256']
    if dest.exists():assert sha(dest.read_bytes())==row['sha256'],dest
    else:dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(data)
  print('All research bundle bytes restored or verified.')
if __name__=='__main__':main()
