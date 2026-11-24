"""Acquire original Kraken molecular metadata and exported conformers without replacement."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib,json,time,requests
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'.stericx/experimental_validation/structures'
IDS=[401,498,723,724,785,906,1057,1058,2062,2063,2064,2067]
API='https://descriptor-libraries.molssi.org/api/kraken'

def fetch(item):
 name,url=item;path=OUT/name
 path.parent.mkdir(parents=True,exist_ok=True)
 if path.exists(): raise FileExistsError(path)
 rcp=dict(url=url,path=str(path.relative_to(ROOT)),retrieved_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
 try:
  r=requests.get(url,timeout=45);path.write_bytes(r.content)
  rcp.update(status=r.status_code,sha256=hashlib.sha256(r.content).hexdigest(),bytes=len(r.content))
 except requests.RequestException as e:rcp['error']=str(e)
 path.with_suffix(path.suffix+'.receipt.json').write_text(json.dumps(rcp,indent=2)+'\n')
 return rcp

def main():
 OUT.mkdir(parents=True,exist_ok=False)
 with ThreadPoolExecutor(max_workers=3) as pool:
  receipts=list(pool.map(fetch,[(f'{i}/molecule.json',f'{API}/molecules/{i}') for i in IDS]))
 print('metadata',[(r['path'],r.get('status')) for r in receipts],flush=True)
 jobs=[]
 for i,r in zip(IDS,receipts,strict=True):
  if r.get('status')!=200:continue
  d=json.loads((OUT/f'{i}/molecule.json').read_text())
  for cid in d.get('conformers_id',[]):
   jobs.extend([(f'{i}/{cid}.sdf',f'{API}/conformers/export/{cid}.sdf'),(f'{i}/{cid}.json',f'{API}/conformers/{cid}')])
 print('conformer downloads',len(jobs),flush=True)
 with ThreadPoolExecutor(max_workers=3) as pool:
  for r in pool.map(fetch,jobs):
   receipts.append(r)
   (OUT/'receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')
 print('completed',len(receipts),'failures',sum(r.get('status')!=200 for r in receipts),flush=True)
 (OUT/'receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')
if __name__=='__main__':main()
