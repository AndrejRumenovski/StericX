"""Cache primary-source bytes and HTTP receipts, including unsuccessful retrievals."""
from pathlib import Path
import hashlib
import json
import time
import requests

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'.stericx/experimental_validation/sources'
URLS={
 'nihda_repo.json':'https://api.github.com/repos/SigmanGroup/Ni-Catalyzed-hDA/commits/main',
 'nihda_tree.json':'https://api.github.com/repos/SigmanGroup/Ni-Catalyzed-hDA/git/trees/main?recursive=1',
 'nihda.csv':'https://raw.githubusercontent.com/SigmanGroup/Ni-Catalyzed-hDA/main/data/kraken.csv',
 'nihda_notebook.json':'https://raw.githubusercontent.com/SigmanGroup/Ni-Catalyzed-hDA/main/Enantioselectivity_Model.ipynb',
 'nihda_si.pdf':'https://authors.library.caltech.edu/records/ef9d0-y2c15/files/ja5c09948_si_001.pdf?download=1',
 'nihda_zenodo.json':'https://zenodo.org/api/records/15338599',
 'bisphosphine_tree.json':'https://api.github.com/repos/SigmanGroup/BisphosphineConformerSelection/git/trees/main?recursive=1',
 'paton_tools.html':'https://patonlab.com/code/',
 'zahrt.html':'https://pmc.ncbi.nlm.nih.gov/articles/PMC6417887/',
 'reid.html':'https://pmc.ncbi.nlm.nih.gov/articles/PMC6641578/',
}

def main():
 OUT.mkdir(parents=True,exist_ok=True)
 receipts=[]
 for name,url in URLS.items():
  path=OUT/name
  if path.exists():
   raise FileExistsError(path)
  receipt={'url':url,'path':str(path.relative_to(ROOT)),'retrieved_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
  try:
   r=requests.get(url,timeout=45)
   path.write_bytes(r.content)
   receipt.update(status=r.status_code,final_url=r.url,content_type=r.headers.get('content-type'),sha256=hashlib.sha256(r.content).hexdigest(),bytes=len(r.content))
  except requests.RequestException as e:
   receipt['error']=str(e)
  receipts.append(receipt)
  (OUT/'receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')
  print(name,receipt.get('status'),receipt.get('bytes'),flush=True)
 (ROOT/'docs/experimental_validation/source_receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')

if __name__=='__main__':main()
