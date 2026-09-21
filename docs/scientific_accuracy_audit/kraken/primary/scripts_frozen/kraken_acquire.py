#!/usr/bin/env python3
"""Acquire immutable primary Kraken API/source responses with provenance."""
from __future__ import annotations
import argparse,datetime,hashlib,json,threading,time
from pathlib import Path
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
ROOT=Path(__file__).resolve().parents[1]/'audit/kraken'
API='https://descriptor-libraries.molssi.org/api/kraken'
_local=threading.local()
def session():
    if not hasattr(_local,'session'):
        s=requests.Session()
        s.headers['User-Agent']='StericX independent scientific-reproducibility audit (public research data)'
        retry=Retry(total=4,backoff_factor=1,status_forcelist=[429,500,502,503,504],allowed_methods=['GET'])
        s.mount('https://',HTTPAdapter(max_retries=retry))
        _local.session=s
    return _local.session

def acquire(url: str,stem: Path):
    stem.parent.mkdir(parents=True,exist_ok=True)
    body=Path(str(stem)+'.body')
    meta=Path(str(stem)+'.metadata.json')
    if meta.exists():
        m=json.loads(meta.read_text())
        if m.get('status')==200 and body.is_file():
            assert hashlib.sha256(body.read_bytes()).hexdigest()==m['sha256']
            assert m['url']==url
            return m
        raise RuntimeError(f'Prior unsuccessful acquisition is immutable: {meta}')
    m={'url':url,'started_utc':datetime.datetime.now(datetime.UTC).isoformat()}
    start=time.monotonic()
    try:
        r=session().get(url,timeout=(20,90))
        body.write_bytes(r.content)
        m.update(status=r.status_code,final_url=r.url,headers=dict(r.headers),bytes=len(r.content),sha256=hashlib.sha256(r.content).hexdigest())
    except Exception as error:
        m['error']=repr(error)
    m.update(finished_utc=datetime.datetime.now(datetime.UTC).isoformat(),elapsed_seconds=time.monotonic()-start)
    meta.write_text(json.dumps(m,indent=2,sort_keys=True)+'\n')
    return m

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('url');parser.add_argument('stem',type=Path)
    args=parser.parse_args()
    print(json.dumps(acquire(args.url,args.stem),indent=2))
if __name__=='__main__':main()
