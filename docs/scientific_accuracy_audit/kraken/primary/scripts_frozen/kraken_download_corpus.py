#!/usr/bin/env python3
"""Download every primary-paper ligand and its available DFT conformers.

The primary publication's ID universe is fixed before geometry eligibility.
No SUT descriptor, geometric cutoff, or published residual filters the inputs.
"""
from __future__ import annotations
import argparse,csv,datetime,hashlib,json,sys
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
from kraken_acquire import ROOT,API,acquire

def save(path,value):path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers',type=int,default=8)
    args=parser.parse_args()
    root=ROOT/'primary';root.mkdir(parents=True,exist_ok=True)
    ref=ROOT/'raw_sources/ni_hda_raw.body'
    with ref.open() as f: ids=sorted(int(r['']) for r in csv.DictReader(f))
    assert len(ids)==1566 and len(set(ids))==1566
    plan={'created_utc':datetime.datetime.now(datetime.UTC).isoformat(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'acquirer_sha256':hashlib.sha256(Path(__file__).with_name('kraken_acquire.py').read_bytes()).hexdigest(),'id_source':str(ref),'id_source_sha256':hashlib.sha256(ref.read_bytes()).hexdigest(),'ids':ids,'eligibility_policy':'Acquire all primary-paper IDs first; independent chemical eligibility is a later recorded stage.','workers':args.workers}
    if not (root/'acquisition_plan.json').exists():save(root/'acquisition_plan.json',plan)
    def ligand(mid):
        row={'molecule_id':mid,'responses':{}}
        for kind,endpoint in [('molecule',f'molecules/{mid}'),('data_types',f'molecules/{mid}/data_types'),('published_dft',f'molecules/data/{mid}?data_type=dft')]:
            stem=root/kind/str(mid)
            try:meta=acquire(f'{API}/{endpoint}',stem)
            except Exception as error:meta={'error':repr(error)}
            row['responses'][kind]={'metadata':str(stem)+'.metadata.json','status':meta.get('status'),'error':meta.get('error')}
        try:
            types=json.loads((root/'data_types'/f'{mid}.body').read_text())
            obj=json.loads((root/'molecule'/f'{mid}.body').read_text())
            row.update(dft_available=bool(types.get('dft_data')),conformer_ids=obj.get('conformers_id',[]),smiles=obj.get('smiles'))
        except Exception as error:row['parse_error']=repr(error)
        return row
    ligands={}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(ligand,mid) for mid in ids]
        for n,future in enumerate(as_completed(futures),1):
            row=future.result();ligands[row['molecule_id']]=row
            if n%100==0:print(f'Metadata {n}/{len(ids)}',flush=True)
    save(root/'ligands_manifest.json',{'ligands':[ligands[k] for k in sorted(ligands)]})
    cids={cid for row in ligands.values() if row.get('dft_available') for cid in row.get('conformer_ids',[])}
    print(f'DFT geometries to acquire: {len(cids)}',flush=True)
    def conformer(cid):
        stem=root/'conformers'/str(cid)
        try:m=acquire(f'{API}/conformers/export/{cid}.sdf',stem)
        except Exception as error:m={'error':repr(error)}
        return {'conformer_id':cid,'metadata':str(stem)+'.metadata.json','status':m.get('status'),'error':m.get('error'),'sha256':m.get('sha256'),'bytes':m.get('bytes')}
    conformers={}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(conformer,cid) for cid in sorted(cids)]
        for n,future in enumerate(as_completed(futures),1):
            row=future.result();conformers[row['conformer_id']]=row
            if n%500==0:print(f'Geometries {n}/{len(cids)}; successful {sum(r["status"]==200 for r in conformers.values())}',flush=True)
    save(root/'conformers_manifest.json',{'conformers':[conformers[k] for k in sorted(conformers)]})
    save(root/'acquisition_summary.json',{'finished_utc':datetime.datetime.now(datetime.UTC).isoformat(),'ligands_attempted':len(ligands),'ligands_dft_available':sum(r.get('dft_available',False) for r in ligands.values()),'conformers_attempted':len(cids),'conformers_successful':sum(r['status']==200 for r in conformers.values()),'failed_conformers':[r for r in conformers.values() if r['status']!=200]})
    print('Acquisition complete',flush=True)
if __name__=='__main__':main()
