#!/usr/bin/env python3
"""One source-bound Python timing block, with scientific construction checks."""
import argparse,hashlib,json,time
from dataclasses import asdict
from pathlib import Path
import numpy as np
import stericx as sx
from stericx._types import ensemble_result

ROOT=Path(__file__).resolve().parents[2]
p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode');p.add_argument('iterations',type=int)
p.add_argument('--threads',type=int,default=1);a=p.parse_args();assert a.iterations>0
ensemble=sx.Ensemble.from_file(ROOT/'python/examples/data/ligand_401_ensemble.xyz')
config=json.loads((ROOT/'python/examples/data/ligand_401_options.json').read_text())
config['volume']=sx.VolumeOptions(**config.get('volume',{}));options=sx.DescriptorOptions(**config)
encoded=json.dumps(asdict(options))
def bridge():return ensemble._engine.analyze(encoded,'{"missing":"require_complete"}','null',None)
raw=bridge();parsed=json.loads(raw);result=ensemble_result(parsed)
assert result.to_dict()==parsed==ensemble.descriptors(options).to_dict()
array=result.to_numpy();assert len(array)==6
for i,row in enumerate(result.to_records()):
    for field in array.dtype.names:
        if array[field].dtype.kind=='f' and row[field] is not None:assert array[field][i]==row[field]
manifest=ROOT/'.stericx/profiling/scientifically_validated_optimization/corrected_workloads_v2/workloads/manifest.json'
paths=[p for p in json.loads(manifest.read_text())['workloads']['conformers_56']['argv'] if p.endswith('.xyz')]
assert len(paths)==56
batch_options=sx.DescriptorOptions(axis='coordination')
calls={'bridge':bridge,'public':lambda:ensemble.descriptors(options),'json_decode':lambda:json.loads(raw),
    'typed_construction':lambda:ensemble_result(parsed),'numpy_export':result.to_numpy,
    'options_encode':lambda:sx._options(options),'record_export':lambda:ensemble.records,
    'per_file_56':lambda:tuple(sx.Molecule.from_file(p).descriptors(batch_options) for p in paths),
    'batch_56':lambda:sx.calculate_files(paths,batch_options,threads=a.threads)}
call=calls[a.mode];call()
start=time.perf_counter_ns()
for _ in range(a.iterations):value=call();del value
elapsed=time.perf_counter_ns()-start
value=call()
if a.mode=='bridge':checked=json.loads(value)
elif a.mode=='public' or a.mode=='typed_construction':checked=value.to_dict()
elif a.mode=='numpy_export':checked={'dtype':value.dtype.descr,'rows':value.tolist()}
elif a.mode=='record_export':checked=value
elif a.mode in ['per_file_56','batch_56']:checked=[v.to_dict() if hasattr(v,'to_dict') else asdict(v) for v in value]
elif a.mode=='options_encode':checked=json.loads(value)
else:checked=value
fingerprint=hashlib.sha256(json.dumps(checked,sort_keys=True,allow_nan=False).encode()).hexdigest()
print(json.dumps({'mode':a.mode,'iterations':a.iterations,'elapsed_ns':elapsed,'result_sha256':fingerprint,
    'module':sx.__file__,'native':sx._native.__file__,
    'native_sha256':hashlib.sha256(Path(sx._native.__file__).read_bytes()).hexdigest(),'numpy_version':np.__version__}))
