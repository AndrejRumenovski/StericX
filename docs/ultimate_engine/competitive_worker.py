#!/usr/bin/env python3
"""Matched finite-grid reference workflows; DBSTEP's default API is not timed."""
import argparse
import importlib.util
import json
import multiprocessing
import os
import time
from pathlib import Path
from types import SimpleNamespace

for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[key]='1'
import numpy as np
from dbstep import sterics
from morfeus.geometry import Sphere

ROOT=Path(__file__).resolve().parents[2]
HISTORICAL=ROOT/'docs/benchmarks/morfeus_current/revision_policy/admitted_campaign/morfeus_batch.py'
spec=importlib.util.spec_from_file_location('historical_morfeus_batch',HISTORICAL)
mb=importlib.util.module_from_spec(spec);spec.loader.exec_module(mb)

def dbstep(path,details=False):
    elements,xyz,donor,neighbors,center,radii=mb.geometry(path)
    keep=np.array([e!='H' for e in elements])
    frames=[]
    for plane in neighbors:
        z=center-xyz[donor];z/=np.linalg.norm(z)
        v=xyz[plane]-center;x=v-np.dot(v,z)*z;x/=np.linalg.norm(x)
        y=np.cross(z,x);aligned=(xyz-center)@np.array([x,y,z]).T
        sphere=Sphere(np.zeros(3),3.5,method='projection',density=.01,filled=True)
        occupied,_,_=sterics.occupied(sphere.points,aligned[keep],radii[keep],
             np.zeros(3),SimpleNamespace(grid=1.,verbose=False,debug=False))
        def counts(points):
            x,y,z=points.T
            q=np.where(y>=0,np.where(x>=0,0,1),np.where(x<0,2,3))
            return np.bincount(q+np.where(z>=0,0,4),minlength=8)
        total=counts(sphere.points);hit=counts(occupied)
        octs=[float(n)/float(d)*sphere.volume/8 for n,d in zip(hit,total,strict=True)]
        qs=[octs[i]+octs[i+4] for i in range(4)]
        row={'plane':plane,'buried_volume':len(occupied)/len(sphere.points)*sphere.volume,
             'octants':octs,'quadrants':qs,'near_vbur':sum(octs[4:]),'far_vbur':sum(octs[:4])}
        if details:row.update(points=sphere.points,occupied=occupied)
        frames.append(row)
    q=np.array([f['quadrants'] for f in frames]);o=np.array([f['octants'] for f in frames])
    volume=float(np.mean([f['buried_volume'] for f in frames]))
    values=dict(file=path,buried_volume=volume,percent_buried_volume=100*volume/(4*np.pi*3.5**3/3),
        qvbur_min=float(q.min()),qvbur_max=float(q.max()),max_delta_qvbur=float(abs(q-np.roll(q,1,axis=1)).max()),
        ovbur_min=float(o.min()),ovbur_max=float(o.max()),
        near_vbur=float(np.mean([f['near_vbur'] for f in frames])),
        far_vbur=float(np.mean([f['far_vbur'] for f in frames])))
    return (frames,values) if details else values

def calculate(pair):
    tool,path=pair
    return mb.calculate(path) if tool=='morfeus' else dbstep(path)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tool',choices=['morfeus','dbstep'],required=True)
    parser.add_argument('--workers',type=int,default=1)
    parser.add_argument('--warm-iterations',type=int,default=0)
    parser.add_argument('paths',nargs='+')
    a=parser.parse_args();jobs=[(a.tool,p) for p in a.paths]
    pool=None if a.workers==1 else multiprocessing.get_context('spawn').Pool(a.workers)
    def run():return list(map(calculate,jobs)) if pool is None else pool.map(calculate,jobs)
    try:
        rows=run();samples=[]
        for _ in range(a.warm_iterations):
            start=time.perf_counter_ns();observed=run();samples.append(time.perf_counter_ns()-start)
            assert rows==observed
        json.dump({'rows':rows,'warm_elapsed_ns':samples},__import__('sys').stdout,allow_nan=False)
        print()
    finally:
        if pool is not None:pool.close();pool.join()

if __name__=='__main__':main()
