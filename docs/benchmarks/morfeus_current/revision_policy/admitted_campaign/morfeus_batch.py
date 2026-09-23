"""Ordinary independent-molecule morfeus loop, optionally Pool.map parallelism.

All input reading, topology/center construction, three reference calculations,
aggregation and output happen inside the measured process tree. No cache.
"""
import argparse
import json
import multiprocessing
import os
import sys

for key in ["OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","BLIS_NUM_THREADS","VECLIB_MAXIMUM_THREADS","NUMEXPR_NUM_THREADS"]:
    os.environ[key]="1"

import numpy as np
from morfeus import BuriedVolume
from morfeus.io import read_xyz

# Explicit shared conventions for all elements present in the frozen corpus.
VDW={"H":1.20,"C":1.70,"N":1.55,"O":1.52,"P":1.80,"S":1.80}
COVALENT={"H":0.31,"C":0.76,"N":0.71,"O":0.66,"P":1.07,"S":1.05}


def geometry(path):
    elements,xyz=read_xyz(path)
    donors=[i for i,e in enumerate(elements) if e=="P"]
    assert len(donors)==1
    donor=donors[0]
    distances=np.sum((xyz-xyz[donor])**2,axis=1)
    neighbors=sorted([i for i,e in enumerate(elements) if i!=donor and distances[i]<=(1.3*(COVALENT[e]+COVALENT["P"]))**2],key=lambda i:(distances[i],i))
    assert len(neighbors)==3
    unit=xyz[neighbors]-xyz[donor]
    unit/=np.linalg.norm(unit,axis=1)[:,None]
    direction=-unit.sum(axis=0)
    assert np.linalg.norm(direction)>0.01
    center=xyz[donor]+2.28*direction/np.linalg.norm(direction)
    radii=np.array([VDW[e]*1.17 for e in elements])
    return elements,xyz,donor,neighbors,center,radii


def calculate(path):
    elements,xyz,donor,neighbors,center,radii=geometry(path)
    volumes=[]
    near=[]
    far=[]
    qs=[]
    os_=[]
    for plane in neighbors:
        z=center-xyz[donor]
        z/=np.linalg.norm(z)
        v=xyz[plane]-center
        x=v-np.dot(v,z)*z
        x/=np.linalg.norm(x)
        y=np.cross(z,x)
        aligned=(xyz-center)@np.array([x,y,z]).T
        b=BuriedVolume(["H",*elements],np.vstack([np.zeros(3),aligned]),1,
                       radii=np.r_[0.,radii],radius=3.5,density=0.01,include_hs=False)
        b.octant_analysis()
        o=[float(b.octants["buried_volume"][i]) for i in [0,1,2,3,7,6,5,4]]
        q=[o[i]+o[i+4] for i in range(4)]
        volumes.append(float(b.buried_volume))
        near.append(sum(o[4:]))
        far.append(sum(o[:4]))
        qs.append(q)
        os_.append(o)
    q=np.array(qs)
    o=np.array(os_)
    volume=float(np.mean(volumes))
    return dict(file=path,buried_volume=volume,percent_buried_volume=100*volume/(4*np.pi*3.5**3/3),
                qvbur_min=float(q.min()),qvbur_max=float(q.max()),max_delta_qvbur=float(abs(q-np.roll(q,1,axis=1)).max()),
                ovbur_min=float(o.min()),ovbur_max=float(o.max()),near_vbur=float(np.mean(near)),far_vbur=float(np.mean(far)))


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--workers",type=int,default=1)
    p.add_argument("paths",nargs="+")
    args=p.parse_args()
    if args.workers==1:
        rows=list(map(calculate,args.paths))
    else:
        with multiprocessing.get_context("spawn").Pool(args.workers) as pool:
            rows=pool.map(calculate,args.paths)
    json.dump(rows,sys.stdout,indent=2,allow_nan=False)
    print()


if __name__=="__main__": main()
