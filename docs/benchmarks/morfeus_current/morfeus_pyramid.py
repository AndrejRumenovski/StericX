"""Separate pyramidalization API driver on shared f32-represented XYZ inputs."""
import argparse
import json
import multiprocessing
import os
import sys
for key in ["OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","BLIS_NUM_THREADS","VECLIB_MAXIMUM_THREADS","NUMEXPR_NUM_THREADS"]:os.environ[key]="1"
import numpy as np
from morfeus import Pyramidalization
from morfeus.io import read_xyz
COVALENT={"H":0.31,"C":0.76,"N":0.71,"O":0.66,"P":1.07,"S":1.05}


def calculate(path):
    elements,xyz=read_xyz(path)
    xyz=xyz.astype(np.float32).astype(float)
    donors=[i for i,e in enumerate(elements) if e=="P"]
    assert len(donors)==1
    d=donors[0]
    ds=np.sum((xyz-xyz[d])**2,axis=1)
    ns=sorted([i for i,e in enumerate(elements) if i!=d and ds[i]<=(1.3*(COVALENT[e]+COVALENT["P"]))**2],key=lambda i:(ds[i],i))
    assert len(ns)==3
    p=Pyramidalization(xyz,d+1,neighbor_indices=[i+1 for i in ns])
    return dict(file=path,pyr_p=float(p.P),pyr_alpha=float(p.alpha))


def main():
    p=argparse.ArgumentParser();p.add_argument("--workers",type=int,default=1);p.add_argument("paths",nargs="+");args=p.parse_args()
    if args.workers==1:rows=list(map(calculate,args.paths))
    else:
        with multiprocessing.get_context("spawn").Pool(args.workers) as pool:rows=pool.map(calculate,args.paths)
    json.dump(rows,sys.stdout,indent=2,allow_nan=False);print()


if __name__=="__main__":main()
