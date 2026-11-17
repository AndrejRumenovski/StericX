"""Matched three-radius external workflows; DBSTEP occupancy kernel only."""
import argparse,importlib.util,json,multiprocessing,os,sys,time
from pathlib import Path
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'docs/ultimate_engine'))
import competitive_worker as old
import numpy as np
from morfeus import BuriedVolume
RADII=[2.5,3.5,4.5]
def calculate(job,details=False):
    tool,path=job;elements,xyz,donor,neighbors,center,radii=old.mb.geometry(path);keep=np.array([e!='H' for e in elements]);profiles=[];all_frames=[]
    for radius in RADII:
        frames=[];density=.01*(radius/3.5)**3
        for plane in neighbors:
            z=center-xyz[donor];z/=np.linalg.norm(z);v=xyz[plane]-center;x=v-np.dot(v,z)*z;x/=np.linalg.norm(x);y=np.cross(z,x);aligned=(xyz-center)@np.array([x,y,z]).T
            if tool=='morfeus':
                b=BuriedVolume(['H',*elements],np.vstack([np.zeros(3),aligned]),1,radii=np.r_[0.,radii],radius=radius,density=density,include_hs=False);b.octant_analysis()
                points=b._sphere.points;occupied=b._buried_points;volume=b._sphere.volume
                octs=[float(b.octants['buried_volume'][i]) for i in [0,1,2,3,7,6,5,4]];vbur=float(b.buried_volume)
            else:
                s=old.Sphere(np.zeros(3),radius,method='projection',density=density,filled=True);points=s.points;volume=s.volume
                occupied,_,_=old.sterics.occupied(points,aligned[keep],radii[keep],np.zeros(3),old.SimpleNamespace(grid=1.,verbose=False,debug=False))
                def counts(p):
                    x,y,z=p.T;q=np.where(y>=0,np.where(x>=0,0,1),np.where(x<0,2,3));return np.bincount(q+np.where(z>=0,0,4),minlength=8)
                octs=[float(n)/float(d)*volume/8 for n,d in zip(counts(occupied),counts(points),strict=True)];vbur=len(occupied)/len(points)*volume
            f={'plane':plane,'buried_volume':vbur,'octants':octs,'quadrants':[octs[i]+octs[i+4] for i in range(4)],'near_vbur':sum(octs[4:]),'far_vbur':sum(octs[:4])}
            if details:f.update(points=points,occupied=occupied)
            frames.append(f)
        q=np.array([f['quadrants'] for f in frames]);o=np.array([f['octants'] for f in frames]);v=float(np.mean([f['buried_volume'] for f in frames]))
        values={'buried_volume':v,'percent_buried_volume':100*v/(4*np.pi*radius**3/3),'qvbur_min':float(q.min()),'qvbur_max':float(q.max()),'max_delta_qvbur':float(abs(q-np.roll(q,1,axis=1)).max()),'ovbur_min':float(o.min()),'ovbur_max':float(o.max()),'near_vbur':float(np.mean([f['near_vbur'] for f in frames])),'far_vbur':float(np.mean([f['far_vbur'] for f in frames]))}
        profiles.append({'radius':radius,'values':values});all_frames.append(frames)
    row={'file':path,'profile':profiles};return (row,all_frames) if details else row
def main():
    p=argparse.ArgumentParser();p.add_argument('--tool',required=True,choices=['morfeus','dbstep']);p.add_argument('--workers',type=int,default=1);p.add_argument('--warm-iterations',type=int,default=0);p.add_argument('paths',nargs='+');a=p.parse_args();jobs=[(a.tool,path) for path in a.paths]
    pool=None if a.workers==1 else multiprocessing.get_context('spawn').Pool(a.workers)
    def run():return list(map(calculate,jobs)) if pool is None else pool.map(calculate,jobs)
    try:
        rows=run();times=[]
        for _ in range(a.warm_iterations):start=time.perf_counter_ns();r=run();times.append(time.perf_counter_ns()-start);assert r==rows
        print(json.dumps({'rows':rows,'warm_elapsed_ns':times}))
    finally:
        if pool:pool.close();pool.join()
if __name__=='__main__':main()
