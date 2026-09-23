"""Independent exact-decimal interval reference and rational IEEE reductions.

No candidate residuals, outputs or observed maxima are inputs to this module.
"""
from fractions import Fraction as Q
import math
import numpy as np


def down(a): return np.nextafter(a, -np.inf)
def up(a): return np.nextafter(a, np.inf)


class I:
    def __init__(self, lo, hi=None):
        self.lo = np.asarray(lo, dtype=float)
        self.hi = np.asarray(lo if hi is None else hi, dtype=float)
        assert np.all(self.lo <= self.hi)
    def __getitem__(self, k): return I(self.lo[k], self.hi[k])
    def __neg__(self): return I(-self.hi, -self.lo)
    def __add__(self, b):
        if not isinstance(b,I): b=I(b)
        return I(down(self.lo+b.lo), up(self.hi+b.hi))
    __radd__=__add__
    def __sub__(self,b):
        if not isinstance(b,I): b=I(b)
        return self+-b
    def __mul__(self,b):
        if not isinstance(b,I): b=I(b)
        p=np.array([self.lo*b.lo,self.lo*b.hi,self.hi*b.lo,self.hi*b.hi])
        return I(down(p.min(axis=0)),up(p.max(axis=0)))
    __rmul__=__mul__
    def __truediv__(self,b):
        if not isinstance(b,I): b=I(b)
        assert np.all((b.lo>0)|(b.hi<0)), 'Unresolved interval division'
        return self*I(down(1/b.hi),up(1/b.lo))
    def sq(self):
        return I(np.where((self.lo<=0)&(self.hi>=0),0,down(np.minimum(self.lo**2,self.hi**2))),up(np.maximum(self.lo**2,self.hi**2)))
    def sum(self,axis=-1):
        # Explicit sequential additions: no unbounded vector reduction.
        lo=np.moveaxis(self.lo,axis,0);hi=np.moveaxis(self.hi,axis,0)
        out=I(np.zeros(lo.shape[1:]))
        for a,b in zip(lo,hi,strict=True): out=out+I(a,b)
        return out
    def sqrt(self):
        assert np.all(self.lo>=0)
        return I(np.maximum(0,down(np.sqrt(self.lo))),up(np.sqrt(self.hi)))
    def unit(self): return self/self.sq().sum().sqrt()


def exact(q):
    q=Q(q);f=float(q)
    if Q(f)==q:return I(f)
    return I(np.nextafter(f,-np.inf) if Q(f)>q else f,np.nextafter(f,np.inf) if Q(f)<q else f)


def array_exact(rows):
    rows=np.asarray(rows,dtype=object)
    cells=[exact(q) for q in rows.flat]
    return I(np.array([x.lo for x in cells]).reshape(rows.shape),np.array([x.hi for x in cells]).reshape(rows.shape))


def cross(a,b):
    c=[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
    return I([x.lo for x in c],[x.hi for x in c])


VDW={k:Q(v) for k,v in dict(H='1.20',C='1.70',N='1.55',O='1.52',P='1.80',S='1.80').items()}
COV={k:Q(v) for k,v in dict(H='.31',C='.76',N='.71',O='.66',P='1.07',S='1.05').items()}
LATTICE=np.array([(i,j,k) for i in range(32) for j in range(32) for k in range(32) if (2*i-31)**2+(2*j-31)**2+(2*k-31)**2<=31**2],dtype=np.int64)
IDS=LATTICE[:,0]*1024+LATTICE[:,1]*32+LATTICE[:,2]
GRID=array_exact([[Q(7*(2*int(t)-31),62) for t in row] for row in LATTICE])
REGIONS=np.array([(0 if i>=16 and j>=16 else 1 if i<16 and j>=16 else 2 if i<16 and j<16 else 3)+(0 if k>=16 else 4) for i,j,k in LATTICE])
assert Q('31.5')**3 < 8*Q('3.5')**3/Q('.01') < Q('32.5')**3
assert len(set(IDS.tolist()))==len(IDS) and np.all((2*LATTICE-31)!=0)


def geometry(path):
    lines=path.read_text().splitlines(); n=int(lines[0]);rows=[s.split() for s in lines[2:2+n]]
    assert len(rows)==n
    elements=[r[0] for r in rows];qxyz=[[Q(x) for x in r[1:4]] for r in rows]
    donors=[i for i,e in enumerate(elements) if e=='P'];assert len(donors)==1;d=donors[0]
    ds=[sum((a-b)**2 for a,b in zip(v,qxyz[d],strict=True)) for v in qxyz]
    ns=sorted([i for i,e in enumerate(elements) if i!=d and ds[i]<=(Q('1.3')*(COV[e]+COV['P']))**2],key=lambda i:(ds[i],i));assert len(ns)==3
    xyz=array_exact(qxyz)
    units=[(xyz[i]-xyz[d]).unit() for i in ns]
    direction=-(units[0]+units[1]+units[2]); assert direction.sq().sum().lo>Q('0.0001')
    center=xyz[d]+exact('2.28')*direction.unit()
    return elements,xyz,d,ns,center


def certify(path):
    elements,xyz,d,ns,center=geometry(path)
    records=[]
    for plane in ns:
        z=(center-xyz[d]).unit();v=xyz[plane]-center
        x=(v-z*(v*z).sum()).unit();y=cross(z,x)
        positions=[];radii=[]
        for i,e in enumerate(elements):
            if e=='H':continue
            r=xyz[i]-center;c=[(r*b).sum() for b in [x,y,z]]
            positions.append(I([a.lo for a in c],[a.hi for a in c]));radii.append(exact(VDW[e]*Q('1.17')).sq())
        # Chunk to keep the independently evaluated full union bounded in RAM.
        lo=np.full(len(IDS),np.inf);hi=np.full(len(IDS),np.inf)
        for p,r in zip(positions,radii,strict=True):
            delta=(GRID-p).sq().sum()-r
            lo=np.minimum(lo,delta.lo);hi=np.minimum(hi,delta.hi)
        occupied=hi<=0;free=lo>0;unknown=~(occupied|free)
        records.append(dict(plane=plane,occupied=occupied,unknown=unknown,clearance_lo=lo,clearance_hi=hi,
                            occupied_counts=[int(np.count_nonzero(occupied&(REGIONS==r))) for r in range(8)],
                            grid_populations=[int(np.count_nonzero(REGIONS==r)) for r in range(8)]))
    return dict(donor=d,neighbors=ns,frames=records)


def ieee(q,p):
    """Correctly round a rational to binary p significant bits, ties to even.

    All reference scalar operations here are normal or exactly zero.
    """
    q=Q(q)
    if not q:return q
    sign=-1 if q<0 else 1;q=abs(q)
    e=q.numerator.bit_length()-q.denominator.bit_length()
    power=lambda n: Q(2**n) if n>=0 else Q(1,2**(-n))
    if q<power(e):e-=1
    assert -126<=e<=127 if p==24 else -1022<=e<=1023
    step=power(e-p+1);t=q/step;n,r=divmod(t.numerator,t.denominator)
    if 2*r>t.denominator or (2*r==t.denominator and n%2):n+=1
    return sign*n*step


def atan_bounds(n,terms=80):
    s=sum((Q(-1 if k%2 else 1,(2*k+1)*n**(2*k+1)) for k in range(terms)),Q(0))
    t=Q(-1 if terms%2 else 1,(2*terms+1)*n**(2*terms+1))
    return min(s,s+t),max(s,s+t)


a,b=atan_bounds(5);c,d=atan_bounds(239);PI_BOUNDS=(16*a-4*d,16*b-4*c)


def reduce_counts(frames,p):
    rnd=lambda x:ieee(x,p)
    pi=rnd(PI_BOUNDS[0]);assert pi==rnd(PI_BOUNDS[1])
    # R^3 and 4*R^3 are exactly binary representable here.
    r3=Q('3.5')**3
    volume=rnd(rnd(rnd(4*pi)*r3)/3) if p==24 else rnd(rnd(4*r3*pi)/3)
    def total(xs):
        t=Q(0)
        for x in xs:t=rnd(t+x)
        return t
    def mean(xs):return total(rnd(x/3) for x in sorted(xs)) if p==24 else rnd(total(xs)/3)
    out=[]
    for f in frames:
        ns=f['occupied_counts'];ds=f['grid_populations']
        frac=lambda n,d:rnd(Q(n,d))
        octs=[rnd(rnd(frac(n,d)*volume)/8) if p==24 else rnd(frac(n,d)*rnd(volume/8)) for n,d in zip(ns,ds,strict=True)]
        qs=[rnd(rnd(frac(ns[i]+ns[i+4],ds[i]+ds[i+4])*volume)/4) if p==24 else rnd(octs[i]+octs[i+4]) for i in range(4)]
        out.append(dict(buried_volume=rnd(frac(sum(ns),sum(ds))*volume),near_vbur=total(octs[4:]),far_vbur=total(octs[:4]),octants=octs,quadrants=qs))
    vals={k:mean([f[k] for f in out]) for k in ['buried_volume','near_vbur','far_vbur']}
    vals.update(qvbur_min=min(x for f in out for x in f['quadrants']),qvbur_max=max(x for f in out for x in f['quadrants']),
                ovbur_min=min(x for f in out for x in f['octants']),ovbur_max=max(x for f in out for x in f['octants']),
                max_delta_qvbur=max(abs(rnd(f['quadrants'][i]-f['quadrants'][(i+3)%4])) for f in out for i in range(4)))
    # Timed morfeus driver expression is 100*V/(4*pi*R**3/3).
    percent_denominator=rnd(rnd(rnd(4*pi)*r3)/3)
    vals['percent_buried_volume']=rnd(100*rnd(vals['buried_volume']/volume)) if p==24 else rnd(rnd(100*vals['buried_volume'])/percent_denominator)
    return [{k:[float(v) for v in x] if isinstance(x,list) else float(x) for k,x in f.items()} for f in out],{k:float(v) for k,v in vals.items()},float(volume)
