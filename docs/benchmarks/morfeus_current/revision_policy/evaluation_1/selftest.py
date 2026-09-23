"""Reference-only arithmetic checks; no molecular/candidate results used."""
from decimal import Decimal, localcontext
from fractions import Fraction as Q
import json, random, struct
import reference as r
import numpy as np
rng=random.Random(220926)
for _ in range(2000):
 a=Q(rng.randint(-10000000,10000000),rng.randint(1,10000));b=Q(rng.randint(1,10000000),rng.randint(1,10000))
 ia,ib=r.exact(a),r.exact(b)
 for got,want in [(ia+ib,a+b),(ia-ib,a-b),(ia*ib,a*b),(ia/ib,a/b),(ia.sq(),a*a)]:
  assert Q(float(got.lo))<=want<=Q(float(got.hi))
 for p in [24,53]:
  got=r.ieee(a,p)
  expected=float(a) if p==53 else struct.unpack('f',struct.pack('f',float(a)))[0]
  assert got==Q(expected)
 with localcontext() as ctx:
  ctx.prec=100
  root=(Decimal(b.numerator)/Decimal(b.denominator)).sqrt();got=ib.sqrt()
  assert Decimal(float(got.lo))<=root<=Decimal(float(got.hi))
# Halfway ties are exact rational inputs, including those lost through double rounding.
for p in [24,53]:
 step=Q(1,2**(p-1))
 assert r.ieee(1+step/2,p)==1
 assert r.ieee(1+3*step/2,p)==1+2*step
 assert r.ieee(1+step/2+Q(1,2**200),p)==1+step
assert r.ieee(r.PI_BOUNDS[0],53)==Q(float.fromhex('0x1.921fb54442d18p+1'))
assert r.ieee(r.PI_BOUNDS[0],24)==r.ieee(r.PI_BOUNDS[1],24)
assert np.all(r.GRID.lo<=r.GRID.hi)
# Exact lattice membership and sign mapping need no geometric floating thresholds.
assert np.all(np.sum((2*r.LATTICE-31)**2,axis=1)<=31**2)
print(json.dumps({'passed':True,'rational_interval_cases':2000,'sqrt_reference_digits':100,'tie_rounding':'exact rational','candidate_results_read':False,'lattice_points':len(r.IDS)}))
