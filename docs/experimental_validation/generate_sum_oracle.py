"""Exact Fraction oracle for tight-add enclosure; no production floating algorithm."""
from fractions import Fraction
from pathlib import Path
import json,math,random,struct
ROOT=Path(__file__).resolve().parents[2]
def bits(x):return struct.unpack('Q',struct.pack('d',x))[0]
def number(x):return struct.unpack('d',struct.pack('Q',x))[0]
def enclosure(a,b):
 exact=Fraction(a)+Fraction(b); nearest=float(exact)
 low=math.nextafter(nearest,-math.inf) if Fraction(nearest)>exact else nearest
 high=math.nextafter(nearest,math.inf) if Fraction(nearest)<exact else nearest
 return [bits(a),bits(b),bits(low),bits(high)]
rng=random.Random(20261124)
pairs=[(0.,0.),(0.,-0.),(1.,-1.),(1.,math.ulp(0.)),(math.ulp(0.),math.ulp(0.))]
for _ in range(6000):
 a=number(rng.getrandbits(64)&0xff6fffffffffffff)
 b=number(rng.getrandbits(64)&0xff6fffffffffffff)
 pairs.append((a,b))
# Signed zero endpoint identity is not a distinct real-number bound.
rows=[enclosure(a,b) for a,b in pairs]
(ROOT/'tests/fixtures/geometric_sum_oracle.json').write_text(json.dumps(rows)+'\n')
print(len(rows))
