#!/usr/bin/env python3
"""Independent high-precision interval search of a transverse support envelope.

This searches angular intervals, not the Rust stationary/pair-event candidate
set. Mathematical lower bounds are evaluated with 80 decimal digits; these are
validation brackets, not a formal interval-arithmetic floating-point certificate.
"""
import heapq
import itertools
import mpmath as mp


def minimum_support(discs, relative_width=1e-10, max_intervals=200000):
    with mp.workdps(80):
        discs = [(mp.mpf(x), mp.mpf(y), mp.mpf(r)) for x,y,r in discs]
        scale = max([abs(v) for d in discs for v in d] + [mp.mpf('1e-300')])
        requested = mp.mpf(relative_width)*scale
        tau = 2*mp.pi
        data = [(x,y,r,mp.sqrt(x*x+y*y),(mp.atan2(y,x)+mp.pi)%tau)
                for x,y,r in discs]
        def support(theta):
            c,s = mp.cos(theta),mp.sin(theta)
            return max([mp.mpf(0)] + [x*c+y*s+r for x,y,r in discs])
        def lower(a,b):
            ca,sa,cb,sb = mp.cos(a),mp.sin(a),mp.cos(b),mp.sin(b)
            bounds = [mp.mpf(0)]
            for x,y,r,rho,stationary in data:
                bounds.append(r-rho if a <= stationary <= b
                              else min(x*ca+y*sa+r,x*cb+y*sb+r))
            # Roundoff margin far larger than 80-digit evaluation error at the
            # tested scales; no claim of a formally proved libmp error bound.
            return max(bounds)-scale*mp.mpf('1e-65')
        sequence = itertools.count()
        heap = []
        best = support(0)
        for i in range(16):
            a,b = tau*i/16,tau*(i+1)/16
            best = min(best,support((a+b)/2),support(b))
            heapq.heappush(heap,(lower(a,b),next(sequence),a,b))
        evaluations = 16
        while best-heap[0][0] > requested and evaluations < max_intervals:
            _,_,a,b = heapq.heappop(heap)
            middle = (a+b)/2
            for lo,hi in [(a,middle),(middle,b)]:
                best = min(best,support((lo+hi)/2))
                heapq.heappush(heap,(lower(lo,hi),next(sequence),lo,hi))
                evaluations += 1
        low = max(mp.mpf(0),heap[0][0])
        return {'lower_decimal':mp.nstr(low,70), 'upper_decimal':mp.nstr(best,70),
                'lower':float(low), 'upper':float(best), 'intervals':evaluations,
                'converged':best-low <= requested, 'requested_width':float(requested)}


if __name__ == '__main__':
    import json,sys
    for line in sys.stdin:
        value=json.loads(line)
        print(json.dumps({'id':value['id'],**minimum_support(value['discs'])}),flush=True)
