# Actual-cell error enclosures for the convergence study

This derivation was fixed before the convergence campaign. It is a conservative
diagnostic error bound, not a fitted accuracy tolerance or a claim that its
endpoints are attainable. The controller is
[`analyze_dbstep_volume_convergence.py`](../../../../scripts/analyze_dbstep_volume_convergence.py).

Let the actual three sorted coordinate axes define a tensor lattice. Give each
node its Cartesian Voronoi box: interior boundaries are exact adjacent-node
midpoints; exterior boundaries extend half the last adjacent gap. Require the
finite box union to contain the entire evaluation sphere. All these coordinates
and widths are exact rational interpretations of the stored binary numbers.
The maximum cell radius `δ` is an outward square-root bound on the sum of the
three largest squared half-spans. Actual nonuniform spacings are retained.

Native full-axis metadata is recovered using the native coordinate-generation
operation order and the observed axis count. Every captured point must match
one actual Cartesian node exactly; ambiguous reconstruction is an error. This
reconstructs coordinate metadata only. It does not implement the native atomic
occupancy predicate. DBSTEP axes follow its documented rounded-linspace
construction over an explicit common cube. Morfeus axes follow its documented
projection grid, with all returned sphere/occupied points verified against the
axes. All full node masks and actual axes are retained in compressed arrays.

The independent binary64 direct union predicate is not presumed exact. With
`u=2^-53`, its subtract-square-add graph has distance relative bound
`γ5=5u/(1-5u)` and the radius square contributes `u*r²`. Exact rational algebra
verifies that its possible radial error is at most `ε(r)=8*u*r`. Finite normal
arithmetic is checked for all supplied point/center differences and radii;
unsupported domains fail. Thus disagreement between the point rule and the
continuous ball can affect its selected-cell union only in the radial shell
between `max(0,r-δ-ε)` and `r+δ+ε`. A rational Machin enclosure of π gives an
outward bound on that shell volume.

Let `H_S` be that shell bound for the evaluation sphere and `H_i` for each
atomic ball. For the sphere-intersected union, the sufficient boundary bound is
`H_S + sum(H_i)` by the union bound. It deliberately overcounts overlapping
boundaries and parts lying outside the evaluation sphere. An empty atomic union
has zero numerator error.

Actual SUT numerator and denominator masks are separately compared with the
independent direct rules on **every full tensor node**. If `D_S` and `D_O` nodes
differ, add `D_S*w_max` and `D_O*w_max` respectively. This accounts for observed
binary32 sphere/occupancy decisions, cKDTree behavior and DBSTEP's distinct
numerator/denominator predicates without assuming their internal arithmetic.
The reference predicate's own rounding remains covered by the shell enlargement.
Every differing node index is retained. Counts alone are insufficient.

DBSTEP's fast denominator has a count-returning API. Its controller mask is
explicitly labelled reconstructed from the documented subtractive expression,
with its count checked against the unmodified routine; it is not an independently
observed pointwise mask. A count-equivalent auxiliary mask is sufficient for the
ratio bound below because the reported estimator depends on its denominator
count. It must not be used as an independent point-decision oracle. The other
actual masks come from returned native/reference points or native decisions.

All new disputed native nodes additionally retain exact rational sphere/atomic
margins and explicit binary32 operation traces. Every trace must reproduce the
actual native sphere and numerator decisions, and every difference from the
exact rational predicate must fall within the separately derived a priori
roundoff bound. Broad cell enclosures do not excuse unexplained native errors.

The engine estimates a count ratio, whereas actual cells can have unequal
volumes. Choose `v0` as the product of the three average adjacent spacings and
let `e_w=max(|w_min-v0|,|w_max-v0|)`, all evaluated rationally. With actual selected
counts `N_S,N_O`, valid quadrature bounds are

```
E_S = H_S + D_S*w_max + N_S*e_w
E_O = H_S + sum(H_i) + D_O*w_max + N_O*e_w
```

For positive `Q_S=N_S*v0` and `Q_O=N_O*v0`, if `V_S` is the exact sphere volume
and `V_O` the exact continuous intersection volume, then

```
|V_S*Q_O/Q_S - V_O| <= E_O + (Q_O/Q_S)*E_S.
```

This follows by adding/subtracting `Q_O` and applying the two quadrature bounds;
it does not assume numerator nodes are a subset of denominator nodes. Finally,
bound the actual reported value's difference from exact `V_S*N_O/N_S` using the
rational count ratio and the independent π enclosure. This explicitly includes
native normalization rounding, reference floating arithmetic and DBSTEP's
percentage clamp. Add that difference to the discretization bound.

These bounds apply to the observed estimator at its actual resolution. They do
not assert that all successively finer numerical errors decrease monotonically.
Analytical lens/union volumes and independently adaptive box intervals provide
separate accuracy information. Broad shell bounds are reported as broad;
matching a broad bound does not justify a fine-accuracy claim.
