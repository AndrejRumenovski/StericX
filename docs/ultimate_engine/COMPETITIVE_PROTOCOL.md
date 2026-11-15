# Prospective competitive comparison scope

The only currently admitted whole-molecule cross-tool speed comparison is the
three-frame, radius-3.5 Å, density-0.01 finite-lattice volume estimator on the
56 frozen Morfeus comparison conformers. Reapply the existing exact-decimal
interval predicate and exact rational IEEE reduction reference to the final
source. Every grid/occupancy/region decision must agree; every tool's output
must agree with its own independently rounded arithmetic. No residual tolerance
is fitted. Any unexplained failure stops that timing lane.

Morfeus uses its unchanged public BuriedVolume/octant APIs and the historical
explicit-radii/frame adapter. DBSTEP uses its unchanged `sterics.occupied`
kernel, with the same Morfeus projection grid supplied explicitly, followed by
the same declared regional reductions. This is a matched finite-grid DBSTEP
kernel workflow, not DBSTEP's default grid/volume API. DBSTEP's internal
`workers=-1` is preserved; process affinity limits available CPUs and any
oversubscription must remain visible. Native StericX uses its public volume API.
All workflows include XYZ parsing, center/frame inference, three calculations
and JSON serialization in end-to-end measurements. Repeated calls may use the
documented StericX bounded grid cache; no molecular result cache is introduced.

Preflight checks all 56 geometries and deterministic order at 1/2/4/6/12 workers.
Single, 56-conformer and repeated large-batch measurements use identical input
bytes, one warmup and eight alternating method orders, preserving all attempts.
Cold process startup includes imports; warm worker measurements are separately
labelled and exclude parsing/output only when every method does so. Process
tree RSS must not be inferred by summing wait4 high-water marks.

Sterimol's finite-angle grids and L conventions differ. The new high-precision
and explicit-angle comparisons establish scientific residuals but do not admit
an unrestricted combined-descriptor speed ratio. Report Sterimol costs against
the independently assessed accuracy separately. No speed ranking is inferred
from density, angle count, or unmatched default APIs. This campaign does not
claim fastest-available status.
