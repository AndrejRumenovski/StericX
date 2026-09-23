# Native occupancy probe and optimized row cache

This is a source-level justification for the already frozen policy's native
predicate probe, not a new numerical acceptance rule. The probe evaluates the
production `glam::Vec3::distance_squared` predicate at every native grid point;
`occupied_volumes` additionally uses a row cache. The probe does not instrument
an optimized-mask output (the public API returns scalar volumes).

For a point and atom, let `dx`, `dy`, `dz` be their f32 component differences.
The frozen glam implementation expands squared distance as
`((dx*dx) + (dy*dy)) + (dz*dz)`, with the same left-associated operations used by
`occupied_volumes`. The recorded build uses ordinary IEEE operations, no
fast-math or reassociation flags, and default x86-64 SSE2 features. The source
identities are retained in `predicate_source_identity.json`; the exact production
source is already in the original frozen source archive.

Within a row, x and y are identical **by their bit patterns**. Therefore an atom's
cached `xy_squared` is bit-identical to recomputing the first two squared terms.
Adding a nonnegative rounded z-square is monotone: if `xy_squared > radius_squared`,
that atom cannot hit any point in the row. Aligned atoms and radii are checked for
finiteness by the unchanged production API. Overflow of a distance square to
positive infinity is still a miss against a finite squared radius.

The cache invariant is that all atoms before `next_atom` have either been retained
with their exact first sum or rejected by that safe row-wide condition. For each
point, every retained candidate is tested before the unexamined suffix is advanced.
An early hit may leave the suffix unexamined for that point; the next point tests
the existing cache and resumes the suffix if necessary. If the suffix is exhausted,
all atoms have been considered. On an x/y change, the cache and suffix index reset.
Consequently the optimized union predicate equals the direct point-by-point union
predicate used by the probe for every point. This argument uses the operation
sequence and cache invariant, not matching observed regional totals.

The retained numerical evidence separately checks optimized per-plane scalars and
public API aggregates against independently reconstructed certified counts. This
is not machine-checked verification of the compiler or all possible program inputs;
it is the source-level bridge between the direct predicate probe and this frozen
optimized implementation, under the recorded ordinary floating-point build.
