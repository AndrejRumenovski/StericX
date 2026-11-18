# Occupancy and precision experiments

All timing in this file is exploratory unless a final controlled receipt is
explicitly cited. No unfavorable sample or failed build was removed.

## Accepted finite-grid changes

O1 (`29bc3313`) replaces scalar tails of sufficiently long rows with an
overlapping final eight-lane AVX2 load. Repeated decisions are combined by OR.
This can **increase** distance evaluations while reducing instructions; it is
not a reduction in required point–atom comparisons. Rows shorter than eight
retain scalar evaluation. Operation grouping, inclusive comparison and original
f32 inputs are preserved; no FMA occurs in the occupancy assembly.

O3 (`34e32d33`) prepares immutable row boundaries, contiguous Z coordinates,
sign masks and regional denominator counts. The original grid points remain
unchanged. Total, quadrant, octant and hemisphere counts use the same occupancy
mask. It removes repeated row traversal and mask construction rather than
changing the sphere predicate. Direct tests reconstruct every original point
and compare every occupied bit, including long rows spanning multiple masks.

O4 (`d732d1cb`) retains up to eight exact radius/density keys in a per-thread
LRU. Retained point, Z and row buffers share an 8 MiB cap. Atomic geometry and
occupied masks are never cached. Adjacent floating-point keys, eviction,
oversized grids and thread independence are tested. This benefits repeated
requests and radius profiles without assuming that conformers share geometry.

## Rejected or restricted spatial strategies

O2 evaluated X-slab atom filtering and squared-X reuse. All focused scientific
tests passed, but eight alternating pairs showed volume and regional kernels
about 12% slower. Its source, patch and complete samples remain in
`.stericx/occupancy_precision/o2_slab`, `o2.patch` and `o2_ab`.

O5 implements research versions of conservative 3D AABB BVH pruning,
uniform XY cell lists, bit-packed row unions and two/four/six-thread row
evaluation. The BVH stores center bounds and a maximum squared radius; its
nearest-box distance is a conservative lower bound under the unchanged
nonnegative f32 distance operations. Cell lists bound actual row coordinates,
then test the original atom predicate. Fixed 8×8 bins avoid hashing overhead
on this bounded integration domain. Leaves retain exact SIMD masks.

All six modes reproduce the brute-force mask on 101 cases: 36 molecular
orientations, six sparse synthetic sizes (1–1,024 atoms), three denser grids,
and 56 real conformers. Preparation and execution were recorded separately.
On those 56 conformers, median baseline/(cell-list preparation + execution)
is 0.807; for the BVH it is 0.749. The cell-list range is 0.660–1.119. Neither
is selected as a general molecular path.

The large sparse synthetic cases do benefit: at 256 atoms the cell-list
execution alone is about 14.8× faster, but preprocessing must also be paid.
On the denser 256-atom grid, preparation plus execution is about 7× faster.
These are narrow exploratory observations, not universal molecular claims.
Internal parallelism helps sufficiently large synthetic work but makes the
small ligands slower. It is not shipped: batch-level Rayon already provides
parallelism, and creating nested local pools would risk oversubscription.
The research pools are created outside the timed execution loop; their startup
cost is therefore an additional cost, not included in those execution ratios.

`o5_measurements.json` binds all raw samples and candidate counts. The first
research build failed on a literal and the raw-grid/prepared-grid distinction;
the corrected version and failure are both retained. The successful spatial
timing runs finished before the subsequent Clippy run; filesystem timestamps
were checked to resolve a possible-overlap concern. They still lack a quiet
host gate and remain exploratory.

Bounding spheres, atom AABBs and influence regions were considered as
conservative coarse rejection devices. The measured BVH/cell-list experiments
demonstrate both their potential on sparse large inputs and their setup cost
on the actual small/medium molecules. No atom-radius padding, approximate
predicate or unsafe geometry cache was introduced to manufacture a gain.

## Continuous method and shared geometry

C1 (`15534fb8`) uses conservative sphere candidates inherited by child XY
rectangles, interval-enclosed chord unions, and a midpoint remainder only on
provably smooth patches. Its proof and explicit Rust interface are in
[CONTINUOUS_MATHEMATICS.md](CONTINUOUS_MATHEMATICS.md).

Two source-bound ablations isolate its mechanisms. Testing every atom afresh
instead of inheriting candidates preserves every scientific bound and cell
decision but uses 1.16–3.98× as many geometric queries on the original twelve
molecules. Removing the smooth remainder leaves all twelve 0.05 Å³ targets
unmet at 262,143 cells, with error bounds of 0.215–0.627 Å³. Full C1 reaches
all twelve. These ablation runtimes were collected under validation load;
the query counts and certificate results, not those timings, support this
comparison. See `chord_ablations.json`.

The useful shared strategy is conservative candidate rejection and geometry
reuse. The arithmetic remains specialized: finite-grid semantics require
the original f32 operations, while continuous certificates require outward
f64 interval bounds for exact supplied spheres. Forcing both into the same
predicate or cache would neither preserve the same mathematical quantity nor
remove their different data-access costs. C1 candidates are local to one
integration and cannot survive a geometry change. Atom order is canonicalized
and duplicates removed without changing a union.

Analytical sphere intersections accelerate empty, contained, single-lens and
eligible two-sphere cases. Higher-order overlaps are not approximated by a
pair-only formula. Independent coaxial references integrate the exact disk
envelope and support genuine triple and higher intersections.

## Other mathematical candidates

Plane clipping was analyzed but not implemented. Around a cell center x₀,
the sphere polynomial is exactly
`F(x)=F(x₀)+2(x₀−c)·(x−x₀)+|x−x₀|²`.
Bounding the final term by [0, ρ²] gives inner and outer halfspaces. Turning
these into tight *union* volumes requires conservative polyhedral clipping and
overlap handling, plus rounding bounds for every generated vertex and volume.
It offers a plausible second-order boundary bound, but adds substantial proof
and implementation work after C1 already meets the molecular target. No speed
or certification claim is made for this unimplemented alternative.

Unverified Richardson extrapolation and agreement of successive quadratures
remain unsuitable stopping certificates. Ordinary adaptive quadrature is
useful as an estimate; C1 admits its derivative-based remainder only where
smoothness and endpoint topology are proved. The retained 3D box reference
and independent rational octree remain useful independent corroboration even
though their broad boundary-volume enclosures converge slowly.

The campaign stops general kernel tuning here because the retained changes
have substantial measured effects, the remaining tested general spatial
changes regress representative molecules, and validation/controlled timing
are more valuable than tuning differences comparable to run variation.


## Later measured continuous arithmetic/storage work

After the general occupancy experiments, allocation diagnostics exposed
0.68–2.98 million allocations per molecular 0.01 calculation. A targeted C1
Callgrind run (asymmetric phosphine, 0.05) models 812,610,878 instructions:
interval multiplication alone accounts for 144,082,719, while malloc/free
and vector-growth functions have smaller exclusive costs. This justified a
separate continuous optimization pass; these are instruction models, not cycles.

C3 specializes finite nonnegative interval multiplication/division while
preserving the general four-product bounds bit for bit. Its 100,000-case
exponent/rounding test passes. Eight alternating pairs across twelve molecules
and two tolerances give a median exploratory ratio of 1.114× (range
1.056–1.174). C4 separately reuses scratch vectors and component index ranges,
with a median 1.076× over C3 (range 1.022–1.162). Both preserve all tested
certificates and work counts. Complete raw pairs remain in `c3_ab` and `c4_ab`;
final controlled gains and allocation reductions must be measured on the
integrated implementation before final performance claims.


Final integrated diagnostics confirm median allocation-call reductions of 82.80%
(0.05) and 83.85% (0.01), while peak live memory changes only by retained scratch.
`final_diagnostics.json` contains every comparison. Final-v2 review/setup failures
are recorded in WORK_LOG.md and raw logs; no acceptance criterion was changed.
