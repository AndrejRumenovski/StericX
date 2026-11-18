# Prospective candidate ledger

No production implementation has changed at creation of this ledger. Each
accepted optimization requires exact finite-grid regression, isolated timing,
and a separate commit. Failed candidates and timing attempts stay on disk.

| ID | Hypothesis | Required evidence | Current decision |
|---|---|---|---|
| O1 | Overlapping final AVX2 vector removes scalar tails without changing a predicate | 23 focused tests including 5,200 boundary masks and all 56 conformers; eight alternating kernel pairs | Retain; final controlled acceptance pending |
| O2 | Exact X-slab filtering and squared-X reuse reduce repeated row–atom work | 23 focused tests; eight alternating equal-output pairs | Rejected general path: volume/regions about 12% slower |
| O3 | Cached row/Z/regional metadata reduces point traversal and mask preparation | 24 focused tests including full prepared-grid point masks; bounded 8 MiB cache; eight alternating pairs | Retain; final controlled acceptance pending |
| O4 | Bounded multi-key grid LRU avoids repeated preparation in radius profiles | 25 focused tests; eight alternating equal-output pairs | Retain: exploratory radius profile 1.893× over O3 |
| O5 | Cell lists/BVH/bitset row spans/internal parallelism help only larger work | 101 exact mask cases; preparation/execution separated | Reject general path; retain large sparse synthetic evidence |
| C1 | Certified chord ranges plus smooth-patch Taylor integration tighten boxes | Proof, 6,000 rational predicates, 167 original and 96 additional analytical observations, molecular nesting | Retain as separate explicit API; full regression/final cost pending |
| C3 | Specialize nonnegative interval arithmetic | 100,000 bit-equivalence cases and full 695-certificate replay | Retain; unchanged bounds, final controlled timing pending |
| C4 | Reuse per-integration scratch vectors | Full 695-certificate replay, separate allocation diagnostics | Retain; about 83% fewer allocations, peak live memory unchanged |
| C2 | Inner/outer plane bounds give certified partial sphere–cell volumes | Mathematical analysis in EXPERIMENTS.md | Not implemented; no performance/certificate claim |
| G1 | Reuse conservative sphere candidates between related cells/rows | Source-bound full/all-atoms/range-only ablations | Inherited candidates save 1.16–3.98× queries; specialize arithmetic |

Prior rejected endpoint searches, unproved Richardson stopping and pairwise-only
many-sphere union approximations remain rejected historical approaches. A new
experiment must identify what differs before revisiting one.

Practical research schedule: prioritize effects large enough to matter on the
measured occupancy hotspot, use bounded work ladders for continuous methods,
and stop tuning differences comparable to run variation. The 0.05 Å³ molecular
target stays fixed even when a budget fails.

O1 exploratory paired results: volume 1.115×, regions 1.135×, six-conformer
ensemble 1.086× and radius profile 1.092×. The unchanged Sterimol control is
0.985× (slower); every sample remains in `o1_measurements.json` and its linked
raw run. All 80 launches have identical scientific outputs. These desktop
measurements support retaining the candidate for full final replay, not a final
controlled speed claim. Overlapping SIMD lanes evaluate the exact same f32
operation sequence and are combined by idempotent OR. Short rows retain the
scalar fallback; no out-of-slice load or FMA is introduced.

O2 source/patch and all measurements are preserved under `o2_slab`, `o2.patch`
and `o2_ab`. Filtering every X slab added work and changed the inner-loop layout;
the smaller number of distance operations did not produce lower runtime on the
declared molecular workload. It is not shipped as a general optimization.

O3 caches immutable row boundaries, contiguous Z values, positive-Z masks and
regional denominator counts. Grid construction and every original f32 point
remain unchanged. The cache counts retained point, Z and metadata allocations
against its 8 MiB cap. It caches no atoms or occupancy and therefore cannot
silently reuse incompatible molecular geometry. The scalar fallback retains
the original point path.

Against O1, O3 exploratory median ratios are volume 1.791×, regions 1.758×,
ensemble 1.528×, radius profile 1.122× and unchanged Sterimol 1.000×. All 80
timed launches preserve exact scientific outputs. New direct checks verify every
prepared-grid point, sign bit and union decision, including rows longer than
64 and inclusive boundaries. Preparation/memory tradeoffs and controlled final
ratios remain to be reported; timing only warmed volume calls would hide the
extra metadata construction and retention cost.

O4 keeps up to eight immutable grids per thread under the same total 8 MiB
buffer cap. Exact binary radius/density keys, eviction, thread independence and
adjacent-float keys are tested. All 80 trials preserve outputs. The main gain is
radius-profile reuse (1.893× over O3); other modes span 1.008–1.025× and are
not claimed as substantive gains. Failed controller invocations (unavailable
`python`, then wrong output-relative path) are retained alongside the successful
run. Molecular geometry and occupancy are never cached.
