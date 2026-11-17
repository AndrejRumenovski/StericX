# Prospective candidate ledger

No production implementation has changed at creation of this ledger. Each
accepted optimization requires exact finite-grid regression, isolated timing,
and a separate commit. Failed candidates and timing attempts stay on disk.

| ID | Hypothesis | Required evidence | Current decision |
|---|---|---|---|
| O1 | Overlapping final AVX2 vector removes scalar tails without changing a predicate | 23 focused tests including 5,200 boundary masks and all 56 conformers; eight alternating kernel pairs | Retain; final controlled acceptance pending |
| O2 | Exact X-slab filtering and squared-X reuse reduce repeated row–atom work | All masks and output bits; small/dense/large overhead; counter and timing reduction | Pending |
| O3 | Cached row/Z/regional metadata reduces point traversal and mask preparation | Exact grid keys/nodes; bounded cache memory; cold/warm/profile A/B | Pending |
| O4 | Cell lists/BVH/bitset row spans/internal parallelism help only larger work | Preprocessing separately timed; mask equality; avoid nested pools | Investigation pending |
| C1 | Certified chord ranges plus smooth-patch Taylor integration tighten boxes | Full proof of range/topology/derivative bounds; rational primitives; independent enclosures | Investigation pending |
| C2 | Inner/outer plane bounds give certified partial sphere–cell volumes | Outward-rounded clipping proof; high-precision and rational adversaries | Investigation pending |
| G1 | Reuse conservative sphere candidates between related cells/rows | Correct precision-specific predicates and invalidation; shared versus specialized A/B | Investigation pending |

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
