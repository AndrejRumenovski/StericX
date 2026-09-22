# C3 implementation admission against accepted C2

This is a **pre-implementation hypothesis**, not acceptance or a native speed
claim. The original [row-rejection proposal](C3_ROW_REJECTION_PROPOSAL.md) remains
unchanged. The present predecessor is accepted C2, build manifest SHA-256
`82b63f16c90e0767c5cb5b60daed61e2ce52f123d4f1a9f3156e791c5a03fd81`,
native executable `b233501640e4d06555a36a278581fbf9cb011962311f205c81f3fc9a83ae89e6`.
Its acceptance is `accepted_c2_v1/acceptance.json` (SHA-256
`071f29b5e52cc4940b9581eb3de90430b3c1d7d4d2d1c4e7d4c9d30750f8cabe`).

The fresh all-40 C2 profile is `accepted_c2_reprofile_v1/manifest.json`, SHA-256
`52e7d287332afdb1edf13348cb73e7bccc243e2c0259d187c5705dc675426187`.
It contains 40 warmups and 200 measured diagnostic launches, with the eight
contemporaneous full-matrix native samples per configuration as timing context.
Serial `descriptors_10000` occupancy consumes **82.4636% of process wall time**.
At 2/4/6 threads, occupancy accounts for 83.3020/83.0824/82.6858% of summed worker
elapsed time; those are not additive process-wall shares. The serial native
predecessor median is 4,909.577 ms and six-thread median is 903.662 ms.

At the measured serial share, eliminating all occupancy work has an Amdahl ceiling
of approximately **5.70×**. The unchanged-source shadow measured a potential net
reduction of 883,481,917 Z predicates in the fixed 10,000-file workload, or 45.4505%
of Z predicates and 33.6009% of combined XY/Z comparisons. It includes the added
86,487,808 endpoint predicates. These are operation counts, not native timings or
equal-cost instructions. The scientific grid still contains 15,408 points in 740
contiguous rows per frame, with the original order retained.

The realistic hypothesis remains **15–25% less serial descriptor wall time**,
approximately 1.18–1.33×. This requires about 18.2–30.3% occupancy improvement and
discounts the comparison-count model for row scans, finite checks, branches and
cache bookkeeping. The candidate may be slower. The predeclared acceptance rule
remains reproducible **at least 10% lower target batch wall time** against actual
C2, beyond paired variation, all scientific gates passing, and no unacceptable
other-workload or memory regression. A 10% whole-process reduction requires about
12.13% occupancy improvement at this measured share.

Implementation will scan each actual contiguous XY-bit row for its minimum and
maximum Z and complete finite status. A private helper may reject a lazily reached
atom only when all required inputs are finite, its Z lies strictly outside that
range, and the nearest actual endpoint misses the original f32 predicate
`xy_squared + dz * dz <= radius_squared`. Rounded subtraction magnitude, square
and addition are monotone on this one-sided interval; a strict endpoint miss
therefore proves every original row-point predicate misses. Equality, interior
atoms and nonfinite cases retain the original path. No square-root/inverted
predicate, precision change, reassociation, point/atom sorting or eager atom
preparation is permitted. Original reductions and all 15 occupied-volume fields
remain unchanged. Scalar row bounds use the stack; the existing candidate vector
is reused, with no new heap allocation required.

Before editing, freeze all current source/test fixtures, Cargo/Python files,
current diff, accepted source/binaries, profile and this addendum under
`candidate_c3_row_rejection/`. Focused tests must force endpoint rejection and
retention, one-ULP boundaries, absorbed addition, signed zero, underflow, overflow,
NaN/Inf fallback, nonfinite later row points, unsorted/reversed/interrupted rows,
lazy prefixes and dense first hits. Use the unchanged direct glam reference for
all 15 fields; that differential oracle supplements the full independent corpus.

Run complete exact observations, independent references, errors/invariance,
prediction bytes and engineering gates **before native timing**. Then run the
mandatory original-baseline all-40 matrix and a predeclared incremental C2→C3
comparison: `descriptors_10000`, `conformers_56`, `ensemble_sdf`,
`db_build_ensembles`, `screen_1000`, `parse_ensembles_1000`, `search_database`, each
at one and six threads, eight alternating measured pairs per configuration with
one warmup per side (252 launches). Keep every sample. If the candidate fails a
gate, retain the patch/evidence and revert its production changes. If accepted,
freeze it, reprofile, and perform the final complete scientific recheck.
