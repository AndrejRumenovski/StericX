# Corrected-baseline row-endpoint shadow hypothesis (pre-edit)

Scope: operation counting only, no production optimization and no native timing claim.
Copy accurate_baseline_v1/source into a new private tree. Retain its whole buried-volume
source as an exact prefix, append an instrumented copy of occupied_volumes and call it
only as additional diagnostics after the original per-conformer calculation. Replay
using the original explicit neighbors chosen by the corrected descriptor pipeline;
no old donor inference is substituted. Actual stdout/stderr come from original code.

For each contiguous row defined by identical x/y bit patterns, scan all actual z
coordinates to obtain min/max; do not assume monotonic or regular point order. For
finite row points, finite atom coordinates, finite r² and finite q=(dx*dx+dy*dy),
an atom strictly outside [minz,maxz] has a nearest actual endpoint. Test that endpoint
with the exact original f32 grouping q + dz*dz <= r². A strict miss can eliminate
this atom for the row. Retain equality, inside-range atoms and all nonfinite cases.
With round-to-nearest IEEE arithmetic, rounded subtraction magnitude on one side,
nonnegative square and addition are monotone. Thus each actual point's rounded
sum is >= the nearest endpoint's rounded sum. Finite subtraction/square overflow
is +infinity and preserves this order; underflow and absorption preserve monotonicity.
No f64 rewrite, FMA, algebraic r²-q inversion or sqrt boundary is allowed. For example,
q=1, dz=1e-4, r²=1 is occupied after absorbed addition, while dz² <= r²-q falsely misses.
Mixed nonfinite rows fall back in full: original NaN and inf<=inf semantics are retained.
Preserve original atom order, lazy prefix extension, point order, row resets,
point denominators, counters, and final reductions. Assert every would-skip ORIGINAL
predicate is false, while still executing it and retaining its original result.

Counts: rows/point visits; actual lazy XY preparations/rejections; cached/new Z tests;
first hits/misses; finite fallback and inside range; endpoint checks/rejections;
would-eliminate cached/new predicates; per-atom counts and a distinct eager row×atom
upper-work model. Compare every frame's all15 f32 bits with unmodified occupancy.
Run exact current manifest arguments on diverse56 and actual10000 repeated files.
Compare complete raw stdout/stderr and scientific fingerprints to frozen baseline
and its current corrected native captures. Actual instrumented process durations
are diagnostic overhead, not evidence for any performance claim.

Necessary extra overhead includes O(points) row scanning, finiteness checks, branches,
cache metadata and endpoint predicates; net Z-test reductions alone are insufficient
for a speedup claim. Candidate implementation would require independent scientific
gates and quiet paired native trials after C1/C2 admission/reprofiling. Stop or defer
sustained counts if a native timing window starts. No native benchmark is run here.
