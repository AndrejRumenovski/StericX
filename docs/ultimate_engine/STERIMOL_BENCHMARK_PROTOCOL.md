# Prospective Sterimol competitive supplement

The final coverage audit requests explicit cross-tool Sterimol timing in addition
to the admitted volume workflow. This supplement uses all twelve `real_*` cases
from the already frozen 39-case continuous-support validation, selected by name
before collecting new outputs or timings. No failed case may be dropped.

The common task is raw forward L, continuous transverse minimum support B1 and
radial maximum B5 on the identical represented coordinates, axis and supplied
radii. The admission target is **0.01 Å maximum absolute error for every value
on every selected case**, fixed before the new comparisons. B1 uses the retained
independent 80-digit search bracket; L/B5 use independent 80-digit coordinate
formulas. These are validation references, not universal floating-point error
certificates. If any method misses the target, do not increase it or claim a
fair speed ratio for that method on the complete task.

StericX's explicit endpoint sampling API uses 361 directions; the separate
continuous-event API is also measured. Morfeus uses its unchanged public
Sterimol constructor with 361 vectors and its uncorrected L output. DBSTEP uses
its unchanged classic kernel (361 directions) after explicitly disclosed
unrounded axis alignment, as in the prior validated reference adapter. Its
absolute-axis L convention must independently satisfy the common target on
these inputs; passing this subset does not make the conventions interchangeable.

Timing covers a batch of all twelve molecules in memory. Parsing and conversion
to each library's ordinary input container occur before timing; alignment,
calculation and result construction occur inside. Native Rust remains the
scientific implementation. JSON output and imports/startup are outside the
internal timer. This is an API-workflow comparison, not an instruction-only
kernel comparison. All numerical residuals accompany the times.

After admission, retain one warmup and eight alternating method orders on CPU 2,
with 100 complete twelve-molecule batches per measurement. Use the same quiet
criteria and artifact flushing as the other final measurements. Preserve every
trial, source/package identity, raw result and failed attempt. Compare methods
only at this shared validated accuracy requirement; do not assert identical
actual errors or extrapolate the result to arbitrary molecules/default APIs.
