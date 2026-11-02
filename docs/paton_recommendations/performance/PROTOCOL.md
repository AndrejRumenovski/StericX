# Paton extension performance protocol

Frozen before timing on 2026-11-02 UTC. This is a new comparison against the
task's source-bound `c9380d2c` baseline, not a rerun or extension of the earlier
corrected-baseline optimization admission.

## Native descriptor workloads

Reuse the immutable `corrected_workloads_v2/workloads/manifest.json` inputs and
commands for `single_small`, `single_large`, `ensemble_sdf`, `conformers_56`,
`descriptors_10000`, `parse_ensembles_1000`, and `db_build_ensembles`. These cover
single-molecule latency, conformer aggregation, batch throughput, packed export
and database construction. The earlier prediction, model-screening and database
search timing workloads are outside this descriptor-extension comparison;
scientific model replay remains a separate validation gate.

Use the existing native `wait4` launcher and measurement/fingerprint functions.
Measure process startup, calculation and regular-file output together. Rehash
the frozen input inventory and both executables before and after the campaign.
For each workload use 1, 2, 4 and 6 Rayon threads with the established affinities
`2`, `2,3`, `0,1,2,3`, and `0,1,2,3,4,5`. Use one warmup per binary and eight
paired repetitions, alternating baseline/candidate then candidate/baseline.
Keep every attempted sample and raw stdout/stderr. Require exact scientific
output fingerprints against the baseline warmup, including packed/table files;
exclude only the already documented volatile timing/path fields.

Report median, MAD, range, paired baseline/final wall-time ratio, per-child peak
RSS, CPU time, and completed records per second. A median slowdown greater than
10% is an investigation trigger, not permission to change scientific constants
or discard samples. Timing must not overlap scientific jobs, builds or other
known heavy task processes. Warm-cache policy; no cache flushing or outlier
removal. Short process timings can remain dominated by startup and OS noise.

Admission requires the complete byte-exact 128,021-record legacy replay, current
engineering checks, and reviewed independent reference results with all known
scientific limitations retained. Exact replay proves numerical preservation,
not universal physical validity.

## Allocations and Python overhead

Build each frozen source tree with its existing `profiling` feature in separate
targets. Run allocation diagnostics separately from native timing. Record the
unchanged `TrackingAllocator` counts, requested bytes and peak live requested
bytes, plus output identity. These counts exclude allocator metadata and are not
RSS. Diagnostic instrumentation is never used for throughput claims.

The new structured ensemble API has no historical Python counterpart. Measure
its in-memory native calculation/statistics and the public Python method on the
same retained six-conformer ligand-401 input and options. Separately measure
Python's private native JSON bridge and public typed-result construction, so
serialization and object-construction costs are visible. Reuse parsed inputs,
one warmup and eight alternating paired blocks, single CPU, fixed iteration
count, no weighting or sensitivity CSV in this comparison. Check identical
result values before timing. Report absolute times and the scope of each layer;
do not call these a baseline speedup. Density performance is supplementary and
does not replace descriptor or ensemble measurements.

All raw samples, build identities, controller source hashes, machine/environment
metadata and failed attempts remain in the local evidence directory; committed
receipts and the report link those artifacts.
