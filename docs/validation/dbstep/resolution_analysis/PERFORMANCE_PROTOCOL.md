# Matched performance after discrepancy resolution

Frozen before new timing. The pre-correction baseline is
`f04b156f7b8b31700e1c6ab978e337b00cc734e0`, freshly built under
`.stericx/dbstep_resolution/baseline/build`. The candidate is the independently
built native executable in `mask_build_v2`. All 53 scientific source/Cargo files
and both native executable hashes are identical. The corrected research
comparison does not enter production execution. This is still a fresh matched
measurement; the prior approximately 10,825 files/s is historical context only.

Reuse `benchmark_paton_extensions.py` and the fixed workload manifest with no
changes to the measurement code: seven workloads, 1/2/4/6 threads, established
CPU affinities, one warmup per executable and eight alternating AB/BA pairs per
configuration. The resulting 504 attempts retain all raw output, process latency,
CPU time, RSS and scientific fingerprints. Apply the existing 10% investigation
trigger to either the ratio of side medians or the median within-pair ratio.
No outlier removal, cache flush or post-hoc selection of individual samples.

Scientific timing admission must bind the approved policy, complete legacy
replay, independently reviewed new volume evidence, advanced Python/conformer
validation and actual three-platform CI. Stop all agent builds/scientific jobs
before timing, record host state, and obtain a fresh quiet-window confirmation.
Record any whole-campaign repeat and its reason; preserve the original attempt.

Allocation diagnostics use separately compiled copies of the frozen source
trees with the existing profiling feature. Their native seven-workload 1/6-thread
results remain separate from throughput measurements. Also run the existing
six-conformer warmed allocation example for each side, recording its complete
stdout and process-resource measurements.

Python overhead uses the same six-conformer input, options and unchanged
`benchmark_paton_python.py`: 100 calls per block, eight blocks, one CPU, parsed
input reuse, and exact native/bridge/typed-result equality checks. Use a fresh
installation of the previous pre-correction wheel as baseline and the newly
tested wheel as candidate. Run whole campaigns in **A/B/B/A** order, retaining
all sixteen block samples per method per side. This adds a matched package
comparison while retaining the established within-campaign measurement method.
Package hashes and scientific source identities must be independently verified.

Report every workload, including short-process OS noise. Distinguish per-process
RSS from allocator requested bytes and from the Python call timing scope.
Because production bytes are unchanged, any measured difference is host/sample
variation rather than an algorithmic speedup or slowdown; a substantial shift
still requires investigation and disclosure.
