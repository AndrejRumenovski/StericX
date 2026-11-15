# Campaign protocol

Baseline: c03630cb; source/datasets frozen before production edits. Existing
`docs/media/` is unrelated user material and is not staged or modified.

Measurements use the actual Rust implementation. Observer code is a measurement
adapter, not an independent mathematical reference. Portable release builds use
thin LTO, one codegen unit and line-table debug information. A native CPU target
is a separately labelled experiment. Never enable fast-math or silently change
FMA, operation grouping, grids, radii, centers, normalization or populations.

First run the existing Rust/Python tests and fresh descriptor profiling. Use
existing 56-conformer/11-ligand and 10,000-file workloads, the six-conformer
public API fixture, density cubes and analytical/adversarial geometries. Hash
every consumed input. Historical comparisons are context, not fresh measurements.

Instrumented timings, Valgrind counts and native timings are separate evidence.
Hardware perf is denied on this host; Cachegrind events are simulated and cannot
provide native cycles, IPC or memory bandwidth. A Callgrind instruction flamegraph
must be labelled as such. Record all missing capabilities and actual CPU affinity.

Exploratory measurements identify candidates; final acceptance requires exact
legacy values/errors and independent scientific comparisons. Each candidate gets
a pre-edit hypothesis, correctness checks, native alternating A/B measurements,
and an explicit retain/reject decision. Rejected source/measurements are preserved
outside production. Accepted changes get separate commits.

New numerical methods must have explicit names, quantity definitions, independent
references, adversarial tests and accuracy/cost measurements. No certified bound
is claimed from an empirical convergence estimate. Existing defaults stay fixed.

Final paired native measurements: one warmup per side, eight alternating AB/BA
pairs, fixed input/affinity/thread settings; report every sample, medians, MAD,
ranges, CPU time, RSS, context switches and output hashes. Include 1/2/4/6/12
threads. Do not run builds, profilers or reference jobs in the timing window.
Do not discard unfavorable samples. If host activity makes final measurement
unreliable, finish all independent work first and request a quiet window.

No releases, publication, GUI work, remote messaging or unrelated changes.
