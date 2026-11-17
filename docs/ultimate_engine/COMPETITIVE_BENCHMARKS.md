# Competitive benchmark admission and measurements

Status: scientific admission, worker equivalence, separate memory sampling and
controlled volume timing are complete. All 432 volume launches pass their
scientific checks. The additional Sterimol cost/accuracy experiment is complete;
two external methods miss its fixed target and are excluded from equal-accuracy
ratios. No whole-library or universal fastest-engine claim is made.

## Matched task

The admitted comparison is total and regional hard-sphere volume on the exact
radius-3.5 Å, density-0.01 lattice, using three donor-plane frames, scale 1.17,
the explicitly matched radius table and hydrogen exclusion. Inputs are the
same 56 immutable XYZ conformers from 11 ligands. The 1,000-calculation batch
repeats those paths; it is not 1,000 distinct molecular structures.

StericX calls its public Rust volume API. Morfeus uses its public BuriedVolume
and octant-analysis APIs with the existing explicit geometry adapter. DBSTEP
uses its unchanged `sterics.occupied` kernel on the common projection grid,
followed by explicitly declared regional reductions. **This is not DBSTEP's
default grid-generation/volume workflow.** Its own default additional sphere
clipping is scientifically different at certain boundaries and is preserved
as a separate failed historical lane.

The unchanged independent reference derives ideal grid membership from integer
lattice indices, encloses geometry/occupancy from exact decimal inputs, and
models each tool's scalar reductions with exact rational IEEE rounding. All
56 conformers, 168 frames and 2,588,544 point decisions pass; no point is
uncertain or differs from the independently certified occupancy. Every native,
Morfeus and DBSTEP-adapter scalar and regional result matches its independently
rounded expression. There is no fitted residual tolerance.

The measured adapters also return exactly the expected ordered per-tool values
at 1, 2, 4, 6 and 12 workers. The final native source is frozen; its competitive
build uses the same dependency versions as the profiling observer. An initial
unpinned adapter build resolved newer transitive dependencies. It was retained
but excluded, then replaced with the explicitly locked build before admission.

## Timing scope and fairness

[COMPETITIVE_PROTOCOL.md](COMPETITIVE_PROTOCOL.md) predates the gate. The
measurement controller preserves every attempted process, output hash, return
code, CPU time, wall time, context switches and wait4 high-water RSS. Final
timing starts only after the host-activity check passes. One warmup and eight
alternating tool orders are retained; no unfavorable trial is discarded.

Fresh-process measurements include imports/startup, XYZ parsing, topology and
center/frame construction, three volume calculations, aggregation and JSON
output. The operating-system file cache is warmed; "fresh process" does not
mean caches were forcibly evicted. Warm workflow timings use an already started
process/pool and repeat the same parsing/calculation/result-construction work,
excluding final JSON serialization. The same scope applies to every tool.

The single-conformer case runs at one worker. The 56-conformer scaling matrix
uses 1/2/4/6/12 workers. The repeated 1,000-calculation case uses one and six
workers for final distributions. StericX's separate controlled internal A/B
also includes the established 10,000-file batch at every thread count.
Python BLAS/OpenMP threads are limited to one. DBSTEP's internal
`query_ball_point(workers=-1)` is retained and disclosed; affinity limits the
CPUs available, while any scheduling/oversubscription cost remains in its data.

Memory is measured in separate launches. Concurrent process RSS is sampled
and summed at each snapshot; independent process high-water marks are never
summed. RSS double-counts shared pages and sampling may miss short peaks, so
the report retains sample spacing and does not call it unique physical memory
or an exact maximum. Instrumented or sampled-memory elapsed time is excluded
from speed claims.

## Scope that remains separate

Legacy Sterimol angular grids, transverse phases, B5 sampling and L conventions
differ between tools. Fresh analytical/high-precision comparisons are reported
in [NUMERICAL_ACCURACY.md](NUMERICAL_ACCURACY.md), including derived angular
discretization envelopes. They do not authorize a blanket combined-descriptor
speed ratio. Density and hard spheres also describe different quantities.

Conformer aggregation, Python bridge/object/NumPy overhead, individual calls
and native batches are measured separately by the campaign's Rust and Python
observers. They are not relabelled equivalent whole-library operations when
the external tools perform different work. Cross-tool speed claims are limited
to the admitted finite-grid volume task.

## Baseline and targets

The fresh original implementation's exploratory 10,000-file medians were
1,977.2 files/s at one thread and 10,264.5 files/s at six threads. These came
from this host and frozen original source, but still include desktop activity.
The historical 10,774.5 files/s value is context only. The campaign's 15,000
and 20,000 files/s targets are aspirations; no accuracy change is justified by
reaching them. The final paired baseline is 2,080.3 files/s at one thread and 11,242.0 at six;
the final engine reaches 4,503.3 and 23,783.8 respectively. The six-thread gain
is 111.56% in throughput, preserving exact scientific outputs. Both targets
are exceeded on this frozen workload; they are not universal throughput claims.

No claim that StericX is universally fastest, universally accurate, or superior
outside the measured inputs and definitions follows from these experiments.

## Reproduction

```bash
.stericx/paton/reference/venv/bin/python docs/ultimate_engine/competitive_gate.py \
  --build .stericx/ultimate_engine/competitive_build_pinned --output NEW_GATE
.venv/bin/python docs/ultimate_engine/competitive_preflight.py \
  --build .stericx/ultimate_engine/competitive_build_pinned \
  --gate NEW_GATE/summary.json --output NEW_PREFLIGHT
.venv/bin/python docs/ultimate_engine/competitive_measure.py timing --output NEW_TIMING
.venv/bin/python docs/ultimate_engine/competitive_measure.py memory --output NEW_MEMORY
```

The measurement controller defaults to the retained campaign gate/preflight
paths; supply `--gate NEW_GATE/summary.json --preflight NEW_PREFLIGHT/summary.json`
to bind the newly reproduced checks. `--native` selects a rebuilt workflow
observer. Frozen build receipts, lockfiles and
all scientific source hashes are under `.stericx/ultimate_engine/competitive_*`.
The retired historical benchmark runners are not used to overwrite old results.

## Separate memory observations

The 1,000-calculation process-tree samples include parent, workers and resource
tracker. These are observed concurrent RSS sums, not unique memory or exact
peak requirements. Maximum observed sample gaps are about 10–26 ms. Timing from
these sampled launches is not used for speed comparisons.

| Workers | StericX MiB | Morfeus workflow MiB | Matched DBSTEP workflow MiB |
|---|---:|---:|---:|
| 1 | 5.24 | 80.46 | 80.23 |
| 2 | 5.47 | 247.22 | 246.58 |
| 4 | 5.95 | 400.32 | 402.04 |
| 6 | 6.48 | 555.68 | 555.28 |
| 12 | 7.91 | 1,017.68 | 1,014.89 |

The reference workflow module imports both reference packages. These process
figures therefore describe the declared benchmark environment, not the minimum
possible standalone installation footprint of either library. All output rows
match their admitted references exactly. Raw per-process snapshots and maximum
sample spacing are retained under `.stericx/ultimate_engine/competitive_memory`.


## Final controlled timing distributions

[FINAL_TABLES.md](FINAL_TABLES.md) and
[final_measurements.json](final_measurements.json) retain all eight independent
trials per configuration/tool/scope, median, unscaled MAD, min/max, process CPU,
RSS and context switches. Warm iterations within a trial are averaged before
summarizing the eight trial means; they are not treated as independent
replicates. Original inner samples and the one declared warmup per setting are
all retained. Both outer and inner quiet gates passed.

| Volume workflow | StericX median ± MAD | Morfeus median ± MAD | Matched DBSTEP median ± MAD |
|---|---:|---:|---:|
| One conformer, fresh process | 1.226 ± 0.051 ms | 291.869 ± 0.859 ms | 290.458 ± 0.377 ms |
| One conformer, warm | 0.162 ± 0.001 ms | 14.118 ± 0.014 ms | 13.553 ± 0.013 ms |
| 56 conformers, one worker, warm | 11.514 ± 0.017 ms | 756.107 ± 1.493 ms | 629.196 ± 2.696 ms |
| 56 conformers, six workers, warm | 2.143 ± 0.008 ms | 166.786 ± 2.280 ms | 155.516 ± 3.865 ms |
| 1,000 calculations, one worker, warm | 204.652 ± 0.188 ms | 13,931.327 ± 158.103 ms | 12,344.950 ± 38.400 ms |
| 1,000 calculations, six workers, warm | 36.521 ± 0.067 ms | 2,769.444 ± 10.631 ms | 2,570.056 ± 19.574 ms |

These are API workflows with equivalent finite-grid scientific work, including
parsing/frame construction, rather than isolated occupancy instructions. The
fresh-process lane includes Python imports and new pools; the warm lane excludes
those costs and final JSON. It is a process-cold start on warmed file caches,
not a controlled disk-cache-cold measurement. Initial launch/page-fault behavior
is retained in warmup records; no cache was evicted to manufacture a cold ratio.
Both reference packages are imported in the same adapter environment.

Warm 56-conformer six-worker efficiencies relative to one worker are
89.53%/75.56%/67.43% for native/Morfeus/DBSTEP workflows. At twelve they fall to
50.73%/45.65%/40.75%. Fresh-process pool startup is substantially less efficient;
all twelve-worker fresh medians are slower than six. Full per-worker results
are published, including configurations that do not benefit from more workers.

## Sterimol supplement: failed accuracy admission retained

[STERIMOL_BENCHMARK_PROTOCOL.md](STERIMOL_BENCHMARK_PROTOCOL.md) fixes all twelve
`real_*` cases, identical represented coordinates/radii/axis, raw forward L,
continuous support B1, radial B5, and a 0.01 Å target on every output before new
results. Independent 80-digit formulas supply L/B5; the retained independent
80-digit angular-search bracket supplies B1. This is empirical validation on
the selected inputs, not a universal floating-point certificate.

| Method | Twelve-molecule batch median ± MAD | Maximum absolute residual | All targets pass |
|---|---:|---:|---|
| StericX explicit 361 directions | 121.073 ± 0.156 µs | 0.0088093 Å | Yes |
| StericX continuous events | 167.353 ± 0.167 µs | 2.93e-12 Å | Yes |
| Morfeus 361 directions | 1,890.307 ± 9.887 µs | 0.0125479 Å | No |
| DBSTEP classic 361 directions | 25,030.612 ± 68.302 µs | 0.0155439 Å | No |

Morfeus fails trimethylphosphine B1; DBSTEP fails triphenylphosphine B1. All
other selected outputs pass. No molecule was dropped, tolerance increased or
sampling grid altered after the result. The external timings are cost/error
observations; no equal-accuracy speed ratio against either failed method is
claimed. The two admitted native methods show a real tradeoff: continuous events
are more accurate but cost more on this batch. Earlier 39-case external
convention validation uses its own declared resolutions and is separate.

The common in-memory scope excludes parsing/input-container preparation and
imports, includes axis alignment/calculation/result construction, and uses 100
complete twelve-molecule batches per trial. One warmup and eight alternating
method orders are retained, with CPU 2 affinity and BLAS/OpenMP limited to one.
DBSTEP receives the same explicitly disclosed unrounded alignment as the prior
reference adapter; its unchanged absolute-axis L kernel happens to meet the
forward-L target on this subset. This does not make the definitions identical
for arbitrary inputs. Morfeus uses `L_value_uncorrected` and 361 vectors.
Versions: Morfeus 0.8.0, DBSTEP 2.2.0, NumPy 2.5.4, SciPy 1.18.1.

Build/input/controller hashes, independent references and all residuals are
retained under `.stericx/ultimate_engine/sterimol_competitive_*`; the final seal
binds them. Reproduction uses `sterimol_competitive.py preflight` followed by
`timing --gate ... --output NEW_PATH`, on a quiet host. The Rust scientific
implementation remains authoritative for StericX.
