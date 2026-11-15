# Competitive benchmark admission and measurements

Status: fresh scientific admission and all-worker output preflight pass.
Controlled competitive timings are pending the final quiet measurement window.
No final cross-tool speedup is claimed in this report yet.

The post-validation host check failed the prospective quiet criterion (0.351
busy core equivalents versus 0.15). A quiet-window confirmation is pending;
no final competitive timing was started. The separate memory observations
below remain valid for their declared sampled-process scope.

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
reaching them. Final original-versus-candidate distributions and throughput
remain unmeasured until the controlled window.

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
