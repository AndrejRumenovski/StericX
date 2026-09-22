# Scientifically exact optimization

**The corrected, independently checked build is faster with exact scientific
outputs preserved.** Against that same corrected baseline, the final six-thread
10,000-file descriptor batch is **5.46× faster** (81.69% lower paired wall time),
and screening takes **16.7–18.6% less time** across 1/2/4/6 threads. C1/C1b and C2
are retained. C3 passed the scientific gates but was slower and was reverted.
The fresh final scientific replay and all 720 native launches passed. Small
parser regressions and increased parallel batch CPU/memory remain disclosed below.

## Scientific baseline

The denominator is the corrected, independently reviewed implementation at
commit `1ecfbd5be8b354711bdadd8a8d148417e03d8e7a`, following the
[scientific remediation](../scientific_remediation/SCIENTIFIC_REMEDIATION.md) and
[scoped admission](../scientific_remediation/ADMISSION.json). It is not the
audit-only `b515c4f` build or the original September 16 executable. The later
`3ea936d` CLI/CI changes are retained and were admitted separately as C1b.

The immutable local evidence root is
`.stericx/profiling/scientifically_validated_optimization/`. Its
`accurate_baseline_v1` manifest binds 76,714 raw evidence files, 211 copies and
22 gate receipts, including Cargo/Rust/scientific Python sources, tools,
executables, inputs, machine identity and outputs/errors. No baseline evidence
was overwritten. SHA-256 identities are:

| Artifact | SHA-256 |
| --- | --- |
| Corrected baseline executable | `102b6883c639bfc9ef211a88abae02391419cbd277c498098fa2bdfc7008206f` |
| Baseline freeze manifest | `feb021566f4a172c8d249184d1e16e39978b59c1784d216f4f50c9261ed27ced` |
| Corrected workload manifest | `a1041cc07fc5583c33d2a7cf584dbc638866c64bdd0b64e3720ef83cc85e1e6a` |
| Retained C2 executable | `b233501640e4d06555a36a278581fbf9cb011962311f205c81f3fc9a83ae89e6` |
| Retained C2 build manifest | `82b63f16c90e0767c5cb5b60daed61e2ce52f123d4f1a9f3156e791c5a03fd81` |

Optimization comparisons require **zero numerical tolerance**, including IEEE
bits, scientific fields, packed records, prediction vectors, ranks, exclusions,
validation and failures. Only previously named timing/resource/creation metadata
are excluded from their specific command comparisons. Independent equations and
reference tools are a separate gate, retaining their original limits. No
integration point, precision, boundary inequality, radius, donor/frame rule,
conformer population, model coefficient or statistical method was changed for
speed.

## Measurement methods and workload scope

Native release processes ran on the same Ryzen 5 5600G Linux machine with fixed
inputs, one warmup per side and alternating AB/BA measurements. The main matrix
uses eight measured pairs for each of ten workloads at 1/2/4/6 threads: 720
launches and 320 measured pairs. Affinity is 1→CPU 2, 2→CPUs 2–3, 4→CPUs 0–3,
6→CPUs 0–5. Warm filesystem cache, process startup, output and I/O are included.
Every sample is retained. Our heavy builds/reference jobs were stopped during
native measurement windows; the desktop and unrelated user applications remained
running. This is controlled application benchmarking, not an isolated host.

Report medians, median absolute deviations (MAD), full ranges, paired ratios,
user/system CPU, utilization, peak RSS, minor/major faults and context switches.
MAD is not a confidence interval. Paired speedup is the median of individual
baseline/candidate ratios and can differ from a ratio of medians. CPU 100% means
one core. Instrumented and modeled counters do not supply native wall timings.

The labels describe these fixed inputs, not newly independent chemistry:

- `single_small` and `single_large`: the frozen individual molecular inputs.
- `conformers_56`: 56 separate XYZ files from 11 Ni-hDA ligands.
- `descriptors_10000`: cyclic repetitions of those 56 files, retaining order.
- `ensemble_sdf`: one original bonded SDF ensemble with 86 conformers; C1 does
  not parallelize conformers inside one file.
- `parse_ensembles_1000` and `screen_1000`: repeated frozen reaction rows with
  complete original ensembles, not 1,000 new chemical structures.
- `db_build_ensembles`: fresh construction from 56 geometries and 11 ligands.
- `predict_1000000`: actual one-million-record inference from the fixed 64 MB
  packed input; every prediction is checked separately by the vector exporter.
- `search_database`: a freshly prepared 1,543-ligand/31,618-conformer database.
  Its 23 whole-ensemble source exclusions were declared before timing; no
  selected runtime geometry was silently dropped.

Coordination axes and supplied-weight aggregation are fixed. Ni-hDA response
metadata at 353.15 K does not reweight the historical 298.15 K populations or
alter targets. Exactness on these repeated throughput inputs does not increase
the independent scientific sample size.

## Profiling and candidate selection

The [new corrected baseline profile](CORRECTED_PROFILE.md) covers all 40 native
and 40 diagnostic configurations. It found descriptor file batches using about
one CPU regardless of requested threads. Occupancy remained the largest
arithmetic cost. Screening's repeated Student-t quantiles, not the old buried-
volume screening hypothesis, provided an independently measured opportunity.

| Measured component | Diagnostic share | Conditional ceiling | Decision |
| --- | ---: | ---: | --- |
| Descriptor per-file work | 99.16–99.27% of serial process wall | 5.76–5.79× with six ideal workers | Implement ordered file parallelism |
| Original descriptor occupancy | 81.77–82.05% of process wall | 5.49–5.57× if eliminated | Revisit after accepted changes |
| Original grid preparation | about 9% of process wall | about 1.10× if eliminated entirely | A realistic partial saving alone is below the major-change target |
| Screening Student-t quantile | about 20.0–20.35% of process wall | about 1.25× if eliminated | Reuse the identical model multiplier |
| C2 occupancy after reprofile | 82.464% of serial process wall | 5.70× if eliminated | Implement and test C3; reject on native result |

Amdahl ceilings assume other work stays fixed. Inclusive scopes overlap; summed
worker elapsed shares are not additive process-wall shares. The full rankings,
parsing, geometry/frame, Sterimol, grids, reductions, models, uncertainty/domain,
ranking, serialization and worker scope records remain in the linked profiles.
Per-worker CPU time, isolated scheduler wait and memory bandwidth were not
measured. The submit/join envelope does not measure scheduler overhead alone.

Fresh [counter checks](CORRECTED_COUNTERS.md) retained denied hardware `perf`
access without changing system permissions. Valgrind modeled approximately
68.08 billion instructions for the baseline 10,000-file workload, with occupancy
about 82.15% and grid preparation about 10.67%. DHAT attributed approximately
3.93 GB of requested allocation traffic to grids. These are baseline modeled
counts/site attribution, not final native hardware-event or allocation-time
measurements. Complete events and measurement limits are preserved separately.

## Candidate decisions

C1, C2 and C3 had a source/output freeze and prediction before editing. C1b
separately re-admitted the already-authored CLI/CI commit. Every stage passed
full scientific/reference/engineering gates before native timing and received
an explicit decision after measurement.

| Candidate | Pre-edit expectation and risk | Outcome |
| --- | --- | --- |
| C1: ordered descriptor file parallelism | Realistic 2–4× at four/six workers; preserve errors, duplicates and deterministic output order; higher worker memory/CPU possible | Accepted: large reproducible batch gain with disclosed resource and small unrelated costs |
| C1b: retain the later CLI/CI commit | Re-admit actual newer executable rather than transfer C1 timings to it | Accepted after full scientific replay and its own 720-launch matrix |
| C2: lazy Student-t multiplier reuse | About 20% removable screening work; retain leverage-first validation and unchanged interval arithmetic | Accepted: 15.5–21.7% direct screening reduction across thread settings |
| C3: exact row-endpoint rejection | Predict 15–25% lower serial descriptor time; require at least 10% target improvement; finite checks and row scanning may erase savings | Rejected and reverted: about 15.6% longer batch time against C2 at one and six threads |

The [pre-remediation checkpoint](PRE_REMEDIATION_OPTIMIZATION_CHECKPOINT.md)
retains earlier withdrawn experiments and hypotheses. Their speedups do not
support this corrected-baseline claim. Historical proposal/profile documents
retain their at-the-time status; this report supplies the final decision chain.

## Accepted implementations

[C1](candidates/C1_ORDERED_FILES.md) schedules independent files in a private
Rayon pool using indexed collection, then emits outcomes in original order.
Each file's scientific code, conformer ordering and arithmetic are unchanged.
Explicit batch-index checks occur before work; one-file/one-thread paths remain
serial; pool-construction failure falls back to the original serial path.
No numeric reduction is reordered. Duplicates and first-error behavior retain
their original meaning on fixed immutable inputs. This reduces serial waiting,
not the number of scientific calculations, and can increase total CPU work.

[C2](candidates/C2_QUANTILE_REUSE_RESULTS.md) immutably borrows one model's
training geometry. Construction performs no numerical work. Each query still
computes its original leverage first; the identical Student-t multiplier,
including an unavailable result, is lazily computed once. Candidate intervals
and validation results are never cached. The public uncached method is unchanged.
Frozen endpoint-bit tests cover invalid geometry, first-query failures, overflow
followed by recovery, signed zero and distinct models. The cache occupies 24
stack bytes and adds no cache heap allocation.

Fresh [C1b](C1B_REPROFILE.md) and [C2](C2_REPROFILE.md) profiles demonstrate the
change: every 1,000-row screen still makes 1,000 interval queries, but quantile
calls fall from 1,000 to one. The quantile share falls from about 20% to about
0.03%. Separate native comparisons establish the gain: 1,320 C2 admission
launches, including 96 direct screening pairs; 94/96 wall and 96/96 CPU pairs
favor C2. Both AB/BA median reductions exceed 10% at every thread setting.

## Rejected buried-volume experiment

[C3](candidates/C3_ROW_REJECTION_RESULTS.md) retained the original XY rows and
atom order. A finite, one-sided nearest-row-endpoint miss used the original f32
predicate to prove every point in that row misses the atom. Equality, interior
atoms and nonfinite inputs retained the original path. Independent source review,
six focused test groups, 128,021 exact responses, private bins and the full
independent reference replay passed.

The [actual operation probe](C3_OPERATION_COUNTS.md) measured 883,481,917 fewer
Z predicates, including the new endpoint checks: 45.4505% fewer Z predicates and
33.6009% fewer combined XY preparations and Z predicates under the declared
count model. Row-discovery comparisons are excluded from that reduction.
All 30,168 frame outputs/452,520 f32
fields matched both the candidate and original shadow. However, 972 native
launches/432 measured pairs showed worse runtime. The direct batch paired ratios
were approximately 0.865× at one and six threads, with 14/16 wall and 16/16 CPU
pairs slower. Both launch orders failed the predeclared threshold.

The operation probe counts row discovery, bounds and finite checks; these
counts help describe the failed hypothesis. Row discovery replaces the original
row-change checks, so not all listed work is net new. Branch behavior and each
check's native cycle cost were not isolated.
Fewer distance predicates did not establish a faster algorithm. Only
`src/geometry/buried_volume.rs` was restored from its pre-edit snapshot, including
removal of the rejected helper's tests. All scientific source bytes then matched
accepted C2. The candidate source, patch, six test groups, failures, scientific
passes and every native sample remain available for reproduction.

## Final scientific verification

The retained executable is the same C2 binary used for admission; no scientific
kernel from rejected C3 remains. The fresh post-performance campaign passed all
128,021 exact responses, 93 CLI cases at each of four thread settings, 24
additional cases and 31 separately frozen command-boundary cases. All 25
independent-reference commands passed with unchanged scientific expectations.
The model replay retained six creation-time-derived provenance differences;
an additive resolution verified their nonnumerical origin and all 7,624 model
comparisons. The raw failed provenance gate remains preserved.

The final engineering gate passed 319 Rust tests, 124 Python tests, Clippy,
rustdoc with warnings denied, Rustfmt and Ruff. Four reference-integrity tests
and the Study 011 verification also passed. The freshly rebuilt delivered
`target/release/stericx` has exactly the retained C2 hash above. Final reference
gate SHA-256 is
`5efb637c764f13a299103509073d7e813e4932d089f5ea2953d182114296d89a`;
final engineering completion is
`3c761d0cc42eb0d38dba2a72a62d9e45627bfbc0a17c6cf2f4c788a44b8da680`.

The [scientific equivalence report](POST_OPTIMIZATION_SCIENTIFIC_EQUIVALENCE.md)
separates optimization equivalence, implementation correctness, methodological
validity and predictive/experimental validity. Full topology coverage retains
31,721 conformers, 95,163 orientations and 1,141,956 bins; the inferred-connectivity
lane retains its eight expected errors. Reference limits, B1 angular sampling,
f32/grid effects, radius/center conventions, source conflicts, failed historical
folds, unfavorable interval coverage and screening outcomes remain. No scientific
classification is promoted merely because an optimized build matches baseline.

## Final native performance

The fresh final matrix compares corrected baseline `102b6883…` directly with
final `b2335016…`. The independent reviewer verified all 720 launches, 320 pairs,
raw streams/artifacts, resources, launch order, source identities and the full
scientific gate. All measured scientific fingerprints match. No sample is excluded.

These are six-thread native end-to-end medians in milliseconds; paired speedup
is calculated from the eight individual pairs, not from the displayed medians.

| Workload | Accurate baseline ms | Accurate optimized ms | Paired speedup | Scientific match |
| --- | ---: | ---: | ---: | --- |
| `descriptors_10000` | 4,969.625 | 911.846 | 5.462× | Exact |
| `conformers_56` | 29.607 | 7.161 | 4.188× | Exact |
| `screen_1000` | 250.585 | 206.862 | 1.206× | Exact |
| `db_build_ensembles` | 8.053 | 8.335 | 1.004× | Exact |
| `predict_1000000` | 9.815 | 8.949 | 1.068× | Exact |
| `search_database` | 7.511 | 10.227 | 0.969× | Exact |

The [complete final table](FINAL_NATIVE_MEASUREMENTS.md) and
[machine-readable review](current_final_measurements.json) retain all 40
configurations, medians/MAD/full ranges, individual ratios, AB/BA results,
user/system CPU, utilization, RSS and faults. The earlier C2 admission and direct
incremental series remain separate; their more favorable samples are not
substituted for this final series. Search is highly dispersed and order-sensitive;
its displayed ratio is not a reliable speedup estimate. Short-workload changes
near measurement noise are not promoted to optimization benefits.

Final descriptor paired speedups at 1/2/4/6 threads are
1.001/1.954/3.784/5.462×; every measured parallel-batch pair is faster.
The 56-file speedups are 1.005/1.817/3.132/4.188×. Final screening paired speedups
are 1.204/1.229/1.201/1.206×, with all 32 wall and CPU pairs faster. Both launch
orders retain more than 10% lower screening time at every thread setting.

Independent final native review SHA-256:
`a428b351d87e4454314230875b02f2696448371e2c61d4c0aa707b4261895757`.

## Memory, allocations and thread scaling

The accepted-C2 reprofile already contains 200 measured allocation captures,
five for each of 40 configurations, with zero allocation failures or dropped
scopes. Its native admission contains 320 measured candidate resource samples.
The fresh final matrix peaks at 26.48 MiB for the six-thread 10,000-file batch
and 69.70 MiB across measured configurations during prediction. For that batch,
the median peak-RSS increase is 5.57 MiB and paired CPU cost is +5.05%; the
56-file batch has +13.73% paired CPU cost and a 2.23 MiB median peak-RSS increase.
Final six-thread batch CPU utilization is 570% of one core. These resources
are bounded on the tested inputs, not a promise for arbitrarily large batches.

Parallel batching adds an ordered outcome buffer and worker-local state. Its
CPU and memory increase is part of the tradeoff, not evidence of less arithmetic.
All screening allocation/reallocation/free call counts match C1b. Tiny requested-
byte differences follow executable/report path lengths and are not cache savings.
Requested heap counters exclude mappings, stacks, allocator bookkeeping and
report construction; they are different from RSS and allocation-time attribution.
Final C2 byte identity permits reusing these complete allocation measurements.

The [C2 profile](C2_REPROFILE.md) records near-even file distribution at two/four
workers and a median max/mean file count of 1.011 at six. Native admission batch
CPU utilization reached about 574% of one core. Single-file ensembles remain
serial, and small jobs have proportionally larger scheduling/startup costs.

## Costs and limits

C2's direct incremental parser regression is reproducible: about 1.35–1.67%
wall and 1.4–1.7% CPU, approximately 2.6–3 ms. Both orders agree; unchanged parser
source does not establish the cause. The screening gain justifies that bounded
cost for this workload. Earlier parallel batch admission also retains increased
CPU/RSS and small startup costs. Search has broad, order-sensitive wall
variation; no precise search speedup or absence of regression is claimed.
The fresh final matrix independently retains 1.46–1.99% paired parser wall costs
and 1.48–1.97% CPU costs, roughly 3.4–3.7 ms higher medians. Search retains
negative rows and CPU costs. Occasional long startup/database samples and the
recorded major page fault remain in the full ranges and raw resource records;
no rerun replaces them. The targeted batch/screening gains justify these
disclosed costs for the project's fixed workloads; this is not an assertion of
zero regression in every command.

All numerical conclusions apply to the frozen inputs, executable, toolchain,
platform and declared conventions. Finite adversarial/reference campaigns are
not proofs for all chemistry or all hardware. Repeated throughput rows do not
validate predictive chemistry. No new prospective experiments, recalibration,
model-selection method, confidence interpretation or applicability guarantee was
introduced. The task establishes targeted performance gains under unchanged
scientific limits, not universal acceleration or experimentally reliable ranking.

## Reproduction and preserved evidence

The local immutable root retains executables and every candidate stage.
Durable, content-addressed packages under `evidence/` retain sources, inputs,
outputs/errors, manifests, controllers, raw measurements and recovery instructions:
[accepted C1](evidence/accepted_c1/README.md),
[accepted C1b](evidence/accepted_c1b/README.md),
[accepted C2](evidence/accepted_c2/README.md),
[rejected C3](evidence/rejected_c3/README.md), and the
[final replay and benchmark](evidence/final_optimization/README.md).
Shared providers avoid duplicating identical multi-gigabyte observation streams.
Archives omit ELF binaries, target directories and virtual environments; native
hashes, exact sources/build commands and the retained local executables identify
what must be rebuilt on a compatible toolchain. Offline restoration/verification
is distinct from rerunning science or native measurements.

Failed tooling attempts remain: initial cross-thread metadata comparison,
C1b documentation formatting, the exact creation-time provenance mismatch,
quantile-probe warmup configuration, overly strict allocation-byte equality,
C3 probe fixture setup, the boundary harness's empty-JSON expectation, and
publication/controller preparation failures. None licenses a changed scientific
expectation, omitted outlier or rewritten baseline. See the individual receipts
and package inventories for the exact failed and successful attempts.
