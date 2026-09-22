# Scientifically exact optimization — incomplete checkpoint

**No production optimization has been accepted. No speedup is claimed.** The
file-parallel experiment was withdrawn before benchmarking when the user explicitly
requested scientific corrections first. Its source was restored and patch retained.
This checkpoint preserves earlier measurements; a corrected build must be validated
and profiled afresh. It is not a completed optimization report.

## Scientific baseline

Current commit: `b515c4f1736a1519be19a3ecdf40c3c1327ac9b4`.
Fresh native release executable SHA-256:
`b577c49d3e98b4745213fce0e70e93a55fc61120c5abe6fdb0658547c43cab94`.
All 42 Cargo/Rust files match the independently audited implementation. This
commit publishes the audit; it does not contain subsequent scientific corrections.

The current audit documents valid-symmetry rejection, automatic atom-order
dependence, a CREST temperature problem and invalid-prediction validation failures,
among other limitations. Preserving these observations cannot establish that
those behaviors are scientifically correct. The requested unrestricted validated
baseline was therefore unavailable at that commit. The user has now authorized
[remediation](../scientific_remediation/REMEDIATION_PLAN.md), in progress in the
working tree. The current status is in [BASELINE_STATUS.md](BASELINE_STATUS.md).

The frozen source, Cargo files, scientific inputs, executable, compiler/package
versions, CPU/OS identity and audit results are under
`.stericx/profiling/scientifically_validated_optimization/current_audited_b515c4f/`.
All 74,298 original sealed audit files were verified. The executable is preserved
independently of the mutable Cargo build directory. Baseline evidence is never
overwritten by a later capture.

## Profiling

Fresh native release measurements cover all ten requested workloads at 1/2/4/6
threads, one warmup and seven measured launches per configuration. Separate
diagnostic, Cachegrind and heap investigations attribute work without substituting
instrumented time for native runtime. See [CURRENT_AUDITED_PROFILE.md](CURRENT_AUDITED_PROFILE.md)
for all 40 native configurations, medians, MAD, full ranges, CPU/RSS, allocation
and counter evidence, affinity and measurement limitations.

The current ranked opportunities are:

| Component/workload | Fresh measured share | Conditional ceiling |
| --- | ---: | ---: |
| Independent per-file descriptor work, batch | 99.2547% diagnostic wall | 5.784× with six workers |
| Buried-volume occupancy, batch | 81.753% diagnostic wall | 5.480× if eliminated completely |
| Repeated Student-t quantile, screen | 15.651% at one thread; 33.60% at six | 1.186× / 1.506× if eliminated |
| Integration-grid preparation, batch | 9.134% diagnostic wall | 1.101× if eliminated |

These scopes are not additive, and the first two overlap. The conditional ceilings
are Amdahl calculations, not performance results. Native batch processing remains
near one-core utilization even with six configured threads.

## Candidate table

| Candidate | Evidence before implementation | Status |
| --- | --- | --- |
| Ordered parallel processing of independent descriptor files | 99.25% file work; kernels already used in database workers | Withdrawn before benchmarking; user selected scientific remediation first |
| Exact row-endpoint buried-volume rejection | Shadow probe finds 45.45% net Z-predicate reduction after endpoint costs | Nonproduction count/proof experiment complete; production candidate pending |
| Cache the invariant model-local Student-t multiplier | 1,000 repeated quantiles and 200,000 bisection iterations per screen | Deferred until corrected baseline is independently checked and profiled |
| Bounded exact grid reuse | Material allocation traffic, smaller measured wall share | Deferred pending higher-impact candidates |

Expected realistic savings, exactness arguments, risks, memory effects and
acceptance thresholds were recorded before implementation in
[CANDIDATE_ADMISSION.md](CANDIDATE_ADMISSION.md).

## Accepted optimizations

None. Neither an exact operation-count reduction nor an Amdahl bound passes the
native performance gate. There is no optimized executable to compare yet.

The isolated occupancy probe executed every original predicate, including every
predicate its proposed rejection rule would skip. It confirmed that all such
predicates miss and preserved all fifteen per-frame f32 values and native
workload output fingerprints. Across 30,000 frames it counted 1,943,834,272 current
Z predicates, 969,971,505 removable predicates and 86,487,986 added endpoint
predicates: a net reduction of 883,483,519. This experiment measures potential
work elimination, not production wall time; row scans, bookkeeping and branches
could offset the arithmetic savings.

## Rejected optimizations

No production candidate has reached an accept/revert gate. Preparatory harness
failures and strengthened guards are retained separately in
[BASELINE_SCIENTIFIC_RECHECK.md](BASELINE_SCIENTIFIC_RECHECK.md). They must not be
misreported as failed performance experiments. Rejected candidate patches,
scientific comparisons and timing evidence will be recorded here when experiments
are actually run.

## Scientific verification

Fresh observations replay the full 31,721-conformer corpus, all available regional
bins and focused geometry, kinetics, aggregation and model-domain witnesses.
The CLI/Python oracle additionally covers errors, packed records, model fitting,
screening, uncertainty, applicability, rankings and exported decks at all four
thread settings. The final oracle versions are Rust v3 and CLI v5; an exact
cross-version bridge binds the stronger Rust capture to all previously completed
independent-reference inventories. New independent reference calculations cover geometry/Morfeus,
Kraken, thermodynamics, Eyring and model/statistical equations.

The individual scientific values and known failures reproduce the audit.
The strict historical reference replay nevertheless records two differences:
traceback source formatting and one ULP in a BLAS-reduced aggregate R². Their
unchanged raw failures, descriptor rows, controlled reduction experiment and
high-precision diagnosis are preserved; no scientific tolerance was widened.
The [scientific recheck](BASELINE_SCIENTIFIC_RECHECK.md) states both its completed
coverage and the supplemental analyses that remain historical evidence.

The current baseline passes 264 Rust tests, the separately run published-screening
regression, Clippy, rustdoc, Rust formatting, the updated 71-test Python suite and
Study 011 verification. Python lint and formatting pass with the user's existing,
untracked `docs/media/` drafts excluded and left untouched. These are engineering
results. They do not erase negative
independent scientific findings. A candidate and its post-optimization independent
recheck are still required.

## Performance comparison

No baseline-to-optimized claim is available. The table below records one-thread
native baseline medians only; the baseline is **audited**, not globally corrected.

| Workload | Current audited baseline (ms) | Optimized (ms) | Speedup | Scientific match |
| --- | ---: | ---: | ---: | --- |
| `descriptors_10000` | 5,044.513 | — | — | No candidate |
| `conformers_56` | 30.084 | — | — | No candidate |
| `screen_1000` | 323.953 | — | — | No candidate |
| `db_build_ensembles` | 31.069 | — | — | No candidate |
| `predict_1000000` | 20.824 | — | — | No candidate |
| `search_database` | 6.240 | — | — | No candidate |

After scientific remediation and the new baseline freeze, each candidate must pass scientific gates first,
then alternating native A/B or AB/BA timing under the identical recorded controls.
An improvement must exceed measurement noise and disclose other-workload
regressions. Final combined measurements must use the same preserved baseline.

## Memory

Current one-thread native RSS medians are 18.13 MiB for the descriptor batch,
8.84 MiB for screening and 69.50 MiB for million-record prediction. Complete
resource distributions and diagnostic allocation counts are retained. DHAT's
56-conformer replay requests 23,340,655 bytes, with a 425,706-byte live peak;
integration grids account for 94.34% of requested allocation traffic. Allocation
traffic is not peak RSS or a measurement of time spent allocating.

No candidate memory delta is measured. File parallelism will require particular
scrutiny of concurrent grids, whole-file parsing, result buffers and worker stacks.

## Thread scaling

Baseline median milliseconds, from the complete native matrix:

| Workload | 1 thread | 2 threads | 4 threads | 6 threads |
| --- | ---: | ---: | ---: | ---: |
| Descriptor batch | 5,044.513 | 5,030.745 | 4,992.289 | 5,008.436 |
| 56 conformers | 30.084 | 29.443 | 29.763 | 30.285 |
| Screening | 323.953 | 218.620 | 155.138 | 138.913 |
| Database build | 31.069 | 17.027 | 10.115 | 7.693 |
| Prediction | 20.824 | 13.861 | 10.096 | 8.813 |

The other five workloads, full ranges and utilization are in the linked complete
profile. Worker timing sums overlap caller waits and must not be interpreted as
disjoint wall-time fractions.

## Limitations and remaining gates

- Scientific remediation is in progress. An optimization against this
  audited build must not be described as starting from a scientifically corrected
  build or as proving that the audit passes.
- The fixed benchmark corpus deliberately repeats geometries/records; it is not
  evidence for all chemistry, all input sizes or experimental prediction.
- Host perf permissions prevent native hardware counters. Cachegrind events are
  modeled; diagnostic wall shares contain instrumentation overhead.
- No accepted candidate, paired candidate timings, final combined corpus replay
  or post-optimization independent scientific audit exists yet.
- Equality with an audited implementation establishes optimization equivalence
  only. Published-algorithm correctness, methodological validity and experimental
  generalization remain separate questions with the original audit classifications.

The requested success criteria are **not met**. This checkpoint must be replaced
with actual accepted/rejected candidate evidence and final comparisons before the
task can be called complete.
