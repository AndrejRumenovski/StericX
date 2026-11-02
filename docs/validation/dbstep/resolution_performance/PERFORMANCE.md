# Controlled performance after DBSTEP discrepancy resolution

The confirmed quiet-window campaign passes the unchanged performance policy:
**no native or Python measurement crosses the 10% investigation threshold**.
All scientific output fingerprints agree. The frozen native baseline and
candidate are byte-identical; observed differences describe host/sample variation,
not an algorithmic speedup. This report admits performance for the measured
workloads; the [final resolution](../FINAL_RESOLUTION.md) separately decides
scientific, packaging and release readiness.

## Frozen builds and measurement scope

The baseline is the fresh build of
`f04b156f7b8b31700e1c6ab978e337b00cc734e0`. The candidate's scientific engine has
the same 53 Rust/Cargo source files and native SHA-256:

`3b9079338386a2532b17ebe844a28d72344a7bdf91d83526d4b3d0fd00438424`

Measurement sources were frozen at
`497de873d503d624d6d47b3100a04b7093bd0218`, before release-version metadata
changes. Both tested Python packages were 0.4.0.dev0, with identical installed
package bytes and extension SHA-256
`2d3c1b407c2cd651880da6c037186cb2a4c74092ffb1b5ef99e4076cc82297fa`.
Preparation receipts bind both wheels, interpreters, dependencies and diagnostic
executables. Subsequent release-version packaging requires its own validation;
these timings must not be relabelled as newly measured versioned wheels.

The [prospective protocol](../resolution_analysis/PERFORMANCE_PROTOCOL.md) was
unchanged. On an AMD Ryzen 5 5600G, Linux x86-64, the controller ran seven fixed
workloads at 1, 2, 4 and 6 Rayon threads. Affinities were CPU 2, CPUs 2–3, CPUs
0–3 and CPUs 0–5 respectively. One warmup per binary and eight alternating AB/BA
pairs per configuration produced **504 native attempts: 56 warmups and 448
measured attempts**. Every sample remains included. There were no campaign
repeats, outlier exclusions or cache flushes. Process timing includes startup,
calculation and regular-file output on the original disk.

Separate diagnostic binaries provided 28 allocation runs: seven workloads at
1/6 threads for both sides. Python campaigns followed A/B/B/A order, with eight
blocks per campaign, alternating native/bridge/public method order, 100 calls per
block and CPU 2 affinity. This gives **96 Python block attempts**, or sixteen
samples per method per side. Input was parsed before timing. Two further runs
measured one warmed six-conformer descriptor/statistics allocation call each.

## Results

The [complete 28-configuration table](final_v1/TIMING_TABLE.md) reports medians,
MAD, full ranges, within-pair ratios, CPU time, RSS and record throughput.

| Measurement | Baseline | Candidate |
|---|---:|---:|
| 10,000-file batch, six threads | 938.394 ms | 928.118 ms |
| Six-thread file throughput | 10,656.50 files/s | 10,774.50 files/s |
| Six-thread batch peak RSS, median | 26.334 MiB | 26.457 MiB |
| Single small molecule, one thread | 1.889 ms | 1.981 ms |
| Single large molecule, one thread | 1.942 ms | 1.933 ms |
| 56-conformer batch, six threads | 7.646 ms | 7.541 ms |
| 86-conformer SDF, one thread | 47.315 ms | 47.353 ms |

The largest native ratio of side medians is **1.04843**, for the approximately
2 ms small-molecule process. All within-pair medians also satisfy the unchanged
criterion. The prior approximately 10,825 files/s is historical context, not
the acceptance baseline.

| Parsed six-conformer call | Baseline median ± MAD, ms | Candidate median ± MAD, ms | Candidate/baseline |
|---|---:|---:|---:|
| Native Rust | 2.967019 ± 0.019373 | 2.960327 ± 0.008337 | 0.997744 |
| Python native bridge | 3.272706 ± 0.017516 | 3.247483 ± 0.007744 | 0.992293 |
| Python typed public API | 4.031144 ± 0.015360 | 4.012544 ± 0.030674 | 0.995386 |

The bridge adds options decoding, the GIL transition and native JSON encoding.
The public API also encodes options, parses JSON and creates typed objects. Its
median-minus-bridge difference is 0.758437 ms baseline and 0.765062 ms candidate;
this descriptive difference isolates Python conversion in the same wheel. The
native executable and portable wheel have different build targets, so their
difference is not a pure FFI overhead measurement. Full per-block values remain
in [the review](final_v1/review.json) and raw archive.

Every native timed and diagnostic output fingerprint matches. All Python
preflight results match exactly; every timed native block's returned result also
matches. The unchanged protocol checks bridge/public equality at preflight,
not every Python return value during timing.

All fourteen native diagnostic pairs have identical allocation, reallocation
and deallocation counts, with zero failures. Whole-process requested bytes
differ by only +2 to +12. These counters include different-length executable
and output paths; they do not isolate descriptor allocations. Their tiny
differences are consistent with that scope, without an allocation trace proving
each byte's cause. Six-thread peak live storage also depends on overlapping
worker lifetimes. All individual differences remain in the review.

The isolated warmed six-conformer call has exactly **871 allocations, 73
reallocations, 871 deallocations and 2,454,897 requested bytes** on each side,
with zero live-byte delta. Process-wide peak requested storage is 446,852 versus
446,853 bytes. RSS is a separate per-child operating-system measurement; the
largest native median RSS increase is 1.77%, for the one-thread database build.

## Host observations and interpretation

The user explicitly confirmed machine availability. Agent builds and scientific
jobs had finished. The recorded window was **2026-11-02 18:52:46.914471–
18:57:08.042674 UTC**, lasting 261.128 seconds. Before/after host records are
included verbatim in the archived `journal.json`.

Aggregate CPU counters show 82.41% idle and 0.36% I/O wait over the window.
Accumulated CPU pressure is 0.214% of elapsed time; I/O pressure is 0.604% some /
0.589% full, and memory pressure does not increase. All 504 native attempts have
zero major page faults. These counters include the benchmark's own activity.
The recorded `ps %CPU` values are lifetime averages and do not establish
instantaneous interference or perfect host isolation. All slow observations
remain included; none warrants changing the frozen acceptance policy or
repeating this campaign.

The [independent review](final_v1/independent/REVIEW.md) separately verifies the
schedule, raw fingerprints, timing reductions, Python block counts, allocation
comparisons and host counters. Its first reviewer-only exception is preserved
alongside the corrected audit; no benchmark measurement or criterion changed.

## Reproduction and complete retained evidence

Use the frozen measurement source commit above and the bound preparation
artifacts. New version metadata intentionally changes some live source hashes;
do not substitute those bytes into a historical replay. The exact executed
measurement and review commands were:

```sh
.venv/bin/python scripts/run_dbstep_resolution_performance.py \
  --baseline .stericx/dbstep_resolution/baseline/build \
  --candidate .stericx/dbstep_resolution/mask_build_v2 \
  --preparation .stericx/dbstep_resolution/performance_preparation_v1 \
  --python-baseline .stericx/dbstep_resolution/performance_preparation_v1/python/baseline_env/bin/python \
  --python-candidate .stericx/dbstep_resolution/advanced_v1/installs/environments/wheel-3.12/bin/python \
  --python-preparation .stericx/dbstep_resolution/performance_preparation_v1/python/receipt.json \
  --review .stericx/dbstep_resolution/admission_v2/review.json \
  --engineering docs/validation/dbstep/resolution_advanced/ci_v2/rust_final.json \
  --quiet-window-receipt .stericx/dbstep_resolution/quiet_window_20261102.json \
  --output .stericx/dbstep_resolution/performance_v1 \
  > .stericx/dbstep_resolution/performance_v1.log 2>&1
.venv/bin/python scripts/report_dbstep_resolution_performance.py \
  --input .stericx/dbstep_resolution/performance_v1 \
  --output .stericx/dbstep_resolution/performance_review_v1
```

The [raw archive](final_v1/raw_evidence.tar.gz) preserves **all 2,449 files and
659,978,708 logical bytes**: stdout/stderr, metrics, profiles, Python results,
journal, launchers and logs. Identical duplicate content uses ordinary tar
hardlinks; standard extraction restores every original path. Every archive
entry and original file was re-read against the complete
[path/hash/size/mode manifest](final_v1/raw_evidence_manifest.json). The archive
is 3,695,157 bytes, SHA-256
`517dba0dcc90f5f707100e4fdbc7c6e6825878e966b12feaf34d28802435f35d`.

The reusable archiver does not calculate or alter scientific results:

```sh
.venv/bin/python scripts/archive_dbstep_performance_evidence.py \
  --input .stericx/dbstep_resolution/performance_v1 \
  --review .stericx/dbstep_resolution/performance_review_v1 \
  --quiet-window .stericx/dbstep_resolution/quiet_window_20261102.json \
  --runner-log .stericx/dbstep_resolution/performance_v1.log \
  --output docs/validation/dbstep/resolution_performance/final_v1
```

It requires a new output directory. The separately archived `independent/`
directory retains all six independent-audit files, including the failed first
review attempt, in its own `raw_evidence.tar.gz` with a complete verified
manifest. Its report, JSON receipt and execution logs are also copied outside
that archive for reading. The immutable audit source snapshots remain inside
the archive with their original bytes. The [final evidence manifest](final_v1/manifest.json)
includes those copies and all derived reports, in addition to the archive
verification receipt. Baseline preparation, scientific admission and referenced
inputs retain their separately recorded manifests and Git identities.
