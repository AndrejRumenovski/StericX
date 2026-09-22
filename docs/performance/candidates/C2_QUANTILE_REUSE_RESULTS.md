# C2: reuse the model's Student-t multiplier

**Accepted for screening, with a measured parsing tradeoff.** The frozen C2
executable reduces screening time by approximately 15.5–21.7% against the actual
accepted C1b executable. Both alternating launch orders exceed the predeclared
10% reduction threshold at every thread setting. All scientific gates passed
before either native campaign.

This is the result of the [pre-edit proposal](C2_QUANTILE_REUSE_PROPOSAL.md) and
[implementation addendum](C2_IMPLEMENTATION_ADDENDUM.md). Those documents remain
unchanged. The root acceptance is `accepted_c2_v1/acceptance.json` under
`.stericx/profiling/scientifically_validated_optimization/`, SHA-256
`071f29b5e52cc4940b9581eb3de90430b3c1d7d4d2d1c4e7d4c9d30750f8cabe`.

## Exactness and eliminated work

The evaluator immutably borrows one model's training geometry and lazily retains
its Student-t multiplier, including an unavailable result. Construction performs
no numerical work. Each candidate still performs the original leverage calculation
first, then the original f64 interval arithmetic and finite-value checks. Neither
candidate intervals nor validation results are cached. The original public
uncached method remains byte-identical. No statistical method, confidence level,
model coefficient, descriptor, error or ranking convention changes.

The accepted C1b profile recorded 1,000 multiplier evaluations for a 1,000-row
screen. All eight C2 diagnostic captures recorded 1,000 interval queries, one
initializer and one multiplier evaluation. This is a measured call-count
reduction; the native timings below independently establish the wall-time gain.

The pre-edit diagnostic multiplier share was about 20.0–20.25%, giving an Amdahl
ceiling of approximately 1.25× if that component alone vanished in the measured
process. This is a profile-based model, not a strict bound on independently
sampled native wall ratios: instrumentation and scheduling affect the two
measurements. Native CPU reductions are approximately 19%.

## Scientific gates

The corrected `1ecfbd5` baseline and its independently established limits remain
the scientific contract. C2 passed:

- Exact comparisons for all 128,021 observations across 13 lanes, including the
  complete available conformer replay, geometry failures and private volume bins.
- All 93 CLI cases at 1/2/4/6 threads and 24 additional CLI cases; error output,
  statuses, packed records, predictions, rankings and exclusions retain their
  established comparison rules.
- The actual 25-command independent replay, including Morfeus, analytical
  thermodynamics/kinetics, full explicit-topology volume references, the fixed
  historical outlier union and 7,624 model checks.
- Eight complete million-value prediction streams with identical bytes.
- 319 Rust tests, 124 Python tests, four additional reference-integrity tests,
  study verification, Clippy, rustdoc, Rust formatting and Ruff checks.

The initial independent model gate retained six timestamp-derived checksum
differences. The additive resolution verified each model's own raw hash and
identical canonical model content after removing only `/created/created_utc`.
No numerical expectation or tolerance changed. Both the failed raw gate and the
successful provenance resolution are preserved. Existing scientific limitations
and unfavorable validation results remain unchanged.

The C2 build manifest SHA-256 is
`82b63f16c90e0767c5cb5b60daed61e2ce52f123d4f1a9f3156e791c5a03fd81`;
its native executable SHA-256 is
`b233501640e4d06555a36a278581fbf9cb011962311f205c81f3fc9a83ae89e6`.

## Native measurements

The cumulative campaign compares the original corrected baseline to C2 across all
ten workloads at four thread settings: 720 launches, including 320 alternating
measured pairs. A separate incremental campaign compares actual C1b to C2 for
screening, parsing and search: 600 launches, including 24 measured pairs per
configuration. Every launch and both orders are retained, with no exclusions.
Affinity, inputs, release configuration and warm-cache policy are fixed.

The following screening timings are from the incremental campaign. Speedup is
the median within-pair ratio; it need not equal the ratio of the two medians.

| Threads | C1b median (ms) | C2 median (ms) | Paired speedup | Pairs faster |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 262.303 | 203.093 | 1.277× | 24/24 |
| 2 | 255.816 | 207.419 | 1.232× | 24/24 |
| 4 | 258.169 | 208.614 | 1.239× | 24/24 |
| 6 | 255.950 | 215.838 | 1.184× | 22/24 |

Both AB and BA median reductions exceed 10% at every setting. These descriptive
statistics support the targeted improvement; they are not a formal confidence
interval or a claim of identical timing on another machine.

The full original-baseline matrix also retains C1's file-parallelism benefits.
For example, six-thread `descriptors_10000` medians are 4.968 s versus 0.904 s.
Those batch gains are cumulative and are **not attributed to the C2 cache**.

## Costs, memory and limitations

Incremental parsing is consistently about 1.35–1.67% slower, with approximately
1.4–1.7% more CPU work and 2.6–3 ms higher median wall time. Both launch orders
agree. This is a real small regression; unchanged parser source does not establish
its cause. The substantial screening gain justifies this bounded cost for the
intended screening workload.

Search wall measurements are strongly dispersed and order-sensitive. Its negative
medians and modest CPU costs remain in the evidence; no search speedup or absence
of regression is claimed. Short-workload and prediction penalties in the full
matrix are also retained. Parallel descriptor batches still cost additional total
CPU and memory. The six-thread 10,000-file peak RSS is 26.55 MiB; the full-matrix
maximum is 69.66 MiB during prediction. All launches recorded zero major faults.

The evaluator occupies 24 stack bytes and adds no heap allocation for its cache.
Matching-phase diagnostic allocation-call counts are identical. Tiny requested-byte
differences track executable/profile-path lengths, so they are not claimed as a
cache memory improvement. Incremental screening median RSS shifts by less than
0.07 MiB in either direction; this does not prove zero memory cost.

The independent combined review is
`candidate_c2_independent_review_v1/review.json`, SHA-256
`bd3527f7dd36c63ac76d57d76681a1bfe7b1c73e8e8914a059e81fb4d380aaec`.
It binds all 52 configuration distributions, medians, MAD, full ranges, resources,
launch orders, negative results and raw-output checks. Scientific agreement is
separate from predictive or experimental validity; this optimization establishes
no new chemical prediction claim.

The next step is a fresh accepted-C2 profile. Any buried-volume change requires
its own pre-edit hypothesis, scientific gates and incremental native comparison.
