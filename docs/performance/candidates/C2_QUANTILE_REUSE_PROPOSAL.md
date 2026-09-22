# C2: reuse the model's Student-t multiplier during one screen invocation

Status: proposal only. C1 must pass all gates, be accepted and receive a fresh
profile before C2 is implemented. This document changes no production code,
tests, scientific expectations or frozen C1 source. Measurements below belong to
the corrected pre-C1 profile and must not be presented as post-C1 measurements.

## Evidence and conditional estimate

In `CORRECTED_PROFILE.md`, the corrected one-thread `screen_1000` diagnostic
attributes 20.035% of external process wall time to 1,000 calls of
`model::student_t_two_sided_quantile`. The one-thread native screen median is
256.825 ms. The ideal complete-removal ceiling is
`1 / (1 - 0.20035) = 1.25055x`. Retaining one of 1,000 calls gives the conditional
model `1 / (1 - 0.20035 + 0.20035 / 1000) = 1.25023x`, or about 20.015% less wall
time before cache overhead. A realistic hypothesis is 15–20% native screen wall
reduction (about 1.18–1.25x speedup), conditional on the repeated work remaining
invariant in the post-C1 profile. These are estimates, not measured gains.

## Current evaluation and error order

`screen_command` first validates temperature, then reads/validates the complete
portable model, hashes model/library inputs, resolves required inputs, loads the
library and applies tested-ligand exclusion. Missing library columns and ranking
option errors are resolved before inference. The candidate loop is serial and
borrows one immutable model/report throughout.

Each candidate is first checked for missing descriptors. Its expanded features
and prediction are calculated; a nonfinite prediction is excluded before any
interval calculation. Applicability, range exceedances, coefficient band and
bootstrap uncertainty follow. Leverage is calculated for the reported diagnostic,
then `TrainingGeometry::prediction_interval` separately performs:

1. `self.leverage(expanded)?`, including `design_vector`/`validate` and the current
   finite/nonnegative leverage guard;
2. `self.t_multiplier()?`, using checked positive `observations - parameters`,
   the identical `usize as f64` conversion and fixed alpha `0.05`;
3. the exact existing left-associated f64 expression
   `multiplier * residual_standard_error * (1.0 + leverage.max(0.0)).sqrt()`;
4. the current finite half-width and finite lower/upper endpoint checks.

Absent training geometry yields no interval. The existing leverage path usually
returns `None` for invalid geometry, but has the unchecked zero-parameter exception
described below; C2 must preserve its actual behavior. Saturated fits return `None`
because residual degrees
of freedom are unavailable. A nonfinite quantile also becomes `None`. Those are
unavailable intervals, not skipped candidates or grounds for an invented fallback.
The public `t_multiplier` itself does not validate all geometry, so computing it
eagerly cannot substitute for candidate-specific leverage validation.

## Proposed structure and exactness argument

Use a short-lived borrowed interval evaluator for this one immutable
`TrainingGeometry`, constructed after existing pre-loop validations. Its only
cache is `OnceCell<Option<f64>>`: uninitialized, initialized with an unavailable
multiplier, or initialized with the exact f64 multiplier. Construction performs
no calculation or validation and allocates no heap storage. A context method
first calls the original leverage path for every candidate, then lazily obtains
`geometry.t_multiplier()` from the cell. Cache initialization occurs only at the
same point where the first eligible candidate previously requested a multiplier.

A borrowed context prevents mutation of model statistics while its cache lives;
two models or two screen invocations receive separate contexts. Do not place a
cache inside the serialized/publicly mutable `TrainingGeometry`, key a global
cache by a partial model identifier, or reuse a multiplier across model lifetimes.
The current serial screen loop needs no synchronized/global cache. If a later
candidate parallelizes screening, its synchronization and lifetimes require a
separate review.

The existing public `TrainingGeometry::prediction_interval`, `confidence_interval`,
`t_multiplier` and serialized model semantics remain unchanged. An additive opaque
borrowed context API may be needed because the CLI binary calls the public library
crate. It must offer the same interval result semantics; it must not accept an
arbitrary caller-supplied multiplier that could belong to another model. If a
private arithmetic helper is shared, retain the exact operations and checks;
keep the original public method's leverage-then-multiplier evaluation order.

The original quantile is deterministic and depends only on `0.05` and the fixed
positive residual degrees of freedom. Reusing its f64 bits therefore preserves
every half-width operand. Keep all subsequent operations, candidate order,
model-validation calls, interval guards, error strings, classifications, ranking,
bootstrap bands and output emission unchanged. Do not also reuse the separately
reported leverage or precompute/reassociate the half-width product in this
candidate. No finite result may replace `None`, and no new tolerance is needed.

Lazy initialization preserves paths that never reach interval estimation. Missing
features/nonfinite predictions still exit at their original checks; absent geometry
constructs no evaluator; invalid leverage leaves its cell uninitialized. If every
candidate is excluded or has invalid leverage, there are zero quantile calls.
A saturated model caches `None` after its first valid-leverage attempt, while still
performing leverage validation on later candidates. The direct public API's
nonfinite-prediction behavior also remains unchanged: it evaluates leverage and
the multiplier before rejecting nonfinite interval bounds.

The current unchecked direct API has a compatibility edge: with zero parameters,
empty feature/scaling/inverse arrays and positive observations, `validate()` rejects
the geometry, but `design_vector()` returns an empty vector that passes the zero-
length checks in `leverage()`. It then returns `Some(0)`, so the interval method can
reach the multiplier and return an interval. Portable model loading rejects this
metadata before CLI inference. C2 must match the existing direct method for this
case without presenting the result as scientifically justified. No new upfront
`validate()` guard or scientific correction belongs in this optimization.

Cache only the multiplier's `Option<f64>`, never the result of an individual
interval attempt. A valid multiplier followed by nonfinite interval bounds must
remain cached as `Some`; later ordinary candidates may have valid intervals.
For the direct API, a nonfinite prediction still reaches the multiplier after
successful leverage. For the CLI, the existing earlier nonfinite-prediction
exclusion still prevents interval evaluation. Student-t intervals and bootstrap
coefficient uncertainty are independent output channels; absent bootstrap data
must not disable a valid Student-t interval.

Extra runtime storage is one borrowed reference and a small stack cell (a few
machine words; record `size_of` for the actual target if implemented). No vectors,
model clones, hash maps, locks or heap allocations are required. All existing
candidate-specific design-vector allocations remain; their elimination would be
a separate candidate. Risks are stale/cross-model cache values, confusing cached
`None` with uninitialized, moving validation, changing f64 grouping and accidentally
sharing prediction-interval and mean-confidence-interval formulas.

## Verification and measured-operation plan

Before any implementation, freeze accepted C1 source, tests, profile, proposal,
workload identities and Git diff in a new C2 evidence directory. Reconfirm that
`screen_1000` still performs 1,000 identical-model quantile evaluations. The
original quantile and independent equations stay unchanged.

Focused differential tests must compare f64 endpoint bits and `None` status
against the original public method, covering:

- Valid models with differing degrees of freedom; multiple predictions/leverage
  levels; zero and positive residual error; large finite predictions and
  overflowing endpoint/half-width cases.
- Invalid feature indices/dimensions, malformed inverse, nonfinite or nonpositive
  scales and invalid residual statistics: preserve the actual public validation,
  leverage and interval results. A `leverage == None` path must not initialize the
  cache; do not infer that result from `validate().is_err()` alone.
- Nonfinite expanded features and finite negative or overflowing leverage; no
  cache initialization until a later valid candidate reaches the original point.
- Saturated/underspecified fits, absent training geometry, and a cached unavailable
  multiplier; preserve the unavailable-interval path without treating a candidate
  as excluded. Cached `None` must not recompute on every valid-leverage attempt.
- Two distinct models evaluated in interleaved order within one process, plus a
  newly constructed context after model changes, proving there is no cross-model
  or stale-state reuse. Existing public methods must retain uncached semantics.
- CLI invalid-model rejection before library evaluation; all-missing/excluded
  candidates; nonfinite-prediction exclusions; JSON/text/CSV full outputs, order,
  interval presence and error messages at 1/2/4/6 thread settings.

A separate diagnostic build records interval attempts, successful leverage paths,
cache initializations, available/unavailable cache values, quantile invocations
and actual bootstrap/domain operations. For the same valid 1,000-candidate model,
the expected quantile count is 1 instead of 1,000 while interval attempts and all
candidate results remain unchanged. Reusing 999 results also avoids their 200
bisection iterations each, plus their bracketing work; record counts directly if
instrumenting these internals rather than inferring native time from iteration
counts. Zero-eligible-candidate and two-model fixtures pin counts separately.
Diagnostic counters and allocation instrumentation are excluded from native timings.

Run the full corrected observation, CLI, independent-reference and engineering
gates before native measurements. Exact output bits and preserved classifications
are mandatory. Then use the same paired AB/BA measurement protocol and frozen
workloads at 1/2/4/6 threads, without concurrent heavy work; retain every sample
and raw stdout/stderr/artifact. Report screen wall/CPU/RSS changes, paired
variation and other-workload regressions. Reconfirm the benefit exceeds noise and
has no material memory or other-workload cost before acceptance. No native speedup
is claimed by this proposal.

## Static review identities

These are the source bytes inspected while C1 validation was running:

- `src/model/domain.rs`: `cdc6b7bbb7fd8e2e25f137be943c7614b03075da91492d46b65204c223f14744`
- `src/commands/screen.rs`: `55a19624daedcbc093283f4521fe8c5cd4a2aadd995e32b469a8334bc5becb2d`
- `src/model/portable.rs`: `ab07d8d060d0cf7d794b56a8d03f64546594d19e306d603f3ee4c1b50faedfec`
- `docs/performance/CORRECTED_PROFILE.md`: `42f94341335a076536d4db55dc0e9ba9becad955f0e452374b2696016a16a3f1`

The C1 `controllers/run_candidate_science.py` was also inspected read-only. Its
build/capture and thirteen-lane observation comparison precede 93-case CLI
comparisons at every thread setting. Live source and tests are checked around
engineering validation, logs are retained, and the completion receipt explicitly
requires independent references separately. No obvious admission blocker was found
in this bounded controller review; it is not the final independent-reference or
performance admission itself.

A later read-only design review and the focused behavior matrix are retained in
[C2_TEST_MATRIX.md](C2_TEST_MATRIX.md). The pre-review proposal is archived at
`.stericx/profiling/scientifically_validated_optimization/c2_static_review_v1/`
with its SHA-256. This update records compatibility/test requirements only;
C1 acceptance and fresh profiling are still required before a C2 source freeze.
