# C2 lazy multiplier reuse: behavior gates and pre-edit freeze

This is a read-only design review and implementation checklist, not a test result
or implementation authorization. C1 acceptance and a fresh profile remain pending.
No build, scientific replay, benchmark, or production edit was performed for this
review. See the [proposal](C2_QUANTILE_REUSE_PROPOSAL.md).

The reviewed source is unchanged from the first C2 proposal:

| File | SHA-256 |
| --- | --- |
| `src/model/domain.rs` | `cdc6b7bbb7fd8e2e25f137be943c7614b03075da91492d46b65204c223f14744` |
| `src/model/portable.rs` | `ab07d8d060d0cf7d794b56a8d03f64546594d19e306d603f3ee4c1b50faedfec` |
| `src/commands/screen.rs` | `55a19624daedcbc093283f4521fe8c5cd4a2aadd995e32b469a8334bc5becb2d` |
| `src/model/mod.rs` | `7c5a419685c90fc576b3c3769d3dc12939deccc60547ff76133af931abc2ce5e` |
| `tests/cli_screen.rs` | `7d6ef0a1b3e785bfbf50d0f5d77d844d1d1da628ad63d46f031dac249e57df43` |

## Exact lazy boundary

Construct one borrowed evaluator after the existing pre-loop validations. Construction
must do no numerical work. Its only state is an initially empty
`OnceCell<Option<f64>>`, tied to the borrowed immutable `TrainingGeometry`.
For every interval attempt, call the existing `leverage(expanded)?` first; only then
initialize/reuse `geometry.t_multiplier()`. Preserve the exact f64 arithmetic and
finite endpoint guards from the original method. Cache only the multiplier, not
leverage, half-width, interval availability, applicability, or the candidate result.
Do not deduplicate the separately reported leverage call in this candidate.

The actual eligibility condition is **the original leverage returns `Some`**.
It is not `validate()` success, a finite prediction, or a successful prior interval.
The direct API currently checks prediction finiteness only through interval endpoints.
The CLI separately excludes nonfinite predictions before this point. Preserve both
paths. Context construction itself must not validate, calculate a quantile, change a
model error, or turn an unavailable interval into a skipped candidate.

A named opaque public evaluator may require a re-export from `model/mod.rs`, because
`domain` is private while the CLI is a separate crate target. Keep its fields private
and its lifetime borrowed; do not introduce serialized cache fields or global state.
Existing `prediction_interval`, `confidence_interval`, `t_multiplier`, and model
validation semantics stay unchanged. Preserve the original arithmetic independently
for the differential comparison; a shared edited helper alone is not an adequate
reference for both sides of a test.

## Direct API matrix

For every row, compare the evaluator with the **existing uncached method** using
`Option` status and both f64 endpoint `to_bits()` values. Do not use numeric tolerances.
Tests of cache initialization distinguish uninitialized, cached `None`, and cached
`Some`; isolated diagnostic counters can additionally establish kernel-call counts.

| Case | Concrete setup / sequence | Required behavior |
| --- | --- | --- |
| Ordinary repeat | Existing one-feature geometry, several predictions and near/far feature rows | Exact endpoints for every row; one multiplier initialization after the first successful leverage |
| Different df | Valid geometries with several positive residual df; interleave two live evaluators | Exact per-model multipliers/endpoints; no cross-model reuse |
| Context lifetime | Drop the first context, change observations on an owned geometry, construct a new context | Recompute for the new model; borrowing must prevent mutation while the old evaluator is live |
| Zero width | Residual error `+0.0` and `-0.0`, predictions `+0.0`/`-0.0` and ordinary finite values | Preserve endpoint bits, including signed zeros; do not substitute a special formula |
| No residual df | Structurally valid saturated `n=p` and underspecified `n<p` geometry, repeated valid leverage | Cached `None`; original intervals stay unavailable. One cache initialization but **zero Student-t kernel calls**, because `t_multiplier` rejects df first |
| Typical invalid metadata | Out-of-bounds/duplicate feature index, malformed inverse, zero/nonfinite scale, invalid residual statistics, inconsistent stored points | Match original validation/leverage/interval behavior; if original leverage is `None`, leave the cache uninitialized |
| Zero-parameter unchecked edge | Empty feature/means/scales/inverse arrays, `parameters=0`, positive observations and finite residual error | `validate()` rejects, but existing `leverage()` returns `Some(0)` and can yield an interval. Match existing behavior; do not label that interval scientifically justified |
| Selected versus unused nonfinite feature | NaN/infinity in a selected column, then in a column absent from `feature_indices` | Match the original design-vector path. A nonfinite unused column is not a blanket reason for direct-API rejection |
| Leverage arithmetic overflow | Finite very small positive scale with a large f32 feature, followed by a zero/ordinary feature | First leverage is unavailable without cache initialization; a later successful leverage initializes normally |
| Negative leverage | Symmetric positive-diagonal but indefinite inverse such as `[[1,2],[2,1]]`; feature `-1`, then `0` | Existing validation can accept the matrix; original negative quadratic form yields `None`, then the ordinary row succeeds. Preserve this path rather than strengthening validation |
| Prediction unavailable | Direct API prediction NaN, positive/negative infinity with valid leverage, then an ordinary prediction | First call still initializes an available multiplier and returns `None`; later ordinary interval succeeds |
| Candidate-specific width overflow | Valid finite residual error/geometry where a distant finite leverage overflows half-width, then a centre row with finite width | Keep cached `Some(multiplier)`; never cache the first interval's `None` for the model |
| Endpoint overflow | Valid finite width with prediction near f64 limit, then prediction zero | First bounds unavailable; second exact interval succeeds using the same available multiplier |
| Invalid row after initialization | Ordinary row, invalid selected feature, ordinary row | Middle failure does not erase/change the cached multiplier or bypass leverage checking |
| Unavailable quantile | If an established baseline fixture yields `t_multiplier()==None` with positive df, repeat it without changing the kernel | Preserve unavailable status and cache it. Do not invent a finite replacement or claim such a fixture exists without evidence; saturated/underspecified cases already exercise cached `None` |
| Existing public methods | Invoke uncached prediction and confidence methods before/after evaluator use | They remain unchanged and independent; mean confidence uses `sqrt(h)`, prediction interval uses `sqrt(1+h)` |

The zero-parameter and indefinite-inverse rows describe compatibility of unchecked
public calls. They are not new scientific validation, endorsements of those metadata,
or scientific repairs to fold into C2. CLI model loading rejects the zero-parameter
geometry before inference. A future scientific correction would be a separate change
with its own evidence and baseline decision.

## Screening and error matrix

Retain complete exit status, stdout/stderr, exclusions, interval fields, ranking,
bootstrap fields, applicability diagnostics and deck/CSV artifacts. Test only gaps
not already covered by the frozen full CLI corpus or existing integration tests.

| Path | Required result / useful fixture |
| --- | --- |
| Invalid temperature plus bad model/library | Existing temperature error wins; no model load or interval evaluation |
| Invalid model plus missing library | Existing `PortableModel::from_json`/model-load error wins before library access; finite invalid geometry metadata avoids JSON-parser ambiguity |
| Library parse, geometry aggregation, missing-column or order error | Preserve the existing pre-loop error and ordering; evaluator construction has no computation or side effects |
| `--exclude-tested` removes everything | Existing dedicated error before inference; zero interval attempts/quantile calls |
| Required column exists but one row lacks its value | Preserve `missing_descriptors`, detail and summary; excluded row does not initialize the cache |
| Finite input expands to a nonfinite prediction | Use an accepted model and finite CSV values whose f32 feature product overflows; preserve `non_finite_prediction`. A later ordinary row may initialize normally |
| Every candidate is missing/nonfinite | Preserve `no library member could be screened with this model`; do not create intervals, alter exclusions or move error handling |
| No training geometry | Every Student-t interval/leverage remains absent; bootstrap uncertainty, if present, is handled independently |
| Training geometry but no bootstrap ensemble | Student-t intervals still appear when available; bootstrap `uncertainty` stays absent |
| Structurally accepted saturated model | Predictions remain screened; Student-t bounds are absent. Ensure the whole portable fixture passes existing validation before testing interval unavailability |
| `--in-domain-only` removes all predictions | Existing error occurs **after** evaluation; do not move domain filtering ahead of the interval call to obtain zero quantile work |
| `--top` truncation | Every eligible candidate is still evaluated before ranking/truncation; output limits do not redefine cache eligibility |
| Invalid diversity weight / no geometry for diversity | These errors currently occur after candidate evaluation; do not hoist them or change which earlier error wins |
| Two screen invocations / two models | Each invocation has its own cache; no values survive from the previous model |
| JSON, text, CSV and review deck | Exact full outputs, interval presence and rankings at 1/2/4/6 thread settings, with only established metadata normalization |

Current unit coverage includes `prediction_interval_widens_with_leverage`,
`degrees_of_freedom_guard_against_saturated_fits`, and
`malformed_geometry_declines_inference_without_panicking`. They do not establish
all lazy cache state transitions above. Existing CLI tests already cover missing
values, tested exclusions, domain filtering, ranking, repeated inference and diverse
selection. The test named
`a_model_without_an_ensemble_reports_no_interval_rather_than_inventing_one` checks
**bootstrap** uncertainty, not absence of Student-t intervals; it cannot substitute
for the two independent-channel cases above.

## Pre-edit freeze plan

After root accepts C1 and records a fresh profile, create a new, overwrite-refusing
C2 directory. Do not reuse a pre-acceptance C1 executable/profile as the C2 baseline.
Before any production edit, retain:

1. The accepted C1 admission receipt, paired native/supplemental decisions, fresh
   profile and raw function-call counts. Reconfirm 1,000 eligible interval attempts
   and repeated identical-model quantiles in the screen workload.
2. Exact full Rust source, Cargo manifests/lock, tests, scripts and relevant model
   documentation; Git HEAD/status/diff and hashes. Scope the expected edits to
   `domain.rs`, the public re-export if needed, the screen call site and focused tests.
   Any additional production change requires an explicit scope review.
3. Native/diagnostic executable identities, toolchain and build flags, workload/input
   manifests, 1/2/4/6 baseline fingerprints and full prediction vectors. Preserve
   the source→build→binary links, not just separately valid hashes.
4. This reviewed proposal/matrix and the exact existing interval/leverage/quantile
   source. Retain differential expected status/bit vectors from the accepted baseline
   before changing a shared arithmetic helper. Record the zero-parameter limitation
   as existing unchecked compatibility, not as scientifically validated behavior.
5. Exact oracle/controller sources, corrected observation baseline, independent-
   reference scripts/results/classifications and allowed metadata normalization.
   Preserve previous evidence; later captures need new versioned directories.

After implementation, run the focused direct-API and screening gates first, then the
full corrected observations/private bins, all CLI cases, independent references and
engineering suite. Keep quantile equations and reference tolerances unchanged.
Separate diagnostics should count interval attempts, successful leverage, cache
initialization and Student-t kernel calls; cached `None` and zero eligible rows must
be reported explicitly. A valid 1,000-candidate model should change kernel calls
from 1,000 to 1 while preserving all scientific fields. Cache-state counters are not
native timings. Only after admission may root schedule quiet paired native trials,
including the parser/screen CPU regressions being investigated for C1.
