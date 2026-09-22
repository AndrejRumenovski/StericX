# Model and statistical corrections

The demonstrated implementation and terminology failures in the pre-fix
[EVIDENCE_LEDGER.md](EVIDENCE_LEDGER.md) have focused regression coverage. This is
not a declaration of general scientific validity or a completed optimization
gate. The original audit, models, predictions, failed folds and references are
unchanged. No performance optimization was made in this work.

The [captured check receipt](verification_20260922/receipt.json) records commands,
return codes, log hashes and reviewed source hashes. All four commands exited 0:
47 model unit tests; 120 integration tests (17 model CLI, 74 screening, four
remediation, six training API, 19 portable-format); two Python classifier tests;
and Ruff for the affected study/test sources. These are focused checks, separate
from the repository-wide engineering and full independent replay gates.
The [text/CSV contract addendum](verification_20260922/text_csv_final.json)
records a further 74/74 screening pass after exposing aggregation in both
formats. The intervening `04_text_csv_contract.log` failure is retained: a
header test required the word `validation`, restored alongside the explicit
fixed-feature scope. It was a display-contract failure, not a numerical one.
[Retained evidence hashes](verification_20260922/retained_evidence.json) verify
all 11 audit references named by the pre-fix ledger remain byte-identical.

| Correction | Verified behavior and limits |
| --- | --- |
| M69 nonfinite evaluation | The original NaN, infinity and finite-mismatch CLI witnesses all fail without writing successful metrics. Nonfinite recomputed predictions and finite-input residual overflow also fail. Ordinary frozen prediction evaluation still passes. The existing finite consistency allowance is unchanged. |
| M22 Student-t clipping | Adaptive finite bracketing replaces the artificial 1,048,576 ceiling. Analytic df=1 Cauchy quantiles pass for alpha 0.05, 1e-8, 1e-20 and 1e-100 (the original alpha=1e-8 expected value is about 63,661,977.23675813). Unsupported floating-point range, exemplified by alpha=1e-300, reports numerical unavailability rather than a clipped interval. This is not an arbitrary-precision special-function implementation. |
| M34 ordinary Mahalanobis distance | Distance uses the unregularized sample covariance of stored training points. The original rank-one/ridge-stabilized case has no ordinary distance; the unavailable reason is exposed. Nonsingular hand calculations still pass. Malformed dimensions, indices, scales, inverse matrices and stored-point shapes are rejected. Missing legacy training points remain unavailable. No ridge or pseudoinverse distance is relabeled ordinary Mahalanobis. |
| M26/M35 uncertainty language | The marginal coefficient envelope calculation is retained, without guaranteed joint-coverage language. Whole-vector bootstrap propagation is unchanged. Range/leverage strings describe those measurements rather than claiming `reliable`. Nominal Student-t intervals remain conditional on the fixed linear model and IID homoscedastic normal errors. |
| M50 descriptor input contract | Schema 3 records `descriptor_aggregation`. New fits default to `supplied_record_values`, accurately leaving population provenance unestablished. Legacy files without the field report `unknown`; both accept precomputed descriptors with caller responsibility and refuse geometry substitution. Explicit `single_geometry` accepts one geometry and honors supplied row axes; partial axes fail. `supplied_weight_mean` uses every listed conformer, explicit validated weights, explicit bond-axis indices and the same weighting kernel as reaction parsing. A test proves exactly matching reported packed-descriptor values for that route. Supplied/computed mixtures are identified per descriptor, with caller provenance retained. No population is invented. |
| M63 yield boundary | All 746 immutable primary-table rows use the inclusive `>=` labels. The seven affected reaction/ligand IDs are I/162, II/88, II/179, IV/310, V/88, VII/84 and VIII/11. Twelve corrected classifiers match the immutable independent inclusive-reference threshold, accuracy, MCC and active count exactly when fed the frozen descriptors. This isolates label correction from geometry changes. It does not assert exact recovery of published summaries. |
| M02/M13/M65/M70 scope | Model description is `fixed_vocabulary_ols`. Regularized `nested_loo` retains its compatibility key with `validation_scope=fixed_feature_nested_alpha_loo`; feature selection is not nested. Study source text identifies repeated ligand identities, resubstitution scores and formula-only composition checks. No historical mapping was demonstrated incorrect. |

Historical empirical failures remain evidence: nominal intervals covered 92/141
outcomes (70/108 labeled interpolation), overlapping panels reused eleven
ligands, seven of ten historical outer CLI folds failed, and prospective
experimental outcomes were unavailable. These corrections neither calibrate
coverage nor replace conditional diagnostics with full-pipeline validation.

The final Science SI was inaccessible during the audit; source/version/subset
and printed-rounding differences remain unresolved. Ligand 2064 has conflicting
published ee/energy information and is not relabeled here. The prepared Ni-hDA
migration retains supplied historical 298.15 K populations while correcting
response-temperature metadata to 353.15 K; this does not establish those
populations at reaction conditions. Direct precomputed inputs cannot establish
their own atom-axis, geometry-preparation or population provenance.

The recorded source hashes are a checkpoint during coordinated scientific
remediation, not a frozen optimization baseline. Full current-source reference
replays and repository engineering checks must be reviewed before freezing one.

A later independent review found a second Student-t numerical defect within the
accepted domain: alpha close to one makes `df/(df+t²)` round to one, erasing the
small central probability. Ten before-fix observations are preserved in
[near_unit_alpha/before_manifest.json](near_unit_alpha/before_manifest.json).
The corrected branch evaluates the complementary central incomplete beta
directly for alpha > 0.5; the ordinary/small-alpha branch is unchanged. All
28 domain tests passed, including exact Cauchy near-unit-alpha values and
analytic central-density checks for df 2, 3, 4 and 10 at the same existing
relative bound ([focused receipt](near_unit_alpha/focused_check.json)).

The full repository run also exposed the old screening regression's implicit
geometry substitution. Its fixture now supplies the immutable packed
training descriptors explicitly. All 22 tests pass, including an explicit
legacy geometry-refusal test. The numerical frozen-prediction assertions and
tolerances are unchanged; the expected ranking is independently derived from
the retained packed values and published coefficients. Its exact range-boundary
calculation uses the declared f32 interaction arithmetic, not a newly widened
tolerance. The [addendum](verification_20260922/screening_regression_addendum.json)
preserves both stale-assumption failures and the final pass.

The [documented-command receipt](documented_commands_v3/README.md) records all
17 tutorial/model-format native commands succeeding on the pinned corrected v3
build, plus seven schema/refusal cases with expected outcomes. Both workflows
write fresh descriptors, models and decks under separate demo output roots.
The [v3 near-unit-alpha replay](near_unit_alpha/after_v3_reference.json) records
all ten observed values and analytic residuals; the maximum relative residual
is 4.64e-15. The complement identity is
[NIST DLMF 8.17.4](https://dlmf.nist.gov/8.17.E4). These receipts identify v3
explicitly and are not a declaration that later source changes have been replayed.

All 17 documented commands have now also succeeded on the final v4 snapshot
([v4 receipt](documented_commands_v4/README.md)). The independent model replay
separately checks their numerical results; runnable examples alone are not a
scientific gate. A subsequent public-API witness exposed an unchecked feature
index in `assess_applicability`, despite the individual geometry methods checking
it. The [original failure](applicability_admission/before_manifest.json) remains
preserved. Early assessment admission now returns an explicit unknown/unavailable
result for malformed metadata, with direct API regression coverage and 29 domain
tests passing. The final v4 independent replay includes all four malformed
geometry witnesses, including the out-of-bounds index, without an observer panic.
