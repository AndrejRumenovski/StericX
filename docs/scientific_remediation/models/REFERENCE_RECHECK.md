# Corrected model reference recheck

The frozen v4 build passes **7,175 numerical and input-contract comparisons**, with
zero failed comparisons, including **4,400 exact comparisons of unchanged
historical numeric values**. This is a numerical-contract result, not a claim of
experimental validity, calibrated chemical uncertainty, or globally validated
science.

The complete receipt is
[reference_recheck_v4/manifest.json](../../../.stericx/scientific_remediation/models/reference_recheck_v4/manifest.json),
SHA-256 `c40d88064a65a72a663a5af512e88adff886d8c6019d4babadab9b19b0460c4f`.
It binds the frozen v4 source/build, native executable and observer, complete
observation receipt, actual CLI captures, corrected input bundle, documented
commands, current Python environment, immutable reference programs, and all raw
results. Inputs, helpers, binaries and receipts were verified before and after.
The original audit remains unchanged.

Coverage includes:

- All **33 original model CLI cases**, including successful fits/screens and the
  seven original unsuccessful outer folds. NaN, infinity and mismatched frozen
  predictions are rejected; an unspecified ranking objective remains rejected.
- All **33 original model-domain observer requests**: Student-t quantiles,
  leverage, intervals, nearest-neighbor distance, ordinary covariance distance,
  and the rank-one covariance witness.
- **20 additional observer requests**: ten preserved near-unit-alpha witnesses,
  six invalid quantile domains, and four malformed training-geometry cases.
- **Eight additional screening contract rejections**: legacy unknown aggregation,
  supplied record values substituted by geometry, multiple geometries under a
  single-geometry declaration, and missing, negative, zero, mismatched or
  nonfinite weights. Weight vectors have matching lengths when the intended
  failure is their numeric domain.
- Both corrected parse-built models and their **three JSON screens** from the
  17-command v4 tutorial receipt. Screen rows returned are 11, 2 and 11; the
  two-row diverse panel selects from three untested candidates. Record target,
  temperature and electronic values match the supplied corrected CSV at the
  declared f32 representation. Fitting and screening use the same supplied
  weighted descriptor means and 353.15 K response metadata.

The newly exposed public applicability panic is preserved in
[attempt1](../../../.stericx/scientific_remediation/models/reference_recheck_v3_attempt1/).
Its out-of-range feature index now produces an identified `unknown` assessment
with an explicit unavailable reason and no leverage, intervals or distance.
[Attempt2](../../../.stericx/scientific_remediation/models/reference_recheck_v3_attempt2/)
and [attempt3](../../../.stericx/scientific_remediation/models/reference_recheck_v3_attempt3/)
retain adapter failures in handling the original CSV label convention and NumPy
boolean serialization. They were not overwritten or presented as successful
validation runs.

The [adapter](recheck.py) copies the sealed `models_audit.py`, `models_domain.py`
and `models_screen.py` unchanged. It calls the original independent
NumPy/SciPy/sklearn equation functions on actual captured data. SUT subprocess
execution is prohibited during the reference phase. A legacy serialized model
is used only when it was the actual input to that historical screen request;
it never substitutes for a newly fitted corrected model.

The original domain reference output is retained verbatim. A separate
[domain input comparison](../../../.stericx/scientific_remediation/models/reference_recheck_v4/reference/models/results/domain_f32_input_comparison.json)
applies the public API's f32 feature conversion before the same covariance,
hat-matrix and Student-t equations. This resolves the original comparison's
approximately 1e-7 input-rounding differences without widening numeric bounds.

The complete [tolerance policy](../../../.stericx/scientific_remediation/models/reference_recheck_v4/tolerance_policy.json)
records its immutable original-audit evidence. Unchanged historical numbers and
hit identity/order require exact equality, except the intentionally corrected
ee and covariance calculations. Independent f64 calculations use 1e-12 absolute
and 2e-12 relative allowances; coefficients use two f32 ULPs; f32 fit metrics use
eight f32 epsilons absolute and relative. Original Student-t references retain
the documented 1.73e-10 residual scale under a 2e-10 absolute plus 2e-12 relative
bound. Near-unit analytic quantiles require zero absolute and 2e-12 relative
error. All actual residuals are retained, and no bound was widened to hide a
failed result.

Eight focused [adapter tests](../../../tests/test_remediation_model_recheck.py)
pass, covering tampered raw files, build substitution, missing observations,
nonfinite values, caught panics, exact historical numbers, missing numeric
fields, and reference-scalar serialization. Ruff passes.

Retained limitations include conditional feature-selection diagnostics, the
seven failed native outer folds, uncalibrated nominal intervals, a marginal
coefficient envelope without universal coverage, source/population uncertainty,
and the absence of new chemical outcome evidence. Bootstrap replicate generation
and permutation RNG are not independently reimplemented here; propagation and
quantiles use the actual retained replicates. Greedy diversity ordering remains
covered by focused tests rather than a new independent implementation. The
training-versus-screen descriptor check is an input-contract check; independent
geometry accuracy is covered by the separate geometry reference gate. The
receipt deliberately records `scientific_global_pass: false`.

A separate [coverage supplement](../../../.stericx/scientific_remediation/models/reference_coverage_v4/manifest.json)
adds **449 passing comparisons**, with zero failures, for **7,624 total checks**.
Its manifest SHA-256 is
`88176550fa4ab8520c70b0bcb72012f2ada718fd10f3f3dede74f71e990ee678`.
The [supplement adapter](supplement_recheck.py) derives means, scales, standardized
training points, row identities and dimensions directly from the packed training
records for all eight successful fitted models. Portable-model geometry must
match the fit report exactly. It also derives every eligible ID and the complete
returned order/count for the three corrected screens, including tested-candidate
exclusion and the two-member diverse panel. That panel's rank-normalized greedy
objective is checked using an independent full pairwise distance matrix. Thus
the earlier diversity-order limitation remains only for the historical screening
campaign; the corrected tutorial panel has an independent order calculation.
Both receipts retain the same scientific scope and limitations.
