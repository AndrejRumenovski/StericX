# Model/statistical scientific correction ledger

Recorded before implementation. The original audit and its expectations remain
immutable. `before_fix_manifest.json` identifies the reports, failing witnesses,
and initial sources. This is a correction task; no performance optimization is
authorized here. The unfinished `scripts/compare_current_scientific_references.py`
is withdrawn: it pins the superseded b515c4f performance baseline and is not a
validation gate for corrected science.

| Claims / discrepancy | Evidence and interpretation | Planned correction or retained limitation |
| --- | --- | --- |
| M69 / DM15 | Frozen evaluation accepts NaN and emits successful NaN/null accuracy; infinity and finite mismatch controls fail. | Reject nonfinite recorded and recomputed predictions before subtraction; reject nonfinite residuals and summaries. Preserve all three original CLI witnesses. |
| M22 / DM06 | df=1, two-sided alpha=1e-8 silently returns 1048576 instead of the analytic Cauchy quantile 63661977.23675813. | Adaptive finite bracketing, finite-argument validation, and documented unavailable result if floating-point range cannot support a bracket. No change to alpha or interval coverage target. |
| M34 / DM07 | Rank-one stored training points with a ridge-stabilized design inverse yield an unlabeled ordinary Mahalanobis value 316227.75293364684. Native selection did not generate this witness. | Compute ordinary sample-covariance distance from the stored training points; reject singular/unsupported covariance and malformed geometry. Do not substitute a ridge or pseudoinverse without a separate named definition. |
| M26 / DM05 | Marginal-coefficient envelope has 90.4% joint empirical coverage in a mathematical counterexample, disproving a universal conservative guarantee. Joint bootstrap propagation already agrees with the independent calculation. | Retain the envelope calculation, remove universal coverage language, identify it as a marginal-coefficient envelope. Do not tune percentiles or fit to held-out outcomes. |
| M35 / DM04 | `trust=reliable` reflects descriptor range and leverage, not chemical calibration. | Replace trust assertions with descriptive range/leverage labels and state the nominal interval assumptions. |
| M50 / DM03 | Representative geometry and fitted ensemble record can disagree by 1.20976108497 kcal/mol. Packed records alone do not establish thermodynamic population provenance. | Add an explicit descriptor-input contract: fitted packed descriptor values retain their supplied aggregation meaning. Geometry substitution requires a declared compatible model; otherwise require precomputed candidate descriptors with matching semantics. No inferred energies or invented Boltzmann weights. Coordinate record semantics with reaction aggregation corrections. |
| M63 / DM08 | Seven printed-yield boundary rows are inactive under study `>` but active under official notebook `>=`. | Use inclusive classification in Studies 007/009; test the seven original boundary witnesses and preserve both old and corrected row labels. No changed cutoffs or tuned class weights. |
| M70 / DM16 | CCCP and CCPC share formula C3H9P but differ in connectivity. Formula equality cannot identify an isomer. No historical mapping was shown wrong. | Describe formula checking as composition consistency only; retain the identity limitation rather than invent connectivity from coordinates or claim known mappings were wrong. |
| M02 / DM14 | Fixed vocabulary/BIC/correlation/term-cap restrictions do not enforce mechanism. | Rename the model description to fixed-feature-vocabulary OLS; preserve fitted arithmetic. |
| M13 / DM01 | `nested_loo` nests alpha/scaling, not feature selection. | Explicitly label fixed-feature nested-alpha diagnostics while retaining compatibility for reading older artifacts. |
| M08, M19 / DM01–02 | Broader audit hypotheses fail: conditional CV/permutations do not validate full selection. Seven of ten historical outer CLI folds fail. The reference intercept fallback is a different workflow. | Preserve the conditional calculations, name their scope, and retain failures. Do not introduce or claim a completed full-pipeline score or repeated-selection permutation test. |
| M17 / historical holdout | Ligand 723 shares the recorded Murcko scaffold with training ligands. Study 001 does not claim scaffold-disjoint validation. | No numerical fix; preserve the narrow historical-holdout interpretation. |
| M27, M48–49 / empirical validity | Historical nominal intervals cover 92/141 outcomes (70/108 interpolation). Overlapping panels reuse eleven ligands. Prospective outcomes were unavailable. | Negative evidence remains; no calibration, new experimental validation, or independent-trial claim follows from code corrections. |
| M52, M54, M56–59 / DM09 | Complete open-preprint rows and reference descriptors already fail to recover some published summaries. Final Science SI unavailable (403). | Retain source/version/subset/rounding uncertainty. Correct only the demonstrated inequality; do not tune thresholds, drop rows, or claim exact publication reproduction. |
| M65–66 / DM10–11 | Reaction holdouts substantially reuse ligand identities. Single-geometry versus ensemble comparison is a sensitivity comparison with uncertain suitability, not an asserted equality. | Clarify the population and comparison scope. No invented new-chemistry validation or unwarranted numerical correction. |
| M67 / DM12 | Primary ligand 2064 lists conflicting 3% ee and energy corresponding to about 2% ee. | Preserve supplied target and source conflict pending author/original-measurement clarification. |
| DM13 | Earlier audit extraction/provenance defects were corrected and intermediates preserved. | No production fix; use the final complete immutable audit evidence. |

The Student-t reference is the inverse survival function at alpha/2, consistent
with [SciPy's t distribution](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.t.html).
An ordinary covariance distance requires an inverse; a singular matrix has none
([NIST](https://www.itl.nist.gov/div898/handbook/pmc/section5/pmc532.htm)).

Focused tests will distinguish corrected failures from retained valid cases,
including nonfinite evaluation controls, analytic Cauchy tails, nonsingular and
singular covariance, malformed metadata, explicit aggregation contracts, and
inclusive yield boundaries. A later results receipt must report actual outcomes
and any unclosed gate; this ledger itself is not evidence of a successful fix.
