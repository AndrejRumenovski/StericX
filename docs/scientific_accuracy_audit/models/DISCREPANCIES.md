# Preserved model and statistical discrepancies

No production algorithm, target or threshold was changed. A finding below can concern an incorrect implementation, misleading interpretation, conflicting primary input, or an audit-harness defect; these are explicitly distinguished. Proposed corrections are recommendations after establishing the stated evidence, not implemented changes.

## DM01 — Conditional diagnostics do not assess the complete selection workflow

1. **Failing input:** `inputs/native_ni_hda.sigpack`, the ten `native_outer_*_labels.csv` files, and the original native model. Every outcome is in `results/outer_fold_observations.json`.
2. **StericX result:** conditional fixed-feature Q² = 0.002041716065248389. The actual CLI fails outer folds 2, 4, 5, 6, 7, 8 and 9 with `no non-constant descriptor improved the intercept model`; three folds succeed. No full-pipeline CLI Q² exists.
3. **Independent result:** feature selection changes across outer folds. A separate, explicitly modified workflow using the training mean when no feature is selected has Q² = -1.5369635574172444 and RMSE = 1.0972319341075567 kcal/mol. This is not a StericX result.
4. **Reference definition:** predictive assessment of a selected workflow must withhold the response from selection and fitting; see [Cawley and Talbot](https://www.jmlr.org/papers/v11/cawley10a.html). Conditional fixed-feature LOO estimates a narrower quantity.
5. **Likely cause:** features are selected once before LOO, group LOO, nested-alpha regularization, bootstrap and permutation diagnostics. The CLI also rejects the intercept-only outcome rather than returning the reference fallback.
6. **Scientific magnitude:** seven of ten full-workflow folds cannot return a model. The available native conditional Q² is already approximately zero, so the evidence does not support this native model as a predictive replacement for the published buried-volume model.
7. **Affected claims:** M07–M13, M16, M45. The current `fixed_feature_loo` name is appropriately explicit; treating it as full-pipeline validation would be incorrect. The `nested_loo` label for regularized baselines needs its fixed-feature scope stated.
8. **Proposed correction:** preserve existing conditional diagnostics, add a separately named fully nested evaluation, and predeclare how no-feature folds are handled. Count failures. Do not retroactively substitute the reference fallback and call its score a completed StericX validation.

## DM02 — Fixed-feature permutations do not test the searched model

1. **Failing input:** the same ten-row native dataset, independent seed 192026 and 2000 paired target permutations; all null statistics are in `results/permutation_selection.json`.
2. **StericX result:** `src/model/fit.rs::permutation_test` retains the features selected from the original responses, shuffles targets and computes an add-one p-value from training R².
3. **Independent result:** with the same permutation for each comparison, fixed-feature p = 0.06746626686656672; repeated-selection p = 0.2513743128435782.
4. **Reference definition:** a permutation test of an entire response-adaptive procedure must repeat that procedure under the null. The add-one arithmetic `(extreme+1)/(B+1)` is itself correct.
5. **Likely cause:** conditioning on the winning features omits the opportunities afforded by searching multiple descriptors.
6. **Scientific magnitude:** the independent p-value changes by about 0.184, or a factor of 3.7. The result does not establish statistical significance of the full search.
7. **Affected claims:** M19–M20. This is a limitation of a broader interpretation, not a falsification of the correctly described fixed-feature calculation. The production and reference RNG sequences differ.
8. **Proposed correction:** report the conditional test as such, and use repeated selection when making a searched-workflow claim. Preserve both distributions and do not choose the more favorable one after inspection.

## DM03 — Training and geometry screening can use different predictors

1. **Failing input:** frozen Study 011 `rankings.csv`; per-row `screen_minus_fit_prediction_kcal_mol`, model artifacts and source paths. Independent aggregates are in `results/study011_recalculated.json`.
2. **StericX result:** geometry screening extracts one representative geometry, while training records may contain Boltzmann-weighted ensemble descriptors. Both flows can use the same nominal model feature names.
3. **Independent result:** maximum absolute prediction disagreement = 1.2097610849678695 kcal/mol on the retained historical rows.
4. **Reference definition:** a fitted predictor must have the same meaning when applied to a new candidate; a representative geometry is not generally equal to a weighted ensemble descriptor.
5. **Likely cause:** separate descriptor-aggregation paths. This mismatch is already disclosed in the current README and Study 011.
6. **Scientific magnitude:** approximately 1.21 kcal/mol can materially change selectivity and ranking; it is much larger than ordinary stored-coefficient rounding. The audit does not assign all historical prediction error to this one cause.
7. **Affected claims:** M50 and any screening interpretation assuming identical ensemble features.
8. **Proposed correction:** make aggregation provenance part of the model/input contract; require matching precomputed candidate descriptors or a verified matching ensemble path. First establish the scientific convention for each model, then implement it in a separate correction task.

## DM04 — Nominal intervals and the reliable label are not chemically calibrated

1. **Failing input:** all 141 frozen predictions from 47 successful Study 011 panels; ten failed panels remain in the all-57 ranking denominator.
2. **StericX result:** nominal 95% Student-t prediction intervals and applicability labels are emitted. The `trust=reliable` branch is based on range and neighbor proximity, not measured calibration.
3. **Independent result:** 92/141 = 65.24822695035462% interval coverage. The interpolation subset covers 70/108 = 64.81481481481481%, with MAE = 0.9332370675934485 kcal/mol. Top-1 accuracy across all 57 panels is 0.15789473684210525, below the random 1/3 baseline.
4. **Reference definition:** an interval's nominal coverage is conditional on assumptions and evaluated over repeated prediction outcomes; [NIST](https://www.itl.nist.gov/div898/handbook/pmd/section5/pmd512.htm) distinguishes future-observation uncertainty from fitted-mean uncertainty. Descriptor proximity is not a probability of correctness.
5. **Likely cause:** a combination of model inadequacy, small-sample selection, aggregation mismatch and possible chemical shift. The audit does not isolate those contributions.
6. **Scientific magnitude:** descriptive coverage is approximately 30 percentage points below nominal. Overlapping panels reuse eleven ligands; no independent-binomial significance or 141-independent-trial claim is made.
7. **Affected claims:** M27, M35 and M49. Current documentation already warns that AD distance is not calibrated reliability; the emitted label is stronger than that documentation.
8. **Proposed correction:** use descriptive proximity labels; report observed coverage with the dependence limitation; resolve predictor consistency and evaluate a locked model on genuinely held-out chemistry before claiming calibrated predictive reliability.

## DM05 — A marginal coefficient band is not universally conservative

1. **Failing input:** `inputs/marginal_ci_counterexample_model.json` and `inputs/marginal_ci_counterexample.csv`; 1000 stored joint coefficient replicates with disjoint 2.4% tails.
2. **StericX result:** the marginal interval-arithmetic band is `[2.63716983795166, 2.63716983795166]`. README and comments describe this construction as conservative.
3. **Independent result:** empirical joint coverage of that band is 0.904. Correct propagation of the same joint replicates yields the 95% interval `[-97.36283016204834, 102.63716983795166]`; StericX's separate joint-bootstrap implementation reproduces it.
4. **Reference definition:** separate marginal 95% intervals do not guarantee a simultaneous 95% coefficient region, hence not 95% coverage of every linear combination.
5. **Likely cause:** confusing an interval-arithmetic box with a simultaneous confidence region; coefficient dependence and marginal tail allocations matter.
6. **Scientific magnitude:** the universal coverage guarantee is mathematically false. The artificial response scale is chosen to expose the failure, not to estimate a chemistry error distribution.
7. **Affected claims:** M26 and the unqualified conservative wording in README and `src/commands/screen.rs`. M25's joint mean-response implementation is supported.
8. **Proposed correction:** call the legacy object a marginal-coefficient envelope without a universal coverage guarantee, or construct a properly justified simultaneous region. Continue preferring the clearly labeled joint bootstrap interval for the quantity it estimates.

## DM06 — Student-t tail bracketing silently truncates an accepted input

1. **Failing input:** `inputs/domain_requests.jsonl`, request `t_df1.0_a1e-08`; frozen output in `raw/domain_observer.jsonl`.
2. **StericX result:** quantile = 1048576.
3. **Independent result:** SciPy and the df=1 analytic Cauchy quantile give 63661977.23675813; absolute error = 62613401.23675813, relative error = 0.9835290067083472.
4. **Reference definition:** for a two-sided alpha with df=1, `t = cot(pi*alpha/2)`; [SciPy's t distribution](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.t.html) supplies an independent implementation.
5. **Likely cause:** the upper bracket growth stops at a fixed ceiling and the routine returns that bound without exposing failure.
6. **Scientific magnitude:** extreme accepted tails are badly wrong. The audited fixed 95% screening intervals at positive integer degrees of freedom are unaffected; their maximum tested error is 1.720721343e-10.
7. **Affected claims:** M21–M22; the failure is scoped to the public tail API, not every ordinary interval.
8. **Proposed correction:** return an explicit numerical-limit error or use a justified adaptive bracket/asymptotic calculation. Preserve and test the exact accepted tail witness before changing production code.

## DM07 — Singular covariance can be reported as ordinary Mahalanobis distance

1. **Failing input:** `singular_covariance_off_manifold` in `inputs/domain_requests.jsonl`. Training standardized points lie on `z1=z2`; the supplied normal matrix uses the fit's slope ridge floor of 1e-10.
2. **StericX result:** Mahalanobis distance = 316227.75293364684; `mahalanobis_unavailable` is null, and the range/neighbor verdict is interpolation.
3. **Independent result:** sample covariance rank is one; its ordinary inverse is undefined. The huge number is a regularization-dependent surrogate, not an ordinary covariance distance.
4. **Reference definition:** ordinary Mahalanobis distance is `sqrt((x-mu)ᵀ S^(-1)(x-mu))` for nonsingular S. A pseudoinverse or ridge alternative requires a specified interpretation.
5. **Likely cause:** deriving covariance distance from a stabilized design inverse without retaining or checking the original rank.
6. **Scientific magnitude:** the accepted public metadata/API path can produce a highly misleading distance scale and missing diagnostic. The witness was injected; native feature selection did not naturally produce this exact fit.
7. **Affected claims:** M33–M34 and any serialization contract promising a usable ordinary covariance distance for every accepted stored geometry.
8. **Proposed correction:** preserve rank/conditioning metadata, reject unsupported ordinary-distance requests, or clearly label a chosen regularized distance. Do not treat the arbitrary ridge scale as chemically calibrated.

## DM08 — Cross-coupling label inequality differs from the official notebook

1. **Failing input:** `results/crosscoupling_boundary_rows.csv`: I/162, II/88, II/179, IV/310, V/88, VII/84 and VIII/11 exactly meet the respective printed yield cutoff.
2. **StericX result:** the frozen Study 007/009 scripts use `yield > cutoff`, making all seven inactive.
3. **Independent result:** the official Threshold notebook uses `0 if yield < cutoff else 1`, making all seven active. All 72 independent classifier fits under both conventions are preserved.
4. **Reference definition:** [official Threshold notebook](https://github.com/SigmanGroup/Threshold/blob/main/Sigman_single_parameter_threshold_analysis.ipynb), explicit inclusive active classification.
5. **Likely cause:** the study reproduction replaced `>=` with `>`.
6. **Scientific magnitude:** seven target labels across six reactions change. Thresholds and MCC can change discontinuously; this is not a floating-point tolerance question. Original unrounded yields may resolve some published-table differences separately.
7. **Affected claims:** M51–M63 and the wording that the published classifiers are reproduced exactly.
8. **Proposed correction:** adopt the verified reference inequality in a separately authorized scientific correction, preserve old results and publish recomputed row-level metrics. First distinguish printed-yield rounding from original source measurements.

## DM09 — Published classifier summary discrepancies are not all StericX descriptor error

1. **Failing input:** complete primary preprint data, official full-precision descriptor CSV and fresh StericX map. See `results/crosscoupling_primary_comparison.json` and all per-row results.
2. **StericX result:** under the official inclusive convention, fresh descriptors give VII threshold 25.1200685501, Right, MCC 0.3687069949; IX threshold 30.2667427063, Right, MCC 0.7825855809.
3. **Independent result:** using the official reference descriptors already gives VII 25.1752935155, Right, MCC 0.4034371921, and IX 30.5333931222, Right, MCC 0.7216878365. Primary summary values are VII 32.43, Left and IX 28.65, Right.
4. **Reference definition:** primary SI Tables S7/S8 with the declared weighted-Gini stump, plus the official notebook. Final Science SI retrieval returned 403; the complete open preprint SI is retained.
5. **Likely cause:** unresolved source/version/subset and yield-precision differences contribute; current evidence cannot uniquely assign the discrepancy. Both independent exhaustive enumeration and sklearn agree on the audited rows.
6. **Scientific magnitude:** VII changes the active side and differs by over seven percentage points in threshold. Therefore a blanket exact-reproduction claim is not supported. No boundary point or outlier was removed to force agreement.
7. **Affected claims:** M57–M59 and the broad Study 009 fidelity interpretation; each of the twelve reactions remains separately reported.
8. **Proposed correction:** obtain the exact unrounded analysis tables and subset masks used for the final publication. Reconcile source identities before attributing the residual to StericX or to the authors, and before changing a scientific kernel.

## DM10 — Reaction-level transfer does not imply unseen-ligand validation

1. **Failing input:** all six Ni reaction tables and the independent ID intersections in `results/crosscoupling_transfer.json`.
2. **StericX result:** Study 007 labels its leave-one-reaction-out exercise out-of-sample transfer. It correctly excludes held-reaction targets from fitting.
3. **Independent result:** for II, III, IV, V and RS1, all 89 test ligand IDs also occur in training on other reactions. For I, 20/34 recur. Recomputed strict-label MCCs remain positive, 0.412445–0.540062.
4. **Reference definition:** withholding a reaction is a different evaluation target from withholding a ligand or scaffold. An evaluation must state the intended prediction population.
5. **Likely cause:** the same ligand panel is measured across reactions. This overlap is not automatically target leakage for the declared reaction-transfer task.
6. **Scientific magnitude:** the experiment supplies evidence of transfer between these measured reactions on substantially known ligand families; it cannot establish broad novelty in ligand chemistry.
7. **Affected claims:** M65 and any broader extrapolation from Study 007. The narrow leave-reaction-out arithmetic is supported.
8. **Proposed correction:** report identity/scaffold overlap and the precise validation target. Add a separately designed reaction-and-ligand holdout if claiming transfer to unseen chemistry; do not relabel the existing exercise as that stronger test.

## DM11 — Reference suitability for the single-geometry comparison remains uncertain

1. **Failing input:** Study 007 section 4, which compares eighteen independent free-ligand geometries to published Boltzmann descriptors.
2. **StericX result:** the prose calls `%Vbur(boltz)` the “fair published reference” and reports imperfect agreement: R² = 0.9735 and MAE = 0.446 percentage points. It does not state that the two quantities are mathematically equal.
3. **Independent result:** the methodological suitability of that reference choice for all eighteen ligands has not been independently established here. Separately, the audit's hypothetical assertion of universal equality is false: `sum_i p_i*d_i` is not generally one selected geometry's descriptor. That hypothetical assertion is not attributed to Study 007.
4. **Reference definition:** normalized statistical-mechanical ensemble averaging; see the independent conformer audit. A useful sensitivity comparison need not compare quantities that are mathematically identical.
5. **Likely cause or evidence gap:** the audit has not assessed the individual geometry choices and ensemble populations needed to establish the reference's suitability in each case. No source-level conflation of equality with comparison is alleged.
6. **Scientific magnitude:** the study's imperfect aggregate agreement is reported, not independently recomputed in this model audit. There is no justified new numerical correction or demonstrated implementation error from this observation alone.
7. **Affected claims:** M66, classified UNCERTAIN for methodological suitability. The earlier audit wording incorrectly attributed exact equivalence to the study and has been corrected.
8. **Proposed follow-up:** examine the eighteen geometries and reference ensembles case by case, reporting the scope of the comparison. Retain geometry provenance and populations. No production correction or retraction of an equality claim is warranted by the evidence presented here, because the study made no such claim.

## DM12 — Conflicting primary ee and ΔΔG for Ni-hDA ligand 2064

1. **Failing input:** corrected Ni-hDA SI Table S3, rendered page S22, and `sources/ni_hda_reaction_data.csv`, row 2064. Both list ee = 3%; target ΔΔG = 0.028072105 kcal/mol (SI prints 0.0281).
2. **StericX result:** the retained model uses the supplied target 0.028072105, subject to f32 record rounding.
3. **Independent result:** at 353.15 K, 3% ee gives 0.042119509923352755 kcal/mol. The stored target corresponds approximately to 2% ee. Difference = 0.014047404923352755 kcal/mol.
4. **Reference definition:** `abs(ΔΔG‡) = RT ln((100+ee)/(100-ee))` under the two-pathway kinetic assumptions. Primary experiment temperature is 80 °C.
5. **Likely cause:** an unresolved inconsistency in the primary table or its conversion. The audit cannot determine the correct experimental ee from tabular evidence alone.
6. **Scientific magnitude:** 0.0140 kcal/mol is small compared with the observed model errors, but it disproves strict internal consistency and matters for exact reproduction claims. Other retained conversion differences are <= 0.0002204 kcal/mol, consistent with the source's rounded constant.
7. **Affected claims:** M44, M67–M68 and exact target-provenance claims. This finding is UNCERTAIN about the experimentally correct value, not an accusation that StericX invented a wrong target.
8. **Proposed correction:** obtain the original measurement or author clarification; retain both primary values and provenance meanwhile. Do not silently change the training target to improve agreement.

## DM13 — Audit-harness extraction and provenance defects, corrected transparently

1. **Failing input:** the first audit parser required a ligand name on the numeric line of every primary SI row. The first acquisition manifest also preceded two later journal entries.
2. **StericX result:** no production failure is implicated; these are audit-process defects.
3. **Independent result:** visual inspection restored VII/401, XI/27 and XII/418,1088,1090,1096,1101,1102, increasing 738 rows to 746. The original journal's exact bytes were recovered with the original SHA-256; current entries have the identical eleven-record prefix.
4. **Reference definition:** the original PDF tables, not an extraction heuristic, define the data; frozen provenance must retain its original hash anchor.
5. **Likely cause:** multiline PDF cell layout, and appending download records after freezing the first phase manifest.
6. **Scientific magnitude:** the omitted rows changed some classifier metrics materially, including XII's apparent threshold. The final comparison now includes all rows. Original negative and erroneous intermediate results remain available.
7. **Affected claims:** earlier draft M51–M62 and the initial acquisition-journal manifest entry; no original source payload or frozen SUT result was overwritten.
8. **Correction performed:** audit-only parser fixed, extraction v2 frozen, all twelve fits rerun, previous results archived, and `provenance_amendment_downloads.json` added. This is not a production scientific code change.

## DM14 — A fixed feature vocabulary is not a mechanistic constraint

1. **Failing input:** every native portable fit labeled `mechanistically_constrained_ols`; for example `raw/native_ni_hda/portable.json`.
2. **StericX result:** the serialized model name implies a mechanistic constraint.
3. **Independent result:** source inspection identifies a fixed descriptor vocabulary, interaction terms, a term cap, correlation screening and forward BIC. It does not impose a reaction-specific mechanistic equation, sign constraint or kinetic consistency condition.
4. **Reference definition:** a mechanistic claim requires a justified relationship between the imposed model restrictions and the reaction mechanism. Interpretable descriptors can motivate hypotheses without enforcing a mechanism.
5. **Likely cause:** model naming conflates chemical interpretability with mathematical mechanistic restriction.
6. **Scientific magnitude:** terminology can overstate the model's scientific basis; this observation alone changes no numerical prediction.
7. **Affected claims:** M02. The observed feature algebra itself remains supported by independent calculations.
8. **Proposed correction:** name the model according to its actual restricted feature vocabulary, or document and validate explicit reaction-specific mechanistic constraints before using the stronger description.

## DM15 — Frozen evaluation accepts NaN and reports successful nonfinite metrics

1. **Failing input:** `inputs/evaluate_nan_predictions.csv`, copied from the finite synthetic holdout prediction CSV with only `Predicted_ddG_kcal_mol` replaced by `NaN`. `inputs/evaluate_edges_manifest_before_execution.json` records hashes before execution.
2. **StericX result:** frozen `evaluate` exits 0 without stderr; stdout prints `mae_kcal_mol=NaN` and `rmse_kcal_mol=NaN`. `raw/evaluate_nan/evaluation.json` writes null for MAE, RMSE, the prediction and its residual.
3. **Independent result:** NaN is not the finite response produced by the supplied model and cannot support finite accuracy metrics. Two controls, positive infinity and finite value 999, are rejected with exit 2. All cases are in `results/evaluate_edge_cases.json`.
4. **Reference definition:** the documented frozen-evaluation contract re-derives each prediction and rejects a mismatch before scoring. A nonfinite value cannot agree with a finite prediction within a numerical tolerance.
5. **Likely cause:** `(recomputed - frozen_prediction).abs() > 1e-4` evaluates false for NaN; CSV parsing accepts the token and the guard lacks an explicit finiteness check. JSON serialization converts the resulting nonfinite numbers to null.
6. **Scientific magnitude:** one malformed accepted row invalidates MAE and RMSE while the command signals success. This is a validation and reporting failure, not evidence that ordinary finite least-squares arithmetic is wrong. R² being unavailable for the one-row target set is independently expected and not the defect.
7. **Affected claims:** M69 and the frozen-evaluation integrity contract. It limits any broader claim that nonfinite scientific results are always rejected.
8. **Proposed correction:** explicitly validate both recorded and recomputed predictions as finite before the consistency comparison, then reject nonfinite scoring results. Preserve the NaN, infinity and finite-mismatch cases as independent validation witnesses in a later production fix.

## DM16 — Formula matching does not reject isomeric geometry mappings

1. **Failing input:** `inputs/formula_identity/n_propylphosphine.xyz` paired with `inputs/formula_identity/methylethylphosphine.sdf`. Both explicit graphs, including hydrogen atoms and all bonds, are frozen in `inputs/formula_identity/explicit_graphs.json`; hashes and executed RDKit version are recorded before helper execution in `inputs/formula_identity_manifest_before_execution.json`.
2. **StericX result:** the frozen Study 007 `_xyz_formula` and `_sdf_formula` definitions both return `((C,3),(H,9),(P,1))`. Thus the formula comparison at lines 503–507 does not reject the cross-isomer pairing. Same-identity positive controls also pass. The source claim says formula matching rejects isomeric structures before scoring.
3. **Independent result:** RDKit independently gives formula C3H9P for both explicit graphs but distinct canonical SMILES, `CCCP` and `CCPC`. The primary n-propylphosphine phosphorus has neighbors C/H/H; the secondary methylethylphosphine phosphorus has C/C/H. These are constitutionally different chemical identities.
4. **Reference definition:** chemical identity includes connectivity, and stereochemistry where relevant; a molecular formula records element counts. Independent graph/formula evidence uses RDKit `CalcMolFormula` and explicit adjacency lists, not StericX bond inference.
5. **Likely cause:** `_xyz_formula` and `_sdf_formula` count atom symbols only. Such a necessary composition check cannot serve as a sufficient identity or isomer check.
6. **Scientific magnitude:** a mapping to a different donor-substitution class can pass this guard. The witness quantifies no historical descriptor error and does not establish that any of the eighteen historical mappings was wrong. Other stem/ID lookup and file-presence conditions were not rerun in this helper-level test.
7. **Affected claims:** M70, the specific universal rejection claim in Study 007 section 4 and its source comment. This is separate from M66's uncertain reference-choice methodology.
8. **Proposed correction:** describe the present guard as a composition check. For identity validation, verify an independently justified graph mapping and relevant stereochemistry with preserved provenance. Do not alter historical mappings without case-specific evidence; no production code was changed during this audit.
