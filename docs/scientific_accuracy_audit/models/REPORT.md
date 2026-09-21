# Independent reaction-model, statistics, uncertainty and screening audit

**Scope:** frozen commit `6393aafe0d983e504baf8abc1e18dc2a0f0d40e7`. No production code was changed. The frozen CLI and observation adapter are the systems under investigation; NumPy SVD, independently enumerated decision-tree cutpoints, SciPy and scikit-learn supply numerical references. This report evaluates model arithmetic and the strength of the associated scientific interpretations separately.

The ordinary model arithmetic is largely reproducible on the audited finite inputs. That does **not** establish predictive validity. The strongest negative evidence is the native Ni-hDA model's failed outer-selection folds, poor historical ranking performance, undercoverage of nominal intervals, and the mismatch between training ensembles and screened representative geometries. Exact published cross-coupling reproduction is also more limited than the studies' broad language suggests.

## Evidence and reproducibility

The original system and dependency versions are recorded in `../manifest_initial.json`; model reference versions are in `results/runtime.json`. The runs used Python 3.12.13, NumPy 2.5.1, SciPy 1.18.0 and scikit-learn 1.9.0. The web documentation read during final review may display newer patch versions; the executed versions above are authoritative.

`inputs/` contains frozen model records, candidate CSVs, metadata variants and adversarial public-API requests. `raw/` contains the unmodified CLI reports, commands, exits and stderr, including all seven failed outer fits. `results/` contains reference calculations and per-row comparisons at available precision. `manifest.json` is the original model snapshot; `manifest_finalize.json` hashes the additional phase artifacts. The latter is a **post-analysis consolidation**, not evidence that every phase hash manifest predated computation. The scripts wrote input bytes before executing each SUT case and used overwrite guards for retained inputs and raw observations. The initial system snapshot predates the audit; the first model manifest was consolidated after its first fits. These timing distinctions matter and are not hidden by rehashing.

The final replay (`results/final_replay_validation.json`) ran eight reference scripts successfully and verified that all frozen raw model observations and the original model manifest remained byte-identical. The audit helper was changed to leave an existing initial model manifest untouched; no production code was involved.

One provenance amendment is explicit: `sources/downloads.json` grew from eleven to thirteen entries after the first manifest. Exact original bytes were recovered as `sources/downloads_initial_manifest_snapshot.json`, matching the original SHA-256. `provenance_amendment_downloads.json` proves that only two acquisition entries were appended. All source payloads were retained.

An **audit parser failure** was also found and corrected. Multiline ligand names caused the first cross-coupling table extraction to omit eight numeric rows. The original extraction and all its derived results remain under `inputs/crosscoupling_preprint_rows.json` and `results/prior_extraction_v1/`. The corrected extraction is `inputs/crosscoupling_preprint_rows_v2.json`; `results/extraction_correction.json` identifies every restored row. Original PDF pages S90, S98, S99 and S100 were rendered and visually checked. The final comparison contains all **746 ligand–reaction rows: 479 Ni and 267 Pd**, with no missing fresh StericX descriptors among these rows.

Reproduce the core analysis from the repository root, preserving the existing frozen observations:

```bash
uv run --extra science python docs/scientific_accuracy_audit/scripts/models_audit.py
uv run --extra science python docs/scientific_accuracy_audit/scripts/models_screen.py
uv run --extra science python docs/scientific_accuracy_audit/scripts/models_domain.py
uv run --extra science python docs/scientific_accuracy_audit/scripts/models_extended.py
uv run --extra science python docs/scientific_accuracy_audit/scripts/models_discrepancies.py
uv run --extra science python docs/scientific_accuracy_audit/scripts/models_outer_folds.py
uv run --extra science python docs/scientific_accuracy_audit/scripts/models_evaluation_edges.py
uv run --extra science python docs/scientific_accuracy_audit/scripts/models_formula_identity.py
uv run --extra science python docs/scientific_accuracy_audit/scripts/models_reactions.py
uv run --extra science python docs/scientific_accuracy_audit/scripts/models_finalize_evidence.py
uv run --extra science python docs/scientific_accuracy_audit/scripts/models_polish_claims.py
```

These scripts reuse frozen SUT observations rather than overwriting them. A new executable comparison should be run in a new audit phase directory. The first draft inventory remains in `results/claims_draft_preserved.json`. Both `models_claims.py` and `models_polish_claims.py` now generate the corrected final inventory.

## Primary and reference sources: what was actually read

| Source | Read scope and use | Limitation |
|---|---|---|
| [Official Ni-hDA repository](https://github.com/SigmanGroup/Ni-Catalyzed-hDA) | Frozen `Enantioselectivity_Model.ipynb`, `reaction_data.csv`, `kraken_features.csv`, repository metadata; identities, feature search and target values | Repository data are primary model inputs, not proof of experimental measurement accuracy |
| [Corrected Ni-hDA supporting information](https://acs.figshare.com/articles/journal_contribution/30933433), associated correction [DOI 10.1021/jacs.5c18763](https://doi.org/10.1021/jacs.5c18763) | Modeling section and Table S3; rendered page S22 checked against eleven retained IDs | The correction letter itself returned HTTP 403. No claim of reviewing all changes, chromatograms or the entire experimental SI |
| [Newman-Stonebraker et al. primary preprint and SI](https://www.cambridge.org/engage/chemrxiv/article-details/60c758aabdbb89d828a3ade9), published as [Science 2021, DOI 10.1126/science.abj4213](https://doi.org/10.1126/science.abj4213) | SI data Tables S1/S2/S9–S12/S15/S16; classifier method and summary Tables S7/S8; key table pages rendered | Final Science SI endpoint returned 403. The open preprint is a primary source, but its identity to every final-paper datum is not established |
| [Official Threshold notebook](https://github.com/SigmanGroup/Threshold/blob/main/Sigman_single_parameter_threshold_analysis.ipynb) | Active-label inequality, class weights, depth-one tree and feature-loop code | Example CSVs are not treated as interchangeable with all twelve SI datasets |
| [Cawley and Talbot, JMLR 2010](https://www.jmlr.org/papers/v11/cawley10a.html) | Accessible abstract and evaluation-bias discussion support separating selection from assessment | This report's numerical falsification comes from explicit fresh fits, not literature authority alone |
| [NIST regression prediction intervals](https://www.itl.nist.gov/div898/handbook/pmd/section5/pmd512.htm) | Single-observation versus fitted-mean uncertainty, residual term, assumptions and coverage interpretation | The equations do not guarantee calibration under chemical distribution shift |
| [SciPy Student-t reference](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.t.html), [Ridge](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.Ridge.html), [LASSO](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.Lasso.html), [VIF](https://www.statsmodels.org/stable/generated/statsmodels.stats.outliers_influence.variance_inflation_factor.html) | Official numerical definitions and executed reference calculations | These are independent implementations, not physical ground truth |

The BIC convention is associated with [Schwarz 1978](https://doi.org/10.1214/aos/1176344136), and the bootstrap with [Efron 1979](https://doi.org/10.1214/aos/1176344552). The inaccessible original full texts were not claimed as fully reviewed. The exact operational equations audited below are stated explicitly and evaluated independently.

## Definitions and independent statistical checks

The declared feature vector is `[1, L, B1, B5, q, B1*q, B5*q, IR]`. Interaction products are rounded once in the packed f32 record contract; reference regressions then use f64 SVD. This verifies the declared algebra, not its chemical sufficiency.

For residuals `e_i = prediction_i - y_i`, the audit computes `MAE = mean(abs(e))`, `RMSE = sqrt(mean(e²))`, and `R² = 1 - sum(e²)/sum((y-mean(y))²)`. Conditional LOO uses held-out residuals in the numerator to obtain Q². The BIC reference is `n ln(RSS/n) + p ln(n)`, followed by the explicitly declared forward-selection improvement cutoff of 2. Forward BIC is not exhaustive global subset selection.

Ridge minimizes `RSS + alpha*sum(beta²)` on standardized predictors. LASSO minimizes `RSS/(2n) + alpha*sum(abs(beta))`; intercepts are unpenalized. Each reference regularization fold relearns means and population standard deviations. VIF is `1/(1-R_j²)`. Production caps a singular VIF at `1e9`, which should not be interpreted as a finite mathematical variance guarantee.

| Calculation | Independent reference | Frozen StericX | Scope |
|---|---:|---:|---|
| Native training R² | 0.362545518711 | 0.362545518711 | 10 observations, selected B5*q |
| Native fixed-feature LOO Q² | 0.002041707210 | 0.002041716065 | Selection held fixed |
| Native leave-group-out Q² | 0.001315520959 | 0.001315540471 | Supplied scaffold groups, selection held fixed |
| Native nested-alpha ridge RMSE | 0.732179928258 | 0.732179919787 | Alpha and scaling nested; selected features fixed |
| Native nested-alpha LASSO RMSE | 0.800454806613 | 0.800454796181 | Same conditioning |
| Synthetic fixed-feature LOO RMSE | 0.153649234398 | 0.153649243499 | 32 observations, two selected predictors |
| Maximum native coefficient difference | 1.163970287e-7 | — | Consistent with f32 coefficients on this design |
| Maximum bootstrap coefficient-quantile difference | 4.440892099e-16 | — | Quantiles recomputed from saved replicates |

The nominal OLS solve includes a standardized ridge floor of `1e-10`; bootstrap replicate fits use `1e-8`. Their effects are negligible on these well-conditioned examples but exact unregularized OLS is not guaranteed on singular designs. The LASSO iteration cap does not provide an independent dual-gap certificate for difficult inputs.

An adversarial frozen-prediction CSV exposes a separate **evaluation validation failure**: `evaluate` accepts `NaN` as a prediction, exits successfully, prints NaN MAE/RMSE and writes null metrics/prediction/residual to JSON. A positive infinity and a large finite mismatch are rejected. The consistency guard's subtraction comparison is bypassed by NaN. The inputs and binary were hashed before this focused test; all three outputs are retained in `results/evaluate_edge_cases.json`. Ordinary finite-input arithmetic agreement does not excuse this accepted nonfinite case.

Perturbing only held-out records' descriptors, electronics, temperature and targets left the entire training report unchanged (`results/holdout_partition.json`). This supports the explicit train/test partition. It does not fix feature selection outside the internal validation folds.

### Selection, permutations and failed folds

**The actual frozen CLI failed 7 of 10 full outer-selection Ni-hDA fits** with `no non-constant descriptor improved the intercept model`. Three fits succeeded; their selected features and predictions are preserved in `results/outer_fold_observations.json`. There is no complete StericX outer-selection LOO Q² to report.

A separate independent workflow explicitly substitutes the training mean in those seven folds. That modified workflow has Q² = **-1.5369635574172444** and RMSE = **1.0972319341075567 kcal/mol**. It is a diagnostic of selection instability, **not a StericX score**. The current code correctly calls its main diagnostic `fixed_feature_loo`; a broader claim of fully nested pipeline validation would be unsupported.

The same issue affects response permutations. With 2000 paired independent permutations, the conditional fixed-feature p-value is **0.0674663**, while repeating feature search yields **0.2513743** (`results/permutation_selection.json`). Different RNG sequences prevent treating these as an exact replay of StericX's p-value. They demonstrate why fixed-feature significance does not test the complete searched workflow. The add-one counter formula itself is correct.

## Published Ni-hDA model and historical variants

Independent fits use the official ten training IDs and the declared `vbur_max_delta_qvbur_min` predictor. Every training and LOO prediction is saved separately.

| Input/model | Training R² | Conditional Q² | LOO RMSE, kcal/mol | Historical 723 absolute error, kcal/mol |
|---|---:|---:|---:|---:|
| Official published-input model | 0.8193060694 | 0.7521012497 | 0.3429876568 | 0.3730081597 |
| Study 002 saved descriptors | 0.7693260539 | 0.6548939768 | 0.4046854518 | 0.1528764892 |
| Study 003 saved descriptors | 0.7162832021 | 0.5940804801 | 0.4388956642 | 0.1107259485 |

For the official model, the independent slope is **0.08920426624752126** and intercept **-0.19793447033801892**. Stored coefficients agree within `3.4e-16`. The conditional Q² values above are the audit/Study validation calculations, not numbers attributed to the original notebook, which reports the training fit. These arithmetic reproductions do not revalidate the original 191-descriptor search, which used the same ten responses, or establish future chemical prediction.

Ligand 723 shares a Murcko scaffold with training ligands 1057 and 1058. No identical isomeric canonical SMILES duplicates occur among the eleven records. A historical holdout is useful evidence with this limitation; it is not a scaffold-disjoint prospective experiment. The prospective Study 003 deck has no new measured outcome dataset available for this audit.

The corrected SI retains the eleven IDs and printed ΔΔG targets. It also contains a conflict for **ID 2064**: both SI and CSV say **3% ee**, while ΔΔG = **0.028072105 kcal/mol** corresponds to approximately **2% ee** at 353.15 K. Independently converting 3% yields **0.042119509923352755 kcal/mol**, a **0.014047404923352755** difference. This is conflicting primary evidence, not proof that StericX chose the wrong value. The retained input is unchanged. `results/corrected_ni_hda_targets.json` records every row. The separate kinetics audit addresses the larger 298.15 K versus 353.15 K model-metadata issue.

## Twelve cross-coupling reproductions, separately

Three predictor inputs are kept distinct: the SI's printed `%Vbur(min)` values, the official full-precision Kraken descriptor table, and the fresh frozen StericX geometry campaign. Both `yield >= cutoff` (official notebook) and the study scripts' `yield > cutoff` are evaluated. No threshold or class weight is tuned against the audit result.

All **72 independent weighted-Gini fits** agree with scikit-learn on every prediction. This verifies the independently stated classifier arithmetic. The table below uses the official `>=` active convention; it does not silently replace the frozen studies' strict-inequality results.

| Reaction | N | Primary summary threshold / direction | Independent reference-descriptor threshold / MCC | Fresh StericX-descriptor threshold / MCC |
|---|---:|---|---|---|
| I | 34 | 32.42 / Left | 32.417944 / 0.624294 | 32.528557 / 0.624294 |
| II | 89 | 32.74 / Left | 32.736501 / 0.541896 | 32.901739 / 0.541896 |
| III | 89 | 31.55 / Left | 31.553997 / 0.494451 | 31.756231 / 0.470566 |
| IV | 89 | 31.89 / Left | 31.892355 / 0.457947 | 31.937955 / 0.457947 |
| V | 89 | 51.53 / Left | 51.526975 / 0.364079 | 50.817759 / 0.364079 |
| RS1 | 89 | 31.89 / Left | 31.892355 / 0.546615 | 31.937955 / 0.546615 |
| VII | 55 | 32.43 / Left | 25.175294 / 0.403437, Right | 25.120069 / 0.368707, Right |
| VIII | 53 | 28.87 / Right | 28.867231 / 0.705276 | 28.936267 / 0.739970 |
| IX | 30 | 28.65 / Right | 30.533393 / 0.721688 | 30.266743 / 0.782586 |
| X | 28 | 30.53 / Right | 30.533393 / 0.864064 | 30.578270 / 0.864064 |
| XI | 30 | 28.82 / Right | 28.819477 / 0.792527 | 29.105011 / 0.782586 |
| XII | 71 | 29.58 / Right | 29.581509 / 0.376270 | 29.611242 / 0.376270 |

For VII and IX, the independent calculation using primary reference descriptors already differs from the primary summary table. Therefore attributing the entire discrepancy to StericX grid error would be unjustified. Potential source-version differences, printed yield rounding and source subset choices remain relevant; access to the final unrounded analysis inputs is needed to resolve them. No outliers or boundary ligands were removed.

The **strict versus inclusive inequality changes seven labels** across I, II, IV, V, VII and VIII. Every affected row is in `results/crosscoupling_boundary_rows.csv`. This is a definite implementation-convention mismatch against the official notebook, separately from uncertainty about rounded versus original yields.

The fresh descriptor comparison against the printed SI gives:

| Family | N, including repeated ligands | MAE, percentage points | RMSE | Maximum error | Slope | Intercept | 1:1 R² |
|---|---:|---:|---:|---:|---:|---:|---:|
| Ni | 479 | 0.142280804 | 0.237909964 | 0.959031677 | 0.978691344 | 0.683025464 | 0.999195775 |
| Pd | 267 | 0.117708696 | 0.170057848 | 0.706700897 | 0.983837234 | 0.538254502 | 0.999353525 |

High R² coexists with nonzero offsets and classifier boundary changes. These are descriptor-reproduction results, not generalization metrics. Median errors and every underlying row are saved in `results/crosscoupling_descriptor_fidelity.json` and `results/crosscoupling_complete_inputs.csv`.

Study 007's six leave-one-reaction-out fits withhold the held reaction's targets correctly. Independent current-descriptor MCCs with the study's strict labels are 0.540062, 0.503173, 0.412445, 0.445671, 0.474522 and 0.491403 for I, II, III, IV, V and RS1. However **all 89 test ligands appear in training in five folds**, and 20 of 34 do for I. This is reaction-level transfer on overlapping ligand families, not unseen-ligand or scaffold-disjoint validation. The exact overlap and both label conventions are recorded in `results/crosscoupling_transfer.json`.

The study bootstrap fits and scores the same resampled rows. Its intervals summarize apparent-fit variability; they do not establish held-out classifier performance. Study 007 calls the Boltzmann descriptor the “fair published reference” for its single-geometry comparison and explicitly reports imperfect agreement (R² = 0.9735, MAE = 0.446 percentage points). It does **not** claim mathematical equality. Whether that reference choice is scientifically appropriate for all eighteen cases remains **UNCERTAIN** in this audit: the historical comparison was inspected, not rerun. The separate hypothetical assertion that a single geometry always equals an ensemble mean is mathematically false; that assertion is the audit's test hypothesis, not a claim attributed to the study.

A different, explicit Study 007 claim is falsified: its formula guard cannot reject all isomeric geometry mappings. Frozen explicit graphs for primary n-propylphosphine (`CCCP`) and secondary methylethylphosphine (`CCPC`) are distinct but both have formula **C3H9P**, independently verified with RDKit. Executing the unmodified frozen `_xyz_formula` and `_sdf_formula` helpers gives identical element-count tuples for the cross-isomer pairing, so the formula comparison does not reject it. The phosphorus atoms have different heavy-neighbor counts (one versus two). This is a **helper-level counterexample**, not evidence that any of the eighteen historical mappings was wrong. Inputs, graphs, pre-execution hashes and observations are retained under `inputs/formula_identity/`, `inputs/formula_identity_manifest_before_execution.json` and `results/formula_identity_counterexample.json` (M70, DM16).

## Uncertainty and applicability domain

For a fixed full-rank model with `p` fitted coefficients, the audited formulas are:

- Residual variance: `s² = RSS/(n-p)`.
- Leverage: `h = xᵀ(XᵀX)^(-1)x`, including the intercept.
- Mean-response confidence interval: `prediction ± t*s*sqrt(h)`.
- New-observation prediction interval: `prediction ± t*s*sqrt(1+h)`.
- Full-rank sample-covariance Mahalanobis distance: `sqrt((n-1)*(h-1/n))`.

On the synthetic screening library, maximum errors are `7.91e-14` for leverage, `5.33e-15` for prediction-interval endpoints, `1.43e-13` for Mahalanobis distance, zero for nearest-neighbor distance, and `1.11e-16` for joint bootstrap endpoints. The main 95% Student-t quantile is accurate to `1.73e-10` on the tested degrees of freedom.

Two preserved edge failures limit broader numerical claims. The public Student-t inverse silently caps an extreme accepted tail: df = 1 and alpha = `1e-8` returns **1048576**, whereas the Cauchy/SciPy reference is **63661977.23675813**. An accepted stored geometry with singular covariance reports an unlabeled ridge-dependent Mahalanobis surrogate of **316227.75293364684**. The latter was injected through the public metadata API; it was not naturally produced by native forward selection.

The bootstrap method is correctly labeled **mean-response coefficient uncertainty** in the current README and serialized method. It excludes residual variation and does not quantify geometry, feature selection, mechanism changes or all predictive uncertainty. In contrast, the older marginal-coefficient band is described as universally conservative. A retained 1000-replicate counterexample gives only **90.4% empirical joint coverage** for that band; its joint propagated 95% interval remains broad. The universal wording is incorrect, though the current joint propagation itself works.

Applicability uses standardized selected-feature ranges and nearest training distances. Training points have distance zero; extreme extrapolations are detected. A point inside every range but outside the training convex hull can still be labeled interpolation. The maximum-neighbor threshold is a policy, not a calibrated probability. Singular covariance and the input precision limits above must be considered separately.

Observed chemical calibration is poor: Study 011 covers **92/141 = 65.248%** of outcomes with nominal 95% intervals; the interpolation subset covers **70/108 = 64.815%**. Those 141 predictions reuse eleven ligands, so they are not independent calibration trials. The evidence rejects a strong reliability interpretation of the emitted `trust=reliable` label. It does not establish which share of the error comes from feature selection, chemical inadequacy or the descriptor-aggregation mismatch.

## Screening and predictive interpretation

Independent synthetic tests confirm ascending, descending and absolute-magnitude objectives, negative predictions, deterministic ties, missing-descriptor exclusions, prior-tested identifier exclusions, explicit OOD filtering and visible unfiltered OOD candidates. An unspecified objective produces an explicit error. Diversity changes ordering but preserves every per-ID prediction exactly; weight zero recovers ordinary ordering. These observations validate the tested mathematical contracts, not experimental effectiveness.

The frozen Study 011 historical panel evidence remains negative. Recomputed all-57 top-1 accuracy is **0.1578947** versus a random **1/3**, and top-2 set overlap is **0.4649123** versus random **2/3**. Ten failed fits remain in the denominator. Mean Spearman rho over successful panels is **-0.3404255**. The largest difference between screening and fitting predictions is **1.20976108497 kcal/mol**, associated with the documented representative-geometry versus ensemble-feature mismatch.

It is safe to claim reproducible conditional linear-model diagnostics and tested screening ordering behavior. It is not safe to claim broad mechanistic constraint from a fixed feature vocabulary, complete nested validation from conditional diagnostics, calibrated chemical reliability from an AD label, exact reproduction of every published cross-coupling classifier, or prospective experimental success without new measurements.

The individual classifications are in `../claims_models.md` and `results/claims.json`. Significant failing witnesses and proposed evidence-backed corrections are in `DISCREPANCIES.md`. Production algorithms and scientific targets remain untouched.
