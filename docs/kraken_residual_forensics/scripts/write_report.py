"""Render the forensic report from preserved measurements, without rescoring."""
from pathlib import Path

import pandas as pd

from freeze import OUT


def csv(name):
    return pd.read_csv(OUT / name, float_precision="round_trip")


def table(df, columns, labels=None):
    def fmt(value):
        if pd.isna(value):
            return "unavailable"
        if isinstance(value, str):
            return value
        return f"{value:.10g}"
    labels = labels or columns
    lines = ["| " + " | ".join(labels) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
    lines += ["| " + " | ".join(fmt(v) for v in row) + " |" for row in df[columns].itertuples(index=False, name=None)]
    return "\n".join(lines)


def main():
    assert (OUT / "historical_target/complete.json").exists()
    final = csv("historical_target/final/metrics.csv")
    baseline = csv("historical_target/baseline/metrics.csv")
    stages = csv("historical_target/stage_metrics.csv")
    native_api = csv("results/final_native/metrics.csv")
    direct = csv("results/exact_available_input_algorithm_metrics.csv")
    ceilings = csv("historical_target/conditional_ceilings.csv")
    rows = csv("historical_target/all_residual_dossiers.csv.gz")
    head = "max_delta_qvbur_min"
    h = final.set_index("metric").loc[head]
    b = baseline.set_index("metric").loc[head]
    a = native_api.set_index("metric").loc[head]
    hm = rows[rows.metric == head].sort_values("absolute_residual", ascending=False)
    ht = hm.head(10).copy()
    ht.insert(0, "rank", range(1, 11))
    direct = direct[direct.comparison == "native_vs_morfeus_same_export"]
    stage_table = table(stages[stages.metric == head],
        ["stage", "R2", "MAE", "RMSE", "maximum_absolute_error", "scientific_reason"],
        ["Stage", "R²", "MAE", "RMSE", "Max AE", "Independent reason"])
    metric_columns = ["metric", "N", "R2", "pearson_r", "MAE", "RMSE", "median_absolute_error", "maximum_absolute_error", "slope", "intercept"]
    metric_labels = ["Descriptor/reduction", "N", "R²", "Pearson r", "MAE", "RMSE", "Median AE", "Max AE", "Slope", "Intercept"]
    variance = csv("historical_target/remaining_residual_variance_components.csv").set_index("metric").loc[head]
    ref = csv("historical_target/morfeus_vs_published_metrics.csv").set_index("metric").loc[head]
    top_fraction = (hm.residual.head(3)**2).sum()/(hm.residual**2).sum()
    report = f'''# StericX ↔ Kraken residual forensics

Investigation completed on the fixed historical cohort: **1,541 ligands, 31,611 conformers, 14 native fields × four ensemble reductions = 56 comparisons**. No ligand, conformer or outlier was removed. Percent buried volume is an explicitly derived reference, not a fourteenth independently published quantity.

**Substantially closer reproduction is possible through documented input and convention changes. Universal R² ≥ 0.99999 against the historical library is not supported by the available inputs.** For the original comparison descriptor, `max_delta_qvbur_min`, historical-reference R² improves from **{b.R2:.12f} to {h.R2:.12f}**, MAE from **{b.MAE:.9g} to {h.MAE:.9g} Å³**, and RMSE from **{b.RMSE:.9g} to {h.RMSE:.9g} Å³**. This is one descriptor, not an aggregate score.

On identical **available** geometries, all 14 direct per-conformer native-versus-Morfeus comparisons exceed 0.99999; the limiting field is sampled B1. This does **not** establish identical historical inputs. The verified exact-historical-geometry subset is **N = 0**, with R² undefined. Some minimum-volume-conformer reductions amplify a one-cell volume difference into a large change in another descriptor.

Every comparison has a machine-readable dossier. Many historical causes remain **UNRESOLVABLE FROM AVAILABLE DATA**. It would be unjustified to claim that every remaining residual has a uniquely identified geometry, conformer or reference-data cause. This report distinguishes measured explanations from those unresolved alternatives.

## 1. Frozen starting point and reference identity

The accepted C2 executable, source, tests, conventions, all input SDFs, request streams, native observations, prior independent observations and reference snapshots are preserved in [frozen/manifest.json](frozen/manifest.json). Source commit: `ce97f6d98b404710eaed38e9af6cad33a0a98e7c`. Native executable SHA-256: `b233501640e4d06555a36a278581fbf9cb011962311f205c81f3fc9a83ae89e6`. Frozen manifest SHA-256: `73ee1518d8ea548662307f3b78fa96435509ef69343d0909e0eae87a6463bc69`.

The baseline is the scientifically validated current implementation, not the older pre-remediation result also preserved in `frozen/historical_1541_predictions.csv`. Its comparison requests already use a 2.28 Å center distance; the general library default remains 2.1 Å. The baseline center direction is the sum of unit donor-bond vectors, its volume density is 0.01 Å³, and its public total/near/far volumes average three planes. Raw replay confirms the frozen outputs exactly.

Two reference versions must be kept separate. The earlier campaigns in `baseline/`, `campaigns/` and `results/` score the **frozen current API**. The authoritative historical comparison in this report uses **`frozen/historical_reference.csv`**, matching the original SI wherever SI rows exist. `historical_target/` rescores every stage against that fixed historical snapshot without changing predictions. API and historical reference differ substantively for **821 and 1036**; the largest difference is 1.34503339 Å³ in ligand 821 far volume. The API-target headline R² is {a.R2:.12f}, MAE {a.MAE:.9g}, RMSE {a.RMSE:.9g}. Both sets remain complete.

The paper SI contains 1,558 ligands, whereas this requested cohort comes from the earlier 1,541-row comparison. Eight cohort IDs—724, 1057, 1058, 2059, 2062, 2063, 2064 and 2067—lack rows in that SI workbook. They remain in the full comparison using the frozen historical CSV; their precise original publication provenance is unresolved. No unavailable SI row is silently replaced or removed.

See [screened reference differences](historical_target/substantive_reference_version_differences.csv), [all snapshot comparisons](provenance_analysis/reference_snapshot_comparison.csv.gz), and [historical-target completion/erratum](historical_target/complete.json). The screening file uses 1e-7 and also flags five SI-absent IDs (2059, 2062, 2063, 2064, 2067) whose smaller differences are explained by CSV decimal precision. An independent [literal-decimal quantization audit](historical_target/reference_precision_v2/complete.json) confirms that only 821 and 1036 differ beyond the sum of half a unit in each source's final serialized decimal place. The older dossier text describing 1e-8 as the maximum reference-snapshot difference was incorrect; its compatibility label is only a descriptive small-error threshold, not a universal uncertainty bound. Also, the older snapshot-audit rows labelled percent volume carry the underlying absolute-volume reference; the historical-target comparison applies the proper percent conversion. Neither issue changes the preserved native outputs.

## 2. Definitions and descriptor-by-descriptor metrics

Residual means StericX minus Kraken. R² is `1 − SSE/SST`, using Kraken variance; it is not Pearson r squared. Slope/intercept describe the ordinary diagnostic regression of StericX on Kraken and are never used to alter predictions. Each descriptor is scored separately in its own units. Min/max span every available conformer; delta is max minus min; `vburminconf` reads the descriptor of the conformer selected by minimum total volume.

All nine requested statistics, plus signed mean error, SSE, SST and the investigative target's SSE budget, are saved for all 56 comparisons. Full baseline and final tables are appended below. [Historical stage metrics](historical_target/stage_metrics.csv), [API stage metrics](results/stage_metrics.csv), [historical final CSV](historical_target/final/metrics.csv), and [baseline CSV](historical_target/baseline/metrics.csv) retain full precision. There are zero missing calculated comparisons: 86,296 rows at every stage.

Boltzmann averages, minimum-energy conformers, electronic descriptors and other Kraken fields outside the implemented comparison are not assigned invented values. Historical corrected energies/weights are unavailable. A complete historical Boltzmann reproduction cannot be established from these data.

## 3. Ranked residuals and distributions

Every metric has a complete rank 1–1541 table with reference, prediction, signed/absolute residual and conformer IDs: [final rankings](historical_target/final/rankings/) and [baseline rankings](historical_target/baseline/rankings/). Top 10, 25, 50 and 100 membership and squared-error concentration are explicitly recorded in [final top-k evidence](historical_target/final/top_10_25_50_100.json). All top-100 rows were traced through geometry, frame, conformer extrema and independent-reference observations; this was an automated full-coverage investigation, supplemented by the focused case diagnoses below, not 5,600 separate manual molecular inspections.

For `max_delta_qvbur_min`, the final top three contribute **{100*top_fraction:.4f}%** of SSE. This concentration is a finding, not grounds for deletion. The complete baseline/API top-k concentration history is in [report_evidence/worst_case_concentration.csv](report_evidence/worst_case_concentration.csv).

{table(ht, ["rank", "molecule_id", "kraken", "stericx", "residual", "absolute_residual"], ["Rank", "Ligand", "Kraken", "StericX", "Residual", "Absolute residual"])}

![Historical-reference residual diagnostics](historical_target/final/plots/max_delta_qvbur_min.png)

Each of the 56 baseline and final plots includes a scatter with y=x, residual versus reference, histogram, and diagnostics versus molecular weight, atom count, heavy-atom count, donor H count, conformer count, bulk radius and donor extent. [Standalone family plots](family_plots_v2/) provide readable labels and every individual residual, plus baseline/final MAE by family, for the frozen API comparison. Its only substantive reference differences are the two ligands identified above. Families are declared structural proxies (donor-neighbor elements and aromatic-neighbor count), not hand-selected classes that optimize correlation.

## 4. Structured chemistry and geometry associations

The baseline error is structured. For the original API comparison, headline absolute error correlates with donor H count (Pearson 0.303) and descriptor magnitude (0.313). Raw-vector and source-grid reproduction reduces those correlations to 0.010 and 0.005. Baseline MAE is 0.781 Å³ for P–H ligands (N=9) and 1.457 Å³ for P–H₂ ligands (N=15); final MAEs are 0.00116 and 0.0181 Å³. One remaining P–H₂ case, ligand 1487, dominates that small group's error. These observations support a donor-frame convention effect; they do not justify class-specific corrections.

The final headline has weak Pearson associations with molecular weight, atom count, conformer count and bulk (absolute values ≤0.022). Spearman associations can differ and are reported, not suppressed. The pattern is now a narrow central distribution plus a few large historical reconstruction discrepancies, rather than a universal proportional bias. There is no scientific basis for rescaling.

[Stratified metrics](historical_target/final/stratified_metrics.csv), [all feature associations](historical_target/final/residual_associations.csv), [stage-by-family effects](report_evidence/stage_family_metrics.csv) and [geometric proxy associations](report_evidence/geometric_proxy_associations.csv) contain the complete analysis. Historical coordinate RMSD cannot be correlated with residuals because it is unavailable. Present-day intraensemble RMSD is a different quantity: it uses one graph isomorphism and is only an upper bound to symmetry-minimized RMSD, not proof of historical mismatch. No causal inference is made from those proxy correlations.

## 5. Geometry provenance for all 1,541 ligands

[Ligand provenance](provenance_analysis/ligand_provenance.csv) and [31,611 conformer records](provenance_analysis/conformer_provenance.csv.gz) identify the source, atom identities, donor/neighbors, API membership, input hashes, available conformer IDs/counts, measured RMSD and missing historical fields.

| Class | Meaning | Verified N | Interpretation |
| --- | --- | ---: | --- |
| A | Exact historical geometry and relevant ensemble identity | 0 | No affirmative historical mapping recovered |
| B | Demonstrably different DFT geometry of the same ligand | 0 | DFT export labels alone do not establish this |
| C | Regenerated geometry | 0 | StericX did not regenerate this cohort |
| D | Demonstrably different historical conformer | 0 | Historical IDs/coordinates unavailable |
| E | Demonstrably different historical ensemble | 0 | Ligand 369 has supporting evidence, not a retention record |
| F | Historical identity unknown/unrecoverable from available data | 1541 | Public DFT exports exist; exact historical identity unverified |

F does not mean all geometries are necessarily wrong. It means no ligand met the affirmative identity standard. Low residuals were never used to promote a ligand to A. Class F metrics equal the full-cohort metrics; A–E R² values are undefined. See [class metrics including empty classes](report_evidence/provenance_class_metrics.csv) for the API reference and the historical `provenance_class` rows in the stratified CSV.

Atom identities and source/native donor neighbor sets agree for every conformer. Native f32 coordinate conversion produces at most 2.13e-7 Å aligned RMSD, far below a new DFT geometry. The 627 initial canonical-SMILES disagreements in 44 ligands are explained by 383 bond-order/formal-valence and 244 metal-coordination representations; no nonmetal adjacency difference was found. These graph representations do not alter this explicit-geometry calculation. [Graph diagnosis](provenance_analysis/graph_diagnosis/complete.json) preserves every case.

Every one of the 31,611 API XYZ exports was retrieved and retained. Only **27 conformers, all ligand 821**, carry more precision than the four-decimal SDF export. Atom/order and rounding consistency were checked independently of descriptors before admission. Their six-decimal coordinates were used in the final profile and recalculated with all three implementations. This is a legitimate input recovery, not evidence that historical retained membership is known. [Recovery receipt](recovered_xyz/complete.json).

## 6. Conformer and source-data forensics

All current API ensemble memberships/counts match the frozen inputs. Historical energies, retention lists, selected IDs and weighting are null and explicitly marked **UNRESOLVABLE FROM AVAILABLE DATA**. The original code reads the last Gaussian optimization geometry, may reuse cached conformer data/virtual-center coordinates, removes failed or imaginary-frequency cases, and removes near duplicates using energy/RMSD and fallback property criteria. Current exports do not reconstruct that historical process.

For ligand **369**, all 84 present-day conformers remain. The lower-ID block contains 44 conformers; a later block contains 40. Fresh independent calculations give maximum total volume 139.60317 Å³ in the lower block and 149.96704 Å³ in all 84, against the historical 139.60212 Å³. Maximum far volume similarly changes from 60.35916 to 70.73035 Å³. This is strong evidence consistent with ensemble-version mismatch, but it does not prove the original retained list. No trimmed score is reported. [Fresh block dossier](report_evidence/ligand_369.json).

Ligand **1290**, conformer 54301, drives both Lmax (+0.95503 Å) and alpha-max (+1.15600°) disagreements. Native, Morfeus and the independent equations agree on the exported geometry, and the residuals exceed conservative coordinate-rounding intervals. The discrepancy therefore cannot be repaired by more angular/grid precision under the documented fixed ensemble. Which historical geometry or retained conformer produced the published extrema is unknown. The same limitation applies to prominent B5 cases such as **1805**, and B1 cases such as **560**.

The source search checked original repository trees/releases, SI archives, the updated repository, API metadata, conformer XYZ and energy endpoints. All 88 selected energy endpoints returned no per-conformer data; the larger XYZ sweep recovered precision but no historical logs. The updated repository's initial `main`-branch 404 is retained; querying its actual default branch succeeded. Public absence is not proof that authors' private archives lack the files. [Source evidence](source_availability/), [historical dependency/source check](historical_dependencies/), and the linked earlier immutable energy search record delimit what was searched.

## 7. Frame, axis and atom-order conventions

The native comparison uses the normalized sum of **unit** P–neighbor directions. Kraken's source constructs its virtual Pd using the normalized sum of **raw** P–neighbor vectors. P–H and unequal donor-bond lengths make these axes differ systematically; the maximum observed angle is 13.97°. The independently sourced raw-vector profile explains most baseline error.

For buried volume, the center is virtual Pd, P lies along negative z, and each of the three donor neighbors defines the positive-x side of the xz plane. Quadrant/octant extrema and maximum **adjacent**, not diagonal, quadrant difference span all three frames. Kraken's loop leaves total/near/far scalars from the **last** neighbor plane in original atom order. The corrected general StericX API averages the three planes to preserve permutation symmetry. Historical reproduction therefore exposes the last-plane result as an explicit audit view; production behavior is not regressed to atom-order dependence.

Native zero-plane bin assignment differs from Morfeus's strict-open regional masks. At source density 0.001, the 70-point Cartesian axis has no zero coordinate, so this distinction does not explain the cohort's source-grid regional residuals. Finite integration remains orientation sensitive; continuous-volume rotational invariance does not imply identical occupancy under arbitrary grid rotation.

Sterimol uses virtual Pd→P, with explicit radii and the source L correction. Native alignment uses its established shortest-arc transverse basis; Morfeus uses a Kabsch alignment to x. Their finite angular scans have different phase even at the same sample count. This is measured and bounded, not treated as exact equality. [All native/source frames](frame_forensics/all_frames.jsonl.gz) record the complete donor, neighbor order, axes, centers, bases, bins and reduction conventions for every conformer.

## 8. Independently sourced constants

The [convention plan](conventions/plan.json) was sealed before measuring effects (SHA-256 `aa50a95f8e3cdfa7f8955c5aaeefabac074f8a85e2f9595a7975343b391f2951`). No parameter sweep against targets was performed.

| Quantity | Reproduction convention | Independent basis |
| --- | --- | --- |
| Virtual Pd distance | 2.28 Å | Original `add_valence` source |
| Center direction | normalized raw P−neighbor sum | Original source and SI |
| Integration sphere | 3.5 Å | SI and Morfeus defaults |
| Volume radii | Bondi table, ×1.17 | SI and captured source table |
| H occupancy | excluded | SI/default; H still defines donor geometry |
| Volume density | 0.001 Å³ per nominal point | Historical Morfeus source and SI |
| Sterimol radii | Bondi; entries equal to 1.20 replaced with 1.09 Å | Original explicit code |
| Sterimol L | +0.40 Å convention | Morfeus implementation/documentation |
| Sterimol sampling | 3600 endpoint-inclusive directions | Original explicit argument and historical dependency snapshots |
| Center/axes/regions | as described above | Original code and saved frame reconstruction |
| Ensemble extrema | min, max, max−min, minimum-volume conformer | Original `read_ligand` |
| Boltzmann | historical corrected gas free energies, 298.15 K, R=0.0019872036 kcal mol⁻¹ K⁻¹ | Original code; required energies unavailable |

The SI's “mesh spacing” wording must not be confused with the implementation's volume-density parameter: at 0.001 the axis has 70 samples and actual Cartesian spacing is about 0.10145 Å. The radius table includes implementation fallbacks (e.g. B/Fe 2.0 Å), not an assertion that every entry is a directly measured Bondi value. Direct full-text access to every original radii/Sterimol paper was not recovered; source compatibility is verified, universal physical truth is not claimed.

Date-selected Morfeus commits preceding the [April 2021 preprint](https://chemrxiv.org/engage/chemrxiv/article-details/60c757f9702a9bdb7018cbd4) and [January 2022 paper](https://doi.org/10.1021/jacs.1c09718) confirm the defaults, endpoint-inclusive scan, +0.4, frame construction and relevant equations. Inspection of the journal snapshot versus 0.8.0 found annotation/index-validation changes in these paths, not a new descriptor definition. These snapshots do **not** establish Kraken's actual dependency pin, NumPy/SciPy versions or cached calculations. [Method diff](historical_dependencies/journal_to_current_methods.diff), [original Kraken source](https://github.com/the-matter-lab/kraken/blob/4eaad505c1343e6083032b4a3fda47e004e19734/conf_selection_and_DFT/PL_dft_library_201027.py), [Morfeus BV documentation](https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html), [Sterimol documentation](https://digital-chemistry-laboratory.github.io/morfeus/sterimol.html).

## 9. Identical-available-input triangulation

All 31,611 conformers were independently evaluated with Morfeus 0.8.0 on the same recovered decimal coordinates and source conventions. A separate, previously sealed mathematical implementation evaluates exact continuous B1 support candidates, analytic L/B5, determinant/signed-angle pyramidalization, and float64 union-of-balls integration. The original consolidated triangulation covers 5,476 adverse conformers, including all final API-reference top-100 extremizers. Fifteen additional cases extend coverage to **5,491**, including every baseline and final historical-reference top-100 extremizer. [Historical coverage receipt](minimal_reference_historical/complete.json).

{table(direct, ["descriptor", "N", "R2", "MAE", "RMSE", "maximum_absolute_error"], ["Direct field", "Conformers", "R² vs Morfeus", "MAE", "RMSE", "Max AE"])}

Morfeus and the minimal implementation give exactly equal finite-grid volumes on those adverse cases. Analytic L/P/alpha differences are near floating-point precision. Native volume differences are generally float representation and rare one-to-three-cell occupancy changes. B1 differs at the expected finite-scan phase scale; native-versus-exact B1 maximum AE is about 0.00552 Å on the adverse subset, and Morfeus is also a sampled approximation. B5's independent analytic value is slightly more exact than Morfeus's sampled maximum. No tool is assumed to be ground truth solely by name.

[Direct metrics](results/exact_available_input_algorithm_metrics.csv) and [ligand-level same-input metrics](results/available_input_algorithm_ligand_metrics.csv) keep these two analyses separate. It would be false to claim every aggregate same-input field exceeds 0.99999: `vburminconf` can amplify selection differences, as the next case shows. These comparisons contain no proof that current exports equal historical inputs.

## 10. Numerical integration and selection floor

The predetermined convergence ladder is 0.1 → 0.01 → 0.001 → 0.000125 → 0.000015625 Å³, with complete ensembles for 55 ligands (636 conformers, 3,180 runs), including representative donor environments and large residuals. Both frame-mean and historical-last views are retained. [All grid values](convergence/analysis/ligand_grid_ladder.csv.gz), [finest estimates and last-step sensitivity](convergence/analysis/finest_estimates.csv).

![Convergence for representative and adverse cases](convergence/analysis/headline_grid_convergence.png)

Ligand 1905's headline descriptor is approximately 3.706 at the source grid and 3.704 at the finest grid, versus published 0.91935. Ligand 163 gives 6.903 and 6.879 versus 8.40698; ligand 873 gives 1.529 and 1.539 versus 2.35119. These discrepancies survive grid refinement. Conversely, 51 of the 55 headline cases are closer to the documented source grid than to the finest estimate. Reproducing a published finite-grid number and estimating the continuous descriptor are different questions. No production grid was changed and no convergence order was fitted to the target. Last-step envelopes are empirical sensitivity estimates, not certified numerical error bars.

**Ligand 390 has an exactly localized numerical selection explanation.** In conformer 38849, one carbon-sphere boundary has squared-distance-minus-radius-squared 2.17618e-6 Å² in float64 and 0 in the native f32 calculation. One grid cell changes occupancy: 50,288 versus 50,289 of 171,712 points, each cell 0.00104590465 Å³. Native volumes for conformers 38847 and 38849 tie, selecting the first; the float64 reference selects 38849. The two conformers' B5 values differ by 2.023125 Å, although the algorithms agree on each fixed conformer's B5 within roughly 1e-6 Å. P, alpha, L and regional quantities are amplified in the same way.

The [minimal point/sphere witness](selection_boundary_v2/complete.json), actual native dump and [selection decomposition for every ligand](results/vburminconf_selection_forensics.csv.gz) identify the operation. It is finite-precision occupancy plus discontinuous argmin, not a wrong B5 formula or permission to choose whichever conformer matches Kraken. A future precision-aware ensemble-selection API could report ambiguity or use a separately specified higher-precision volume evaluator; this investigation does not silently substitute Morfeus-selected IDs into native results.

## 11. Accepted changes and independently established expected behavior

The raw-bond center, source radii, documented grid and historical scalar reduction are explicit reproduction inputs/views. They do not overwrite general scientific defaults. Recovering 27 more precise XYZ inputs is a separate source-data correction. Each full-library stage retains every member and has all 56 metrics.

The sole production change is optional `SterimolCalculator::compute_with_dummy_with_sampling`. The old 360-direction methods retain their exact arithmetic and outputs. The new method permits the independently documented 3600-direction convention. A frozen three-equal-sphere geometry has continuous B1 exactly equal to the sphere radius (1.7 Å), whereas the old scan yields 1.72617936 Å. The expected bound, established before implementation, is maximum transverse atom distance × π/(sample count−1), apart from roundoff. The adverse-set maximum angular error changes from 0.05655 to 0.00552 Å, with no bound failures, no L/B5 changes and no acceptance based on Kraken R².

This is an optional accuracy/convention capability, **not a newly confirmed defect in the documented coarse-scan API**. No new erroneous core mathematical definition was established. [Pre-edit candidate plan](candidates/angular_sampling/plan_before_edit.json), [frozen witness](candidates/angular_sampling/witness_before.json), [full-library acceptance checks](campaigns/primary_sampling_3600/acceptance_checks.json), [production implementation](../../src/geometry/sterimol.rs), [focused tests](../../tests/sterimol_sampling.rs).

## 12. R² tracked only as an outcome

For the original headline descriptor, against the frozen historical reference:

{stage_table}

Stages can affect different fields differently; unchanged headline values do not imply the other 55 comparisons are unchanged. The full stage CSV and family table report those effects. The optional 3600 scan was accepted from its mathematical witness and error bound, not a best-fit scan count. Some individual target residuals worsen after legitimate changes and remain visible.

## 13. Rejected or unsupported explanations

| Hypothesis/action | Finding |
| --- | --- |
| Universal multiplicative or additive correction | Rejected; structured frame differences and case-specific historical gaps do not justify calibration |
| Donor H atoms can be omitted from frame construction | Rejected; occupancy exclusion is distinct from donor geometry |
| Connectivity errors explain all canonical-SMILES differences | Rejected; representation differences explain the audited cases |
| More precise grid resolves the largest residuals | Rejected by the convergence ladder |
| Matching Morfeus proves exact historical geometry | Rejected; it verifies algorithms on the available export only |
| Every current API reference equals the original SI | Rejected for 821 and 1036; both versions retained |
| All differences are SDF rounding | Rejected conditionally by the analytic intervals below |
| All large discrepancies are implementation bugs | Rejected; independent same-input agreement and selection sensitivity distinguish them |
| Ligand 369's higher-ID block can be deleted | Not authorized or scientifically established; all 84 remain |
| A favorable R² proves a convention is correct | Rejected; constants and expected behavior precede outcome calculations |
| Historical software/energy/ensemble identity can be inferred from targets | Unsupported; fields remain unavailable |

## 14. Remaining error and variance attribution

For each metric the exact identity is `native − historical = (native − same-input reference) + (same-input reference − historical)`. The corresponding SSE contains a cross term. For the headline: native SSE **{variance.native_vs_published_SSE:.10g}**, same-input algorithm SSE **{variance.same_input_algorithm_SSE:.10g}**, historical-reconstruction SSE **{variance.published_reconstruction_SSE:.10g}**, cross term **{variance.cross_term:.10g}**. The target R² budget is **{h.target_0_99999_SSE_budget:.10g}**. Even the independent reference on these inputs gives only **{ref.R2:.12f}** against the historical target.

[All metric decompositions](historical_target/remaining_residual_variance_components.csv), [exclusive evidence-status SSE partitions](historical_target/remaining_evidence_class_SSE.csv), and [full API-stage covariance terms](results/all_stage_SSE_cross_terms.csv) quantify what is identifiable. Do not add component SSEs without covariance, or equate an evidence-status partition with proven causal fractions.

| Remaining category | What is established | Quantification/limit |
| --- | --- | --- |
| StericX implementation error | No new confirmed mathematical kernel bug | Cannot assign unexplained errors to a bug by elimination alone |
| Recoverable convention mismatch | Raw/unit center, radii, angular/grid and frame reduction differences measured | Per-stage prediction changes and cross terms retained |
| Numerical discretization/precision | B1 phase/error bound, finite-grid convergence, explicit ligand 390 cell | Same-input deltas; selection component separated |
| Geometry mismatch | Same-input independent agreement; some values beyond rounding bounds | Historical RMSD unavailable, exact causal SSE not identifiable |
| Conformer/ensemble mismatch | 369 ID-block evidence and selection sensitivity | Original retention list absent; no definitive historical SSE fraction |
| Unavailable historical data | Original logs/energies/retention/dependency pins missing | All 1,541 have unverified historical identity |
| Reference-data uncertainty | SI/API differences for 821/1036; eight SI-absent IDs | Explicit version-difference rows and two complete score sets |
| Unresolved | Historical alternatives cannot be uniquely separated | Every residual retains its numeric decomposition and missing-data flags |

The dossier's small-error precision label is descriptive, not a proof that the geometry is exact. Nor does being outside a rounding interval identify which missing historical input is responsible. Those distinctions are necessary for an honest attribution.

## 15. Evidence-supported ceiling and exact-geometry result

There is **no identifiable universal theoretical R² ceiling** without the missing historical inputs. The observed source-profile result is an achieved reconstruction, not a proof that no future source recovery can improve it. Exact-historical-input algorithm reproduction is **N=0, R² undefined**, recorded in [exact subset receipt](results/exact_historical_input_subset.json).

A useful **conditional** ceiling can be calculated. Assume the same complete ensemble, atom identities, donor/neighbors, documented raw-vector center/radii and nearest rounding to four-decimal SDF coordinates (±0.00005 Å per component). Conservative analytic support/determinant/angle perturbation bounds give each ligand an allowed descriptor interval; B1 additionally permits the documented finite-scan phase bound. The minimum possible residual outside that interval yields `R² ≤ 1 − Σ(distance to interval)²/SST`. This does not fit coordinates or constants to targets. The inequalities are evaluated in float64 with an outward guard; they are not a formal interval-arithmetic certificate. Full derivations/assumptions are in [geometry_bounds/plan.json](geometry_bounds/plan.json).

{table(ceilings[~ceilings.metric.str.endswith("vburminconf")], ["metric", "N", "outside_interval", "conditional_R2_upper_bound", "RMSE_lower_bound"], ["Descriptor", "N", "Outside rounding interval", "Conditional upper R²", "Lower RMSE"])}

These optimistic bounds already fall below 0.99999 for several analytic descriptors. They exclude a universal 0.99999 result from **rounding-only recovery under the documented fixed-ensemble model**. They do not exclude improvement from genuinely recovered historical geometries, membership or reference corrections. Bounds for `vburminconf` intentionally allow any ensemble member and are therefore noninformative (upper R²=1); pretending the historical selected ID is known would create a false bound. No rigorous buried-volume ceiling is claimed from the empirical grid ladder.

The answer to “is ≥0.99999 achievable with scientifically justified corrections and the available original inputs?” is therefore: **not across all 56 historical comparisons on the presently documented inputs and fixed ensemble; no evidence-supported correction closes the remaining gaps.** The exact historical inputs needed to decide a broader claim are unavailable. The final historical descriptor R² values range from {final.R2.min():.9f} to {final.R2.max():.9f}; {int((final.R2 >= .99999).sum())} of 56 exceed the investigative target. That count is not an acceptance criterion. On identical available per-conformer inputs, the algorithms already agree beyond that target.

## 16. Generality, validation and preserved failures

After the production API edit, all **322 Rust tests** and **124 Python tests** passed. Formatting, Clippy with warnings denied and rustdoc passed. The scientific oracle was rebuilt and all **13 lanes / 128,021 records** replayed with identical raw outputs, IEEE values and private frame/bin observations for the existing APIs. Coverage includes Morfeus equivalence, atom-order, rotation and translation cases, symmetric geometries, non-phosphorus donors, descriptor regression and the broader scientific audit. Option-specific tests add exact-envelope, count/input validation, non-P spherical geometry, atom permutations and rigid transforms with the mathematically justified finite-scan error bound. [Replay comparison](candidates/angular_sampling/full_replay_comparison.json).

Failures are preserved: the initial baseline CSV-parser equality assertion, the interrupted Morfeus stream, the wrong-branch 404, and the first boundary probe's missing-radius error. Their corrected/recovered successors are named separately. The CSV parser issue was checked against raw JSON/IEEE data; the interrupted stream contributed only its complete ordered prefix. No missing calculation was replaced with a fabricated value.

[REPRODUCE.md](REPRODUCE.md) describes execution, immutable outputs, dependency/source receipts and limitations. [REQUIREMENTS_AUDIT.md](REQUIREMENTS_AUDIT.md), `verification.json`, and `BUNDLE_MANIFEST.json` provide the final coverage and integrity checks. Earlier validation evidence is unchanged. Production defaults are unchanged.

## Appendix A. Final historical-reference metrics

Every row has N=1541. Volume fields use Å³; Sterimol Å; alpha degrees; P is dimensionless; percent volume uses percentage points. Printed values are rounded for readability; CSVs retain full precision.

{table(final, metric_columns, metric_labels)}

## Appendix B. Baseline historical-reference metrics

{table(baseline, metric_columns, metric_labels)}
'''
    (OUT / "KRAKEN_RESIDUAL_FORENSICS.md").write_text(report)


if __name__ == "__main__":
    main()
