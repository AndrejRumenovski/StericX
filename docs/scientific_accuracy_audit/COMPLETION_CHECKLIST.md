# Audit coverage and evidence checklist

This checklist assesses whether the requested **investigation** was performed. It does not convert an incorrect, uncertain or experimentally untested scientific claim into a verified claim. Evidence unavailable after documented primary-source searches remains explicitly limited in the report and inventory.

| User requirement | Authoritative evidence | Scope / negative result retained |
|---|---|---|
| 1. Freeze system, executable, inputs, versions and outputs | `manifest_initial.json`, `observer/manifest.json`, component input/SUT manifests, `manifest_final.json` | Original commit/binary/source and subsequent experiments are separately frozen; full raw observations preserved |
| 2. Inventory every important claim with eight fields | `CLAIMS.md`, component claim inventories and JSON | Claims have individual statuses, sources, conventions, implementation locations, methods, confidence and limitations; no unrelated collective PASS |
| 3. Verify definitions against primary/reference sources | `geometry/METHODS.md`, `kinetics/KINETICS_CONFORMERS.md`, component source archives/manifests | Exact equations and conventions; original Verloop/Radhakrishnan full-text and final Science SI access limits disclosed |
| 4. Independent reference harness | `scripts/geometry_reference.py`, `scripts/kinetics_audit.py`, `scripts/models_*.py`, `kraken/scripts/` | Observation adapter is explicitly SUT; independent equations use NumPy/Decimal/SciPy/Morfeus rather than StericX kernels |
| 5. Sterimol validation | `geometry/descriptor_comparisons.csv`, `geometry/metrics.json`, `geometry/alignment/`, Kraken analytic/reference outputs | Small/bulky/asymmetric/H-containing/non-P/synthetic and permuted cases; L/B1/B5 separate; worst discrepancies minimized |
| 6. All buried-volume descriptors and convergence | `geometry/bin_comparisons.csv`, `geometry/convergence.csv`, `geometry/focused/reference_results.json` | Every quadrant/octant and orientation, near/far, extrema and total; default preserved; coarse/default/fine/very-fine and analytic lens |
| 7. Adversarial geometry | `geometry/inputs.json`, `supplement_inputs.json`, `focused/`, `numerical/`, `alignment/` | Transforms, atom order, representable-float boundaries, donor H, close/collinear/finite extremes; errors and NaNs/infinities retained |
| 8. Historical frame bug | `geometry/historical/`, `geometry/focused/reference_results.json`, `geometry/REPORT.md` | Six original zero-result ligands/twenty conformers, plus explicit-graph controls and primary/secondary phosphines with nonbonded contacts; independent centers and descriptor consequences |
| 9. Donor/bond inference | `geometry/connectivity.json`, focused CLI outputs, `kraken/analysis/connectivity_disagreements.json` | Explicit graphs and SDF bonds are independent of the SUT cutoff; eight real conformer disagreements and ambiguity/scope restrictions retained |
| 10. Conformers / Boltzmann | `kinetics/inputs/ensembles.json`, frozen Python observations and native `.sigpack` records, `kinetics/results/` | Manual one/two/three/shifted/extreme states, temperature, normalization, invalid/missing states and values; CREST temperature failure |
| 11. Kinetics / Eyring | `kinetics/results/kinetics.csv`, `same_selectivity_different_rates.json`, native simulate captures | 195 requests, high-precision equations/constants, zero/sign symmetry, units/ranges; absolute-vs-difference interface failure |
| 12. Full Kraken reproduction | `kraken/analysis/coverage.json`, `metrics.json`, `comparisons.csv`, `kraken/delta_interpretation/all_top20_outlier_dossiers.json`, `visuals/kraken_all_metrics/INDEX.md`, `kraken/REPORT.md` | All 1,566 IDs attempted; 31,721 available conformers; 56 metric records/plots and 1,120 outlier entries. All 280 range entries map both extrema; 77 implicated BV ligands have all 2,009 conformers independently recalculated. Missing historical selections and Boltzmann energies remain explicit |
| 13. Morfeus equivalence | Geometry reference results; `kraken/morfeus_outliers/`, `kraken/delta_outliers/`; final Kraken interpretations | Descriptor-specific controlled geometry/radii/center/density experiments; enlarged matched-default BV sample N=2,372; minimized B1/alignment/boundary cases; Morfeus limitations also reported |
| 14. Every reaction reproduction | `models/REPORT.md`, corrected primary-table extraction and reaction results | Ni-hDA and each of twelve cross-coupling reactions separately; targets, temperature, fitting, thresholds, splits and leakage limitations |
| 15. Statistical methods | Model per-claim inventory and independent math/CV/resampling outputs | R²/Q²/MAE/RMSE, OLS, LOO, bootstrap, permutation, VIF, ridge, LASSO, BIC and group validation independently assessed; selection conditioning explicit |
| 16. Uncertainty meaning/calibration | Model interval comparisons, marginal-CI counterexample, Study011 row recalculation | Parameter confidence vs observation prediction vs bootstrap boxes distinguished; 65.25% recorded coverage of nominal 95% intervals; dependence of repeated splits disclosed |
| 17. Applicability domain | `models/results/domain_comparison.json`, synthetic geometry and frozen model observer outputs | Standardization/ranges/NN/calibration/Mahalanobis/leverage; singular witness; interpolation is not proof of accuracy |
| 18. Screening | `models/results/screen_*.json`, objective/diversity/exclusion cases | Ascending/descending/negative/tied predictions, missing features, OOD and diverse selection; predictions kept separate from rank and experimental truth |
| 19. Large invariance campaign | `geometry/invariance.csv`, `geometry/inputs.json`, alignment rotation supplement | Six representatives ×100 random rotations plus X/Y/Z and translations; per-descriptor moments and maxima; scientific scale quantified |
| 20. Precision / stability | Geometry numerical/alignment/analytic outputs; kinetic IEEE-bit captures | Matched f32 inputs promoted to f64 vs full coordinates, cancellation, short axes, near-Z quaternion approximation, overflow/underflow/signed zero and boundary behavior |
| 21. Classify all claims | `CLAIMS.md` | Only the seven requested status categories; unsupported conclusions remain uncertain/out of scope |
| 22. Separate A/B/C/D correctness | Master report section11 | Prior optimization parity is not evidence for scientific implementation, methodology or experimental prediction |
| 23. Final report/table/14 sections | `SCIENTIFIC_ACCURACY_AUDIT.md` and linked component reports | Opening component table, all14 requested sections, negative findings, defensible claims and recommendations |
| 24. Audit rules / preserved failures | Initial/final manifest verification, `git diff HEAD`, component failure archives/dossiers | No production changes, threshold tuning or outlier removal; all identified discrepancies documented before any proposed production correction |

## Reproducibility gates

`scripts/audit_package.py check-original` checks the original freeze and later immutable captures against their original hashes. It also verifies that tracked production files and the commit have not changed. The model download-journal append is handled through an explicit provenance amendment with a recovered original byte-identical snapshot, not by replacing its original digest.

`scripts/replay_frozen_observations.py` re-executes complete SUT request streams and compares output bytes. The full 31,721-conformer Kraken stream and the main geometry/focused/alignment streams reproduced exactly on 21 September. This proves replay consistency, **not** scientific correctness. Independent references supply the scientific evidence separately.

The final package verification covers every listed artifact byte and fails on altered/missing files. It intentionally excludes Python bytecode caches and its own self-referential verification log. New verification/replay logs can be added without changing historical observations.

## Evidence limits that remain scientifically unresolved

- Original historical full-text access limits and the primary/preprint-final SI distinction remain explicit.
- Published historical Kraken per-conformer populations/energies were not recovered; aggregate published Boltzmann outputs cannot serve as their own independent inputs.
- Conflicting upstream Ni-hDA target/ee data require experimental provenance or author clarification. They are not silently corrected.
- Experimental conformer populations, universal chemical donor assignments, prospective predictive validity and broad calibration require additional empirical evidence. Computational checking cannot establish them by itself.

These are limitations of the conclusions, not permission to assert successful numerical or experimental validation where evidence is missing.
