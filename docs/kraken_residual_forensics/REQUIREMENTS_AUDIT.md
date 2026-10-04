# Requirement audit

“Complete with a data limit” means the requested investigation was performed and its unavailable historical input is explicitly recorded. It does not mean that an unknown cause was inferred or a missing experiment succeeded.

| User requirement | Evidence | Status / scientific limit |
| --- | --- | --- |
| 1. Freeze executable, inputs, references, outputs, conventions, hashes and commit | `frozen/manifest.json`, archives, raw predictions, baseline replay | Complete; earlier evidence retained |
| 2. Each descriptor: N, R², r, MAE, RMSE, median/max AE, slope/intercept | `historical_target/{baseline,final}/metrics.csv`, `stage_metrics.csv`; API equivalents | Complete, 56 comparisons, each N=1541 |
| 3. Rank every ligand; investigate top 10/25/50/100 | Both reference versions' `rankings/`, `top_10_25_50_100.json`, top100 dossiers, minimal-reference coverage plans | Complete automated forensic coverage plus focused diagnoses; no exclusions |
| 4. Residual scatter, y=x, reference magnitude, histogram, family and structural features | 56 baseline + 56 final historical plots, API plots, standalone `family_plots/`, association/stratification CSVs | Complete; historical RMSD unavailable, available proxies clearly labelled |
| 5. Provenance classes A–F and class metrics | `provenance_analysis/`, class-F historical stratified metrics, `report_evidence/provenance_class_metrics.csv` | Complete with data limit: verified A–E N=0, F N=1541; no residual-based classification |
| 6. Exact historical geometry test | `results/exact_historical_input_subset.json`; separate full identical-export native/Morfeus campaign | Historical N=0, not an unperformed test on recoverable exact inputs; full available-input triangulation completed |
| 7. Atom identity/connectivity, coordinates, energies, IDs/counts/membership, weights, RMSD | All-conformer provenance, graph diagnosis, source searches, XYZ recovery, 369 and 390 dossiers | Available comparisons complete; historical energies, coordinates, retained IDs and weights unavailable |
| 8. Donor/neighbors, centers, axes, ordering, quadrants/octants, atom-order invariance | `frame_forensics/all_frames.jsonl.gz`, source plan, scientific replay | Complete; historical dependency phase not asserted as exactly recovered |
| 9. Independently source constants before measuring effects | `conventions/plan.json`, frozen source, captured dependency files, SI, later date-selected source verification | Complete source-compatibility audit; original dependency pin and some original literature full texts unavailable |
| 10. Volume convergence study | `convergence/` and `convergence/analysis/` | Complete: 55 full ensembles, 636 conformers, five grids; empirical estimates not rigorous intervals |
| 11. Native/Morfeus/minimal/published triangulation | Full 31,611 Morfeus observations; 5,491 minimal cases; direct and aggregate metrics | Complete for all direct fields; all adverse extremizers covered; no automatic ground truth |
| 12. Minimize apparent implementation failures | Three-sphere B1 witness; ligand390 one-point/one-sphere occupancy witness | Complete; no newly confirmed erroneous mathematical kernel; numerical limitations distinguished |
| 13. Independently justified candidate fixes, full re-evaluation | Pre-edit sampling plan, optional API, complete profile stages and XYZ recovery | Complete; only optional API changes production code; default behavior preserved |
| 14. Track outcomes after changes, including chemical classes | Both `stage_metrics.csv`, per-family stage metrics, raw row decomposition | Complete; R² not used as acceptance criterion |
| 15. Separate full-library and exact-geometry results | Report §§5, 9, 15; exact subset receipt | Complete with data limit: exact historical N=0; available-input algorithm result separate |
| 16. Remaining categories, variance and achievable ceiling | Historical variance/cross-term tables, evidence-status SSE, conditional analytic bounds, report §§14–15 | Complete to identifiable extent; geometry vs membership vs reference causes cannot be uniquely partitioned; universal ceiling not identifiable |
| 17. Generality after accepted production change | 322 Rust/124 Python tests, Clippy/fmt/rustdoc, 13-lane/128,021-record exact scientific replay, option-specific tests | Complete; finite angular/grid invariance interpreted with justified numerical bounds |
| 18. Final scientific report and preserved machine evidence | `KRAKEN_RESIDUAL_FORENSICS.md`, `REPRODUCE.md`, `verification.json`, `BUNDLE_MANIFEST.json` | Complete after final integrity check |

The reference audit found substantive historical/API revisions for 821 and 1036. A fixed 1e-7 screening file additionally flags decimal-precision differences for five SI-absent IDs; `historical_target/reference_precision_v2/` distinguishes these using serialized decimal places. All versions and earlier inaccurate precision-note text remain visible with an explicit erratum. The final report does not silently substitute one reference for the other.

To resolve the remaining historical causes uniquely requires original optimized coordinates (or matching calculation logs), cached virtual-center/conformer data, retained conformer IDs and corrected energies, the actual dependency/environment versions, and versioned publication-to-API mappings. Those cannot be reconstructed from aggregate descriptor targets without making the forbidden inverse fit.
