# Completion audit

All fourteen requested phases have supporting evidence. Scientific and controlled performance completion are separate verified claims. Raw source/binary/input/measurement hashes checked by `audit_completion.py` are recorded in [completion_audit.json](completion_audit.json). Local evidence paths below are relative to `.stericx/occupancy_precision`; tracked report filenames are relative to this directory.

| Phase | Verified requirement | Evidence |
|---|---|---|
| 0 | Frozen production source, input/compiler identities and fresh scientific baseline | BASELINE.md; baseline_manifest.json; baseline/build.json; baseline_checks.json |
| 1 | Deep occupancy profiling, assembly, counters and before flamegraph; hardware limitations explicit | OCCUPANCY_PROFILE.md; baseline_profile_events.json; baseline_counters.json; baseline_flamegraph.svg |
| 2 | Exact occupancy candidates implemented, compared, retained/rejected; masks and boundaries validated | CANDIDATES.md; EXPERIMENTS.md; o1–o5_measurements.json; final_exact_v2/summary.json |
| 3 | Continuous target defined; analytical, interval, subdivision, clipping and quadrature methods analyzed | CONTINUOUS_MATHEMATICS.md; EXPERIMENTS.md; continuous_volume/certified.rs |
| 4 | Explicit certified API, failure status and fixed tolerance results verified | c1_validation.json; molecular_validation.json; final_certified_exact/summary.json |
| 5 | Shared conservative candidate strategy compared against specialized arithmetic and all-atom ablation | chord_ablations.json; o5_measurements.json; EXPERIMENTS.md |
| 6 | Independent rational/analytical/coaxial references and extreme predicates verified | c1_validation.json; c1_coaxial/summary.json; final_predicates/summary.json; final_certified_exact/summary.json |
| 7 | Frozen real conformers, multiple radii/centers/policies, budgets and nesting evaluated | c1_molecular_broad_v3/summary.json; molecular_validation.json; final_science_review_v2/review.json |
| 8 | Isolated accepted commits, measured iterations, rejected evidence and final reprofile retained | CANDIDATES.md; EXPERIMENTS.md; final_diagnostics.json; WORK_LOG.md; Git commits in FINAL_REPORT.md |
| 9 | Controlled accuracy ladder, certified batch scaling, CPU/wall/memory/evaluations, explicit unmet targets and separate Pareto classes | final_measurements.json; certified_scaling.json; accuracy_pareto.json; accuracy_cost_pareto.svg; FINAL_TABLES.md; CERTIFIED_SCALING.md |
| 10 | Matched Morfeus/DBSTEP finite-grid and radius workflows, warm/fresh trials and scaling completed | competitive_gate_v2/summary.json; competitive_radius_gate_v2/summary.json; final_window_4/competitive; final_window_4/radius_competitive |
| 11 | Full regression, byte-exact historical replay, external references, Python and available platform checks | final_admission.json; final_source_transition.json; final_engineering_v2; final_arm_v2; final_science_review_v2 |
| 12 | Final controlled A/B stages complete under unchanged quiet gates; all failures and samples retained | final_window_3/receipt.json; final_window_4/receipt.json; timing_admission.json; final_measurements.json |
| 13 | Final scientific/engineering report, distributions, figures, commits and reproduction instructions complete | FINAL_REPORT.md; FINAL_TABLES.md; REPRODUCIBILITY.md; completion_audit.json |

The audit verifies all 13 raw historical response streams, final scientific source hashes, 695 certified outputs through their sealed admission, all eight timing-stage gates and trial counts, the additional 45-process certified-scaling supplement, linked raw metrics/stdout/stderr hashes, unchanged accuracy responses across 1,152 observations, and Warp restoration. Every accepted stage is the first qualifying stage. The rejected 630-trial CLI attempt and preflight failures remain recorded.

Documented limits remain: no native hardware-counter access; no physical ARM/macOS/Windows timing; eighteen analytical stress targets unmet; eight broad molecular targets require a larger budget; no certified regional continuous API; no equivalent external certified-continuous comparison. Quiet endpoints cannot prove zero interruption inside a trial, and all outliers remain in the distributions.

No historical report was edited. No scientific constant, acceptance tolerance, finite-grid default or release state was changed. Separate meaningful production commits and complete reproduction commands are listed in FINAL_REPORT.md and REPRODUCIBILITY.md.
