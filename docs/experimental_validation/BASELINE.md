# Frozen scientific baseline

Campaign opened 2026-11-24. Production HEAD is
`e0f1dc79d23bf6638614e85cfc0c9f2282b6b6eb`. The last implementation identified by
the preceding campaign is `aeb3ca747ffbefa5003854e5709f2efa0ba0905e`.
No release, publication, GUI change, or alteration of historical evidence is authorized
by this campaign. The pre-existing untracked `docs/media/` directory is excluded.

`baseline_manifest.json` records SHA-256 and byte counts for source, tests, data,
reference receipts and reports; toolchain, installed Python distributions, CPU,
kernel and relevant environment variables. `.stericx/experimental_validation/baseline/`
contains the complete Git index identities and a Git archive of production source.
Historical response streams and analytical requests have separate hashes. Git
objects preserve all other tracked historical artifacts. `freeze_baseline.py`
reproduces the capture in a new empty directory. An initial attempt failed because
Valgrind is absent; its controller, source archive and traceback remain under
`baseline_attempt1`. The corrected capture records unavailable tools explicitly.

## Scientific inventory

| Family | Implemented baseline | Important convention |
|---|---|---|
| Sterimol | `geometry/sterimol.rs`: sampled L/B1/B5, atom/bond/dummy axes; `continuous_sterimol.rs`: finite support-envelope events | Ordinary coordinates/results are f32; explicit continuous method uses f64 arithmetic. Raw L differs from the optional +0.40 Å convention. Continuous B1 has no universal floating-point certificate. |
| Hard-sphere volume | `buried_volume.rs`: prepared lattice, scalar/AVX2 membership, total and three donor-plane orientations | Inclusive squared distance in preserved f32 operation order. Atomic union, no additive-overlap approximation. Default R=3.5 Å, density=0.01 Å³, Bondi scale=1.17, H excluded. |
| Continuous volume | `continuous_volume.rs`: midpoint slices and outward-rounded adaptive 3D boxes; `continuous_volume/certified.rs`: adaptive XY chord integration with proven smooth-patch remainders | Explicit binary64 sphere API; represented inputs are exact. Volume certificates exclude chemical radius and coordinate uncertainty. |
| Regional descriptors | Four quadrants, eight octants, near/far and anisotropy extrema in finite-grid API | Quadrants: (+,+),(-,+),(-,-),(+,-); octants 0–3 have z>=0, 4–7 z<0. Near is negative local Z. Separate grid denominators mean regional sums need not equal total. No certified regional API yet. |
| Radius profiles | `profiles.rs` and `volume_profiles` tests | Fixed center and orientation policy; requested order, duplicates and per-radius failures retained. No certified profile API yet. |
| Pyramidalization | `geometry/pyramidalization.rs` | Scalar triple product and alpha angle; explicit neighbor mapping and degeneracy checks. |
| Ensembles | `ensemble/{input,statistics,sensitivity}.rs` | Per-conformer identity/failure, min/max including ties, arithmetic mean, median, population SD, range/sensitivity and weighted summaries. |
| Populations | `ensemble/population.rs` | Explicit weights or complete energy/unit/model metadata; g exp(-ΔE/RT), stable normalization, temperature and missing-mass policy. Approximate electronic/MM energies do not establish observed populations. |
| Molecular features | `descriptors.rs`, `model/features.rs`, input/topology/radius modules | Steric/electronic terms and interactions; published and generated geometries must remain distinct. |
| Selectivity models | `model/{fit,regress,training,evaluation,domain,portable}.rs`, `reaction.rs`, `kinetics` | OLS/BIC, ridge/LASSO, grouped/fixed-feature validation, uncertainty and domain diagnostics. A legacy field named nested_loo does not imply nested feature selection. Eyring conversion requires justified T/sign/kinetic assumptions. |

Gas constant is 0.0019872042586408316 kcal mol⁻¹ K⁻¹; kB is
1.380649e-23 J K⁻¹; h is 6.62607015e-34 J s. The existing thermochemical
kcal conversion is 4184 J. No constant, radius or scale is to be fitted to
reference agreement. π enclosure and IEEE round-to-nearest/gradual-underflow
assumptions are specified in the preceding continuous mathematics reports.

## Historical evidence and limits

Read: `ultimate_engine/{FINAL_REPORT,SCIENTIFIC_AUDIT,NUMERICAL_ACCURACY}.md`,
`occupancy_precision/{FINAL_REPORT,CONTINUOUS_MATHEMATICS}.md`,
`validation/dbstep/FINAL_RESOLUTION.md`, and relevant study/audit sources.
The preceding campaign preserved 128,021 historical responses across 13 lanes;
372 native Rust tests and 71 Python tests passed there. Those are historical
counts, not new campaign results. Morfeus/Kraken validation includes 31,721
conformers, 95,163 frames and 1,141,956 bins; historical topology, geometry,
energy-provenance and DBSTEP default-grid disagreements remain visible.

The exact eighteen unmet requests are locked by
`occupancy_precision/c1_validation.json:analytical.unmet`, with original inputs
in `.stericx/occupancy_precision/analytical_chord_requests.json` and outputs in
`c1_public_analytical.jsonl`. Ten are repeated budgets/rigid transformations of
two three-sphere geometries at 1e-5 Å³; eight use 1e-200 or 1e-150 Å radii and
minimum-positive-subnormal absolute tolerance. They are eighteen requests, not
eighteen chemically distinct molecules. New evidence will retain this denominator.

Ni-hDA has only eleven experimentally labeled rows among 1,566 source ligands.
Earlier study outcomes are already known, so no new split is prospectively blind.
Historical 298.15 K metadata was corrected to 353.15 K in the audit, and ligand
2064 has inconsistent source ee/ΔΔG metadata. Published conformer energies have
unresolved availability/identity. New analysis must use direct primary evidence,
lock its response/partitions/model comparisons before computing new predictions,
and report the retrospective nature and source sensitivity analyses.
