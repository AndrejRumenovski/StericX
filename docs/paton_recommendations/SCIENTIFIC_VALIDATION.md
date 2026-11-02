# Scientific validation after the Paton extensions

The extensions preserve the frozen engine's default numerical outputs. Fresh
independent comparisons reproduce the admitted legacy reference results, including
their negative findings. New conformer, population, density and profile methods
have separately stated validation scopes. This report does not assert a global
scientific or predictive pass.

## Source and numerical preservation

The baseline is `c9380d2c24ba5d48cd52ec13c139ee6e57f04991`. The source-bound
baseline and final native/observer builds, toolchains and executable hashes are
recorded in `baseline/manifest.json` and `validation/final_build_manifest.json`.
The final build contains the implemented scientific kernels; subsequent
documentation, orchestration, benchmark adapters and CI corrections do not
change those kernels.

All **128,021 records in 13 lanes** have identical request, response and error
streams against the freshly frozen baseline. Both sets of raw files were
rehashed. The separate four-thread-count CLI campaign retains all **93 cases**
at 1, 2, 4 and 6 threads, with zero differences under the existing explicit
timing/path/timestamp exclusions. Packed artifacts, scientific values, ordered
results, failed commands and error text remain covered.

No production constant, radius table, default integration density, coordinate
precision or angular sampling default was changed. The new per-conformer API
reuses those kernels. The optional 3,600-direction angular API already existed
and was reused.

## Fresh independent legacy replay

The unmodified admitted reference equations and Morfeus 0.8.0 were rerun on both
inferred-connectivity and explicit-topology input lanes. Exact reference
comparisons found **zero changed scientific values or structures**. The review
includes all original errors and retained residual classifications, rather than
testing only correlations or favorable cases.

| Scope | Coverage | Finding |
|---|---:|---|
| Full topology volume reference | 31,721 conformers | No reference execution errors; complete source stream validated |
| Regional integration | 95,163 frames; 1,141,956 quadrant/octant bins | All public and private results compared; original radius/grid strata retained |
| Fixed variant union | 2,517 cases in each of two lanes | Membership, status, raw reference and recombined means preserved |
| Model equations | 7,175 + 449 comparisons | Zero failed numerical/domain-contract checks |
| Thermodynamics | 744 kinetic scalars; 51 aggregation scalars/rejections; 14 CLI cases | Scoped numerical/domain pass |
| Preparation and selectivity | 6 CREST parser mocks, 11 MMFF cases, 24 ee cases | Existing provenance and invalid-input contracts preserved |

The original geometry campaign contains atom permutations, rigid rotations,
symmetry, planar/axis ambiguities, precision and integration-density variants.
Its unchanged results remain a finite validation campaign. Default finite-angle
B1 is not exactly rotationally invariant. Larger full-corpus residuals associated
with B/Br and other reference-radius conventions remain in the complete tables;
they are not recast as numerical agreement.

The model checks establish implementation fidelity, not successful chemical
prediction. Failed outer folds, poor interval coverage, unfavorable ligand
ranking, missing prospective outcomes and the ligand-2064 source conflict remain.
All 1,566 Ni-hDA records retain their published targets and 353.15 K temperature.

## New conformer and population methods

The native API retains every source record, identity and failure. Per-descriptor
outcomes permit partial descriptor availability without silently dropping a
conformer. Means, extrema, population standard deviation, median and explicit
weighted means follow the definitions in [CONFORMER_API.md](CONFORMER_API.md).
Incomplete coverage and inconsistent regional correspondence are visible.

The independent 80-digit Decimal population reference covers **42 valid and 11
invalid cases**, including all four energy units, two-state and multistate
ensembles, temperature, degeneracy, extreme separations and 4,096 states. Maximum
absolute population error is **6.05e-17** and maximum relative error is
**1.59e-14**, within the prospectively fixed limits. Electronic, force-field and
Gibbs free-energy interpretations remain separate; missing or incomparable
energies cannot silently become equilibrium populations.

The fixed **11-ligand / 56-conformer** study compares 784 direct conformer fields
and 770 ensemble reductions with fresh Morfeus calculations. No reference call
failed. Representative, arithmetic, recorded-MMFF-weighted, minimum and maximum
reductions are all retained. Historical-target MAE/RMSE are not uniformly improved
by ensemble averaging. The study fits no predictive model and invents no missing
historical energies. See [the study](conformers/ENSEMBLE_STUDY.md) for the full
design, residuals and unsuccessful comparisons.

## Independent DBSTEP comparison

DBSTEP 2.2.0 is pinned to official commit
`4365353e0a2f24d05336e619759e350bc630a3cb`. Its source was not modified. The
environment, dependencies, geometries, commands, settings and prospective
acceptance rules are in [the DBSTEP evidence](../validation/dbstep/).

* **Sterimol:** all 61 cases complete. Common-radii L/B5 agree at the f32
  arithmetic scale; B1's maximum difference is 0.01891206 Å under the fixed
  one-degree angular bound. Default H radii and signed versus absolute axial L
  definitions produce separately labelled differences. High R² does not remove
  these offsets.
* **Hard-sphere volume:** 68 of 69 matched-lattice cases pass the original
  arithmetic bounds. The R=2.5 Å ligand-401 witness fails in all three frames:
  native f32 inclusion retains 5,565 sphere nodes, whereas DBSTEP's f64 reclip
  retains 5,533. All occupied-point decisions on the common supplied grid agree;
  six occupied boundary nodes per frame are removed by the reclip. The strict
  volume gate remains **false**. The complete final summary is exactly equal to
  the original `volume_v2` summary, including its three denominator errors and
  failed residuals. This is a numerical integration difference, not evidence for
  changing production constants.
* **Independent lattices:** DBSTEP spacings 0.20, 0.10 and 0.05 Å are retained.
  Grid-count volume, analytic-sphere-normalized volume and whole-grid `occ_vol`
  are distinguished. The finest independent-grid maximum volume difference is
  1.19544539 Å³ on the common-surface subset; no fitted tolerance hides it.

The [discrepancy study](../validation/dbstep/DISCREPANCIES.md) freezes the failing
geometry and verifies mapping, axes, radii, boundary arithmetic and Morfeus
triangulation at the earliest mathematical divergence.

## Density and profiles

All **12 density comparisons** pass against unmodified DBSTEP, including Ne,
benzene, a larger organic cube and analytical Gaussian grids. Occupied and sphere
node counts agree exactly; unit-converted grid volumes pass the fixed f64
arithmetic bounds. The Gaussian resolution sequence reduces absolute volume
error from 0.03955 to 0.01206 to 0.00312 Å³. Analytic boundary-layer bounds are
independent of the observed errors.

This validates scalar, positive, axis-aligned isotropic grids, covered evaluation
spheres and an explicitly declared electron-density threshold. It does not
validate a universal physical threshold, orbital fields, oblique-grid analysis,
interpolation, padding, quantum-chemical density quality or interchangeability
with hard-sphere volume. The reader preserves full grid metadata and rejects
unsupported field/unit conventions; analysis rejects unsupported grids.

The radius-profile study contains **20 ordered positions**, **16 independent
valid-position comparisons**, and **four retained expected errors**. Each valid
position reproduces its individual call. Radius changes the sampling sphere at
a fixed center. The known hard-sphere R=2.5 Å strict DBSTEP failure remains;
density profile comparisons pass their separate criteria.

## Engineering, portability and retained attempts

Local checks pass **348 Rust tests**, **128 research-tool Python tests**, Clippy,
Rust formatting and rustdoc. The installed Python distribution passes **71 tests
per environment**, including all 56 frozen conformers, table adapters, RDKit,
errors, populations, density and profiles. The notebook executes all 11 code
cells and its plotted conformer dependence was visually inspected.

The actual GitHub wheel matrix now passes all 12 Linux/Windows/macOS × Python
3.11–3.14 jobs, plus source-distribution installation. The separate Rust CI passes
on all three operating systems. See [distribution validation](python/DISTRIBUTION.md)
for exact runner/interpreter/artifact identities and retained initial failures.

Interrupted reference jobs retain their partial outputs; fresh complete runs use
`_v2` directories. The historical supplement controller's directory-name
assumption failed before calculation; its adapter binds explicit input receipts
and reuses unchanged equations. An initial review incorrectly required no DBSTEP
denominator errors despite exact equality to the known failed witness; the
corrected review requires preservation of the entire original summary and keeps
the strict gate false. Neither correction changes a numerical acceptance limit.

Raw local campaigns live under `.stericx/paton/`. Committed manifests and summaries
bind them by SHA-256; replay commands in the logs/controllers create new output
directories. `scripts/review_paton_science.py` compares all legacy reference
artifacts against the frozen admitted evidence and admits only descriptor
performance measurements. It explicitly does **not** admit a release or set a
global scientific pass. Performance and release decisions have separate reports.
