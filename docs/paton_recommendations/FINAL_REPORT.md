# Paton recommendations: implementation outcome

StericX now has a public conformer-aware Rust API, a tested Python distribution,
independent DBSTEP comparisons, and a separately validated density-volume mode.
Existing default scientific outputs remain unchanged. The final confirmed
benchmark supports performance preservation for the frozen workloads. A release
candidate is **withheld** because the original strict DBSTEP hard-sphere volume
gate remains false; this report does not claim every requested gate passed.

## What already existed

The baseline contained native Sterimol, hard-sphere total/regional buried volume,
pyramidalization, geometry parsers, arithmetic conformer aggregation, ordered
parallel file processing, preparation/thermodynamic helpers, and extensive
Morfeus/Kraken audits and performance machinery. Optional denser angular scans
also existed. These kernels and tools were reused. There were no PyO3 bindings,
installable `stericx` Python package, density-grid reader or public density-volume
analysis. [The initial inventory](IMPLEMENTATION_PLAN.md) records the distinction
for every requested step before production changes.

## What was implemented

* Ordered per-conformer records preserve source geometry, identity, metadata and
  descriptor-specific failures. Sterimol, total/percent volume, quadrants,
  octants and pyramidalization are available through the shared engine.
* Applicable fields have arithmetic mean, extrema, population standard deviation,
  median and explicitly justified weighted means. Sensitivity reports identify
  all tied extreme conformers. Missing results and regional-frame incompatibility
  remain visible.
* Stable population calculations support kcal/mol, kJ/mol, Hartree and eV,
  temperature, degeneracy and explicit energy provenance. Electronic, force-field
  and Gibbs-energy interpretations remain separate. Missing energy information
  does not silently produce uniform populations in the new API.
* A fixed 11-ligand/56-conformer study compares representative, mean, energy-weighted,
  minimum and maximum reductions. Arithmetic averaging improves one historical
  RMSE but worsens MAE; energy weighting worsens both. No universal ensemble
  superiority or held-out predictive success is claimed.
* Official, unmodified DBSTEP 2.2.0 is pinned and independently exercised for
  Sterimol, hard-sphere volume, density volume and radius profiles. Criteria,
  geometries, settings, residuals and discrepancies are retained.
* A strict cube reader preserves grid metadata. Separate density occupation uses
  an explicit threshold and documented node-volume integration. Validated analysis
  supports covered spheres on positive, axis-aligned isotropic scalar grids.
* Profiles vary the **sampling-sphere radius at a fixed center**. They preserve
  input order, duplicate radii and errors. They do not silently move the center.
* Thin PyO3 bindings expose typed molecule/ensemble/density results, structured
  errors, NumPy/pandas/RDKit adapters and ordered batches. Maturin builds wheels
  and a source archive. The tutorial includes an executed 11-cell notebook.

See [the API contract](CONFORMER_API.md), [Python examples](../../python/examples/README.md),
[density methodology](../methods/electron_density_buried_volume.md), and
[the ensemble study](conformers/ENSEMBLE_STUDY.md).

## Scientific results and independent agreement

**No legacy scientific result changed:** all 128,021 request/response/error
records in 13 lanes match the freshly frozen baseline byte for byte. Fresh full
reference replay found zero changed legacy scientific values. The audit includes
31,721 conformers, 1,141,956 regional bins, invariance/symmetry/precision studies,
7,624 model-equation checks and the original thermodynamic coverage. This proves
preservation of both successful and unsuccessful results, not universal validity.

Independent population comparisons cover 42 valid and 11 invalid cases; maximum
absolute population error is 6.05e-17. Fresh Morfeus comparisons cover 784 direct
conformer fields and 770 ensemble reductions, with no reference failures or
predeclared-bound exceedances. DBSTEP Sterimol completes 61 cases with separately
identified radii, signed-axis and finite-angle conventions. All 12 density cases
pass their independent criteria, including exact node counts and analytic
resolution convergence.

DBSTEP hard-sphere matched-grid agreement passes 68 of 69 cases. The R=2.5 Å
ligand-401 witness fails in all three frames because f32 and f64 sphere-boundary
inclusion differ. The strict gate remains false. Common-grid occupied-point
decisions agree, and the earliest divergence is frozen and classified. Numerical
constants were not adjusted to improve agreement. Historical Kraken reproduction
is kept separate from direct same-geometry descriptor agreement.

The [complete scientific report](SCIENTIFIC_VALIDATION.md) and
[final machine-readable review](validation/final_scientific_review.json) retain
coverage, errors, reference definitions and exact evidence bindings.

## Installation and performance

Clean wheel installations pass 71 package tests per environment. Actual GitHub
Actions runs pass all 12 Linux/Windows/macOS × Python 3.11–3.14 wheel jobs and the
source-install job. Tested architectures are Linux/Windows x86_64 and macOS
arm64; untested platforms are not implied. Rust CI also passes on all three OSes.
The initial Linux interpreter-selection failure is retained and fixed without
chemistry changes. The development version is **0.4.0.dev0**. No PyPI upload has
occurred, so registry `pip install stericx` is not claimed to work.

The confirmed quiet-window campaign completes 504 samples, preserves every
scientific fingerprint, and finds no >10% median regression across 28 native
configurations. The largest side-median slowdown is 2.25%. The 10,000-file batch
processes about 10,825 files/s at six threads; the public six-conformer Python
call takes 4.096 ms. Existing allocation/reallocation counts are unchanged.
All four earlier, noisier campaigns remain separate and auditable; no outliers
were discarded. [The benchmark report](performance/PERFORMANCE.md) provides
complete timing, memory, allocation, overhead, host and source identities.

## Unresolved scope

The strict hard-sphere DBSTEP witness blocks an all-gates-passed release claim.
Native f32 arithmetic and the default finite-angle B1 scan retain their precision
and rotational-invariance limits. Different reference radii, axial definitions,
coordination-center defaults and integration conventions remain explicit.

Historical Kraken geometries/energies are incomplete, and the ligand-2064 source
conflict remains. Current MMFF weighting does not reconstruct historical or true
free-energy populations. Ensemble superiority and catalyst prediction have not
been established. Density threshold choice is not a universal physical constant;
orbital fields, oblique/anisotropic analysis and uncovered spheres are unsupported.
The reader does not turn these into silent approximate calculations.

[Release readiness](RELEASE_READINESS.md) records each gate and explains why no
release candidate was created. Large raw evidence remains in the retained local
`.stericx/` trees; committed receipts provide hashes and reproduction commands.
They do not imply every raw reference artifact ships in a fresh Git clone.

## Git history

The work starts after `c9380d2c24ba5d48cd52ec13c139ee6e57f04991` and uses separate
implementation, validation, study, packaging, performance and documentation
commits. [The commit ledger](COMMITS.txt) records full identities through the
performance commit. This report and ledger are added by the following documentation
commit, identifiable without a self-referential hash using:

```bash
git log --reverse --format='%H %s' c9380d2c24ba5d48cd52ec13c139ee6e57f04991..HEAD
```

The validation branch carries the implementation used by remote CI. No remote
main merge, release tag or package publication is part of this work.
