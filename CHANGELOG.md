# Changelog

All notable changes to StericX are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- Internal refactor of the `stericx` CLI binary: the single 2,343-line
  `src/main.rs` was split into focused modules — `cli` (argument types),
  `reaction` (reaction-CSV ingestion and Sterimol ensemble aggregation),
  `descriptors` (geometry-only featurization and the `descriptors` command),
  `output` (atomic writes and run metrics), and a `commands/` module with one
  file per subcommand — leaving `main.rs` a ~110-line dispatch shell. Purely
  structural: no logic changed, each test moved alongside the code it exercises,
  and the full 57-test suite plus `fmt` / `clippy` / `ruff` gates stay green.
  User-facing behavior and command output are identical.
- README restructured as a scannable landing page (~1,030 → ~305 lines): a
  one-minute overview, a curated six-number **Key Results** block, and a quick
  start stay inline; the nine-study detail, CLI reference, internal architecture,
  benchmarks, and roadmap move into collapsed sections and links to the existing
  per-study docs. Statistical density was deliberately cut — the reading path now
  carries only the handful of numbers that matter, with every RMSE / Pearson /
  slope / intercept / median-AE table left in `docs/`. The opening now leads with a
  plain-language description before any chemistry, and the ordering puts validation
  ahead of speed (benchmarks collapsed low; speed is the last Key Results row). No
  content was removed. Added a **Technical challenges** section documenting the
  engineering decisions behind the reproduction — donor detection from raw geometry,
  isolating the kernel from the input geometry, matching Kraken's coordinate
  conventions, the phosphine frame bug found at scale, the unsafe-SIMD/scalar-fallback
  discipline, and the honestly-scoped reproducibility story (machine-native format
  guarded by an endian marker; not claimed as cross-platform).

### Added

- `stericx search` — a ligand similarity and constraint search built on the
  validated descriptors, turning the featurizer into a discovery tool. Given a
  query geometry and a library (a directory of structures, or a CSV previously
  written by `stericx descriptors --format csv`), it ranks candidates by
  Euclidean distance in **standardized** descriptor space: each descriptor is
  z-scored against the library first, so an Å-scale Sterimol `L` and a
  percent-scale `%Vbur` contribute comparably rather than the largest-unit
  feature dominating the ranking. The default similarity space is the shape
  envelope (`L`/`B1`/`B5`), `%Vbur`, the quadrant asymmetry
  `max_delta_qvbur`, and the donor pyramidalization `pyr_P`; `buried_volume` is
  deliberately excluded from it because it is `percent_buried_volume` rescaled
  by a constant sphere volume and would otherwise double-weight the same
  quantity. Repeatable `--filter` constraints (`vbur=30..35`, `b5<7`, `l>=8`,
  with short aliases) narrow the candidate set before ranking, and
  `--less-bulky` / `--more-bulky` express a `%Vbur` bound relative to the query.
  Descriptors that are constant across a library are reported and dropped from
  the distance so a ranking is never silently computed on fewer axes than
  requested. Reported honestly as a **steric-similarity ranking, not a
  prediction of reactivity**. Library featurization is parallelized; the
  `descriptors` command it shares kernels with is untouched, so Study 008's
  single-core benchmark still measures exactly what it claims.
- Cross-platform CI: a `cargo test` matrix on `macos-latest` (Apple silicon,
  aarch64) and `windows-latest` (x86_64) runs alongside the existing Linux
  quality job, which keeps the `fmt` / `clippy` / Python / `ruff` gates. This
  exercises the endian-guarded binary format on a second OS and both sides of
  the SIMD split — the AVX2 inference kernel is `x86_64`-gated with a scalar
  fallback, so the Apple-silicon runner builds and tests the scalar path while
  the Windows x86_64 runner exercises the AVX2 path. All three jobs are green,
  so the build and test suite are now demonstrated on Linux, macOS, and Windows.
- Study 010 — grid convergence of the buried-volume integrator
  (`studies/study_010_grid_convergence.py`, `docs/study_010/`). Studies 002–004
  validated the descriptor at a single grid resolution (Kraken's 0.01 Å³/point);
  this sweeps the integration grid from 0.5 down to 0.001 Å³/point on 60 diverse
  ligands and shows %Vbur **converges** — coarse grids (≥0.05) deviate by ~0.4
  %Vbur, the default 0.01 sits within ~2× the ~0.03 %Vbur discretization floor, and
  refining further costs ~10× the wall time for no meaningful gain. Reported
  honestly: the default's grid uncertainty is comparable in size to the Study 004
  reproduction residual, but that residual is not grid-limited (Kraken uses the same
  grid, so discretization is largely common-mode, and Study 006 localizes the
  residual to structure a finer grid would not touch).
- Non-phosphorus donor validation (`scripts/validate_nonp_donor.py`,
  `docs/validation/NONP_DONOR.md`). The `descriptors` command advertises
  `--donor-element` for any trivalent donor, but every study validated only
  phosphorus. This checks the nitrogen path the same way Studies 002/005 check
  phosphorus: on ten tertiary amines, StericX's `--donor-element N` reproduces
  morfeus's pyramidalization and buried volume to machine precision and Sterimol
  L/B5 essentially exactly (B1 carries the same ~0.01 Å 1°-scan difference seen for
  phosphorus), with donor detection finding the nitrogen and its three substituents
  from geometry alone in every case. Scoped honestly to trivalent nitrogen — oxygen
  and other donors remain unexercised.
- Test asserting the unsafe AVX2 inference kernel agrees with the portable scalar
  fallback across 5,000 random inputs (to f32 tolerance, since the two sum in
  different orders), so the SIMD path cannot silently diverge from the reference
  arithmetic. Runs only where AVX2 is present.
- Study 009 — the opposite direction of the buried-volume ligation cliff
  (`studies/study_009_pd_crosscoupling.py`, `docs/study_009/`). Study 007
  reproduced the Newman-Stonebraker *Science* 2021 classifier for the **nickel**
  reactions, where the *small* ligands are active (`Left` of the ~32 % %Vbur(min)
  cliff). Study 009 reproduces the **palladium** reactions VII–XII, where the
  relationship inverts and the *bulky* ligands are active (`Right` of the cliff),
  using the same descriptor and the paper's mechanistically-preferred *balanced*
  class weight (Tables S12/S14). StericX's %Vbur(min) matches the published value
  at R² = 0.9994 (MAE 0.116 %, n = 267) and its single-node classifier
  independently recovers the `Right` direction for all six reactions (mean MCC
  0.64 vs the paper's 0.67), **including two reactions taken from other groups'
  published data** (Zhao *Science* 2018; Stambuli *et al.*). Honesty is preserved:
  Reaction VII — which the paper's own SI flags as pathological — is where the
  independent fit also disagrees, and the weak Heck dataset (XII) is weak for
  everyone. The copyrighted SI (AAAS) is read locally and never redistributed;
  only StericX's own values and the aggregate comparison are committed.

## [0.1.0] - 2026-07-31

First public release: a complete, independent, honest reproduction of the Kraken
physical-organic ligand descriptors and two published reaction models, as a
single dependency-free Rust binary. Highlights: descriptor kernels validated
against `morfeus` to numerical precision and against Kraken's published values
across 1,541 ligands (buried-volume family, Sterimol, and pyramidalization); a
second reaction model (the Newman-Stonebraker cross-coupling classifier); a
head-to-head speed benchmark (~14x faster than `morfeus` at library scale); and a
pre-registered, falsifiable prospective prediction. Negative results are led with,
not hidden — most notably, StericX's compact native descriptors do **not** model
Ni-hDA enantioselectivity (leave-one-out Q² ≈ 0.002).

### Added

- Pre-registered prospective prediction (`scripts/preregister_prediction.py`,
  `docs/study_003/PREREGISTRATION.md`): the frozen Study 003 Ni-hDA ligand deck is
  elevated to a falsifiable pre-registration. Without touching the frozen deck (it
  is committed to by SHA-256 and verified byte-for-byte), each point prediction of
  ΔΔG‡ gains a 95% OLS prediction interval and an applicability-domain judgement
  (leverage vs h* = 3p/n = 0.60; 9 of 10 candidates inside). The document fixes,
  before any measurement, an exact experimental protocol (asymmetric Ni-hDA,
  cited to the Sigman public repo; ee → ΔΔG‡ via 2RT·atanh(ee)) and a
  pre-registered falsification rule (measured ΔΔG‡ within the 95% interval for ≥6
  of 8 primary candidates AND positive predicted-vs-measured rank correlation,
  else falsified; two boundary ligands scored separately as an extrapolation
  stress-test). It keeps the honest limit in view: for 7 of 10 candidates the 95%
  ee interval spans both enantiomers, so those near-racemate predictions are not
  meaningfully falsifiable on ee, and the test is run on ΔΔG‡. Candidates are
  DFT-characterized Kraken ligands, not synthesis or safety instructions.
- Study 008 (`studies/study_008_speed_benchmark.py`, `docs/study_008/STUDY_008.md`): a
  head-to-head throughput benchmark against morfeus, the reference Python
  implementation. Computing the flagship buried-volume descriptor on the same
  1,546-ligand library, on the same single CPU core, StericX is ~14x faster
  (~1,110 vs ~81 structures/s) — a conservative figure, since StericX also
  computes Sterimol and pyramidalization in the same timed pass. The speedup is
  not from a cheaper computation: of the 1,534 phosphines morfeus can frame,
  1,518 agree to R2 = 0.999999 (max |diff| 0.42 %Vbur); the 16 that differ are
  frame-topology cases (morfeus's nearest-3-heavy vs StericX's covalent bonding,
  Study 006), reported separately rather than averaged away. Only StericX's own
  timings and the aggregate agreement are committed; the Kraken SDF cache stays
  local and gitignored.
- Study 007 (`studies/study_007_crosscoupling.py`, `docs/study_007/STUDY_007.md`): an
  independent second reaction model. StericX reproduces the %Vbur(min)
  cross-coupling reactivity classifier of Newman-Stonebraker et al. (Science
  2021, 374, 301) on the authors' Reactions I-V and RS1 — matching the published
  %Vbur(min) across 479 ligand-reaction points at R2 = 0.9992 (MAE 0.144%), and
  recovering the paper's single-node threshold classifier (Table S11): same
  thresholds (~32%), same direction, matching mean accuracy (0.69 vs 0.69). The
  study now also (A) tests out-of-sample transferability of the ~32% ligation
  cliff — a pooled universal threshold (32.8%, MCC 0.49) and leave-one-reaction-
  out cross-validation (OOS MCC 0.41-0.54), with Reaction V documented as a
  threshold-space outlier that nonetheless still transfers; (B) removes the
  shared-geometry circularity by running StericX on the authors' own DFT
  free-ligand structures (18 ligands, matched to Kraken IDs by molecular formula),
  reproducing published %Vbur(boltz) at R2 = 0.9735 with 16/18 inside StericX's
  own Kraken-conformer range; and (C) reports bootstrap 95% CIs on the per-
  reaction accuracy/MCC. The paper's copyrighted supplementary data is read
  locally (data/external/, gitignored) and never redistributed; only StericX's
  computed values and the comparison are committed.
- Study 006 (`studies/study_006_residual_localization.py`, `docs/study_006/STUDY_006.md`): a
  controlled test that localizes the buried-volume residual to the coordination
  centre. On the same 1,541 ligands, the four descriptors anchored on the centre
  (buried volume, Sterimol L/B1/B5) shift by a mean 1.54 residual-sigma per P-H
  bond, while the two centre-free pyramidalization descriptors stay flat at 0.04
  sigma. Since pyramidalization shares the donor, geometries, frame, and kernel
  and differs only in never placing the centre, this rules out the kernel, the
  frame, and the geometries — confirming the residual is the geometric lone-pair
  centre diverging from Kraken's xTB centre, not a bug.
- Pyramidalization descriptors (`src/geometry/pyramidalization.rs`): a native
  Rust kernel for Kraken's `pyr_P` (Radhakrishnan) and `pyr_alpha` (mean
  out-of-plane angle), reduced to order-invariant closed forms on the donor's
  three unit bond vectors — `pyr_P = |det[â, b̂, ĉ]|` with morfeus's `2 − P`
  acute correction, and `pyr_alpha` as the mean signed out-of-plane angle. Both
  were verified against morfeus to machine precision (4.4×10⁻¹⁶ and 2.8×10⁻¹⁴)
  before implementation. Surfaced through `stericx descriptors` (text/json/csv).
- Study 005 (`studies/study_005_pyramidalization.py`, `docs/study_005/STUDY_005.md`):
  the native kernel reproduces Kraken's published `pyr_P` and `pyr_alpha` (min
  and max conformer reductions) across the same 1,541 ligands at mean
  R² = 0.99998. The small residual (RMSE ~2×10⁻⁴ for `pyr_P`, ~0.03° for
  `pyr_alpha`) reflects the cached DFT SDFs' 4-decimal coordinate precision, not
  a method difference.
- Sterimol reproduced against Kraken at scale (`studies/study_004_sterimol.py`,
  `STUDY_004_STERIMOL.md`). Kraken measures Sterimol along the coordination axis
  (a virtual metal 2.28 Å from the donor on the lone pair, +0.40 Å Verloop `L`
  correction) — a convention recovered by a distance sweep, not assumed. With it
  matched, StericX reproduces Kraken's published `sterimol_L/B1/B5` across the
  1,541 ligands at mean R² = 0.9887. Exposed as a new kernel method
  (`SterimolCalculator::compute_with_dummy`), the public `coordination_center`
  helper, and `descriptors --sterimol-axis coordination`.
- Full buried-volume descriptor-family validation (`studies/study_004_vbur_family.py`,
  `STUDY_004_FAMILY.md`): StericX reproduces Kraken's *entire* published `vbur`
  family — buried volume, quadrant and octant extrema, and near/far hemispheres,
  eight descriptors — across all 1,541 ligands at mean R² = 0.9925. The headline
  `max_delta_qvbur` lands at 0.9852 by an independent path, matching the scaled
  Study 004 as an internal consistency check.
- `stericx descriptors <file>...` — compute Sterimol (L, B₁, B₅) and
  buried-volume descriptors for any ligand geometry directly. The phosphorus
  donor and its substituents are detected from the geometry (covalent-radius
  bonding), so no reaction CSV or manual atom indices are needed. Accepts
  `.xyz` and `.sdf`/`.mol`, treats a multi-model file as a conformer ensemble,
  and reports Kraken's `max_delta_qvbur_min` as the headline descriptor.
  Supports `--format text|json|csv`, batch runs over many files (unparseable or
  non-trivalent inputs are skipped on stderr), and `--donor-element` /
  `--donor-index` for non-phosphine or ambiguous donors.
- Residual-anatomy study (`studies/study_004_frame_residual.py`, `STUDY_004_RESIDUAL.md`,
  `residual_by_phosphine_class.png`): the full-set residual is dissected by
  donor class, showing tertiary phosphines are unbiased and the remaining bias
  is confined to primary/secondary phosphines and grows ~0.7 Å³ per P–H bond.
- Robust statistics for the full-set validation (mean/median absolute error,
  abs-residual percentiles, and trimmed R² of 0.9897 at 1 % / 0.9936 at 5 %).
- Rust unit tests for donor detection, CSV quoting, and the descriptors path;
  SDF tests for blank-title and multi-record files. Continuous integration now
  also syntax-checks the Kraken-reproduction study drivers.

### Fixed

- **Buried-volume frame for primary and secondary phosphines.** The quadrant
  frame identified a donor's substituents as its three nearest *heavy* atoms,
  which mis-framed R–PH₂ and R₂P–H donors by discarding their bonded hydrogens
  (producing a spurious `max_delta_qvbur = 0` or a gross overestimate). The
  frame now uses the donor's covalently bonded atoms, hydrogens included. This
  is identical to the old behaviour for trisubstituted donors, and raised the
  full-set R² from 0.9649 to 0.9852 without discarding any ligand (the validated
  count rose from 1,535 to 1,541).
- **SDF parsing of blank title lines.** A blank first (title) line — which
  OpenBabel and RDKit both emit — was stripped together with the `$$$$` record
  separator, misaligning the fixed four-line header so the counts line was
  misread. Only the inter-record separator newline is now removed.

### Foundation (initial reimplementation)

- Foundational reimplementation captured in this release: the Sterimol and
  buried-volume kernels (validated against `morfeus`, Sterimol R² ≥ 0.9999 and
  buried-volume geometry to numerical precision), the Ni-catalyzed
  homo-Diels–Alder reproduction, the SIMD `.sigpack` storage format, the Eyring
  kinetic link, quantum-geometry tooling (CREST/xTB), and the reproduction report.
