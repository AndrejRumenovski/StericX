# Changelog

All notable changes to StericX are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

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

## [0.1.0]

Initial reproduction study: from-scratch Rust reimplementation of the Kraken
Sterimol and buried-volume ligand descriptors, validated against `morfeus`
(Sterimol R² ≥ 0.9999; buried-volume geometry to numerical precision) and
against Kraken's published `vbur_max_delta_qvbur_min` on the authors' DFT
geometries (11 Ni-hDA ligands, R² = 0.9986). Includes the Ni-catalyzed
homo-Diels–Alder reproduction, the SIMD `.sigpack` storage format, the Eyring
kinetic link, quantum-geometry tooling (CREST/xTB), and the reproduction report.
