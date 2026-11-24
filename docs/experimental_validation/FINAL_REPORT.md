# StericX: experimental chemistry and certified geometry campaign

The campaign establishes a reproducible association with a published catalytic
selectivity dataset and implements independently validated continuous regional and
radius-profile APIs. The ensemble-minimum descriptor improves the primary
cross-validation point estimate, but the evidence **does not confirm the broader
conformer-aware predictive hypothesis**. All **eighteen original numerical target
failures** are resolved by a separately selectable method without changing their
tolerances or default behavior. Four different, deliberately stringent regional
budget tests remain unresolved and explicitly reported.

## 1. Experimental chemistry

### Dataset and provenance

Seven candidate literature/library resources were reviewed, spanning Sigman/Reid,
Paton, Zahrt, phosphine cross-coupling and bisphosphine studies. Counts, structures,
conditions, missing records, energy availability and usage limitations are recorded
in [DATASET_REVIEW.md](DATASET_REVIEW.md). A descriptor library alone was not counted
as experimental validation.

The selected study is Cadge et al., **A Data Science-Guided Approach for the Development
of Nickel-Catalyzed Homo-Diels–Alder Reactions**, [JACS 2025, DOI
10.1021/jacs.5c09948](https://doi.org/10.1021/jacs.5c09948). SI Table S3 reports the
norbornadiene/methyl-vinyl-ketone ligand screen: Ni(cod)2 6.5 mol%, ligand 13 mol%,
1,2-DCE 0.1 M, 80 °C, 16 h. There are twelve quantified ee observations and one
unquantified result. Eleven quantified ligands have verified-source structure exports
and author descriptor rows; 906 is missing and its structure API fails. Yield and
unquantified ee are preserved, not imputed or substituted as endpoints.

The later [published correction](https://doi.org/10.1021/jacs.5c18763) and
[corrected SI](https://doi.org/10.1021/jacs.5c18763.s001) were retrieved. Table S3 is
identical after whitespace normalization and was visually checked. Corrections to
other conditions, products and characterization were not ignored or merged into the
selected endpoint. The 2064 3% ee/printed ΔΔG discrepancy remains in both SI versions;
primary response uses ee and a separate sensitivity uses author ΔΔG.

### Structures and calculation

**135 original library conformers** were retrieved for the eleven ligands. Neutrality,
connectivity, explicit hydrogens, original atom order, P coordination and encoded
stereochemistry were audited. Native Rust calculates L/B1/B5, percent buried volume,
all quadrants/octants, pyramidalization and ensemble statistics. The standard
Kraken-matched convention uses a 2.28-Å virtual metal, R=3.5 Å, radius scale 1.17,
H excluded from volume and included in Sterimol. All three donor-plane frames and
all failures are retained. [DESCRIPTOR_REPRODUCTION.md](DESCRIPTOR_REPRODUCTION.md)
records the complete numerical and molecular conventions.

The data-quality audit finds a source/export P-configuration disagreement for all
15 records of 401 and a mixture in 1058 (eight matching, six opposite). A relative
axial audit finds one opposite sector among 2063's twelve records. Absolute
atropisomer identity is not fully established by the source SMILES. These limitations
prevent treating every distributed record as a verified thermally accessible state
of the experimentally used ligand. The original library analysis remains visible,
with exclusion/filtering sensitivities; no geometry is silently replaced.

No usable published per-conformer energies were obtained. Generated ETKDG/MMFF
structures are a distinct preparation experiment. The initial incomplete axial
checks and poor models are preserved. The corrected preparation returns 218
embeddings from 440 requests and retains 69 states across two fixed seeds, with
actual computed energies and explicit rejected candidates. Its Boltzmann-weighted
means at 298.15/353.15/400 K are approximate force-field summaries, not measured
populations or solution free energies. [CONFORMER_STUDY.md](CONFORMER_STUDY.md)
reports energy distributions, search/optimization sensitivity and coverage limitations.

### Held-out models and hypothesis

The protocol fixed the primary dataset, apparent |ΔΔG| response, family splits,
feature and metric before new model evaluation. This remains an outcome-aware
retrospective design: the published feature choice and prior repository studies
already used these outcomes. All current fitting/tuning excludes held-out responses,
and multivariate ridge tuning/scaling uses nested family splits. Eleven ligands,
not 135 conformers or the 1,566-row largely unlabeled library, define sample size.

| Family-held-out model | MAE kcal/mol | RMSE kcal/mol | R² |
|---|---:|---:|---:|
| Training-response mean | 0.956 | 1.055 | -1.025 |
| Published electronic baseline | 0.947 | 1.177 | -1.520 |
| First documented conformer, one steric predictor | 0.685 | 0.850 | -0.313 |
| Ensemble minimum, same predictor/capacity | 0.323 | 0.389 | 0.724 |
| Published steric reference | 0.311 | 0.378 | 0.740 |
| Published reference plus native ensemble | 0.336 | 0.404 | 0.703 |

Primary ΔRMSE(ensemble−single) is **−0.460 kcal/mol**, with descriptive paired
family-bootstrap 95% interval **[−0.728,+0.044]** and exact four-family sign-flip
p=**0.25**. Four clusters, dependent CV fits, historical feature selection and
structural uncertainty prohibit a strong confirmation claim. Conditional prediction
intervals, all predictions/folds, bias, fold variability, correlations and model
complexity are retained in [STATISTICAL_RESULTS.md](STATISTICAL_RESULTS.md).

The point advantage persists in conservative stereo exclusion sensitivities, but
uncertainty remains large. Eleven-feature nested ridge gives RMSE 1.039 for a single
state, 1.053 for ensemble means and 1.093 for 33 ensemble extrema/mean features.
Corrected generated MMFF-weighted descriptors give RMSE 1.040/1.157 across two seeds.
Ensemble averaging, extra features and approximate weighting therefore do not
reproduce a general predictive benefit. The combined model adds no demonstrated
information beyond the published steric reference.

Conformer variation is larger than observed numerical discrepancies for many
ligands, but not for the three nearly duplicate two-record ensembles. Replacing
finite-grid regional features with certified continuous midpoints changes primary
single/ensemble RMSE to 0.84988/0.39548, much less than the statistical uncertainty.
This directly connects the mathematical and experimental parts of the campaign.
Predictive association here is not a mechanistic or causal conclusion.

## 2. Certified regional volume

### Definitions and implementation

For integration ball S and atomic union A, each regional occupied volume is
∫1S1A1region dV. Quadrants use the established XY ordering; octants add positive
then negative Z. Region percentages divide by the corresponding quarter/eighth
sphere volume. Boundary planes have zero continuous volume. Existing finite-grid
normalizations and defaults are unchanged.

The new Rust `sphere_union_volume_certified_regions` shares adaptive XY cells,
spatial candidates and chord geometry between positive/negative Z halves. Plane
clipping uses the exact zero endpoint. Conservative chord-range bounds remain valid
at tangencies, uncertain endpoint order, degeneracies and singular projected
boundaries. Smooth midpoint remainders are admitted only after interval proofs of
endpoint identity, connectivity, separation and positive curvature denominator.

Every region returns estimate, lower/upper bounds, absolute/relative error where
meaningful, percentage enclosure and convergence status. The enclosing result
reports the common budget, evaluations, subdivisions, heap cells and runtime.
Outward accumulation and nested interval intersection preserve global bounds;
quadrants and total are consistent sums of the octant quantities. Independent
estimates need not sum bitwise. The certificate covers represented local balls,
not exact source rotations, radii uncertainty or chemistry.
[REGIONAL_MATHEMATICS.md](REGIONAL_MATHEMATICS.md) gives the derivation and assumptions.

### Radius profiles and validation

`sphere_union_volume_certified_radius_profile` reuses canonical atomic geometry and
relative bounds. Each distinct radius has its own candidate decisions, integration
tree and budget. Exact duplicate radii reuse a certificate; no unproved cross-radius
bound is reused. Input order and per-radius errors remain explicit. Occupied volume
is nondecreasing at fixed center; percent buried volume need not be.
[CERTIFIED_RADIUS_PROFILES.md](CERTIFIED_RADIUS_PROFILES.md) documents these contracts.

Independent rational cross-section/cap/union integration with rational π bounds
contains **1,248/1,248** total/regional quantities in 96 cases. Four 1e-8-Å³ regional
stress requests at small budgets remain unmet with valid intervals. Radius testing
covers 24 profiles, 120 requested entries and 96 distinct integrations; all **936**
analytical quantities are enclosed and all ordinary requested targets are met.
Molecular profile/total agreement is labeled algorithm-path consistency, not an
exact molecular reference. All 405 experimental frames meet 0.05 Å³.

Exploratory median profile time is 65.43 ms, maximum 442.32 ms on the tested cases.
Geometry preparation is a negligible fraction of total cost; only exact duplicate
reuse demonstrably avoids a whole repeated integration. Claims of arbitrary-radius
speedup or universal convergence are not supported. [PERFORMANCE.md](PERFORMANCE.md)
reports CPU, allocation, memory, simulated instructions and parallel diagnostics.

## 3. Numerical stress cases

The eighteen original failed requests comprise ten three-ball union requests and
eight underflow/contact requests at radii 1e-150 or 1e-200. Their original input
geometries, tolerances, budgets, errors and outcomes remain frozen. They are repeated
variants of four mathematical families, not eighteen distinct chemical systems.
[NUMERICAL_FAILURE_ANALYSIS.md](NUMERICAL_FAILURE_ANALYSIS.md) classifies each mechanism.

The separately selectable `sphere_union_volume_certified_geometric` adds conservative
exact-contact separation and a triangle-free overlap-graph proof. When all atomic
balls are contained and no positive-volume triple intersection can exist, exact
pairwise inclusion–exclusion applies. Otherwise it falls back to the existing
certified algorithm. This is not a general pair-overlap approximation. Error-free
TwoSum residuals bound exact endpoint distances without gratuitous subnormal widening.
This isolated change resolves **14/18** requests.

A second isolated proof uses the integration-ball upper bound: for positive
R≤2^-359, V<8R³≤2^-1074. Returning [0,minimum-positive-f64] preserves a possible
positive real volume and avoids falsely identifying underflow with exact zero.
It resolves the remaining **4/18**. Relative accuracy is not claimed when the
certified lower bound is zero.

Combined, the explicit method meets **18/18 original targets within original
budgets**, with **167/167** original analytical references enclosed and 6,005 exact
arithmetic witnesses passing. An extended 180-observation suite has no invalid
comparisons. The largest new bound among the eighteen is 3.20e-14 Å³. Original
baseline success remains 0/18 even at four times its budgets; new success stays
18/18. Analytical proofs eliminate subdivision on these cases. General clipped,
triangle-containing or stringent regional geometries can still exhaust budgets;
this finite result is not universal convergence.

## 4. Scientific integrity and reproducibility

The production baseline is `e0f1dc79d23bf6638614e85cfc0c9f2282b6b6eb`, with source,
constants, toolchains, datasets and historical output hashes frozen in
[BASELINE.md](BASELINE.md). **128,021 historical outputs in thirteen lanes match
exactly**, including regional bins, failures and ordering. All 385 Rust tests pass
on both Rust 1.97 and actual Rust 1.85; 71 fresh Python-binding and 156 research tests
pass. Available macOS/Windows library cross-checks pass. Independent geometry/Kraken
references were freshly replayed, preserving known topology and convention differences.
New experimental Morfeus/DBSTEP matched-input checks pass on all 135 conformers.

A fresh source-bound rebuild reproduces every analytical/regional/profile result,
published/generated descriptor and locked statistical result exactly, excluding only
named timing/allocation fields. Deterministic data bundles preserve original inputs,
actual energies, attempted/rejected structures, source failures and both favorable
and unfavorable models. [REPRODUCIBILITY.md](REPRODUCIBILITY.md) gives commands;
[VALIDATION_REPORT.md](VALIDATION_REPORT.md) states every gate and limitation.

Confirmed issues include source ee/ΔΔG inconsistency, source/export stereochemical
conflicts, incomplete first-generation axial checks and conservative numerical
budget/rounding failures. Failed builds, invalid profiler loading, corrected test
assumptions and unsuccessful models remain archived. No unsupported population,
chemical-causation or blanket data-license claim is made. Computational reference
agreement is not substituted for experimental evidence.

Profiling uses ordinary workstation load. Hardware counters are unavailable;
Callgrind instructions are simulated. Valid Memcheck reruns report zero errors/lost
bytes, with small reachable runtime storage. No controlled wall-time speedup is
asserted, so no quiet-window approval was required. Accepted changes are separated
into Git commits, listed in [WORK_LOG.md](WORK_LOG.md); historical reports and unrelated user files are preserved. No
release, publication, deployment or GUI work was performed.

## 5. Future research

1. Resolve ligand identities and stereochemistry against authoritative experimental
   structures, especially 401/1058 and absolute axial assignments. Recover matched
   source energies/calculation logs before interpreting conformer populations.
2. Test the locked single-versus-ensemble comparison on a genuinely independent,
   larger experimental dataset with more ligand families and measured ee uncertainty.
   The reviewed bisphosphine datasets are promising only after metal/atom mapping,
   conditions and original geometries are verified.
3. Compare physically relevant bound/free-ligand ensembles and energy methods with
   explicit solvent, entropy and temperature assumptions. The present generated
   force-field results show why generic Boltzmann weighting is insufficient evidence.
4. Use certified continuous features to quantify numerical uncertainty when it is
   material. Further boundary refinement is valuable only for applications requiring
   tighter regional bounds than current geometry/experimental uncertainty; general
   performance optimization is not the dominant scientific need exposed here.
