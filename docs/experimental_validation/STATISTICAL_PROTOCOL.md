# Outcome-aware retrospective design, fixed before new descriptor/model evaluation

This is not a prospective registration or blind discovery: this repository and the
investigator have seen previous Ni-hDA outcomes and models. The purpose of this lock
is to prevent adaptive choices during the **new** single-versus-ensemble comparison.

## Primary dataset and endpoint

Cadge et al., JACS 2025, DOI 10.1021/jacs.5c09948, SI Table S3 (PDF page S22),
norbornadiene + methyl vinyl ketone, Ni(cod)2/monophosphine, 353.15 K.
There are twelve quantified ee results and one unquantified result (L9). The author
CSV includes eleven of these twelve quantified ligands. Ligand 906 is absent from
that CSV and currently fails the primary structure API (HTTP 500); it remains an
explicit excluded observation unless independently verified structures become
available before prediction. Do not substitute an invented geometry or identity.
The intended common primary set is the eleven IDs 401,498,723,724,785,1057,1058,
2062,2063,2064,2067 with retrieved published-library structures.

Response: magnitude RT ln((100+ee)/(100-ee)), kcal/mol, R from the frozen engine,
T=353.15 K, using integer ee from SI S3. This is an *apparent* selectivity free-energy
difference under the competing-pathway kinetic interpretation; it does not establish
measured activation barriers or equilibration. The absolute response does not predict
which product enantiomer is favored. Ligand 2064's 3% ee and printed 0.0281 kcal/mol
are inconsistent; use measured ee for the primary endpoint and preserve both fields.
All eleven authors' original ddG_abs values are a predeclared sensitivity endpoint.
No yield endpoint replaces this primary endpoint after results are calculated.

## Primary estimand and partitions

Primary metric: pooled held-out RMSE(C) minus RMSE(B), lower favors ensemble.
Every model uses identical leave-one-family-out (LOFO) outer predictions. Families
are fixed from chemical structures, not labels or outcomes:

- cage phosphine: 401;
- spiro phosphoramidite: 498;
- benzoxaphosphole family: 723,724,1057,1058;
- BINOL phosphoramidite family: 785,2062,2063,2064,2067.

These broad families intentionally test a demanding scaffold transfer. Eleven
leave-one-ligand-out predictions provide a secondary interpolation assessment.
Repeated CV folds are not independent experiments. Exact memberships will be saved.

## Prespecified model comparison

Primary models are equal-complexity univariate OLS fits (intercept + one slope):

- A0: training-response mean (no StericX).
- A: published donor NBO charge `nbo_P_boltz` (electronic baseline).
- B: native maximum adjacent-quadrant occupied-volume difference from the first
  published-library conformer (lowest numeric conformer ID; **not** called lowest energy).
- C: minimum of the same native descriptor over every available conformer.
- D: published Kraken `vbur_max_delta_qvbur_min` reference descriptor.
- E: published reference plus native ensemble descriptor, two slopes; its extra
  capacity and collinearity are explicitly reported.

B/C use identical atomic/radius/center/frame conventions; only aggregation differs.
The single feature is chosen from the published scientific model before new
predictions, not selected using the new held-out errors. D's feature was originally
selected with these experimental outcomes; D is therefore a retrospective reference,
not a clean independently selected prospective benchmark.

Secondary descriptor panels: L, B1, B5, total %Vbur, adjacent-quadrant difference,
quadrant min/max, octant min/max, pyramidalization P and alpha. Single conformer,
ensemble arithmetic mean (equal dimension), and min/max/mean panels are compared
using ridge regression. Fit scaling inside each inner training fold; tune alpha
from [0.01,0.1,1,10,100,1000] using inner leave-one-family-out CV. No global feature
selection. Use the same inner splits for B/C/D/E. Simple one-feature Pearson and
Spearman associations are descriptive and include all values. Random forests and
gradient boosting are exploratory capacity checks only, fixed small settings and
no test-set tuning; omit them if training folds cannot support their complexity
and document the reason. OLS and ridge are the main evidence.

## Uncertainty, failures and sensitivities

Report N, feature/parameter count, MAE, RMSE, R², bias, all held-out predictions,
per-family errors, fold dispersion, and the paired C-B metric difference.
Use 10,000 paired family-cluster bootstrap draws (seed 20261124) for a descriptive
95% interval; only four clusters make it unstable. Exact sign-flip analysis of
family-level loss differences has only 16 permutations; its one-sided minimum
p=1/16 cannot by itself establish 5% significance. Do not claim support merely
from a favorable point estimate. Include ligand-bootstrap sensitivity, clearly
not family-independent. OLS t prediction intervals are conditional on a fixed
linear iid Gaussian-error model and must be labeled as such; report empirical
held-out coverage and width, not a guaranteed domain-general coverage.

Repeat the complete primary comparison with original author ΔΔG, excluding 2064,
and the eleven ligand-wise splits. Preserve all failures; never drop a failed fold
from headline scoring. Common complete-descriptor eligibility is fixed across
models, and exclusions and missing conformers are explicit.

Published-library conformers have no demonstrated historical energy provenance;
no Boltzmann weighting is inferred from descriptor aggregates or conformer IDs.
Generated ETKDG/MMFF94 or archived CREST/xTB structures form separate analyses,
using documented energies and approximate weights only when complete. Temperature
sensitivity uses 298.15/353.15/400 K for those approximate weights, while the
experimental response remains at its reported reaction temperature. Geometry,
optimization and incomplete coverage sensitivity must not be interpreted as
experimentally measured populations.

## Pre-prediction data-quality amendment (2026-11-24)

After native descriptors were computed, but **before any new experimental model was
fit or scored**, a mapped geometric audit found a phosphorus-handedness disagreement
for all 15 ligand-401 exports and a mixture for ligand 1058 (6 opposite, 8 matching)
relative to the encoded source SMILES. Ordinary RDKit 3D assignment had left these
centers unassigned. Independent stereo-constrained source embeddings and an explicit
virtual-point diagnostic agree after accounting for implicit-neighbor parity.
An initial diagnostic without that parity conversion misleadingly implicated 1057;
that exploratory interpretation was rejected. 1057's five structures match.

The original eleven-ligand comparison will still be reported as the **as-distributed
library association**, with this limitation, rather than silently replacing its
structures or dropping unfavorable results. Add a common nine-ligand sensitivity
excluding 401 and 1058, and a filtered-published sensitivity retaining only the eight
matching 1058 conformers and excluding 401. Neither is a new primary endpoint.
Unknown absolute axial/spiro assignments remain limitations. A statistical improvement
on the as-distributed lane alone cannot establish a conformer-population hypothesis
for the experimentally used stereoisomers.

Generated structures will preserve source-encoded stereochemistry and separately
record relative handedness filters for source-unencoded P/axial geometry. This is a
separate structure-preparation experiment, including a deliberate stereoisomer
correction for 401; it is not a silent replacement of the primary input library.

## Post-primary preparation sensitivity (2026-11-24)

After primary modeling, inspection of the separately generated lane found that
standard RDKit aromaticity marks P-containing BINOL rings aromatic and hides the
single biaryl bond from the initial relative-handedness filter. Preserve generated
v1 and all its outcomes as a preparation experiment with incomplete stereo control.
Generated v2 uses the original SDF bond orders and MDL aromaticity/hybridization
before ETKDG, then the same MMFF94 optimization, seeds, budgets and RMS criteria.
It transfers source-encoded R/S centers explicitly, filters indexed biaryl/spiro
handedness relative to the first original geometry, and evaluates serialized-geometry
energies using the same preparation. MMFF's own aromaticity typing is recorded as
an implementation convention. This is a disclosed post-primary quality correction,
not a newly registered confirmatory comparison. No generated result replaces the
published primary endpoint, feature or split.
