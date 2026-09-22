# Reaction screening: experimental data to candidate deck

This tutorial walks through the complete StericX v0.3 reaction-screening workflow: fit a
model to experimental observations, inspect what was learned, rank a ligand library, assess
uncertainty and model applicability, remove ligands that have already been tested, choose a
diverse subset, and export an auditable candidate deck.

The commands use the repository's Ni-catalyzed homo-Diels–Alder (Ni-hDA) example. They run
offline after the repository has been cloned; no data are downloaded. Run them from the
repository root.

> **This is a workflow tutorial, not a claim that this small model is ready to direct an
> experimental campaign.** The dataset contains only 11 outcomes. The retrospective ranking
> study found below-random top-1 recovery, negative mean rank correlation, and failed fits;
> see [Study 011](study_011/STUDY_011.md). Those results remain part of the interpretation.

## The example data and split

The tutorial preserves the checked-in historical evidence and creates fresh records:

- [`data/reactions_raw.csv`](../data/reactions_raw.csv) is the historical source.
  The preparation helper creates
  `.stericx/scientific_remediation/ni_hda_response_metadata_v2/reactions.csv`,
  correcting response-temperature metadata to **353.15 K**. It preserves the
  published targets and supplied historical **298.15 K** conformer weights; it
  does not re-equilibrate those populations at reaction conditions. The ligand
  2064 ee/energy source conflict remains unresolved.
- `stericx parse` recomputes descriptors from every listed conformer with the
  supplied weights and explicit atom axes into a new tutorial sigpack.
- [`docs/examples/ni_hda_screening_split.csv`](examples/ni_hda_screening_split.csv)
  is the row-aligned `S001` split predefined by Study 011. Eight observations are
  marked `train`; three singleton scaffold groups are marked `blind`.

`stericx fit` learns descriptor selection, scaling, and coefficients from `train` rows only.
It writes predictions for the three `blind` rows without evaluating their known outcomes.
This is a retrospective exercise, so the responses exist in the packed source data, but the
split prevents them from entering the fit.

Build the corrected implementation and prepare a new destination. The helper
requires a new input destination; if `ni_hda_response_metadata_v2` already exists,
verify its receipt and reuse it rather than overwriting it. Use a new demo suffix
for a repeated tutorial run.

```bash
cargo build --release
python3 scripts/prepare_remediation_inputs.py \
  --output .stericx/scientific_remediation/ni_hda_response_metadata_v2
mkdir -p .stericx/scientific_remediation/demo
mkdir .stericx/scientific_remediation/demo/screening_v1

./target/release/stericx parse \
  --csv .stericx/scientific_remediation/ni_hda_response_metadata_v2/reactions.csv \
  --xyz-dir data \
  --output .stericx/scientific_remediation/demo/screening_v1/reactions.sigpack
```

The fresh sigpack is deliberately separate from the historical
`data/reactions.sigpack`. Scientific descriptor corrections can change fitted
coefficients and predictions; old study numbers are retained as historical
results, not substituted for this run.

## 1. Fit a model

```bash
./target/release/stericx fit \
  --data .stericx/scientific_remediation/demo/screening_v1/reactions.sigpack \
  --metadata docs/examples/ni_hda_screening_split.csv \
  --output .stericx/scientific_remediation/demo/screening_v1/fit-report.json \
  --predictions .stericx/scientific_remediation/demo/screening_v1/frozen-predictions.csv \
  --portable-model .stericx/scientific_remediation/demo/screening_v1/model.json \
  --model-id ni-hda-screening-tutorial \
  --reaction-family "Ni-catalyzed homo-Diels-Alder" \
  --catalyst-metal Ni \
  --ligand-class "monodentate phosphorus(III)" \
  --source-url "https://raw.githubusercontent.com/SigmanGroup/Ni-Catalyzed-hDA/main/data/kraken.csv" \
  --descriptor-aggregation supplied_weight_mean \
  --response-temp-k 353.15 \
  --response-sign-convention "Magnitude |ddG| from ddG_abs; larger values mean greater enantioselectivity; no R/S assignment." \
  --optimize maximize \
  --bootstrap 1000 \
  --permutations 500 \
  --seed 20260904
```

The important outputs are:

- `fit-report.json`: feature selection, coefficients, validation metrics, baselines, and
  diagnostics;
- `frozen-predictions.csv`: predictions for the three rows withheld by the split;
- `model.json`: the portable model consumed by `model inspect` and `screen`.

The explicit Study 011 seed makes bootstrap and permutation results repeatable. This is the
first predefined split, not a split selected after comparing its performance. `--optimize maximize`
records the chemical objective: a larger predicted |ΔΔG‡| ranks higher. The explicit
response convention records that `ddG_abs` predicts magnitude, not the favored enantiomer.
When the flag is omitted, new models state that the convention is unspecified. Change the optimization
direction only when the response definition for your reaction requires it.

## 2. Inspect and validate the fitted model

First check that the saved document is internally consistent, then read its scientific
summary:

```bash
./target/release/stericx model validate \
  .stericx/scientific_remediation/demo/screening_v1/model.json

./target/release/stericx model inspect \
  .stericx/scientific_remediation/demo/screening_v1/model.json
```

Do not skip the validation block in `model inspect`. In particular, compare training R² with
leave-one-out Q²/RMSE and leave-group-out Q². Training fit alone answers how well the model
describes observations it already saw; held-out diagnostics are more relevant to screening.
Read the values from this fresh model. The earlier tutorial's numerical metrics
came from the historical descriptor implementation and are not a reference for
this corrected run. The reported LOO and group-LOO diagnostics condition on the
selected feature set; feature selection is not repeated within these diagnostics.
They are not a full-pipeline cross-validation estimate.
Also inspect:

- the selected descriptors and coefficient signs;
- the number of training observations and scaffold groups;
- descriptor ranges used by the applicability-domain check;
- the reaction, target, units, optimization direction, and dataset digests;
- warnings or fields reported as `not_recorded`.

For machine-readable inspection and validation:

```bash
./target/release/stericx model inspect \
  .stericx/scientific_remediation/demo/screening_v1/model.json --format json \
  > .stericx/scientific_remediation/demo/screening_v1/model-summary.json

./target/release/stericx model validate \
  .stericx/scientific_remediation/demo/screening_v1/model.json --format json \
  > .stericx/scientific_remediation/demo/screening_v1/model-validation.json
```

## 3. Screen the ligand library

Use the reaction CSV as the example library:

```bash
./target/release/stericx screen \
  .stericx/scientific_remediation/demo/screening_v1/model.json \
  --library .stericx/scientific_remediation/ni_hda_response_metadata_v2/reactions.csv \
  --temperature 353.15 \
  --top 11
```

This first pass intentionally shows every row, including training ligands. It is useful for
checking descriptor mapping and the shape of the report before applying campaign filters.
The model declares `supplied_weight_mean`. When selected terms require missing
steric descriptors, geometry-driven screening uses all `Conformer_XYZ_Paths`,
explicit `Conformer_Boltzmann_Weights`, and the row
`Attach_Atom_Idx`/`Primary_Bond_Vector_Idx` through the same aggregation routine as
parsing. A model selecting electronics alone consumes those supplied columns
and need not recompute unused geometry descriptors. The separate full-record
example in [MODEL_FORMAT.md](MODEL_FORMAT.md#producing-a-document) exercises
the declared ensemble route when steric terms are selected. Screening does not
replace the ensemble with `Ligand_XYZ_Path`. A directory of
single geometries is incompatible with this contract. Precomputed descriptor
columns are also accepted, with source consistency remaining the caller's responsibility; mixed supplied/computed values are identified per descriptor.

The example library contains an `Exp_ddG_kcal_mol` column because it began as an experimental
dataset. `screen` does not read that response when making predictions, and candidate decks
use an explicit output allow-list that cannot copy it. For a prospective screen, prefer a
target-free library anyway: it makes the separation between prediction and evaluation
obvious to reviewers.

Save the complete, unrounded result rather than scraping the terminal table:

```bash
./target/release/stericx screen \
  .stericx/scientific_remediation/demo/screening_v1/model.json \
  --library .stericx/scientific_remediation/ni_hda_response_metadata_v2/reactions.csv \
  --temperature 353.15 \
  --format json \
  > .stericx/scientific_remediation/demo/screening_v1/ranking-all.json

./target/release/stericx screen \
  .stericx/scientific_remediation/demo/screening_v1/model.json \
  --library .stericx/scientific_remediation/ni_hda_response_metadata_v2/reactions.csv \
  --temperature 353.15 \
  --format csv \
  > .stericx/scientific_remediation/demo/screening_v1/ranking-all.csv
```

## 4. Understand the rankings

`rank = 1` means the candidate is most desirable under the model's recorded optimization
direction. It does not mean that ligand has the largest measured response. In this tutorial,
`maximize` gives a descending ranking by predicted ΔΔG‡. StericX also reports the implied ee
at the explicitly requested 353.15 K, but the raw predicted ΔΔG‡ remains the model output.

Read a ranking row in this order:

1. `predicted_ddg_kcal_mol` is the point prediction used for ranking.
2. `predicted_ee_percent` is a temperature-dependent interpretation of that prediction.
3. `prediction_interval_*` and `uncertainty` describe different sources of uncertainty.
4. `domain_verdict`, `trust`, leverage, and nearest-neighbour fields describe how far the
   candidate is from the training information.

`--top N` is applied after every eligible ligand has been predicted and ranked. The report
therefore distinguishes `library_size`, `screened`, and `returned`. Equal scores are resolved
deterministically by stable identifier and then library position.

The model's stored direction is the default. `--ascending` or `--descending` can override it,
but the report records `ranking_overridden` and warns when the override conflicts with what
the model calls better. Treat an override as a change in the scientific question, not a
display preference.

## 5. Interpret uncertainty

StericX reports two useful but different intervals:

- `prediction_interval_low` / `prediction_interval_high` is a nominal 95% Student-t interval for a
  new observation, conditional on the fixed linear model and IID homoscedastic
  normal errors. It includes residual scatter and widens with leverage; empirical
  chemical coverage is not guaranteed.
- `uncertainty.lower` / `uncertainty.upper` is a percentile-bootstrap interval for the mean
  response caused by coefficient uncertainty. It does not include residual scatter and does
  not diagnose extrapolation.

The older `coefficient_band` propagates marginal coefficient intervals independently.
It discards coefficient correlation, has no guaranteed joint coverage and need not
be conservative. Use the joint bootstrap distribution for the stated coefficient
uncertainty calculation.

Wide or strongly overlapping intervals mean that the order of nearby candidates is not
well resolved. A precise-looking bootstrap interval does not rescue an extrapolative
candidate: coefficient uncertainty and applicability answer separate questions.

## 6. Interpret applicability-domain warnings

Applicability is computed without looking at the prediction. StericX combines descriptor
range checks with distance and leverage relative to the training set:

- `interpolation`: inside every selected-descriptor range and near sampled training space;
- `sparse_interpolation`: inside individual ranges but farther from a training neighbour
  than the selected training-derived spacing rule allows;
- `extrapolation`: outside at least one selected-descriptor range;
- `unknown`: the model lacks enough recorded training geometry for the full assessment.

Read the actual verdicts from the corrected run. `interpolation` means the
chosen descriptor-space checks find support; it does not establish an accurate
ranking or chemically calibrated uncertainty.

The compatibility field `trust` now uses descriptive labels such as
`inside_range:ordinary_leverage` and `outside_range:high_leverage`, with
`range_only:*` when leverage is unavailable. It no longer claims `reliable`.
`nearest_training_ligand`, `nearest_training_distance`, `nearest_training_ratio`,
`leverage`, and `leverage_ratio` show the measured evidence. Ordinary Mahalanobis
distance is unavailable for singular covariance or missing stored training
points; `mahalanobis_unavailable` explains why. JSON `descriptors[].source` and
CSV `descriptor_sources` distinguish supplied values from computed values.

The default neighbour boundary is the largest nearest-neighbour spacing in the training set.
You can inspect a stricter training-derived rule without refitting:

```bash
./target/release/stericx screen \
  .stericx/scientific_remediation/demo/screening_v1/model.json \
  --library .stericx/scientific_remediation/ni_hda_response_metadata_v2/reactions.csv \
  --domain-rule mean-plus-2sd \
  --temperature 353.15 \
  --top 11
```

`--in-domain-only` keeps only `interpolation` candidates. Use it only after examining the
unfiltered screen: filtering first can conceal that the nominally best predictions are all
unsupported. The rule derivation and model fields are documented in
[`MODEL_FORMAT.md`](MODEL_FORMAT.md).

## 7. Exclude ligands that have already been tested

The tutorial exclusion list contains the exact eight `S001` training identifiers:
[`docs/examples/ni_hda_screening_tested.csv`](examples/ni_hda_screening_tested.csv).
Applying it leaves the three held-out singleton-scaffold ligands as the candidate set.

```bash
./target/release/stericx screen \
  .stericx/scientific_remediation/demo/screening_v1/model.json \
  --library .stericx/scientific_remediation/ni_hda_response_metadata_v2/reactions.csv \
  --exclude-tested docs/examples/ni_hda_screening_tested.csv \
  --temperature 353.15 \
  --top 3
```

Matching uses a stable identifier such as `Reaction_ID` first. Exact SMILES text is a
secondary fallback, not chemical structure matching. The report counts matched and excluded
rows and prints every unresolved identifier. The exclusion file is hashed into screening
provenance because it changes the candidate set.

Maintain this file as an experimental ledger: it should include the model's training
ligands and every ligand tested since the model was frozen.

## 8. Select a diverse candidate set

Plain top-k can spend a limited experimental budget on nearly identical candidates. Add
`--diverse` to choose a subset by greedy maximal marginal relevance in the same standardized
descriptor space used for applicability:

```bash
./target/release/stericx screen \
  .stericx/scientific_remediation/demo/screening_v1/model.json \
  --library .stericx/scientific_remediation/ni_hda_response_metadata_v2/reactions.csv \
  --exclude-tested docs/examples/ni_hda_screening_tested.csv \
  --temperature 353.15 \
  --top 2 \
  --diverse \
  --diversity-weight 0.5
```

The first selection is the ordinary rank-1 candidate. Later selections balance rank-based
desirability with distance from candidates already selected. A weight of `0` reproduces the
ordinary top-k; a weight of `1` ignores desirability after the first pick. The output retains
each ligand's original `rank` and adds the selection step, separation, and objective, so the
diversity tradeoff remains visible.

Diversity is not an applicability waiver. A candidate can be different from the rest of the
deck and still be outside the training domain.

## 9. Export a candidate deck

First capture the final, machine-readable screen report:

```bash
./target/release/stericx screen \
  .stericx/scientific_remediation/demo/screening_v1/model.json \
  --library .stericx/scientific_remediation/ni_hda_response_metadata_v2/reactions.csv \
  --exclude-tested docs/examples/ni_hda_screening_tested.csv \
  --temperature 353.15 \
  --top 2 \
  --diverse \
  --diversity-weight 0.5 \
  --format json \
  > .stericx/scientific_remediation/demo/screening_v1/final-screen-report.json
```

Then run the same deterministic selection with deck export enabled:

```bash
./target/release/stericx screen \
  .stericx/scientific_remediation/demo/screening_v1/model.json \
  --library .stericx/scientific_remediation/ni_hda_response_metadata_v2/reactions.csv \
  --exclude-tested docs/examples/ni_hda_screening_tested.csv \
  --temperature 353.15 \
  --top 2 \
  --diverse \
  --diversity-weight 0.5 \
  --export-deck .stericx/scientific_remediation/demo/screening_v1/candidate-deck.csv \
  > .stericx/scientific_remediation/demo/screening_v1/deck-export.txt
```

This writes:

- `candidate-deck.csv`, one review-ready row per selected ligand;
- `candidate-deck.meta.json`, containing the configuration, model and library hashes,
  exclusion accounting, model context, uncertainty method, and deck hash;
- `final-screen-report.json`, the full screening report captured from standard output.
- `deck-export.txt`, the human-readable record of the export run.

The deck contains identifiers, predictions, intervals, descriptors, applicability evidence,
selection rationale, and provenance hashes. It contains no experimental response. Identical
inputs produce identical deck and sidecar bytes.

Before sending a deck to the laboratory, archive the model, target-free library, exclusion
list, command, JSON report, CSV deck, and metadata sidecar together. Freeze them before
measuring outcomes.

## What StericX rankings do not mean

- **A prediction is not experimental evidence.** A ranked deck proposes experiments; it does
  not establish yield, selectivity, activity, robustness, or mechanism.
- **Extrapolation is riskier.** A high predicted value outside the training domain may be an
  artifact of extending a fitted relationship beyond observed chemistry. Domain warnings
  qualify the prediction; they are not optional decoration.
- **Performance depends on the reaction and training data.** A model fitted for one catalyst,
  substrate family, response definition, or condition set does not automatically transfer to
  another. Small, narrow, noisy, or confounded datasets support weaker rankings.
- **A high rank does not guarantee experimental success.** Rank is relative to the screened
  candidates under one fitted model. It is not a calibrated probability of success and does
  not account for unmodeled chemistry, availability, handling, or experimental failure.

The repository's own evidence demonstrates these limits. In the historical audit, across all 57 predefined screens
in [Study 011](study_011/STUDY_011.md), StericX recovered the best held-out ligand at rank 1
in 0.158 of screens versus 0.333 for random selection. Among the 47 successful rankings,
mean Spearman correlation was negative; the other 10 predefined splits failed to produce a
usable fit. A transparent screening workflow preserves such failures rather than treating
the ranking interface as evidence that a model is predictive.

## Adapting the workflow to your reaction

Replace the example files while preserving the separation of roles:

1. Pack row-aligned structures and experimental responses with `stericx parse`.
2. Put only allowed observations in the `train` partition; keep evaluation rows `blind`,
   `test`, or `external` until predictions are frozen.
3. Record `Ligand_Group` so leave-group-out validation reflects scaffold transfer.
4. Record the justified descriptor aggregation, response temperature, reaction context
   and optimization direction when fitting the portable model; retain population provenance.
5. Screen a target-free library that supplies every descriptor selected by the model.
6. Review uncertainty and applicability before excluding or filtering candidates.
7. Exclude every previously tested ligand, then apply diversity within the eligible pool.
8. Export and archive the deck plus sidecar before revealing or generating new outcomes.

If the fitted model selected `nbo_charge`, `ir_frequency`, or an interaction involving one of
them, a geometry-only library is insufficient. StericX reports the missing descriptor rather
than substituting a value. If a saved model lacks recorded training geometry, leverage-based
intervals and full applicability assessment are unavailable; refit the model to add them.

The corrected tutorial is checked independently of historical study artifacts;
see [model remediation receipts](scientific_remediation/models/RESULTS.md).
Legacy schema 1/2 models without aggregation metadata report `unknown` and reject
geometry substitution. Use matching precomputed inputs or refit from a documented
method; changing the metadata alone does not validate old descriptor values.
