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

The tutorial uses three checked-in inputs:

- [`data/reactions.sigpack`](../data/reactions.sigpack) contains the 11 packed descriptor and
  response records.
- [`data/reactions_raw.csv`](../data/reactions_raw.csv) provides reaction identifiers,
  structures, electronic descriptors, groups, outcomes, and source provenance.
- [`docs/examples/ni_hda_screening_split.csv`](examples/ni_hda_screening_split.csv) is the
  row-aligned `S001` split predefined by Study 011. Eight observations are marked `train`;
  three singleton scaffold groups are marked `blind`.

`stericx fit` learns descriptor selection, scaling, and coefficients from `train` rows only.
It writes predictions for the three `blind` rows without evaluating their known outcomes.
This is a retrospective exercise, so the responses exist in the packed source data, but the
split prevents them from entering the fit.

Build StericX and create an ignored working directory for the generated files:

```bash
cargo build --release
mkdir -p .stericx/reaction-screening-tutorial
```

## 1. Fit a model

```bash
./target/release/stericx fit \
  --data data/reactions.sigpack \
  --metadata docs/examples/ni_hda_screening_split.csv \
  --output .stericx/reaction-screening-tutorial/fit-report.json \
  --predictions .stericx/reaction-screening-tutorial/frozen-predictions.csv \
  --portable-model .stericx/reaction-screening-tutorial/model.json \
  --model-id ni-hda-screening-tutorial \
  --reaction-family "Ni-catalyzed homo-Diels-Alder" \
  --catalyst-metal Ni \
  --ligand-class "monodentate phosphorus(III)" \
  --source-url "https://raw.githubusercontent.com/SigmanGroup/Ni-Catalyzed-hDA/main/data/kraken.csv" \
  --response-temp-k 298.15 \
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
records the chemical objective: a larger predicted ΔΔG‡ ranks higher. Change the optimization
direction only when the response definition for your reaction requires it.

## 2. Inspect and validate the fitted model

First check that the saved document is internally consistent, then read its scientific
summary:

```bash
./target/release/stericx model validate \
  .stericx/reaction-screening-tutorial/model.json

./target/release/stericx model inspect \
  .stericx/reaction-screening-tutorial/model.json
```

Do not skip the validation block in `model inspect`. In particular, compare training R² with
leave-one-out Q²/RMSE and leave-group-out Q². Training fit alone answers how well the model
describes observations it already saw; held-out diagnostics are more relevant to screening.
The exact tutorial run reports training R² `0.6929` and ordinary LOO Q² `0.4625`, but
group-LOO Q² `-192.5603`. That extreme group-validation failure is a deployment-stopping
warning: the apparently reasonable row-wise fit does not transfer across the training
scaffold groups.
Also inspect:

- the selected descriptors and coefficient signs;
- the number of training observations and scaffold groups;
- descriptor ranges used by the applicability-domain check;
- the reaction, target, units, optimization direction, and dataset digests;
- warnings or fields reported as `not_recorded`.

For machine-readable inspection and validation:

```bash
./target/release/stericx model inspect \
  .stericx/reaction-screening-tutorial/model.json --format json \
  > .stericx/reaction-screening-tutorial/model-summary.json

./target/release/stericx model validate \
  .stericx/reaction-screening-tutorial/model.json --format json \
  > .stericx/reaction-screening-tutorial/model-validation.json
```

## 3. Screen the ligand library

Use the reaction CSV as the example library:

```bash
./target/release/stericx screen \
  .stericx/reaction-screening-tutorial/model.json \
  --library data/reactions_raw.csv \
  --top 11
```

This first pass intentionally shows every row, including training ligands. It is useful for
checking descriptor mapping and the shape of the report before applying campaign filters.
In a real campaign, the library can instead be a reaction CSV, a descriptor CSV, or a
directory of ligand geometries.

The example library contains an `Exp_ddG_kcal_mol` column because it began as an experimental
dataset. `screen` does not read that response when making predictions, and candidate decks
use an explicit output allow-list that cannot copy it. For a prospective screen, prefer a
target-free library anyway: it makes the separation between prediction and evaluation
obvious to reviewers.

Save the complete, unrounded result rather than scraping the terminal table:

```bash
./target/release/stericx screen \
  .stericx/reaction-screening-tutorial/model.json \
  --library data/reactions_raw.csv \
  --format json \
  > .stericx/reaction-screening-tutorial/ranking-all.json

./target/release/stericx screen \
  .stericx/reaction-screening-tutorial/model.json \
  --library data/reactions_raw.csv \
  --format csv \
  > .stericx/reaction-screening-tutorial/ranking-all.csv
```

## 4. Understand the rankings

`rank = 1` means the candidate is most desirable under the model's recorded optimization
direction. It does not mean that ligand has the largest measured response. In this tutorial,
`maximize` gives a descending ranking by predicted ΔΔG‡. StericX also reports the implied ee
at 298.15 K, but the raw predicted ΔΔG‡ remains the model output.

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

- `prediction_interval_low` / `prediction_interval_high` is a 95% Student-t interval for a
  new observation. It includes residual scatter and widens with leverage. Use this as the
  primary interval when asking how variable an experimental outcome may be.
- `uncertainty.lower` / `uncertainty.upper` is a percentile-bootstrap interval for the mean
  response caused by coefficient uncertainty. It does not include residual scatter and does
  not diagnose extrapolation.

The older `coefficient_band` propagates marginal coefficient intervals independently. It is
conservative and discards coefficient correlation, so it is a weaker signal than the joint
bootstrap interval.

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

All three `S001` held-out candidates happen to report `interpolation` and `reliable` in this
tutorial. That means the model is not extrapolating for those descriptor values; it does not
mean the ranking is accurate. Applicability is a support diagnostic, not a validation score.

The `trust` field combines range and leverage information. Read
`do_not_trust:extrapolation` literally; it is not a lower-confidence version of an otherwise
supported estimate. `nearest_training_ligand`, `nearest_training_distance`,
`nearest_training_ratio`, `leverage`, and `leverage_ratio` show why a warning was assigned.

The default neighbour boundary is the largest nearest-neighbour spacing in the training set.
You can inspect a stricter training-derived rule without refitting:

```bash
./target/release/stericx screen \
  .stericx/reaction-screening-tutorial/model.json \
  --library data/reactions_raw.csv \
  --domain-rule mean-plus-2sd \
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
  .stericx/reaction-screening-tutorial/model.json \
  --library data/reactions_raw.csv \
  --exclude-tested docs/examples/ni_hda_screening_tested.csv \
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
  .stericx/reaction-screening-tutorial/model.json \
  --library data/reactions_raw.csv \
  --exclude-tested docs/examples/ni_hda_screening_tested.csv \
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
  .stericx/reaction-screening-tutorial/model.json \
  --library data/reactions_raw.csv \
  --exclude-tested docs/examples/ni_hda_screening_tested.csv \
  --top 2 \
  --diverse \
  --diversity-weight 0.5 \
  --format json \
  > .stericx/reaction-screening-tutorial/final-screen-report.json
```

Then run the same deterministic selection with deck export enabled:

```bash
./target/release/stericx screen \
  .stericx/reaction-screening-tutorial/model.json \
  --library data/reactions_raw.csv \
  --exclude-tested docs/examples/ni_hda_screening_tested.csv \
  --top 2 \
  --diverse \
  --diversity-weight 0.5 \
  --export-deck .stericx/reaction-screening-tutorial/candidate-deck.csv \
  > .stericx/reaction-screening-tutorial/deck-export.txt
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

The repository's own evidence demonstrates these limits. Across all 57 predefined screens
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
4. Supply the reaction context and correct optimization direction when fitting the portable
   model.
5. Screen a target-free library that supplies every descriptor selected by the model.
6. Review uncertainty and applicability before excluding or filtering candidates.
7. Exclude every previously tested ligand, then apply diversity within the eligible pool.
8. Export and archive the deck plus sidecar before revealing or generating new outcomes.

If the fitted model selected `nbo_charge`, `ir_frequency`, or an interaction involving one of
them, a geometry-only library is insufficient. StericX reports the missing descriptor rather
than substituting a value. If a saved model lacks recorded training geometry, leverage-based
intervals and full applicability assessment are unavailable; refit the model to add them.
