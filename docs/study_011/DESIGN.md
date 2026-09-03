# Study 011 design lock: retrospective ligand ranking

Status: **outcome-independent design, frozen before Study 011 evaluation**.

## Question and estimand

Can a StericX model fitted only to observed ligand outcomes rank a small panel of
previously unseen ligands so that the more enantioselective Ni-hDA ligands are
prioritized? The response is the published absolute transition-state free-energy
difference, `ddG_abs` (kcal/mol); larger is defined as better. The primary estimand
is top-1 recovery across all planned screens, including failed screens.

## Dataset and eligibility

Use the 11 experimentally labeled Ni-catalyzed homo-Diels-Alder ligands already
cached in `data/official/ni_hda_kraken.csv` and aligned to StericX records in
`data/reactions.sigpack`. The authors publicly distribute the source table at the
URL recorded in `data/official/provenance.json`. Unlabeled Kraken rows are ineligible
because they have no outcome against which a ranking can be evaluated.

## Splits fixed without outcomes

`Ligand_Group` is the grouping variable. Enumerate every combination of complete
groups whose combined size is exactly three ligands. This produces 57 exhaustive,
non-random splits: one holdout containing the sole three-member scaffold and 56
holdouts containing three singleton scaffolds. Each split has eight training
observations and a three-ligand candidate panel. No group may cross train/test.
Exact membership is in `split_definitions.json` and row-expanded `splits.csv`.

The exhaustive construction prevents choosing a favorable split. It does not make
the overlapping splits statistically independent. The three-member scaffold enters
one split, whereas each singleton enters 21 splits, so the two strata must also be
reported separately.

## Model and screening protocol

For each split, in this order:

1. Write row-aligned split metadata. Only the eight `train` records may contribute
   outcomes to fitting.
2. Run `stericx fit` with the native physical-organic feature space, BIC forward
   selection, at most two non-intercept terms, the built-in `|r| <= 0.95` guard,
   1,000 bootstrap samples, 1,000 response permutations, and the split seed in
   `seeds.json`. Scaling and selection occur on training rows only.
3. Create a three-row, target-free candidate CSV using an explicit predictor
   allow-list. Run `stericx screen` with `--optimize maximize`, bond-axis Sterimol,
   and the training-derived `max-neighbor` applicability rule. Do not filter
   extrapolations. Rank descending; StericX's identifier tie-break is retained.
4. Finish all 57 fits/screens, write `frozen_predictions.csv`, and hash it.
5. Only then read `ddG_abs`, join outcomes for scoring, and write revealed results.

A split that cannot fit or rank all three candidates remains a failure. It is not
dropped from primary intention-to-screen metrics.

## Metrics fixed in advance

Primary (all 57 planned splits):

- top-1 recovery: whether predicted rank 1 is the observed best ligand; a failed
  screen scores 0;
- top-2 recovery: recall of the observed top two among the predicted top two; a
  failed screen scores 0;
- rank of the observed best ligand (1--3); a failed screen scores 4.

Secondary (successfully ranked splits unless stated otherwise):

- Spearman correlation between predicted and observed response when neither vector
  is constant;
- nDCG over all three non-negative outcomes;
- MAE and RMSE in kcal/mol, both per-split and over all repeated predictions;
- observed outcome of the predicted top 1/top 2 and its gain and ratio relative to
  the candidate-panel mean;
- 95% prediction-interval empirical coverage;
- absolute error by StericX domain verdict and range status, plus descriptive
  association with nearest-training-distance ratio and leverage ratio.

Observed ranks use descending outcome then numeric ligand ID to resolve exact ties.
Pooled errors are descriptive because the same ligand appears in multiple splits.
No hypothesis-test threshold is used to select or suppress a result.

## Random baseline

Enumerate all `3! = 6` possible candidate rankings in every panel. This exact
baseline needs no random seed. Chance expectations are top-1 recovery `1/3`, top-2
recall `2/3`, mean best-ligand rank `2`, and mean Spearman `0`. Enrichment is the
StericX mean divided by the corresponding exact random mean. Outcome-selection
gain is compared with the exact panel-wise random expectation.

## Interpretation constraints

This is a small retrospective test of ranking three measured ligands, not a claim
about recovery from the full 1,566-ligand Kraken library, yield, synthetic
accessibility, prospective performance, or causal ligand effects. Failed and poor
results must remain visible. Because splits overlap, no naive confidence interval
or p-value treating 57 folds as independent will be reported.
