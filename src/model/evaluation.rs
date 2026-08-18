//! Scoring of frozen predictions after their experimental targets are revealed.
//!
//! Evaluation is deliberately separate from training: predictions are produced
//! and hashed first, and only a later call to [`score_frozen_predictions`] joins
//! them to measured values. The scorer re-derives every prediction from the
//! stored model weights and refuses inputs that disagree, so a frozen artifact
//! cannot be silently replaced by a refitted one.

use super::{RegressXPredictor, ScientificFitReport};
use crate::model::dataset::ReactionLabel;
use crate::model::training::FrozenPrediction;
use crate::storage::PackedReactionRecord;
use serde::Serialize;

/// One revealed prediction paired with its experimental target.
#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct ScoredPrediction {
    pub reaction_id: String,
    pub ligand_group: String,
    pub dataset_split: String,
    pub predicted_ddg_kcal_mol: f32,
    pub experimental_ddg_kcal_mol: f32,
    pub residual_kcal_mol: f32,
    pub applicability_domain: String,
}

/// Aggregate accuracy of a frozen prediction set against revealed targets.
#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct EvaluationSummary {
    pub evaluated_records: usize,
    pub mae_kcal_mol: f64,
    pub rmse_kcal_mol: f64,
    /// `None` when the revealed targets carry no variance to explain.
    pub r2: Option<f64>,
    pub applicability_warnings: usize,
    pub scored_predictions: Vec<ScoredPrediction>,
}

/// Joins frozen predictions to revealed targets and scores them.
///
/// `labels` must be positionally aligned with `records`. Each frozen row is
/// matched by `Reaction_ID`, rejected if it belongs to the training split, and
/// re-predicted from `model` to confirm it was produced by that exact model.
pub fn score_frozen_predictions(
    records: &[PackedReactionRecord],
    labels: &[ReactionLabel],
    model: &ScientificFitReport,
    frozen: &[FrozenPrediction],
) -> Result<EvaluationSummary, String> {
    if labels.len() != records.len() {
        return Err(format!(
            "metadata has {} rows but sigpack contains {} records",
            labels.len(),
            records.len()
        ));
    }
    if frozen.is_empty() {
        return Err("frozen prediction file contains no rows".into());
    }

    let predictor = RegressXPredictor::new(model.weights);
    let mut scored = Vec::with_capacity(frozen.len());
    let mut actual = Vec::with_capacity(frozen.len());
    let mut predicted = Vec::with_capacity(frozen.len());
    for frozen_row in frozen {
        let index = labels
            .iter()
            .position(|label| label.reaction_id == frozen_row.reaction_id)
            .ok_or_else(|| {
                format!(
                    "frozen reaction {} is absent from metadata",
                    frozen_row.reaction_id
                )
            })?;
        if labels[index].is_training() {
            return Err(format!(
                "frozen reaction {} is marked as training data",
                frozen_row.reaction_id
            ));
        }
        let recomputed = predictor.predict(&records[index]);
        if (recomputed - frozen_row.predicted_ddg).abs() > 1.0e-4 {
            return Err(format!(
                "frozen prediction for {} does not match the supplied model",
                frozen_row.reaction_id
            ));
        }
        let experimental = records[index].exp_ddg;
        if !experimental.is_finite() {
            return Err(format!(
                "reaction {} has no finite revealed target",
                frozen_row.reaction_id
            ));
        }
        actual.push(f64::from(experimental));
        predicted.push(f64::from(frozen_row.predicted_ddg));
        scored.push(ScoredPrediction {
            reaction_id: frozen_row.reaction_id.clone(),
            ligand_group: frozen_row.ligand_group.clone(),
            dataset_split: frozen_row.dataset_split.clone(),
            predicted_ddg_kcal_mol: frozen_row.predicted_ddg,
            experimental_ddg_kcal_mol: experimental,
            residual_kcal_mol: frozen_row.predicted_ddg - experimental,
            applicability_domain: frozen_row.applicability_domain.clone(),
        });
    }

    let residual_sum = actual
        .iter()
        .zip(&predicted)
        .map(|(actual, predicted)| (actual - predicted).powi(2))
        .sum::<f64>();
    let mean = actual.iter().sum::<f64>() / actual.len() as f64;
    let total_sum = actual
        .iter()
        .map(|actual| (actual - mean).powi(2))
        .sum::<f64>();
    Ok(EvaluationSummary {
        evaluated_records: actual.len(),
        mae_kcal_mol: actual
            .iter()
            .zip(&predicted)
            .map(|(actual, predicted)| (actual - predicted).abs())
            .sum::<f64>()
            / actual.len() as f64,
        rmse_kcal_mol: (residual_sum / actual.len() as f64).sqrt(),
        r2: (total_sum > f64::EPSILON).then_some(1.0 - residual_sum / total_sum),
        applicability_warnings: frozen
            .iter()
            .filter(|row| row.applicability_domain != "inside_training_range")
            .count(),
        scored_predictions: scored,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::training::train_scientific_model;
    use crate::model::{FitOptions, TrainedModel};

    fn fixture() -> (Vec<PackedReactionRecord>, Vec<ReactionLabel>, TrainedModel) {
        let records = (0..12)
            .map(|index| {
                let l = 1.0 + index as f32 * 0.25;
                PackedReactionRecord {
                    l,
                    b1: 1.8,
                    b5: 3.4,
                    nbo_charge: -0.35,
                    ir_freq: 1_650.0,
                    temp_k: 298.15,
                    exp_ddg: 0.5 + 1.25 * l,
                    ..PackedReactionRecord::default()
                }
            })
            .collect::<Vec<_>>();
        let labels = (0..records.len())
            .map(|index| {
                let split = if index < records.len() - 2 {
                    "train"
                } else {
                    "blind"
                };
                ReactionLabel::new(format!("R{index}"), split, format!("group_{}", index % 3))
            })
            .collect::<Vec<_>>();
        let options = FitOptions {
            bootstrap_samples: 25,
            permutation_samples: 25,
            ..FitOptions::default()
        };
        let trained = train_scientific_model(&records, &labels, options).unwrap();
        (records, labels, trained)
    }

    #[test]
    fn scores_revealed_targets() {
        let (records, labels, trained) = fixture();

        let summary = score_frozen_predictions(
            &records,
            &labels,
            &trained.report,
            &trained.frozen_predictions,
        )
        .unwrap();

        assert_eq!(summary.evaluated_records, 2);
        assert_eq!(summary.scored_predictions.len(), 2);
        // Both held-out rows extend past the largest training `L`, so the
        // domain check must flag them even though the fit extrapolates well.
        assert_eq!(summary.applicability_warnings, 2);
        assert!(
            summary
                .scored_predictions
                .iter()
                .all(|row| row.applicability_domain.starts_with("outside:"))
        );
        assert!(summary.mae_kcal_mol < 0.05);
        assert!(summary.rmse_kcal_mol >= summary.mae_kcal_mol);
        for row in &summary.scored_predictions {
            assert!(
                (row.residual_kcal_mol
                    - (row.predicted_ddg_kcal_mol - row.experimental_ddg_kcal_mol))
                    .abs()
                    < 1.0e-9
            );
        }
    }

    #[test]
    fn rejects_predictions_that_do_not_match_the_model() {
        let (records, labels, trained) = fixture();
        let mut tampered = trained.frozen_predictions.clone();
        tampered[0].predicted_ddg += 1.0;

        assert_eq!(
            score_frozen_predictions(&records, &labels, &trained.report, &tampered).unwrap_err(),
            "frozen prediction for R10 does not match the supplied model"
        );
    }

    #[test]
    fn rejects_training_rows_and_unknown_reactions() {
        let (records, labels, trained) = fixture();

        let mut training_row = trained.frozen_predictions.clone();
        training_row[0].reaction_id = "R0".into();
        assert_eq!(
            score_frozen_predictions(&records, &labels, &trained.report, &training_row)
                .unwrap_err(),
            "frozen reaction R0 is marked as training data"
        );

        let mut unknown = trained.frozen_predictions.clone();
        unknown[0].reaction_id = "R999".into();
        assert_eq!(
            score_frozen_predictions(&records, &labels, &trained.report, &unknown).unwrap_err(),
            "frozen reaction R999 is absent from metadata"
        );
    }

    #[test]
    fn rejects_an_empty_frozen_set() {
        let (records, labels, trained) = fixture();

        assert_eq!(
            score_frozen_predictions(&records, &labels, &trained.report, &[]).unwrap_err(),
            "frozen prediction file contains no rows"
        );
    }
}
