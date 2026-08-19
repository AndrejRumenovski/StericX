//! Reusable end-to-end model training.
//!
//! [`train_scientific_model`] is the single entry point callers need: it takes a
//! packed record matrix plus row labels, fits the mechanistically constrained
//! model on the training partition, and returns both the diagnostic report and
//! the frozen predictions for every non-training row. All scientific
//! methodology lives in [`fit_scientific_model_grouped`]; this module only
//! orchestrates it so the work is not duplicated by each caller.

use super::{
    FitOptions, RegressXPredictor, ScientificFitReport, expand_features,
    fit_scientific_model_grouped,
};
use crate::model::dataset::{ReactionLabel, TrainingSplit};
use crate::storage::PackedReactionRecord;
use serde::{Deserialize, Serialize};

/// One prediction recorded before its experimental target is revealed.
///
/// Field names match the frozen-prediction CSV schema so the row round-trips
/// through `stericx fit` and `stericx evaluate` unchanged.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
pub struct FrozenPrediction {
    #[serde(rename = "Reaction_ID")]
    pub reaction_id: String,
    #[serde(rename = "Ligand_Group")]
    pub ligand_group: String,
    #[serde(rename = "Dataset_Split")]
    pub dataset_split: String,
    #[serde(rename = "Predicted_ddG_kcal_mol")]
    pub predicted_ddg: f32,
    #[serde(rename = "Applicability_Domain")]
    pub applicability_domain: String,
}

/// A fitted model together with its frozen non-training predictions.
#[derive(Clone, Debug)]
pub struct TrainedModel {
    /// Model artifact, diagnostics, and validation statistics.
    pub report: ScientificFitReport,
    /// Predictions for every non-training row, in ascending record order.
    pub frozen_predictions: Vec<FrozenPrediction>,
}

/// Fits the scientific model and freezes predictions for all held-out rows.
///
/// `labels` must be positionally aligned with `records`. Training rows are
/// selected by [`TrainingSplit`], descriptor scaling and selection use those
/// rows only, and every remaining row receives a prediction plus an
/// applicability-domain verdict.
pub fn train_scientific_model(
    records: &[PackedReactionRecord],
    labels: &[ReactionLabel],
    options: FitOptions,
) -> Result<TrainedModel, String> {
    if labels.len() != records.len() {
        return Err(format!(
            "metadata has {} rows but sigpack contains {} records",
            labels.len(),
            records.len()
        ));
    }
    let split = TrainingSplit::from_labels(labels)?;
    let mut report = fit_scientific_model_grouped(
        records,
        split.training_indices(),
        split.training_groups(),
        options,
    )?;
    // Name the training observations. `standardized_training_points` is built
    // from the training rows in `training_indices` order, so the identifiers
    // are positionally aligned by construction.
    if let Some(geometry) = report.training_geometry.as_mut() {
        let identifiers = split
            .training_indices()
            .iter()
            .map(|&index| labels[index].reaction_id.clone())
            .collect::<Vec<_>>();
        if identifiers.len() == geometry.standardized_training_points.len() {
            geometry.training_labels = identifiers;
        }
    }

    let predictor = RegressXPredictor::new(report.weights);
    let frozen_records = split
        .frozen_indices()
        .iter()
        .map(|&index| records[index])
        .collect::<Vec<_>>();
    let predicted = predictor.predict_batch(&frozen_records);
    let frozen_predictions = split
        .frozen_indices()
        .iter()
        .zip(predicted)
        .map(|(&index, predicted_ddg)| FrozenPrediction {
            reaction_id: labels[index].reaction_id.clone(),
            ligand_group: labels[index].ligand_group.clone(),
            dataset_split: labels[index].dataset_split.clone(),
            predicted_ddg,
            applicability_domain: applicability_status(&records[index], &report),
        })
        .collect();

    Ok(TrainedModel {
        report,
        frozen_predictions,
    })
}

/// Reports whether a record lies inside the fitted applicability domain.
///
/// Returns `inside_training_range`, or `outside:` followed by the
/// pipe-separated selected descriptors whose values leave their training range.
#[must_use]
pub fn applicability_status(record: &PackedReactionRecord, report: &ScientificFitReport) -> String {
    let features = expand_features(record);
    let outside = report
        .selected_feature_indices
        .iter()
        .zip(&report.applicability_domain)
        .filter(|(column, domain)| {
            let value = f64::from(features[**column]);
            value < domain.minimum || value > domain.maximum
        })
        .map(|(_, domain)| domain.feature.as_str())
        .collect::<Vec<_>>();
    if outside.is_empty() {
        "inside_training_range".into()
    } else {
        format!("outside:{}", outside.join("|"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn records() -> Vec<PackedReactionRecord> {
        (0..12)
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
            .collect()
    }

    fn labels(count: usize) -> Vec<ReactionLabel> {
        (0..count)
            .map(|index| {
                let split = if index < count - 2 { "train" } else { "blind" };
                ReactionLabel::new(format!("R{index}"), split, format!("group_{}", index % 3))
            })
            .collect()
    }

    fn options() -> FitOptions {
        FitOptions {
            bootstrap_samples: 25,
            permutation_samples: 25,
            ..FitOptions::default()
        }
    }

    #[test]
    fn trains_and_freezes_every_non_training_row() {
        let records = records();
        let labels = labels(records.len());

        let trained = train_scientific_model(&records, &labels, options()).unwrap();

        assert_eq!(trained.report.training_count, 10);
        assert_eq!(trained.frozen_predictions.len(), 2);
        assert_eq!(trained.frozen_predictions[0].reaction_id, "R10");
        assert_eq!(trained.frozen_predictions[1].dataset_split, "blind");
        assert!(
            trained
                .frozen_predictions
                .iter()
                .all(|row| row.predicted_ddg.is_finite())
        );
    }

    #[test]
    fn frozen_predictions_match_the_reported_weights() {
        let records = records();
        let labels = labels(records.len());

        let trained = train_scientific_model(&records, &labels, options()).unwrap();

        let predictor = RegressXPredictor::new(trained.report.weights);
        for (row, &index) in trained.frozen_predictions.iter().zip([10_usize, 11].iter()) {
            assert_eq!(
                row.predicted_ddg.to_bits(),
                predictor.predict(&records[index]).to_bits()
            );
        }
    }

    #[test]
    fn rejects_misaligned_labels() {
        let records = records();
        let labels = labels(records.len() - 1);

        assert_eq!(
            train_scientific_model(&records, &labels, options()).unwrap_err(),
            "metadata has 11 rows but sigpack contains 12 records"
        );
    }

    #[test]
    fn flags_records_outside_the_training_descriptor_range() {
        let records = records();
        let labels = labels(records.len());
        let trained = train_scientific_model(&records, &labels, options()).unwrap();
        let selected = trained.report.selected_feature_indices[0];
        let domain = &trained.report.applicability_domain[0];

        let inside = records[0];
        assert_eq!(
            applicability_status(&inside, &trained.report),
            "inside_training_range"
        );

        let mut outside = records[0];
        let far = (domain.maximum + 10.0) as f32;
        match selected {
            1 => outside.l = far,
            2 => outside.b1 = far,
            3 => outside.b5 = far,
            _ => panic!("unexpected selected descriptor {selected}"),
        }
        assert!(applicability_status(&outside, &trained.report).starts_with("outside:"));
    }
}
