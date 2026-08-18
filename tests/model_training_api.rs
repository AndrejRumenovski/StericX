//! Regression tests pinning the reusable training API to the numerical
//! behavior of the pre-refactor implementation.
//!
//! `tests/data/golden_training_report.json` was produced by the previous
//! CLI-coupled code path (`fit_scientific_model_grouped` called directly from
//! `stericx fit`) on the fixture built below. Every field of that artifact is
//! re-checked here through the new [`train_scientific_model`] entry point, so
//! any change to descriptor selection, coefficients, validation statistics, or
//! diagnostics fails the suite instead of silently rewriting a study result.

use steric_x::model::expand_features;
use steric_x::model::{ModelMetrics, ScientificFitReport};
use steric_x::{
    FitOptions, PackedReactionRecord, ReactionLabel, RegressXPredictor, TrainingSplit,
    fit_scientific_model_grouped, score_frozen_predictions, train_scientific_model,
};

/// Absolute/relative slack allowed on f64 diagnostics. The refactor is pure
/// code motion, so agreement is expected to be far tighter than this.
const TOLERANCE: f64 = 1.0e-12;

/// Deterministic 24-row fixture: 20 training rows across 6 scaffold groups and
/// 4 held-out rows. The response depends on `L` and on the `B5 * NBO`
/// interaction, so BIC forward selection has a genuine multi-term choice to
/// make and the correlated-pair and term-count guards are both exercised.
fn fixture_records() -> Vec<PackedReactionRecord> {
    (0..24_usize)
        .map(|index| {
            let l = 1.0 + 0.13 * index as f32;
            let b1 = 2.0 + 0.07 * ((index * 7) % 5) as f32;
            let b5 = 3.5 + 0.11 * ((index * 3) % 7) as f32;
            let nbo_charge = -0.6 + 0.02 * ((index * 5) % 9) as f32;
            let ir_freq = 1_600.0 + ((index * 11) % 13) as f32;
            let noise = ((index * 37) % 11) as f32 * 0.01 - 0.05;
            PackedReactionRecord {
                l,
                b1,
                b5,
                nbo_charge,
                ir_freq,
                temp_k: 298.15,
                exp_ddg: 0.4 + 0.8 * l - 0.35 * b5 * nbo_charge + noise,
                ..PackedReactionRecord::default()
            }
        })
        .collect()
}

fn fixture_labels() -> Vec<ReactionLabel> {
    (0..24_usize)
        .map(|index| {
            let split = if index < 20 { "train" } else { "blind" };
            ReactionLabel::new(
                format!("R{index:02}"),
                split,
                format!("group_{}", index % 6),
            )
        })
        .collect()
}

fn fixture_options() -> FitOptions {
    FitOptions {
        max_terms: 3,
        bootstrap_samples: 250,
        permutation_samples: 200,
        seed: 20_260_725,
    }
}

fn golden_report() -> ScientificFitReport {
    let json = include_str!("data/golden_training_report.json");
    serde_json::from_str(json).expect("golden report is valid ScientificFitReport JSON")
}

#[track_caller]
fn assert_close(actual: f64, expected: f64, label: &str) {
    let slack = TOLERANCE * expected.abs().max(1.0);
    assert!(
        (actual - expected).abs() <= slack,
        "{label}: got {actual:.17e}, expected {expected:.17e}"
    );
}

#[track_caller]
fn assert_optional_close(actual: Option<f64>, expected: Option<f64>, label: &str) {
    match (actual, expected) {
        (Some(actual), Some(expected)) => assert_close(actual, expected, label),
        (None, None) => {}
        (actual, expected) => panic!("{label}: got {actual:?}, expected {expected:?}"),
    }
}

#[track_caller]
fn assert_metrics(actual: &ModelMetrics, expected: &ModelMetrics, label: &str) {
    assert_eq!(actual.count, expected.count, "{label}.count");
    assert_optional_close(actual.r2, expected.r2, &format!("{label}.r2"));
    assert_close(actual.mae, expected.mae, &format!("{label}.mae"));
    assert_close(actual.rmse, expected.rmse, &format!("{label}.rmse"));
}

#[test]
fn training_api_reproduces_the_pre_refactor_report() {
    let expected = golden_report();
    let trained = train_scientific_model(&fixture_records(), &fixture_labels(), fixture_options())
        .expect("fixture trains");
    let actual = &trained.report;

    // Artifact identity and descriptor selection.
    assert_eq!(actual.schema_version, expected.schema_version);
    assert_eq!(actual.model, expected.model);
    assert_eq!(actual.training_count, expected.training_count);
    assert_eq!(actual.training_group_count, expected.training_group_count);
    assert_eq!(actual.feature_names, expected.feature_names);
    assert_eq!(
        actual.selected_feature_indices,
        expected.selected_feature_indices
    );
    assert_eq!(actual.selected_features, expected.selected_features);
    assert_eq!(actual.notes, expected.notes);

    // Coefficients are f32 and must match exactly, not merely closely.
    assert_eq!(actual.weights, expected.weights, "raw-scale weights");
    for column in 0..actual.standardized_means.len() {
        assert_close(
            actual.standardized_means[column],
            expected.standardized_means[column],
            &format!("standardized_means[{column}]"),
        );
        assert_close(
            actual.standardized_scales[column],
            expected.standardized_scales[column],
            &format!("standardized_scales[{column}]"),
        );
    }

    // Training fit, leave-one-out, and leave-one-scaffold-group-out.
    assert_metrics(&actual.training, &expected.training, "training");
    assert_metrics(
        &actual.fixed_feature_loo,
        &expected.fixed_feature_loo,
        "fixed_feature_loo",
    );
    assert_metrics(
        &actual.fixed_feature_group_loo,
        &expected.fixed_feature_group_loo,
        "fixed_feature_group_loo",
    );

    // Nested ridge and LASSO baselines, including the tuned penalties.
    for (actual, expected, label) in [
        (&actual.ridge_baseline, &expected.ridge_baseline, "ridge"),
        (&actual.lasso_baseline, &expected.lasso_baseline, "lasso"),
    ] {
        assert_eq!(actual.model, expected.model, "{label}.model");
        assert_close(
            actual.regularization,
            expected.regularization,
            &format!("{label}.regularization"),
        );
        assert_eq!(actual.weights, expected.weights, "{label}.weights");
        assert_metrics(
            &actual.training,
            &expected.training,
            &format!("{label}.training"),
        );
        assert_metrics(
            &actual.nested_loo,
            &expected.nested_loo,
            &format!("{label}.nested_loo"),
        );
    }

    // Bootstrap intervals, Y-scrambling, VIF, correlations, and the domain.
    assert_eq!(
        actual.coefficient_intervals.len(),
        expected.coefficient_intervals.len()
    );
    for (actual, expected) in actual
        .coefficient_intervals
        .iter()
        .zip(&expected.coefficient_intervals)
    {
        assert_eq!(actual.feature, expected.feature);
        assert_close(
            actual.estimate,
            expected.estimate,
            &format!("{}.estimate", actual.feature),
        );
        assert_close(
            actual.lower_95,
            expected.lower_95,
            &format!("{}.lower_95", actual.feature),
        );
        assert_close(
            actual.upper_95,
            expected.upper_95,
            &format!("{}.upper_95", actual.feature),
        );
    }
    assert_close(
        actual.response_permutation_p_value,
        expected.response_permutation_p_value,
        "response_permutation_p_value",
    );
    for (row, (actual_row, expected_row)) in actual
        .correlation_matrix
        .iter()
        .zip(&expected.correlation_matrix)
        .enumerate()
    {
        assert_eq!(actual_row.len(), expected_row.len());
        for (column, (actual, expected)) in actual_row.iter().zip(expected_row).enumerate() {
            assert_close(*actual, *expected, &format!("correlation[{row}][{column}]"));
        }
    }
    for (column, (actual, expected)) in actual
        .variance_inflation_factors
        .iter()
        .zip(&expected.variance_inflation_factors)
        .enumerate()
    {
        assert_optional_close(*actual, *expected, &format!("vif[{column}]"));
    }
    assert_eq!(
        actual.applicability_domain.len(),
        expected.applicability_domain.len()
    );
    for (actual, expected) in actual
        .applicability_domain
        .iter()
        .zip(&expected.applicability_domain)
    {
        assert_eq!(actual.feature, expected.feature);
        assert_close(
            actual.minimum,
            expected.minimum,
            &format!("{}.minimum", actual.feature),
        );
        assert_close(
            actual.maximum,
            expected.maximum,
            &format!("{}.maximum", actual.feature),
        );
    }
}

#[test]
fn training_api_matches_the_previous_cli_composition_bit_for_bit() {
    let records = fixture_records();
    let labels = fixture_labels();
    let options = fixture_options();

    // Exactly what `fit_command` used to do inline, reproduced here.
    let training_indices = labels
        .iter()
        .enumerate()
        .filter(|(_, label)| label.dataset_split.eq_ignore_ascii_case("train"))
        .map(|(index, _)| index)
        .collect::<Vec<_>>();
    let frozen_indices = labels
        .iter()
        .enumerate()
        .filter(|(_, label)| !label.dataset_split.eq_ignore_ascii_case("train"))
        .map(|(index, _)| index)
        .collect::<Vec<_>>();
    let training_groups = training_indices
        .iter()
        .map(|&index| {
            let group = labels[index].ligand_group.trim();
            if group.is_empty() {
                labels[index].reaction_id.clone()
            } else {
                group.to_owned()
            }
        })
        .collect::<Vec<_>>();
    let legacy_report =
        fit_scientific_model_grouped(&records, &training_indices, &training_groups, options)
            .expect("legacy path fits");
    let legacy_predictor = RegressXPredictor::new(legacy_report.weights);
    let legacy_predictions = legacy_predictor.predict_batch(
        &frozen_indices
            .iter()
            .map(|&index| records[index])
            .collect::<Vec<_>>(),
    );

    let trained = train_scientific_model(&records, &labels, options).expect("new path fits");

    assert_eq!(trained.report.weights, legacy_report.weights);
    assert_eq!(
        trained.report.selected_feature_indices,
        legacy_report.selected_feature_indices
    );
    assert_eq!(
        serde_json::to_string(&trained.report).unwrap(),
        serde_json::to_string(&legacy_report).unwrap(),
        "serialized model artifacts must be identical"
    );
    assert_eq!(trained.frozen_predictions.len(), legacy_predictions.len());
    for (row, expected) in trained.frozen_predictions.iter().zip(&legacy_predictions) {
        assert_eq!(
            row.predicted_ddg.to_bits(),
            expected.to_bits(),
            "frozen prediction for {} changed",
            row.reaction_id
        );
    }
}

#[test]
fn training_split_reproduces_the_previous_inline_partition() {
    let labels = fixture_labels();

    let split = TrainingSplit::from_labels(&labels).expect("fixture partitions");

    assert_eq!(split.training_indices(), (0..20).collect::<Vec<_>>());
    assert_eq!(split.frozen_indices(), [20, 21, 22, 23]);
    assert_eq!(split.training_group_count(), 6);
    for (position, &index) in split.training_indices().iter().enumerate() {
        assert_eq!(
            split.training_groups()[position],
            format!("group_{}", index % 6)
        );
    }
}

#[test]
fn frozen_predictions_keep_the_study_001_csv_schema() {
    let trained = train_scientific_model(&fixture_records(), &fixture_labels(), fixture_options())
        .expect("fixture trains");

    let mut writer = csv::Writer::from_writer(Vec::new());
    for row in &trained.frozen_predictions {
        writer.serialize(row).unwrap();
    }
    let csv = String::from_utf8(writer.into_inner().unwrap()).unwrap();
    let header = csv.lines().next().unwrap();

    assert_eq!(
        header,
        "Reaction_ID,Ligand_Group,Dataset_Split,Predicted_ddG_kcal_mol,Applicability_Domain"
    );
    assert_eq!(csv.lines().count(), 1 + trained.frozen_predictions.len());
}

#[test]
fn frozen_predictions_round_trip_into_evaluation() {
    let records = fixture_records();
    let labels = fixture_labels();
    let trained =
        train_scientific_model(&records, &labels, fixture_options()).expect("fixture trains");

    let summary = score_frozen_predictions(
        &records,
        &labels,
        &trained.report,
        &trained.frozen_predictions,
    )
    .expect("frozen predictions score");

    assert_eq!(summary.evaluated_records, trained.frozen_predictions.len());
    for (scored, frozen) in summary
        .scored_predictions
        .iter()
        .zip(&trained.frozen_predictions)
    {
        assert_eq!(scored.reaction_id, frozen.reaction_id);
        assert_eq!(
            scored.predicted_ddg_kcal_mol.to_bits(),
            frozen.predicted_ddg.to_bits()
        );
        assert_eq!(scored.applicability_domain, frozen.applicability_domain);
    }
}

/// Second fixture, tuned so the correlated-pair gate is the deciding factor.
///
/// `B5` is `B1` plus a small independent term, giving `|r| ≈ 0.971` — above the
/// 0.95 rejection threshold but below 1.0. The response depends on both, so
/// admitting the second descriptor would improve BIC by roughly 96 points. Only
/// the correlation gate keeps it out, which is what makes this fixture able to
/// detect a change to that threshold; the first fixture never reaches it.
fn collinear_records() -> Vec<PackedReactionRecord> {
    (0..18_usize)
        .map(|index| {
            let b1 = 1.5 + 0.05 * index as f32;
            let delta = 0.046 * (((index * 7) % 5) as f32 - 2.0);
            PackedReactionRecord {
                l: 4.0,
                b1,
                b5: b1 + delta,
                nbo_charge: -0.4,
                ir_freq: 1_650.0,
                temp_k: 298.15,
                exp_ddg: 0.3 + 0.9 * b1 + 3.0 * delta + 0.004 * (((index * 13) % 7) as f32 - 3.0),
                ..PackedReactionRecord::default()
            }
        })
        .collect()
}

#[test]
fn forward_selection_still_rejects_strongly_correlated_descriptors() {
    let records = collinear_records();
    let labels = (0..records.len())
        .map(|index| {
            let split = if index < records.len() - 3 {
                "train"
            } else {
                "blind"
            };
            ReactionLabel::new(
                format!("C{index:02}"),
                split,
                format!("group_{}", index % 4),
            )
        })
        .collect::<Vec<_>>();

    let trained = train_scientific_model(&records, &labels, fixture_options()).expect("fits");
    let selected = &trained.report.selected_feature_indices;

    assert_eq!(
        selected.len(),
        1,
        "expected the collinear family to contribute one term, got {:?}",
        trained.report.selected_features
    );
    // The guard is defined on the training rows, so check it there.
    let training_rows = labels
        .iter()
        .enumerate()
        .filter(|(_, label)| label.is_training())
        .map(|(index, _)| expand_features(&records[index]).map(f64::from))
        .collect::<Vec<_>>();
    for (position, &left) in selected.iter().enumerate() {
        for &right in selected.iter().skip(position + 1) {
            let r = pearson(&training_rows, left, right);
            assert!(
                r.abs() <= 0.95,
                "selected descriptors {left} and {right} correlate at r = {r}"
            );
        }
    }
    // Every rejected member of the family really was correlated past the gate.
    for column in [2_usize, 3, 5, 6] {
        if column == selected[0] {
            continue;
        }
        let r = pearson(&training_rows, selected[0], column);
        assert!(
            r.abs() > 0.95,
            "column {column} should have been rejected as collinear, r = {r}"
        );
    }
}

fn pearson(rows: &[[f64; 8]], left: usize, right: usize) -> f64 {
    let count = rows.len() as f64;
    let left_mean = rows.iter().map(|row| row[left]).sum::<f64>() / count;
    let right_mean = rows.iter().map(|row| row[right]).sum::<f64>() / count;
    let covariance = rows
        .iter()
        .map(|row| (row[left] - left_mean) * (row[right] - right_mean))
        .sum::<f64>();
    let left_norm = rows
        .iter()
        .map(|row| (row[left] - left_mean).powi(2))
        .sum::<f64>()
        .sqrt();
    let right_norm = rows
        .iter()
        .map(|row| (row[right] - right_mean).powi(2))
        .sum::<f64>()
        .sqrt();
    covariance / (left_norm * right_norm)
}
