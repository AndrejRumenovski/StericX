//! `stericx fit`: fit and freeze an interpretable physical-organic model.

use crate::commands::FrozenPredictionRow;
use crate::output::{
    atomic_write_csv_rows, atomic_write_json, display_optional_metric, millis,
    print_memory_metrics, resident_memory_bytes,
};
use crate::reaction::load_reaction_metadata;
use std::error::Error;
use std::path::Path;
use std::time::Instant;
use steric_x::model::expand_features;
use steric_x::{
    FitOptions, PackedReactionRecord, RegressXPredictor, ScientificFitReport, SigPackReader,
    fit_scientific_model_grouped,
};

pub(crate) fn fit_command(
    data: &Path,
    metadata_path: &Path,
    output: &Path,
    predictions_path: &Path,
    options: FitOptions,
) -> Result<(), Box<dyn Error>> {
    let total_started = Instant::now();
    let rss_start = resident_memory_bytes();
    let reader = SigPackReader::open(data)?;
    let records = reader.records();
    if records.is_empty() {
        return Err("sigpack matrix contains no records".into());
    }
    let metadata = load_reaction_metadata(metadata_path)?;
    if metadata.len() != records.len() {
        return Err(format!(
            "metadata has {} rows but sigpack contains {} records",
            metadata.len(),
            records.len()
        )
        .into());
    }
    let training_indices = metadata
        .iter()
        .enumerate()
        .filter(|(_, row)| row.dataset_split.eq_ignore_ascii_case("train"))
        .map(|(index, _)| index)
        .collect::<Vec<_>>();
    let frozen_indices = metadata
        .iter()
        .enumerate()
        .filter(|(_, row)| !row.dataset_split.eq_ignore_ascii_case("train"))
        .map(|(index, _)| index)
        .collect::<Vec<_>>();
    if frozen_indices.is_empty() {
        return Err("metadata contains no non-training rows to freeze".into());
    }
    let training_groups = training_indices
        .iter()
        .map(|&index| {
            let group = metadata[index].ligand_group.trim();
            if group.is_empty() {
                metadata[index].reaction_id.clone()
            } else {
                group.to_owned()
            }
        })
        .collect::<Vec<_>>();

    let fit_started = Instant::now();
    let report =
        fit_scientific_model_grouped(records, &training_indices, &training_groups, options)?;
    let fit_time = fit_started.elapsed();
    atomic_write_json(&report, output)?;

    let predictor = RegressXPredictor::new(report.weights);
    let frozen_records = frozen_indices
        .iter()
        .map(|&index| records[index])
        .collect::<Vec<_>>();
    let predicted = predictor.predict_batch(&frozen_records);
    let frozen_rows = frozen_indices
        .iter()
        .zip(predicted)
        .map(|(&index, predicted_ddg)| FrozenPredictionRow {
            reaction_id: metadata[index].reaction_id.clone(),
            ligand_group: metadata[index].ligand_group.clone(),
            dataset_split: metadata[index].dataset_split.clone(),
            predicted_ddg,
            applicability_domain: applicability_status(&records[index], &report),
        })
        .collect::<Vec<_>>();
    atomic_write_csv_rows(&frozen_rows, predictions_path)?;

    let total_time = total_started.elapsed();
    println!("command=fit");
    println!("records_total={}", records.len());
    println!("training_records={}", training_indices.len());
    println!("training_groups={}", report.training_group_count);
    println!("frozen_prediction_records={}", frozen_rows.len());
    println!("selected_features={}", report.selected_features.join(","));
    println!(
        "training_r2={}",
        display_optional_metric(report.training.r2)
    );
    println!(
        "fixed_feature_loo_r2={}",
        display_optional_metric(report.fixed_feature_loo.r2)
    );
    println!(
        "fixed_feature_group_loo_r2={}",
        display_optional_metric(report.fixed_feature_group_loo.r2)
    );
    println!("training_rmse={:.8}", report.training.rmse);
    println!(
        "fixed_feature_loo_rmse={:.8}",
        report.fixed_feature_loo.rmse
    );
    println!(
        "fixed_feature_group_loo_rmse={:.8}",
        report.fixed_feature_group_loo.rmse
    );
    println!(
        "response_permutation_p_value={:.8}",
        report.response_permutation_p_value
    );
    println!("model_output={}", output.display());
    println!("frozen_predictions_output={}", predictions_path.display());
    println!("fit_ms={:.3}", millis(fit_time));
    println!("total_ms={:.3}", millis(total_time));
    print_memory_metrics(rss_start, resident_memory_bytes());
    Ok(())
}

fn applicability_status(record: &PackedReactionRecord, report: &ScientificFitReport) -> String {
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
