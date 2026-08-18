//! `stericx evaluate`: reveal and score previously frozen non-training predictions.
//!
//! Scoring lives in [`steric_x::model::evaluation`]. This module reads the
//! frozen artifacts, calls it, and reports.

use crate::output::{
    atomic_write_json, display_optional_metric, fnv1a64, millis, print_memory_metrics,
    resident_memory_bytes,
};
use crate::reaction::load_reaction_metadata;
use serde::Serialize;
use std::error::Error;
use std::fs;
use std::path::Path;
use std::time::Instant;
use steric_x::{
    FrozenPrediction, ScientificFitReport, ScoredPrediction, SigPackReader,
    score_frozen_predictions,
};

#[derive(Debug, Serialize)]
struct EvaluationReport {
    schema_version: u32,
    model_path: String,
    frozen_predictions_path: String,
    frozen_predictions_fnv1a64: String,
    evaluated_records: usize,
    mae_kcal_mol: f64,
    rmse_kcal_mol: f64,
    r2: Option<f64>,
    applicability_warnings: usize,
    scored_predictions: Vec<ScoredPrediction>,
}

pub(crate) fn evaluate_command(
    data: &Path,
    metadata_path: &Path,
    model_path: &Path,
    predictions_path: &Path,
    output: &Path,
) -> Result<(), Box<dyn Error>> {
    let total_started = Instant::now();
    let rss_start = resident_memory_bytes();
    let reader = SigPackReader::open(data)?;
    let records = reader.records();
    let labels = load_reaction_metadata(metadata_path)?;
    let model_contents = fs::read_to_string(model_path)?;
    let model: ScientificFitReport = serde_json::from_str(&model_contents)?;
    let frozen_bytes = fs::read(predictions_path)?;
    let mut frozen_reader = csv::ReaderBuilder::new()
        .trim(csv::Trim::All)
        .from_reader(frozen_bytes.as_slice());
    let frozen = frozen_reader
        .deserialize::<FrozenPrediction>()
        .collect::<Result<Vec<_>, _>>()?;

    let summary = score_frozen_predictions(records, &labels, &model, &frozen)?;

    let evaluation = EvaluationReport {
        schema_version: 1,
        model_path: model_path.display().to_string(),
        frozen_predictions_path: predictions_path.display().to_string(),
        frozen_predictions_fnv1a64: format!("{:016x}", fnv1a64(&frozen_bytes)),
        evaluated_records: summary.evaluated_records,
        mae_kcal_mol: summary.mae_kcal_mol,
        rmse_kcal_mol: summary.rmse_kcal_mol,
        r2: summary.r2,
        applicability_warnings: summary.applicability_warnings,
        scored_predictions: summary.scored_predictions,
    };
    atomic_write_json(&evaluation, output)?;

    println!("command=evaluate");
    println!("evaluated_records={}", evaluation.evaluated_records);
    println!("mae_kcal_mol={:.8}", evaluation.mae_kcal_mol);
    println!("rmse_kcal_mol={:.8}", evaluation.rmse_kcal_mol);
    println!("r2={}", display_optional_metric(evaluation.r2));
    println!(
        "applicability_warnings={}",
        evaluation.applicability_warnings
    );
    println!(
        "frozen_predictions_fnv1a64={}",
        evaluation.frozen_predictions_fnv1a64
    );
    println!("evaluation_output={}", output.display());
    println!("total_ms={:.3}", millis(total_started.elapsed()));
    print_memory_metrics(rss_start, resident_memory_bytes());
    Ok(())
}
