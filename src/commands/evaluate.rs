//! `stericx evaluate`: reveal and score previously frozen non-training predictions.

use crate::commands::FrozenPredictionRow;
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
use steric_x::{RegressXPredictor, ScientificFitReport, SigPackReader};

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

#[derive(Debug, Serialize)]
struct ScoredPrediction {
    reaction_id: String,
    ligand_group: String,
    dataset_split: String,
    predicted_ddg_kcal_mol: f32,
    experimental_ddg_kcal_mol: f32,
    residual_kcal_mol: f32,
    applicability_domain: String,
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
    let metadata = load_reaction_metadata(metadata_path)?;
    if metadata.len() != records.len() {
        return Err(format!(
            "metadata has {} rows but sigpack contains {} records",
            metadata.len(),
            records.len()
        )
        .into());
    }
    let model_contents = fs::read_to_string(model_path)?;
    let model: ScientificFitReport = serde_json::from_str(&model_contents)?;
    let frozen_bytes = fs::read(predictions_path)?;
    let mut frozen_reader = csv::ReaderBuilder::new()
        .trim(csv::Trim::All)
        .from_reader(frozen_bytes.as_slice());
    let frozen = frozen_reader
        .deserialize::<FrozenPredictionRow>()
        .collect::<Result<Vec<_>, _>>()?;
    if frozen.is_empty() {
        return Err("frozen prediction file contains no rows".into());
    }

    let predictor = RegressXPredictor::new(model.weights);
    let mut scored = Vec::with_capacity(frozen.len());
    let mut actual = Vec::with_capacity(frozen.len());
    let mut predicted = Vec::with_capacity(frozen.len());
    for frozen_row in &frozen {
        let index = metadata
            .iter()
            .position(|row| row.reaction_id == frozen_row.reaction_id)
            .ok_or_else(|| {
                format!(
                    "frozen reaction {} is absent from metadata",
                    frozen_row.reaction_id
                )
            })?;
        if metadata[index].dataset_split.eq_ignore_ascii_case("train") {
            return Err(format!(
                "frozen reaction {} is marked as training data",
                frozen_row.reaction_id
            )
            .into());
        }
        let recomputed = predictor.predict(&records[index]);
        if (recomputed - frozen_row.predicted_ddg).abs() > 1.0e-4 {
            return Err(format!(
                "frozen prediction for {} does not match the supplied model",
                frozen_row.reaction_id
            )
            .into());
        }
        let experimental = records[index].exp_ddg;
        if !experimental.is_finite() {
            return Err(format!(
                "reaction {} has no finite revealed target",
                frozen_row.reaction_id
            )
            .into());
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
    let evaluation = EvaluationReport {
        schema_version: 1,
        model_path: model_path.display().to_string(),
        frozen_predictions_path: predictions_path.display().to_string(),
        frozen_predictions_fnv1a64: format!("{:016x}", fnv1a64(&frozen_bytes)),
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
