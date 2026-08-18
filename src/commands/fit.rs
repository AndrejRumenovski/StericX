//! `stericx fit`: fit and freeze an interpretable physical-organic model.
//!
//! Training itself lives in [`steric_x::model::training`]. This module only
//! reads inputs, calls that API, and reports.

use crate::output::{
    atomic_write_csv_rows, atomic_write_json, display_optional_metric, fnv1a64, millis,
    print_memory_metrics, resident_memory_bytes, write_atomic_text,
};
use crate::reaction::load_reaction_metadata;
use std::error::Error;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;
use steric_x::model::{
    CreationMetadata, DatasetDigest, ModelProvenance, Optimization, PortableModel,
    ReactionProvenance, ResponseSpec, TrainingProvenance,
};
use steric_x::{FitOptions, ScientificFitReport, SigPackReader, train_scientific_model};

/// Optional portable-model output requested on the command line.
///
/// Chemistry context arrives here exactly as the operator typed it. Anything
/// left unset stays `None` in the document rather than being guessed.
pub(crate) struct PortableModelRequest {
    pub(crate) path: Option<PathBuf>,
    pub(crate) model_id: Option<String>,
    pub(crate) reaction: ReactionProvenance,
    pub(crate) response_temp_k: Option<f32>,
    pub(crate) optimization: Optimization,
}

pub(crate) fn fit_command(
    data: &Path,
    metadata_path: &Path,
    output: &Path,
    predictions_path: &Path,
    options: FitOptions,
    portable: PortableModelRequest,
) -> Result<(), Box<dyn Error>> {
    let total_started = Instant::now();
    let rss_start = resident_memory_bytes();
    let reader = SigPackReader::open(data)?;
    let records = reader.records();
    if records.is_empty() {
        return Err("sigpack matrix contains no records".into());
    }
    let labels = load_reaction_metadata(metadata_path)?;

    let fit_started = Instant::now();
    let trained = train_scientific_model(records, &labels, options)?;
    let fit_time = fit_started.elapsed();
    let report = &trained.report;

    atomic_write_json(report, output)?;
    atomic_write_csv_rows(&trained.frozen_predictions, predictions_path)?;

    let portable_output = match portable.path.as_deref() {
        Some(path) => {
            let document =
                build_portable_model(report.clone(), data, metadata_path, options, &portable)?;
            write_atomic_text(path, &document.to_json()?)?;
            Some((path.to_owned(), document))
        }
        None => None,
    };

    let total_time = total_started.elapsed();
    println!("command=fit");
    println!("records_total={}", records.len());
    println!("training_records={}", report.training_count);
    println!("training_groups={}", report.training_group_count);
    println!(
        "frozen_prediction_records={}",
        trained.frozen_predictions.len()
    );
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
    if let Some((path, document)) = &portable_output {
        println!("portable_model_output={}", path.display());
        println!(
            "portable_model_schema_version={}",
            document.schema_version()
        );
        let missing = document.missing_provenance();
        println!(
            "portable_model_missing_provenance={}",
            if missing.is_empty() {
                "none".to_owned()
            } else {
                missing.join(",")
            }
        );
    }
    println!("fit_ms={:.3}", millis(fit_time));
    println!("total_ms={:.3}", millis(total_time));
    print_memory_metrics(rss_start, resident_memory_bytes());
    Ok(())
}

/// Assembles a portable model from the fit plus command-line provenance.
fn build_portable_model(
    report: ScientificFitReport,
    data: &Path,
    metadata_path: &Path,
    options: FitOptions,
    request: &PortableModelRequest,
) -> Result<PortableModel, Box<dyn Error>> {
    let dataset_digests = vec![file_digest(data)?, file_digest(metadata_path)?];
    let model_id = match request.model_id.as_deref() {
        Some(id) if !id.trim().is_empty() => id.trim().to_owned(),
        _ => derive_model_id(&report, &dataset_digests),
    };
    let provenance = ModelProvenance {
        model_id,
        stericx_version: env!("CARGO_PKG_VERSION").to_owned(),
        record_format: "sigpack_v1".to_owned(),
        training: TrainingProvenance {
            record_count: report.training_count,
            group_count: report.training_group_count,
            dataset_digests,
            fit_options: options,
        },
        reaction: request.reaction.clone(),
    };
    Ok(PortableModel::from_fit_report(
        report,
        ResponseSpec::transition_state_energy_difference(
            request.response_temp_k,
            request.optimization,
        ),
        provenance,
        CreationMetadata::now("stericx fit"),
    )?)
}

/// Digests one training input.
///
/// FNV-1a is not a cryptographic hash. The algorithm is recorded alongside the
/// value so a consumer can see exactly what the digest does and does not prove.
fn file_digest(path: &Path) -> Result<DatasetDigest, Box<dyn Error>> {
    let bytes = fs::read(path)?;
    Ok(DatasetDigest {
        artifact: path.file_name().map_or_else(
            || path.display().to_string(),
            |name| name.to_string_lossy().into_owned(),
        ),
        algorithm: "fnv1a64".to_owned(),
        digest: format!("{:016x}", fnv1a64(&bytes)),
        byte_count: bytes.len() as u64,
    })
}

/// Derives a deterministic identifier from the model and its training inputs.
fn derive_model_id(report: &ScientificFitReport, digests: &[DatasetDigest]) -> String {
    let mut seed = format!("{}|{}", report.model, report.selected_features.join(","));
    for digest in digests {
        seed.push_str(&format!("|{}:{}", digest.artifact, digest.digest));
    }
    format!("stericx-{:016x}", fnv1a64(seed.as_bytes()))
}
