//! `stericx model`: inspect and validate saved model documents.
//!
//! Both subcommands read the document through [`steric_x::model::diagnose`], so
//! a file that cannot be parsed is still reported field by field rather than as
//! a decoder message about the whole file.

use crate::cli::ReportFormat;
use serde::Serialize;
use std::error::Error;
use std::fs;
use std::path::Path;
use steric_x::model::{Diagnosis, ModelIssue, ModelSummary, Severity, diagnose};

/// Machine-readable result of `stericx model validate`.
#[derive(Debug, Serialize)]
struct ValidationReport<'a> {
    model_path: String,
    valid: bool,
    errors: usize,
    warnings: usize,
    issues: &'a [ModelIssue],
}

/// Machine-readable result of `stericx model inspect`.
#[derive(Debug, Serialize)]
struct InspectionReport<'a> {
    model_path: String,
    #[serde(flatten)]
    summary: &'a ModelSummary,
}

pub(crate) fn inspect_command(
    model_path: &Path,
    format: ReportFormat,
) -> Result<(), Box<dyn Error>> {
    let diagnosis = read_model(model_path)?;
    let Some(model) = diagnosis.model.as_ref() else {
        report_unreadable(model_path, &diagnosis, format)?;
        return Err(format!("{} could not be read", model_path.display()).into());
    };
    let summary = model.summary();

    match format {
        ReportFormat::Json => {
            let report = InspectionReport {
                model_path: model_path.display().to_string(),
                summary: &summary,
            };
            println!("{}", serde_json::to_string_pretty(&report)?);
        }
        ReportFormat::Text => print_summary(model_path, &summary, &diagnosis),
    }
    Ok(())
}

pub(crate) fn validate_command(
    model_path: &Path,
    format: ReportFormat,
    strict: bool,
) -> Result<(), Box<dyn Error>> {
    let diagnosis = read_model(model_path)?;
    let errors = diagnosis.count(Severity::Error);
    let warnings = diagnosis.count(Severity::Warning);
    let failed = errors > 0 || (strict && warnings > 0);

    match format {
        ReportFormat::Json => {
            let report = ValidationReport {
                model_path: model_path.display().to_string(),
                valid: !failed,
                errors,
                warnings,
                issues: &diagnosis.issues,
            };
            println!("{}", serde_json::to_string_pretty(&report)?);
        }
        ReportFormat::Text => {
            println!("command=model-validate");
            println!("model={}", model_path.display());
            for issue in &diagnosis.issues {
                println!(
                    "{}: {} [{}] {}",
                    severity_label(issue.severity),
                    issue.location,
                    issue.code,
                    issue.message
                );
            }
            println!("errors={errors}");
            println!("warnings={warnings}");
            println!("valid={}", !failed);
        }
    }

    if failed {
        return Err(format!(
            "{} failed validation with {errors} error(s) and {warnings} warning(s)",
            model_path.display()
        )
        .into());
    }
    Ok(())
}

fn read_model(model_path: &Path) -> Result<Diagnosis, Box<dyn Error>> {
    let text = fs::read_to_string(model_path)
        .map_err(|error| format!("could not read {}: {error}", model_path.display()))?;
    Ok(diagnose(&text))
}

/// Prints the findings for a document that could not be parsed into a model.
fn report_unreadable(
    model_path: &Path,
    diagnosis: &Diagnosis,
    format: ReportFormat,
) -> Result<(), Box<dyn Error>> {
    match format {
        ReportFormat::Json => {
            let report = ValidationReport {
                model_path: model_path.display().to_string(),
                valid: false,
                errors: diagnosis.count(Severity::Error),
                warnings: diagnosis.count(Severity::Warning),
                issues: &diagnosis.issues,
            };
            println!("{}", serde_json::to_string_pretty(&report)?);
        }
        ReportFormat::Text => {
            for issue in &diagnosis.issues {
                eprintln!(
                    "{}: {} [{}] {}",
                    severity_label(issue.severity),
                    issue.location,
                    issue.code,
                    issue.message
                );
            }
        }
    }
    Ok(())
}

fn print_summary(model_path: &Path, summary: &ModelSummary, diagnosis: &Diagnosis) {
    println!("command=model-inspect");
    println!("model={}", model_path.display());
    println!("schema_version={}", summary.schema_version);
    println!("portable={}", summary.portable);
    println!("model_name={}", summary.model);
    println!("reaction_family={}", optional(&summary.reaction_family));
    match summary.target.as_ref() {
        Some(response) => {
            println!("target={} [{}]", response.name, response.units);
            println!("target_description={}", response.description);
            println!("target_sign_convention={}", response.sign_convention);
            println!(
                "target_temperature_k={}",
                response
                    .temperature_k
                    .map_or_else(|| NOT_RECORDED.to_owned(), |value| format!("{value}"))
            );
        }
        None => println!("target={NOT_RECORDED}"),
    }
    println!("training_observations={}", summary.training_observations);
    println!("training_groups={}", summary.training_groups);
    println!("intercept={:.8}", summary.intercept);

    println!("descriptors={}", summary.descriptors.len());
    for descriptor in &summary.descriptors {
        println!(
            "  {name}: coefficient={coefficient:.8} mean={mean:.8} sd={sd:.8} \
             range=[{minimum:.8}, {maximum:.8}]",
            name = descriptor.name,
            coefficient = descriptor.coefficient,
            mean = descriptor.training_mean,
            sd = descriptor.training_standard_deviation,
            minimum = descriptor.training_minimum,
            maximum = descriptor.training_maximum,
        );
    }

    println!("training_r2={}", optional_metric(summary.training_r2));
    println!("training_rmse={:.8}", summary.training_rmse);
    println!("loo_q2={}", optional_metric(summary.loo_q2));
    println!("loo_rmse={:.8}", summary.loo_rmse);
    println!("group_loo_q2={}", optional_metric(summary.group_loo_q2));
    println!("group_loo_rmse={:.8}", summary.group_loo_rmse);

    if summary.dataset_digests.is_empty() {
        println!("dataset_digests={NOT_RECORDED}");
    } else {
        println!("dataset_digests={}", summary.dataset_digests.len());
        for digest in &summary.dataset_digests {
            println!(
                "  {artifact}: {algorithm}={digest} ({bytes} bytes)",
                artifact = digest.artifact,
                algorithm = digest.algorithm,
                digest = digest.digest,
                bytes = digest.byte_count,
            );
        }
    }

    println!("model_id={}", optional(&summary.model_id));
    println!("stericx_version={}", optional(&summary.stericx_version));
    println!("created_utc={}", optional(&summary.created_utc));
    println!(
        "missing_provenance={}",
        if summary.missing_provenance.is_empty() {
            "none".to_owned()
        } else {
            summary.missing_provenance.join(",")
        }
    );

    let errors = diagnosis.count(Severity::Error);
    let warnings = diagnosis.count(Severity::Warning);
    println!("validation_errors={errors}");
    println!("validation_warnings={warnings}");
    if errors > 0 {
        println!("note=run `stericx model validate` for the full list");
    }
}

/// Marker for a value the document genuinely does not record.
const NOT_RECORDED: &str = "not_recorded";

fn optional(value: &Option<String>) -> String {
    value
        .clone()
        .filter(|text| !text.trim().is_empty())
        .unwrap_or_else(|| NOT_RECORDED.to_owned())
}

fn optional_metric(value: Option<f64>) -> String {
    value.map_or_else(|| "unavailable".to_owned(), |value| format!("{value:.8}"))
}

fn severity_label(severity: Severity) -> &'static str {
    match severity {
        Severity::Error => "error",
        Severity::Warning => "warning",
    }
}
