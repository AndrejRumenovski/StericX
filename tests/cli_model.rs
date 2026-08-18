//! End-to-end tests for `stericx model inspect` and `stericx model validate`.
//!
//! These drive the built binary rather than the library, so they cover argument
//! parsing, exit codes, and both output formats. Each corruption case edits one
//! field of a valid document, so a failure points at one check.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

const EXE: &str = env!("CARGO_BIN_EXE_stericx");

/// A valid schema-2 document, built once and edited per test.
fn portable_document() -> serde_json::Value {
    let model: serde_json::Value = serde_json::from_str(
        &std::fs::read_to_string(fixture_path("portable_model.json")).unwrap(),
    )
    .unwrap();
    model
}

fn fixture_path(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

/// Writes a document to a unique temporary file and returns its path.
fn write_temp(value: &serde_json::Value) -> PathBuf {
    static COUNTER: AtomicUsize = AtomicUsize::new(0);
    let path = std::env::temp_dir().join(format!(
        "stericx_cli_model_{}_{}.json",
        std::process::id(),
        COUNTER.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::write(&path, serde_json::to_string_pretty(value).unwrap()).unwrap();
    path
}

struct Output {
    status: i32,
    stdout: String,
    stderr: String,
}

fn run(args: &[&str]) -> Output {
    let output = Command::new(EXE).args(args).output().expect("binary runs");
    Output {
        status: output.status.code().unwrap_or(-1),
        stdout: String::from_utf8_lossy(&output.stdout).into_owned(),
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
    }
}

fn validate(document: &serde_json::Value) -> Output {
    let path = write_temp(document);
    let output = run(&["model", "validate", path.to_str().unwrap()]);
    let _ = std::fs::remove_file(path);
    output
}

/// Asserts validation fails and names the expected issue code and message.
#[track_caller]
fn assert_rejected(document: &serde_json::Value, code: &str, fragment: &str) {
    let output = validate(document);
    assert_ne!(output.status, 0, "expected a non-zero exit for {code}");
    assert!(
        output.stdout.contains(code),
        "expected issue code {code}, got:\n{}",
        output.stdout
    );
    assert!(
        output.stdout.contains(fragment),
        "expected message containing {fragment:?}, got:\n{}",
        output.stdout
    );
    assert!(
        output.stdout.contains("valid=false"),
        "expected valid=false, got:\n{}",
        output.stdout
    );
}

#[test]
fn inspect_reports_the_full_scientific_summary() {
    let path = write_temp(&portable_document());

    let output = run(&["model", "inspect", path.to_str().unwrap()]);
    let _ = std::fs::remove_file(&path);

    assert_eq!(output.status, 0, "stderr: {}", output.stderr);
    let stdout = &output.stdout;
    for expected in [
        "command=model-inspect",
        "schema_version=2",
        "portable=true",
        "model_name=mechanistically_constrained_ols",
        "reaction_family=Ni-catalyzed homo-Diels-Alder",
        "target=ddg_double_dagger [kcal/mol]",
        "target_sign_convention=",
        "training_observations=10",
        "training_groups=9",
        "intercept=",
        "B5_x_nbo_charge: coefficient=",
        "mean=",
        "sd=",
        "range=[",
        "training_r2=",
        "loo_q2=",
        "loo_rmse=",
        "reactions.sigpack: fnv1a64=",
        "model_id=",
        "stericx_version=",
        "created_utc=",
        "validation_errors=0",
    ] {
        assert!(
            stdout.contains(expected),
            "inspect output is missing {expected:?}:\n{stdout}"
        );
    }
}

#[test]
fn inspect_emits_machine_readable_json() {
    let path = write_temp(&portable_document());

    let output = run(&[
        "model",
        "inspect",
        path.to_str().unwrap(),
        "--format",
        "json",
    ]);
    let _ = std::fs::remove_file(&path);

    assert_eq!(output.status, 0, "stderr: {}", output.stderr);
    let report: serde_json::Value = serde_json::from_str(&output.stdout).expect("stdout is JSON");
    assert_eq!(report["schema_version"], 2);
    assert_eq!(report["portable"], true);
    assert_eq!(report["training_observations"], 10);
    assert_eq!(report["target"]["units"], "kcal/mol");
    assert_eq!(report["descriptors"][0]["name"], "B5_x_nbo_charge");
    assert!(report["descriptors"][0]["coefficient"].is_number());
    assert!(report["descriptors"][0]["training_standard_deviation"].is_number());
    assert!(report["loo_q2"].is_number());
    assert!(report["loo_rmse"].is_number());
    assert_eq!(
        report["dataset_digests"][0]["algorithm"], "fnv1a64",
        "the digest algorithm must be stated, not assumed"
    );
    assert!(report["stericx_version"].is_string());
}

#[test]
fn validate_accepts_a_well_formed_model() {
    let output = validate(&portable_document());

    assert_eq!(output.status, 0, "stderr: {}", output.stderr);
    assert!(output.stdout.contains("errors=0"));
    assert!(output.stdout.contains("valid=true"));
}

#[test]
fn validate_emits_machine_readable_json() {
    let mut document = portable_document();
    document["provenance"]["training"]["dataset_digests"][0]["digest"] =
        serde_json::json!("nothexadecimal");
    let path = write_temp(&document);

    let output = run(&[
        "model",
        "validate",
        path.to_str().unwrap(),
        "--format",
        "json",
    ]);
    let _ = std::fs::remove_file(&path);

    assert_ne!(output.status, 0);
    let report: serde_json::Value = serde_json::from_str(&output.stdout).expect("stdout is JSON");
    assert_eq!(report["valid"], false);
    assert_eq!(report["errors"], 1);
    let issue = &report["issues"][0];
    assert_eq!(issue["severity"], "error");
    assert_eq!(issue["code"], "malformed_digest");
    assert_eq!(
        issue["location"], "provenance.training.dataset_digests[0]",
        "an issue must point at the field that carries it"
    );
}

#[test]
fn validate_rejects_an_unsupported_schema_version() {
    let mut document = portable_document();
    document["schema_version"] = serde_json::json!(99);

    assert_rejected(
        &document,
        "unsupported_schema_version",
        "newer than the supported maximum",
    );
}

#[test]
fn validate_rejects_a_missing_descriptor() {
    let mut document = portable_document();
    document["applicability_domain"] = serde_json::json!([]);

    assert_rejected(
        &document,
        "domain_length_mismatch",
        "does not cover every selected descriptor",
    );
}

#[test]
fn validate_rejects_dimension_mismatches() {
    let mut document = portable_document();
    let truncated = document["weights"].as_array().unwrap()[..7].to_vec();
    document["weights"] = serde_json::Value::Array(truncated);

    assert_rejected(&document, "dimension_mismatch", "`weights` has 7 entries");

    // A term count that does not match the selected descriptors.
    let mut document = portable_document();
    let extra = document["inference"]["terms"][0].clone();
    document["inference"]["terms"]
        .as_array_mut()
        .unwrap()
        .push(extra);

    assert_rejected(
        &document,
        "term_count_mismatch",
        "does not cover every selected descriptor",
    );
}

#[test]
fn validate_rejects_invalid_numerical_values() {
    let mut document = portable_document();
    document["weights"][0] = serde_json::Value::Null;

    assert_rejected(&document, "non_numeric_value", "`weights[0]` is null");
}

#[test]
fn validate_rejects_zero_and_negative_scaling() {
    for scale in [0.0_f64, -2.5] {
        let mut document = portable_document();
        let column = document["selected_feature_indices"][0].as_u64().unwrap() as usize;
        document["standardized_scales"][column] = serde_json::json!(scale);
        document["inference"]["terms"][0]["training_standard_deviation"] = serde_json::json!(scale);

        assert_rejected(&document, "invalid_scale", "not a positive finite number");
    }
}

#[test]
fn validate_rejects_missing_provenance() {
    let mut document = portable_document();
    document.as_object_mut().unwrap().remove("provenance");

    assert_rejected(&document, "missing_section", "missing the required");

    // A portable model must say which data trained it.
    let mut document = portable_document();
    document["provenance"]["training"]["dataset_digests"] = serde_json::json!([]);

    assert_rejected(
        &document,
        "missing_dataset_digest",
        "must identify its training data",
    );
}

#[test]
fn validate_rejects_malformed_hashes() {
    let mut document = portable_document();
    document["provenance"]["training"]["dataset_digests"][0]["digest"] =
        serde_json::json!("zzz-not-a-hash");

    assert_rejected(&document, "malformed_digest", "is not hexadecimal");
}

#[test]
fn validate_rejects_inconsistent_metadata() {
    // Provenance that disagrees with the fit about the training set.
    let mut document = portable_document();
    document["provenance"]["training"]["record_count"] = serde_json::json!(999);

    assert_rejected(
        &document,
        "training_count_mismatch",
        "provenance records 999 training rows but the fit reports 10",
    );

    // A coefficient edited in one copy but not the other.
    let mut document = portable_document();
    document["inference"]["terms"][0]["coefficient"] = serde_json::json!(42.0);

    assert_rejected(
        &document,
        "inference_disagrees_with_fit",
        "but the fit report says",
    );

    // A descriptor name that contradicts its column.
    let mut document = portable_document();
    document["selected_features"][0] = serde_json::json!("ir_frequency");

    assert_rejected(
        &document,
        "descriptor_name_mismatch",
        "does not match column",
    );
}

#[test]
fn validate_names_a_missing_required_field() {
    let mut document = portable_document();
    document.as_object_mut().unwrap().remove("training");

    assert_rejected(
        &document,
        "missing_field",
        "required field `training` is absent",
    );
}

#[test]
fn validate_explains_a_syntax_error_instead_of_a_decoder_dump() {
    let path = std::env::temp_dir().join(format!("stericx_cli_bad_{}.json", std::process::id()));
    std::fs::write(&path, "{ \"schema_version\": 2, oops }").unwrap();

    let output = run(&["model", "validate", path.to_str().unwrap()]);
    let _ = std::fs::remove_file(&path);

    assert_ne!(output.status, 0);
    assert!(output.stdout.contains("invalid_json"), "{}", output.stdout);
    assert!(
        output.stdout.contains("line 1 column"),
        "a syntax error must locate itself:\n{}",
        output.stdout
    );
}

#[test]
fn validate_reports_every_problem_at_once() {
    let mut document = portable_document();
    let column = document["selected_feature_indices"][0].as_u64().unwrap() as usize;
    document["standardized_scales"][column] = serde_json::json!(0.0);
    document["provenance"]["training"]["record_count"] = serde_json::json!(999);
    document["provenance"]["training"]["dataset_digests"][0]["digest"] = serde_json::json!("nope!");

    let output = validate(&document);

    assert_ne!(output.status, 0);
    for code in [
        "invalid_scale",
        "training_count_mismatch",
        "malformed_digest",
    ] {
        assert!(
            output.stdout.contains(code),
            "validate must list {code} alongside the others:\n{}",
            output.stdout
        );
    }
    assert!(
        output.stdout.contains("errors=4"),
        "expected four errors (two scaling, one count, one digest):\n{}",
        output.stdout
    );
}

#[test]
fn a_legacy_model_is_readable_but_reported_as_incomplete() {
    let legacy = Path::new(env!("CARGO_MANIFEST_DIR")).join("docs/study_001/stericx_model.json");

    let inspect = run(&["model", "inspect", legacy.to_str().unwrap()]);
    assert_eq!(inspect.status, 0, "stderr: {}", inspect.stderr);
    assert!(inspect.stdout.contains("schema_version=1"));
    assert!(inspect.stdout.contains("portable=false"));
    // Absent metadata is reported as absent, never invented.
    assert!(inspect.stdout.contains("target=not_recorded"));
    assert!(inspect.stdout.contains("dataset_digests=not_recorded"));
    assert!(inspect.stdout.contains("stericx_version=not_recorded"));
    // The science it does record is still shown.
    assert!(inspect.stdout.contains("training_observations=10"));
    assert!(inspect.stdout.contains("B5_x_nbo_charge: coefficient="));

    // A legacy model is valid, but warns that it cannot describe itself.
    let validate = run(&["model", "validate", legacy.to_str().unwrap()]);
    assert_eq!(validate.status, 0, "stderr: {}", validate.stderr);
    assert!(validate.stdout.contains("legacy_schema"));
    assert!(validate.stdout.contains("valid=true"));

    // Under --strict that warning becomes a failure.
    let strict = run(&["model", "validate", legacy.to_str().unwrap(), "--strict"]);
    assert_ne!(strict.status, 0);
    assert!(strict.stdout.contains("valid=false"));
}

#[test]
fn a_missing_file_is_reported_clearly() {
    let output = run(&["model", "inspect", "does/not/exist.json"]);

    assert_ne!(output.status, 0);
    assert!(
        output.stderr.contains("could not read"),
        "stderr: {}",
        output.stderr
    );
}
