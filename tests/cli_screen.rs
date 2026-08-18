//! End-to-end tests for `stericx screen`.
//!
//! The central guarantee is that screening a library reproduces the ordinary
//! inference path: whatever route a prediction takes through the library
//! loader, it must equal `RegressXPredictor` on the same descriptors. The rest
//! covers the reporting contract — identifiers, names, the descriptors a
//! prediction consumed, and the accounting for candidates that were excluded.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};
use steric_x::{
    FitOptions, PackedReactionRecord, ReactionLabel, RegressXPredictor, train_scientific_model,
};

const EXE: &str = env!("CARGO_BIN_EXE_stericx");

/// `stericx evaluate` accepts a frozen prediction within 1e-4 of the model's
/// own output; screening is held to the same standard. The gap is the f32
/// engine kernel against the f64 accumulation the screen loop uses.
const INFERENCE_TOLERANCE: f64 = 1.0e-4;

fn repo(path: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(path)
}

fn temp_path(extension: &str) -> PathBuf {
    static COUNTER: AtomicUsize = AtomicUsize::new(0);
    std::env::temp_dir().join(format!(
        "stericx_cli_screen_{}_{}.{extension}",
        std::process::id(),
        COUNTER.fetch_add(1, Ordering::Relaxed)
    ))
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

fn screen_json(model: &Path, library: &Path, extra: &[&str]) -> serde_json::Value {
    let mut args = vec![
        "screen",
        model.to_str().unwrap(),
        "--library",
        library.to_str().unwrap(),
        "--format",
        "json",
    ];
    args.extend_from_slice(extra);
    let output = run(&args);
    assert_eq!(output.status, 0, "screen failed: {}", output.stderr);
    serde_json::from_str(&output.stdout).expect("screen --format json emits JSON")
}

/// The Study 001 native model: selects `B5_x_nbo_charge`, so it needs both a
/// Sterimol term and donor electronics.
fn study_001_model() -> PathBuf {
    repo("docs/study_001/stericx_model.json")
}

/// The reaction table Study 001 was fitted from, usable directly as a library.
fn study_001_library() -> PathBuf {
    repo("data/reactions_raw.csv")
}

/// Fits a geometry-only model so the shipped Kraken library can be screened.
fn geometry_only_model() -> PathBuf {
    let records = (0..24_usize)
        .map(|index| {
            let l = 5.0 + 0.15 * index as f32;
            // Only `L` varies: every other column has zero training variance,
            // so forward selection cannot reach for an electronic term the
            // Kraken library could not supply.
            PackedReactionRecord {
                l,
                b1: 3.0,
                b5: 7.0,
                nbo_charge: -0.4,
                ir_freq: 1_650.0,
                temp_k: 298.15,
                exp_ddg: 0.5 + 0.6 * l + 0.01 * ((index * 11) % 5) as f32,
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
            ReactionLabel::new(
                format!("G{index:02}"),
                split,
                format!("group_{}", index % 4),
            )
        })
        .collect::<Vec<_>>();
    let trained = train_scientific_model(
        &records,
        &labels,
        FitOptions {
            bootstrap_samples: 30,
            permutation_samples: 30,
            ..FitOptions::default()
        },
    )
    .expect("geometry fixture trains");
    assert_eq!(
        trained.report.selected_features,
        ["L_boltz"],
        "the fixture must stay geometry-only for the library to satisfy it"
    );
    let path = temp_path("json");
    std::fs::write(
        &path,
        serde_json::to_string_pretty(&trained.report).unwrap(),
    )
    .unwrap();
    path
}

/// Recomputes a hit through the ordinary predictor using only the descriptors
/// the report says were consumed.
fn predict_from_reported_descriptors(weights: [f32; 8], hit: &serde_json::Value) -> f32 {
    let mut record = PackedReactionRecord::default();
    for descriptor in hit["descriptors"].as_array().expect("descriptors array") {
        let value = descriptor["value"].as_f64().unwrap() as f32;
        match descriptor["name"].as_str().unwrap() {
            "sterimol_l" => record.l = value,
            "sterimol_b1" => record.b1 = value,
            "sterimol_b5" => record.b5 = value,
            "nbo_charge" => record.nbo_charge = value,
            "ir_frequency" => record.ir_freq = value,
            other => panic!("unexpected descriptor {other}"),
        }
    }
    RegressXPredictor::new(weights).predict(&record)
}

fn model_weights(path: &Path) -> [f32; 8] {
    let model: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    let mut weights = [0.0_f32; 8];
    for (index, value) in model["weights"].as_array().unwrap().iter().enumerate() {
        weights[index] = value.as_f64().unwrap() as f32;
    }
    weights
}

#[test]
fn library_accepts_both_the_flag_and_positional_forms() {
    let model = study_001_model();
    let library = study_001_library();

    let flagged = run(&[
        "screen",
        model.to_str().unwrap(),
        "--library",
        library.to_str().unwrap(),
        "--format",
        "csv",
    ]);
    let positional = run(&[
        "screen",
        model.to_str().unwrap(),
        library.to_str().unwrap(),
        "--format",
        "csv",
    ]);

    assert_eq!(flagged.status, 0, "{}", flagged.stderr);
    assert_eq!(positional.status, 0, "{}", positional.stderr);
    assert_eq!(
        flagged.stdout, positional.stdout,
        "both spellings must screen identically"
    );

    // Supplying neither is an error, not an empty screen.
    let neither = run(&["screen", model.to_str().unwrap()]);
    assert_ne!(neither.status, 0);
    assert!(
        neither.stderr.contains("library is required"),
        "{}",
        neither.stderr
    );
}

#[test]
fn screening_agrees_with_the_ordinary_inference_path() {
    let model = study_001_model();
    let weights = model_weights(&model);
    let report = screen_json(&model, &study_001_library(), &[]);

    let hits = report["hits"].as_array().expect("hits");
    assert_eq!(hits.len(), 11, "every Study 001 reaction should screen");
    for hit in hits {
        let screened = hit["predicted_ddg_kcal_mol"].as_f64().unwrap();
        let engine = f64::from(predict_from_reported_descriptors(weights, hit));
        assert!(
            (screened - engine).abs() <= INFERENCE_TOLERANCE,
            "{}: screen said {screened} but the predictor said {engine}",
            hit["ligand"]
        );
    }
}

#[test]
fn study_001_blind_ligand_matches_its_frozen_prediction() {
    let report = screen_json(&study_001_model(), &study_001_library(), &[]);

    // The value frozen in docs/study_001/stericx_frozen_predictions.csv before
    // the experimental target was revealed.
    const FROZEN_723: f64 = 1.2330213;
    let hit = report["hits"]
        .as_array()
        .unwrap()
        .iter()
        .find(|hit| hit["ligand"] == "SIG-NIHDA-723")
        .expect("ligand 723 is in the library");
    let screened = hit["predicted_ddg_kcal_mol"].as_f64().unwrap();

    assert!(
        (screened - FROZEN_723).abs() <= INFERENCE_TOLERANCE,
        "screen predicted {screened} for ligand 723 but the frozen artifact says {FROZEN_723}"
    );
}

#[test]
fn screening_the_shipped_kraken_library_agrees_with_the_predictor() {
    // Studies 007 and 009 publish aggregate metrics rather than a fitted model
    // artifact, so their ligand space is covered here instead: the full shipped
    // Kraken phosphine database, screened end to end.
    let model = geometry_only_model();
    let weights = model_weights(&model);
    let report = screen_json(&model, &repo("data/ligand_db/kraken_phosphines.csv"), &[]);
    let _ = std::fs::remove_file(&model);

    let hits = report["hits"].as_array().expect("hits");
    assert!(
        hits.len() > 1_000,
        "expected the full library, screened {} ligands",
        hits.len()
    );
    for hit in hits {
        let screened = hit["predicted_ddg_kcal_mol"].as_f64().unwrap();
        let engine = f64::from(predict_from_reported_descriptors(weights, hit));
        assert!(
            (screened - engine).abs() <= INFERENCE_TOLERANCE,
            "{}: screen said {screened} but the predictor said {engine}",
            hit["ligand"]
        );
    }
}

#[test]
fn hits_carry_identifier_name_raw_prediction_and_descriptors() {
    let report = screen_json(&study_001_model(), &study_001_library(), &[]);
    let hit = &report["hits"][0];

    assert!(hit["ligand"].as_str().unwrap().starts_with("SIG-NIHDA-"));
    // The reaction table supplies SMILES, so a name is available here.
    assert!(
        hit["ligand_name"]
            .as_str()
            .is_some_and(|name| !name.is_empty()),
        "a library that carries names must report them"
    );
    assert!(hit["predicted_ddg_kcal_mol"].is_number(), "raw prediction");
    // The interpreted value is reported alongside, never instead of, the raw one.
    assert!(hit["predicted_ee_percent"].is_number());

    let descriptors = hit["descriptors"].as_array().unwrap();
    let names = descriptors
        .iter()
        .map(|descriptor| descriptor["name"].as_str().unwrap())
        .collect::<Vec<_>>();
    assert_eq!(
        names,
        ["sterimol_b5", "nbo_charge"],
        "only the descriptors the model consumes should be reported"
    );
    assert!(descriptors.iter().all(|d| d["value"].is_number()));
}

#[test]
fn a_library_without_a_required_descriptor_is_refused_not_guessed() {
    // The Kraken database has Sterimol terms but no donor electronics, and the
    // Study 001 model needs `nbo_charge`.
    let output = run(&[
        "screen",
        study_001_model().to_str().unwrap(),
        "--library",
        repo("data/ligand_db/kraken_phosphines.csv")
            .to_str()
            .unwrap(),
    ]);

    assert_ne!(output.status, 0, "screening must refuse, not approximate");
    assert!(
        output.stderr.contains("does not provide nbo_charge"),
        "the refusal must name the missing descriptor: {}",
        output.stderr
    );
    assert!(
        output.stderr.contains("will not guess"),
        "the refusal must say why: {}",
        output.stderr
    );
}

#[test]
fn ligands_missing_a_required_value_are_excluded_and_summarized() {
    // Blank one row's NBO charge; that ligand can no longer be screened.
    let source = std::fs::read_to_string(study_001_library()).unwrap();
    let base = repo("data");
    let mut lines = source.lines().map(str::to_owned).collect::<Vec<_>>();
    let header = lines[0].split(',').collect::<Vec<_>>();
    let nbo = header.iter().position(|c| *c == "NBO_Charge").unwrap();
    let geometry = header.iter().position(|c| *c == "Ligand_XYZ_Path").unwrap();
    // The copy lives outside the repo, so its geometry references have to be
    // absolute for the Sterimol terms to still resolve.
    for line in lines.iter_mut().skip(1) {
        let mut fields = line.split(',').map(str::to_owned).collect::<Vec<_>>();
        fields[geometry] = base.join(&fields[geometry]).display().to_string();
        *line = fields.join(",");
    }
    let mut fields = lines[1].split(',').map(str::to_owned).collect::<Vec<_>>();
    let excluded_id = fields[0].clone();
    fields[nbo] = String::new();
    lines[1] = fields.join(",");
    let library = temp_path("csv");
    std::fs::write(&library, lines.join("\n")).unwrap();

    let report = screen_json(&study_001_model(), &library, &[]);
    let _ = std::fs::remove_file(&library);

    assert_eq!(report["screened"], 10);
    assert_eq!(report["skipped"], 1);
    let excluded = report["excluded"].as_array().unwrap();
    assert_eq!(excluded.len(), 1);
    assert_eq!(excluded[0]["ligand"], serde_json::json!(excluded_id));
    assert_eq!(excluded[0]["reason"], "missing_descriptors");
    assert_eq!(
        excluded[0]["missing_descriptors"],
        serde_json::json!(["nbo_charge"]),
        "the exclusion must name the descriptor that was absent"
    );
    assert_eq!(
        report["exclusion_summary"],
        serde_json::json!([["missing_descriptors", 1]]),
        "exclusions must be summarized by reason"
    );
}

#[test]
fn every_output_format_carries_the_reported_fields() {
    let model = study_001_model();
    let library = study_001_library();

    let text = run(&[
        "screen",
        model.to_str().unwrap(),
        "--library",
        library.to_str().unwrap(),
        "--top",
        "3",
    ]);
    assert_eq!(text.status, 0, "{}", text.stderr);
    assert!(text.stdout.contains("rank"), "{}", text.stdout);
    assert!(text.stdout.contains("pred ddG"));
    assert!(text.stdout.contains("SIG-NIHDA-723"));
    // Terminal output pairs the identifier with the name and lists descriptors.
    assert!(text.stdout.contains("sterimol_b5="), "{}", text.stdout);
    assert!(text.stdout.contains("nbo_charge="));

    let csv = run(&[
        "screen",
        model.to_str().unwrap(),
        "--library",
        library.to_str().unwrap(),
        "--format",
        "csv",
        "--top",
        "3",
    ]);
    assert_eq!(csv.status, 0, "{}", csv.stderr);
    let header = csv.stdout.lines().next().unwrap();
    for column in [
        "rank",
        "ligand",
        "ligand_name",
        "predicted_ddg_kcal_mol",
        "predicted_ee_percent",
        "descriptors",
    ] {
        assert!(
            header.contains(column),
            "CSV header lacks {column}: {header}"
        );
    }
    assert_eq!(csv.stdout.lines().count(), 4, "header plus three rows");

    let json = screen_json(&model, &library, &["--top", "3"]);
    assert_eq!(json["hits"].as_array().unwrap().len(), 3);
    assert!(json["library_size"].is_number());
    assert!(json["screened"].is_number());
    assert!(json["excluded"].is_array());
    assert!(json["exclusion_summary"].is_array());
}

#[test]
fn a_model_declaring_an_unknown_feature_space_is_refused() {
    // A portable model whose feature space this build does not implement must
    // not be screened by assuming the built-in descriptor layout.
    let source = std::fs::read_to_string(repo("tests/data/portable_model.json")).unwrap();
    let mut document: serde_json::Value = serde_json::from_str(&source).unwrap();
    document["inference"]["feature_space"]["definition"] =
        serde_json::json!("someone.elses.space.v9");
    let path = temp_path("json");
    std::fs::write(&path, serde_json::to_string_pretty(&document).unwrap()).unwrap();

    let output = run(&[
        "screen",
        path.to_str().unwrap(),
        "--library",
        study_001_library().to_str().unwrap(),
    ]);
    let _ = std::fs::remove_file(&path);

    assert_ne!(output.status, 0);
    assert!(
        output.stderr.contains("someone.elses.space.v9"),
        "the refusal must name the unknown feature space: {}",
        output.stderr
    );
    assert!(
        output.stderr.contains("will not guess"),
        "{}",
        output.stderr
    );
}

#[test]
fn a_portable_model_screens_through_its_stored_feature_space() {
    let report = screen_json(
        &repo("tests/data/portable_model.json"),
        &study_001_library(),
        &[],
    );

    // The stored transformations name `sterimol_b5` and `nbo_charge`, so those
    // are what the screen requires and reports.
    assert_eq!(
        report["required_inputs"],
        serde_json::json!(["sterimol_b5", "nbo_charge"])
    );
    assert_eq!(report["hits"].as_array().unwrap().len(), 11);
}
