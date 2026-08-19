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
        "--descending",
    ]);
    let positional = run(&[
        "screen",
        model.to_str().unwrap(),
        library.to_str().unwrap(),
        "--format",
        "csv",
        "--descending",
    ]);

    assert_eq!(flagged.status, 0, "{}", flagged.stderr);
    assert_eq!(positional.status, 0, "{}", positional.stderr);
    assert_eq!(
        flagged.stdout, positional.stdout,
        "both spellings must screen identically"
    );

    // Supplying neither is an error, not an empty screen.
    let neither = run(&["screen", model.to_str().unwrap(), "--descending"]);
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
    let report = screen_json(&model, &study_001_library(), &["--descending"]);

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
    let report = screen_json(&study_001_model(), &study_001_library(), &["--descending"]);

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
    let report = screen_json(
        &model,
        &repo("data/ligand_db/kraken_phosphines.csv"),
        &["--descending"],
    );
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
    let report = screen_json(&study_001_model(), &study_001_library(), &["--descending"]);
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
        "--descending",
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

    let report = screen_json(&study_001_model(), &library, &["--descending"]);
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
        "--descending",
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
        "--descending",
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

    let json = screen_json(&model, &library, &["--top", "3", "--descending"]);
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
        "--descending",
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

// ---------------------------------------------------------------------------
// Ranking
// ---------------------------------------------------------------------------

/// Writes a portable model with a chosen optimization direction.
fn model_with_direction(direction: &str) -> PathBuf {
    let source = std::fs::read_to_string(repo("tests/data/portable_model.json")).unwrap();
    let mut document: serde_json::Value = serde_json::from_str(&source).unwrap();
    document["inference"]["response"]["optimization"] = serde_json::json!(direction);
    let path = temp_path("json");
    std::fs::write(&path, serde_json::to_string_pretty(&document).unwrap()).unwrap();
    path
}

/// A library of explicit descriptor values, so predictions are controllable.
///
/// The Study 001 model is `intercept + coefficient * (B5 * NBO)`, and its
/// coefficient is negative, so a larger `B5 * NBO` product gives a smaller
/// prediction.
fn descriptor_library(rows: &[(&str, f64, f64)]) -> PathBuf {
    let mut csv = String::from("ligand,sterimol_b5,nbo_charge\n");
    for (ligand, b5, nbo) in rows {
        csv.push_str(&format!("{ligand},{b5},{nbo}\n"));
    }
    let path = temp_path("csv");
    std::fs::write(&path, csv).unwrap();
    path
}

fn ranked_ligands(report: &serde_json::Value) -> Vec<String> {
    report["hits"]
        .as_array()
        .unwrap()
        .iter()
        .map(|hit| hit["ligand"].as_str().unwrap().to_owned())
        .collect()
}

fn predictions(report: &serde_json::Value) -> Vec<f64> {
    report["hits"]
        .as_array()
        .unwrap()
        .iter()
        .map(|hit| hit["predicted_ddg_kcal_mol"].as_f64().unwrap())
        .collect()
}

/// Spans a wide `B5 * NBO` range so predictions are clearly separated.
fn spread_library() -> PathBuf {
    descriptor_library(&[("low", 6.0, 0.9), ("mid", 9.0, 1.0), ("high", 12.0, 1.2)])
}

#[test]
fn a_minimize_model_ranks_ascending() {
    let model = model_with_direction("minimize");
    let library = spread_library();

    let report = screen_json(&model, &library, &[]);
    let _ = std::fs::remove_file(&model);
    let _ = std::fs::remove_file(&library);

    assert_eq!(report["ranking_order"], "ascending");
    assert_eq!(report["model_optimization"], "minimize");
    assert_eq!(report["ranking_overridden"], false);
    let values = predictions(&report);
    assert!(
        values.windows(2).all(|pair| pair[0] <= pair[1]),
        "a minimize model must list the smallest prediction first: {values:?}"
    );
}

#[test]
fn a_maximize_model_ranks_descending() {
    let model = model_with_direction("maximize");
    let library = spread_library();

    let report = screen_json(&model, &library, &[]);
    let _ = std::fs::remove_file(&model);
    let _ = std::fs::remove_file(&library);

    assert_eq!(report["ranking_order"], "descending");
    assert_eq!(report["model_optimization"], "maximize");
    let values = predictions(&report);
    assert!(
        values.windows(2).all(|pair| pair[0] >= pair[1]),
        "a maximize model must list the largest prediction first: {values:?}"
    );
}

#[test]
fn negative_predictions_rank_by_the_stated_direction() {
    // A large enough B5 * NBO product drives the prediction below zero, so this
    // library straddles the sign change.
    let library = descriptor_library(&[
        ("positive", 5.0, 1.0),
        ("near_zero", 16.3, 1.0),
        ("negative", 30.0, 1.0),
    ]);

    let maximize = model_with_direction("maximize");
    let report = screen_json(&maximize, &library, &[]);
    let values = predictions(&report);
    assert!(values.iter().any(|value| *value < 0.0), "{values:?}");
    assert!(
        values.windows(2).all(|pair| pair[0] >= pair[1]),
        "maximize must place the least negative first: {values:?}"
    );
    assert_eq!(ranked_ligands(&report).first().unwrap(), "positive");

    // Magnitude ranking treats a large negative as just as good as a large
    // positive — the selectivity is equal, the favoured enantiomer is not.
    let magnitude = model_with_direction("maximize_magnitude");
    let report = screen_json(&magnitude, &library, &[]);
    let values = predictions(&report);
    assert!(
        values.windows(2).all(|pair| pair[0].abs() >= pair[1].abs()),
        "magnitude ranking must order by absolute value: {values:?}"
    );
    assert_eq!(
        ranked_ligands(&report).last().unwrap(),
        "near_zero",
        "the least selective ligand ranks last under magnitude"
    );

    let _ = std::fs::remove_file(&maximize);
    let _ = std::fs::remove_file(&magnitude);
    let _ = std::fs::remove_file(&library);
}

#[test]
fn ties_break_deterministically_by_identifier() {
    // Identical descriptors, so identical predictions; only the tiebreak can
    // decide the order, and it must not depend on input order.
    let forward = descriptor_library(&[
        ("charlie", 8.0, 1.0),
        ("alpha", 8.0, 1.0),
        ("bravo", 8.0, 1.0),
    ]);
    let reversed = descriptor_library(&[
        ("bravo", 8.0, 1.0),
        ("alpha", 8.0, 1.0),
        ("charlie", 8.0, 1.0),
    ]);
    let model = model_with_direction("maximize");

    let first = screen_json(&model, &forward, &[]);
    let second = screen_json(&model, &reversed, &[]);
    let _ = std::fs::remove_file(&model);
    let _ = std::fs::remove_file(&forward);
    let _ = std::fs::remove_file(&reversed);

    let values = predictions(&first);
    assert!(
        values.windows(2).all(|pair| pair[0] == pair[1]),
        "the fixture must actually tie: {values:?}"
    );
    assert_eq!(
        ranked_ligands(&first),
        ["alpha", "bravo", "charlie"],
        "ties resolve by identifier"
    );
    assert_eq!(
        ranked_ligands(&first),
        ranked_ligands(&second),
        "the tiebreak must not depend on the order the library was written in"
    );
}

#[test]
fn repeated_runs_are_byte_identical() {
    let model = model_with_direction("maximize_magnitude");
    let library = study_001_library();

    let first = run(&[
        "screen",
        model.to_str().unwrap(),
        "--library",
        library.to_str().unwrap(),
        "--format",
        "csv",
    ]);
    let second = run(&[
        "screen",
        model.to_str().unwrap(),
        "--library",
        library.to_str().unwrap(),
        "--format",
        "csv",
    ]);
    let _ = std::fs::remove_file(&model);

    assert_eq!(first.status, 0, "{}", first.stderr);
    assert_eq!(
        first.stdout, second.stdout,
        "screening the same library twice must produce the same table"
    );
}

#[test]
fn top_limits_output_after_ranking_and_counts_are_reported() {
    let model = model_with_direction("maximize");

    let full = screen_json(&model, &study_001_library(), &[]);
    let limited = screen_json(&model, &study_001_library(), &["--top", "4"]);
    let _ = std::fs::remove_file(&model);

    assert_eq!(full["screened"], 11);
    assert_eq!(full["returned"], 11);
    // `--top` narrows what is reported, not what is predicted.
    assert_eq!(limited["screened"], 11, "every candidate is still screened");
    assert_eq!(limited["returned"], 4);
    assert_eq!(limited["hits"].as_array().unwrap().len(), 4);
    assert_eq!(
        ranked_ligands(&limited),
        ranked_ligands(&full)[..4],
        "the top N must be the first N of the full ranking"
    );
    // Ranks are positions in the returned table, starting at one.
    for (index, hit) in limited["hits"].as_array().unwrap().iter().enumerate() {
        assert_eq!(hit["rank"], serde_json::json!(index + 1));
    }
}

#[test]
fn an_override_is_reported_and_only_applies_when_requested() {
    let model = model_with_direction("maximize");

    let natural = screen_json(&model, &study_001_library(), &[]);
    assert_eq!(natural["ranking_overridden"], false);

    let overridden = screen_json(&model, &study_001_library(), &["--ascending"]);
    assert_eq!(overridden["ranking_order"], "ascending");
    assert_eq!(overridden["model_optimization"], "maximize");
    assert_eq!(overridden["ranking_overridden"], true);
    let mut reversed = ranked_ligands(&natural);
    reversed.reverse();
    assert_eq!(ranked_ligands(&overridden), reversed);

    // The terminal report says so in words, not just in a JSON field.
    let text = run(&[
        "screen",
        model.to_str().unwrap(),
        "--library",
        study_001_library().to_str().unwrap(),
        "--ascending",
        "--top",
        "1",
    ]);
    let _ = std::fs::remove_file(&model);
    assert!(
        text.stdout
            .contains("reverses the direction the model records"),
        "an override must be called out: {}",
        text.stdout
    );

    // Asking for the direction the model already states is not an override.
    let model = model_with_direction("maximize");
    let agreeing = screen_json(&model, &study_001_library(), &["--descending"]);
    let _ = std::fs::remove_file(&model);
    assert_eq!(agreeing["ranking_overridden"], false);
}

#[test]
fn a_model_without_a_direction_is_not_ranked_on_a_guess() {
    let model = model_with_direction("unspecified");

    let output = run(&[
        "screen",
        model.to_str().unwrap(),
        "--library",
        study_001_library().to_str().unwrap(),
    ]);

    assert_ne!(output.status, 0, "screening must refuse to invent an order");
    assert!(
        output.stderr.contains("will not rank on a guess"),
        "{}",
        output.stderr
    );
    assert!(
        output.stderr.contains("--optimize"),
        "the refusal must say how to record the direction: {}",
        output.stderr
    );

    // An explicit flag is enough to proceed.
    let explicit = screen_json(&model, &study_001_library(), &["--ascending"]);
    let _ = std::fs::remove_file(&model);
    assert_eq!(explicit["ranking_order"], "ascending");
    assert_eq!(explicit["model_optimization"], "unspecified");
    assert_eq!(
        explicit["ranking_overridden"], false,
        "a model with no stated direction cannot be contradicted"
    );
}

#[test]
fn csv_and_json_keep_full_precision() {
    let model = model_with_direction("maximize");
    let json = screen_json(&model, &study_001_library(), &["--top", "1"]);
    let csv = run(&[
        "screen",
        model.to_str().unwrap(),
        "--library",
        study_001_library().to_str().unwrap(),
        "--format",
        "csv",
        "--top",
        "1",
    ]);
    let _ = std::fs::remove_file(&model);

    let expected = json["hits"][0]["predicted_ddg_kcal_mol"].as_f64().unwrap();
    let row = csv.stdout.lines().nth(1).expect("one data row");
    let reported: f64 = row.split(',').nth(3).unwrap().parse().unwrap();

    assert_eq!(
        reported.to_bits(),
        expected.to_bits(),
        "the CSV value must read back as the same f64, not a rounded copy"
    );
    // A rounded field would have lost digits; a round-tripping one has not.
    assert!(
        (reported
            - f64::from(predict_from_reported_descriptors(
                model_weights(&repo("tests/data/portable_model.json")),
                &json["hits"][0]
            )))
        .abs()
            <= INFERENCE_TOLERANCE
    );
}

// ---------------------------------------------------------------------------
// Applicability domain
// ---------------------------------------------------------------------------

#[test]
fn screened_candidates_carry_structured_applicability_information() {
    // The Study 001 training set spans B5 * NBO from about 5.05 to 14.48. The
    // first row sits inside it; the second is far outside on both counts.
    let library = descriptor_library(&[("inside", 8.0, 1.0), ("far_outside", 25.0, 2.0)]);
    let model = model_with_direction("maximize_magnitude");

    let report = screen_json(&model, &library, &[]);
    let _ = std::fs::remove_file(&model);
    let _ = std::fs::remove_file(&library);

    let hits = report["hits"].as_array().unwrap();
    let find = |name: &str| {
        hits.iter()
            .find(|hit| hit["ligand"] == name)
            .unwrap_or_else(|| panic!("{name} was screened"))
            .clone()
    };

    let inside = find("inside");
    assert_eq!(inside["domain_verdict"], "interpolation");
    assert_eq!(inside["maximum_extrapolation"], 0.0);
    assert!(inside["nearest_training_distance"].as_f64().unwrap() >= 0.0);
    assert!(inside["nearest_training_ratio"].as_f64().unwrap() <= 1.0);
    assert!(inside["outside_domain"].as_array().unwrap().is_empty());

    let outside = find("far_outside");
    assert_eq!(outside["domain_verdict"], "extrapolation");
    assert!(outside["maximum_extrapolation"].as_f64().unwrap() > 0.0);
    assert!(!outside["outside_domain"].as_array().unwrap().is_empty());
    assert!(
        outside["nearest_training_distance"].as_f64().unwrap()
            > inside["nearest_training_distance"].as_f64().unwrap(),
        "the distant ligand must be farther from the training set"
    );
    assert!(
        outside["nearest_training_ratio"].as_f64().unwrap() > 1.0,
        "beyond the calibrated neighbour boundary"
    );

    // Mahalanobis is reported here because the covariance is estimable.
    assert!(outside["mahalanobis_distance"].as_f64().unwrap() > 0.0);

    // The report carries the derivation of the boundary, not just the verdict.
    let rule = report["neighbor_rule"].as_str().expect("boundary recorded");
    assert!(rule.contains("nearest other training point"), "{rule}");
    assert!(report["domain_summary"].is_array());
}

#[test]
fn applicability_does_not_depend_on_how_good_the_prediction_looks() {
    // Same two ligands, screened under opposite optimization directions. The
    // ranking flips; every applicability figure must be identical, because the
    // assessment never sees a prediction.
    let library = descriptor_library(&[("a", 7.0, 0.9), ("b", 13.0, 1.1)]);
    let maximize = model_with_direction("maximize");
    let minimize = model_with_direction("minimize");

    let up = screen_json(&maximize, &library, &[]);
    let down = screen_json(&minimize, &library, &[]);
    let _ = std::fs::remove_file(&maximize);
    let _ = std::fs::remove_file(&minimize);
    let _ = std::fs::remove_file(&library);

    assert_ne!(
        ranked_ligands(&up),
        ranked_ligands(&down),
        "the fixture must actually reorder"
    );

    let domain_of = |report: &serde_json::Value, name: &str| {
        let hit = report["hits"]
            .as_array()
            .unwrap()
            .iter()
            .find(|hit| hit["ligand"] == name)
            .unwrap();
        serde_json::json!({
            "verdict": hit["domain_verdict"],
            "distance": hit["nearest_training_distance"],
            "threshold": hit["nearest_training_threshold"],
            "ratio": hit["nearest_training_ratio"],
            "mahalanobis": hit["mahalanobis_distance"],
            "extrapolation": hit["maximum_extrapolation"],
        })
    };
    for ligand in ["a", "b"] {
        assert_eq!(
            domain_of(&up, ligand),
            domain_of(&down, ligand),
            "{ligand}: applicability changed with the ranking direction"
        );
    }
}
