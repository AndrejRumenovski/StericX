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

/// The Study 001 schema-2 portable artifact, emitted alongside the legacy
/// schema-1 model so the flagship model can be screened without direction flags.
fn study_001_portable_model() -> PathBuf {
    repo("docs/study_001/stericx_portable_model.json")
}

#[test]
fn the_study_001_portable_artifact_is_a_strict_superset_of_the_legacy_model() {
    let legacy: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(study_001_model()).unwrap()).unwrap();
    let portable: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(study_001_portable_model()).unwrap())
            .unwrap();

    assert_eq!(legacy["schema_version"], 1, "legacy artifact must stay v1");
    assert_eq!(portable["schema_version"], 2);

    // Every key the legacy document publishes must survive unchanged. A moved
    // validation statistic would silently republish a different scientific claim.
    // Supersetting is recursive: a nested object may gain keys (`training_geometry`
    // gains the applicability calibration) but never alter one it already had.
    fn assert_superset(legacy: &serde_json::Value, portable: &serde_json::Value, path: &str) {
        let Some(fields) = legacy.as_object() else {
            assert_eq!(legacy, portable, "portable model changed `{path}`");
            return;
        };
        let nested = portable
            .as_object()
            .unwrap_or_else(|| panic!("portable model replaced object `{path}`"));
        for (key, value) in fields {
            let child = nested
                .get(key)
                .unwrap_or_else(|| panic!("portable model dropped `{path}{key}`"));
            assert_superset(value, child, &format!("{path}{key}."));
        }
    }

    for (key, value) in legacy.as_object().unwrap() {
        if key == "schema_version" {
            continue;
        }
        let child = portable
            .get(key)
            .unwrap_or_else(|| panic!("portable model dropped legacy field `{key}`"));
        assert_superset(value, child, &format!("{key}."));
    }

    for added in ["inference", "provenance", "created"] {
        assert!(
            portable.get(added).is_some(),
            "portable model is missing `{added}`"
        );
    }
}

#[test]
fn the_study_001_portable_artifact_validates_without_findings() {
    let output = run(&[
        "model",
        "validate",
        study_001_portable_model().to_str().unwrap(),
    ]);
    assert_eq!(output.status, 0, "validate failed: {}", output.stderr);
    assert!(output.stdout.contains("errors=0"), "{}", output.stdout);
    assert!(output.stdout.contains("warnings=0"), "{}", output.stdout);
    assert!(output.stdout.contains("valid=true"), "{}", output.stdout);
}

#[test]
fn the_study_001_portable_artifact_ranks_without_a_direction_flag() {
    // The whole point of the artifact: the legacy model refuses to rank, this
    // one carries its own direction.
    let refused = run(&[
        "screen",
        study_001_model().to_str().unwrap(),
        "--library",
        study_001_library().to_str().unwrap(),
    ]);
    assert_ne!(
        refused.status, 0,
        "the legacy model must still refuse to rank on a guess"
    );

    let report = screen_json(&study_001_portable_model(), &study_001_library(), &[]);
    assert_eq!(report["model_optimization"], "maximize");
    assert_eq!(report["ranking_order"], "descending");
    assert_eq!(report["ranking_overridden"], false);

    // Recorded direction means descending raw value, ties aside.
    let values = predictions(&report);
    let mut sorted = values.clone();
    sorted.sort_by(|left, right| right.total_cmp(left));
    assert_eq!(values, sorted, "hits are not in the model's stated order");
}

#[test]
fn the_study_001_portable_artifact_reports_applicability_not_unknown() {
    // The legacy artifact predates `neighbor_calibration`, so it can only say
    // `unknown`. The portable one carries the calibration and can say more.
    let report = screen_json(&study_001_portable_model(), &study_001_library(), &[]);
    for hit in report["hits"].as_array().unwrap() {
        assert_ne!(
            hit["domain_verdict"], "unknown",
            "{} has no applicability verdict",
            hit["ligand"]
        );
    }
}

#[test]
fn the_study_001_portable_artifact_reproduces_the_frozen_blind_prediction() {
    let report = screen_json(&study_001_portable_model(), &study_001_library(), &[]);
    let blind = report["hits"]
        .as_array()
        .unwrap()
        .iter()
        .find(|hit| hit["ligand"] == "SIG-NIHDA-723")
        .expect("the blind ligand is in the library");
    let predicted = blind["predicted_ddg_kcal_mol"].as_f64().unwrap();
    // `stericx evaluate`'s own tolerance for a frozen prediction: the screen
    // loop accumulates in f64 while the engine kernel is f32.
    assert!(
        (predicted - 1.233_021_3_f64).abs() < 1e-4,
        "blind prediction drifted: {predicted}"
    );
}

/// The Study 001 model's own training geometry: standardized 1-D points, the
/// standardization constants, and the recorded neighbour calibration.
fn study_001_geometry() -> (Vec<f64>, f64, f64, serde_json::Value) {
    let document: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(study_001_portable_model()).unwrap())
            .unwrap();
    let geometry = &document["training_geometry"];
    let mut points = geometry["standardized_training_points"]
        .as_array()
        .expect("the schema-2 artifact records its training points")
        .iter()
        .map(|row| row[0].as_f64().unwrap())
        .collect::<Vec<_>>();
    points.sort_by(f64::total_cmp);
    (
        points,
        geometry["means"][0].as_f64().unwrap(),
        geometry["scales"][0].as_f64().unwrap(),
        geometry["neighbor_calibration"].clone(),
    )
}

/// Turns a standardized coordinate back into a `B5 x NBO` pair the screen
/// library can carry. `nbo_charge` is 1.0, so the product is `sterimol_b5`.
fn at_standardized(z: f64, mean: f64, scale: f64) -> f64 {
    z * scale + mean
}

#[test]
fn a_real_gap_in_the_training_set_reports_sparse_interpolation() {
    let (points, mean, scale, calibration) = study_001_geometry();
    let threshold = calibration["threshold"].as_f64().unwrap();

    // The widest gap the published training set actually leaves.
    let (width, low, high) = points
        .windows(2)
        .map(|pair| (pair[1] - pair[0], pair[0], pair[1]))
        .max_by(|left, right| left.0.total_cmp(&right.0))
        .unwrap();
    assert!(
        width / 2.0 > threshold,
        "the fixture assumes the widest training gap outruns the calibrated \
         spacing: half-gap {} vs threshold {threshold}",
        width / 2.0
    );

    let library = descriptor_library(&[
        // Centre of that gap: in range, but farther from every training point
        // than the training set's own sparsest spacing.
        (
            "in_gap",
            at_standardized((low + high) / 2.0, mean, scale),
            1.0,
        ),
        // An actual training coordinate: distance zero.
        ("on_training_point", at_standardized(high, mean, scale), 1.0),
    ]);

    let report = screen_json(&study_001_portable_model(), &library, &[]);
    let hit = |name: &str| {
        report["hits"]
            .as_array()
            .unwrap()
            .iter()
            .find(|hit| hit["ligand"] == name)
            .unwrap_or_else(|| panic!("{name} was screened"))
            .clone()
    };

    let gap = hit("in_gap");
    assert_eq!(gap["domain_verdict"], "sparse_interpolation");
    // Sparse is emphatically not extrapolation: every descriptor is in range.
    assert_eq!(gap["maximum_extrapolation"].as_f64().unwrap(), 0.0);
    assert!(
        gap["outside_domain"]
            .as_array()
            .is_none_or(std::vec::Vec::is_empty)
    );
    assert!(gap["nearest_training_ratio"].as_f64().unwrap() > 1.0);

    let on_point = hit("on_training_point");
    assert_eq!(on_point["domain_verdict"], "interpolation");
    assert!(on_point["nearest_training_distance"].as_f64().unwrap() < 1e-9);
}

#[test]
fn a_stricter_domain_rule_shrinks_the_boundary_without_moving_the_measurement() {
    let (points, mean, scale, calibration) = study_001_geometry();
    let permissive = calibration["threshold"].as_f64().unwrap();
    let strict = calibration["mean"].as_f64().unwrap()
        + 2.0 * calibration["standard_deviation"].as_f64().unwrap();
    assert!(
        strict < permissive,
        "Study 001's spacing must actually be uneven for this test to mean anything"
    );

    // A coordinate between the two boundaries: accepted by the default rule,
    // rejected by mean + 2σ. It has to sit *inside* the training range, so
    // place it in the widest interior gap rather than past an end point --
    // outside the range it would be extrapolation and the rule would not matter.
    let (width, low, _high) = points
        .windows(2)
        .map(|pair| (pair[1] - pair[0], pair[0], pair[1]))
        .max_by(|left, right| left.0.total_cmp(&right.0))
        .unwrap();
    let offset = f64::midpoint(strict, permissive);
    assert!(
        offset < width / 2.0,
        "the gap must be wide enough to hold a point at {offset} from its edge"
    );
    let library = descriptor_library(&[(
        "between_boundaries",
        at_standardized(low + offset, mean, scale),
        1.0,
    )]);

    let verdict_and_distance = |extra: &[&str]| {
        let report = screen_json(&study_001_portable_model(), &library, extra);
        let hit = report["hits"][0].clone();
        (
            report["domain_rule"].as_str().unwrap().to_owned(),
            report["domain_threshold"].as_f64().unwrap(),
            hit["domain_verdict"].as_str().unwrap().to_owned(),
            hit["nearest_training_distance"].as_f64().unwrap(),
        )
    };

    let (default_rule, default_threshold, default_verdict, default_distance) =
        verdict_and_distance(&[]);
    let (strict_rule, strict_threshold, strict_verdict, strict_distance) =
        verdict_and_distance(&["--domain-rule", "mean-plus-2sd"]);

    assert_eq!(
        default_rule, "max_neighbor",
        "the default must stay permissive"
    );
    assert_eq!(strict_rule, "mean_plus_2sd");
    assert_eq!(default_verdict, "interpolation");
    assert_eq!(strict_verdict, "sparse_interpolation");
    assert!(strict_threshold < default_threshold);

    // Only the boundary moved. The measured distance is a property of the
    // candidate and the training set, so it must be identical under both rules.
    assert_eq!(default_distance, strict_distance);
}

#[test]
fn the_domain_rule_is_reported_so_a_stricter_run_is_never_ambiguous() {
    let report = screen_json(
        &study_001_portable_model(),
        &study_001_library(),
        &["--domain-rule", "mean-plus-sd"],
    );
    assert_eq!(report["domain_rule"], "mean_plus_sd");
    assert!(
        report["domain_rule_description"]
            .as_str()
            .unwrap()
            .contains("mean + 1 standard deviation"),
        "the report must carry the derivation, not just a name"
    );

    let text = run(&[
        "screen",
        study_001_portable_model().to_str().unwrap(),
        "--library",
        study_001_library().to_str().unwrap(),
        "--domain-rule",
        "mean-plus-sd",
    ]);
    assert_eq!(text.status, 0, "{}", text.stderr);
    assert!(
        text.stdout.contains("domain rule    mean_plus_sd"),
        "{}",
        text.stdout
    );
}

#[test]
fn the_domain_rule_cannot_change_a_prediction() {
    // Applicability is scored from descriptors alone; changing the boundary
    // must not perturb a single predicted value or the ranking.
    let library = spread_library();
    let default = screen_json(&study_001_portable_model(), &library, &[]);
    let strict = screen_json(
        &study_001_portable_model(),
        &library,
        &["--domain-rule", "mean-plus-2sd"],
    );
    assert_eq!(predictions(&default), predictions(&strict));
    assert_eq!(ranked_ligands(&default), ranked_ligands(&strict));
}

/// A library whose best-predicted member is deliberately out of domain.
///
/// The Study 001 coefficient on `B5 x NBO` is negative, so the *smallest*
/// product gives the largest predicted ddG. Putting that smallest product below
/// the training minimum makes the top-ranked candidate an extrapolation.
fn best_candidate_is_out_of_domain_library() -> PathBuf {
    descriptor_library(&[
        ("best_but_ood", 3.0, 1.0),
        ("safe_mid", 10.0, 1.0),
        ("safe_high", 12.0, 1.0),
    ])
}

#[test]
fn the_highest_predicted_candidate_can_still_be_flagged_out_of_domain() {
    let report = screen_json(
        &study_001_portable_model(),
        &best_candidate_is_out_of_domain_library(),
        &[],
    );

    let top = &report["hits"][0];
    assert_eq!(top["ligand"], "best_but_ood");
    assert_eq!(top["rank"], 1);
    assert_eq!(
        top["domain_verdict"], "extrapolation",
        "the best-looking prediction must not buy its way into the domain"
    );
    assert!(
        top["maximum_extrapolation"].as_f64().unwrap() > 0.0,
        "an extrapolating hit must quantify how far outside it is"
    );
    assert!(
        !top["outside_domain"].as_array().unwrap().is_empty(),
        "the offending descriptor must be named"
    );

    // Ranking is by prediction alone: being out of domain neither promotes nor
    // demotes it.
    let values = predictions(&report);
    assert!(
        values[0] > values[1],
        "the out-of-domain candidate really is the top prediction"
    );
}

#[test]
fn an_out_of_domain_candidate_is_never_dropped_or_altered_by_default() {
    let library = best_candidate_is_out_of_domain_library();
    let report = screen_json(&study_001_portable_model(), &library, &[]);

    assert_eq!(report["screened"], 3, "nothing is dropped without the flag");
    assert_eq!(report["returned"], 3);
    assert_eq!(report["domain_filter_applied"], false);
    assert!(report["domain_filtered"].as_array().unwrap().is_empty());
    assert_eq!(ranked_ligands(&report).len(), 3);

    // The prediction must be the raw model output, not softened toward the
    // training mean because the candidate is outside the domain.
    let top = &report["hits"][0];
    let predicted = top["predicted_ddg_kcal_mol"].as_f64().unwrap();
    let recomputed = predict_from_reported_descriptors(
        model_weights(&study_001_portable_model()),
        &serde_json::json!({ "descriptors": top["descriptors"].clone() }),
    );
    assert!(
        (predicted - f64::from(recomputed)).abs() < 1e-4,
        "an out-of-domain prediction was modified: {predicted} vs {recomputed}"
    );
}

#[test]
fn in_domain_only_is_opt_in_and_reports_what_it_removed() {
    let library = best_candidate_is_out_of_domain_library();
    let filtered = screen_json(&study_001_portable_model(), &library, &["--in-domain-only"]);

    assert_eq!(filtered["domain_filter_applied"], true);
    assert_eq!(
        filtered["screened"], 2,
        "the extrapolating candidate is gone"
    );
    let removed = filtered["domain_filtered"].as_array().unwrap();
    assert_eq!(removed.len(), 1);
    assert_eq!(removed[0][0], "best_but_ood");
    assert_eq!(removed[0][1], "extrapolation");

    // The census still describes the whole library, so filtering never hides
    // how many candidates existed.
    let census = filtered["domain_summary"].as_array().unwrap();
    let total: u64 = census.iter().map(|entry| entry[1].as_u64().unwrap()).sum();
    assert_eq!(total, 3, "the census must cover every screened candidate");

    // Surviving predictions are untouched by the filter.
    let unfiltered = screen_json(&study_001_portable_model(), &library, &[]);
    assert_eq!(
        predictions(&filtered),
        predictions(&unfiltered)[1..].to_vec()
    );
}

#[test]
fn the_terminal_table_shows_the_domain_beside_the_prediction() {
    let output = run(&[
        "screen",
        study_001_portable_model().to_str().unwrap(),
        "--library",
        best_candidate_is_out_of_domain_library().to_str().unwrap(),
    ]);
    assert_eq!(output.status, 0, "{}", output.stderr);

    assert!(
        output.stdout.contains("domain") && output.stdout.contains("pred ddG"),
        "the table must carry a domain column beside the prediction:\n{}",
        output.stdout
    );
    // The warning is visible, and the prediction it applies to is still printed.
    assert!(
        output.stdout.contains("extrapolation!"),
        "an extrapolating row must be marked:\n{}",
        output.stdout
    );
    assert!(
        output.stdout.contains("best_but_ood"),
        "a flagged ligand must still appear:\n{}",
        output.stdout
    );
    assert!(
        output.stdout.contains("--in-domain-only to exclude them"),
        "the report must say the filter exists:\n{}",
        output.stdout
    );
}

#[test]
fn filtering_every_candidate_explains_itself_rather_than_reporting_nothing() {
    // A library that is entirely out of domain: the error must name the filter
    // instead of implying the model could not screen anything.
    let library = descriptor_library(&[("way_out", 2.0, 1.0), ("also_out", 2.5, 1.0)]);
    let output = run(&[
        "screen",
        study_001_portable_model().to_str().unwrap(),
        "--library",
        library.to_str().unwrap(),
        "--in-domain-only",
    ]);
    assert_ne!(output.status, 0);
    assert!(
        output.stderr.contains("--in-domain-only removed all"),
        "{}",
        output.stderr
    );
}

/// A library spanning the training centre, an in-domain extreme, and a point
/// just past the training maximum.
fn uncertainty_probe_library() -> PathBuf {
    descriptor_library(&[
        ("centre", 9.97, 1.0),
        ("in_domain_low", 5.10, 1.0),
        ("ood_above_max", 14.70, 1.0),
    ])
}

fn uncertainty_of(report: &serde_json::Value, ligand: &str) -> serde_json::Value {
    report["hits"]
        .as_array()
        .unwrap()
        .iter()
        .find(|hit| hit["ligand"] == ligand)
        .unwrap_or_else(|| panic!("{ligand} was screened"))["uncertainty"]
        .clone()
}

#[test]
fn every_bootstrap_interval_is_ordered_and_brackets_its_point_estimate() {
    let report = screen_json(&study_001_portable_model(), &study_001_library(), &[]);
    assert_eq!(
        report["uncertainty_method"],
        "percentile_bootstrap_mean_response"
    );
    assert_eq!(report["uncertainty_level"], 0.95);
    assert!(report["uncertainty_replicates"].as_u64().unwrap() > 0);

    for hit in report["hits"].as_array().unwrap() {
        let uncertainty = &hit["uncertainty"];
        let low = uncertainty["lower"].as_f64().unwrap();
        let high = uncertainty["upper"].as_f64().unwrap();
        let central = hit["predicted_ddg_kcal_mol"].as_f64().unwrap();
        assert!(low <= high, "{}: inverted interval", hit["ligand"]);
        assert!(
            low <= central && central <= high,
            "{}: point estimate {central} outside [{low}, {high}]",
            hit["ligand"]
        );
        assert_eq!(uncertainty["level"], 0.95);
        assert_eq!(
            uncertainty["replicates"], report["uncertainty_replicates"],
            "per-hit replicate count must match the ensemble"
        );
        // The name must not claim more than the interval delivers.
        let method = uncertainty["method"].as_str().unwrap();
        assert!(method.contains("mean_response"), "{method}");
        assert!(
            !method.contains("prediction_interval"),
            "a coefficient-only interval must not be called predictive: {method}"
        );
    }
}

#[test]
fn a_narrow_bootstrap_interval_does_not_make_a_candidate_in_domain() {
    // The whole point of keeping the two axes separate. Just past the training
    // maximum the coefficient band is *narrower* than it is at an in-domain
    // point near the other end of the range, because interval width tracks the
    // bootstrap coefficient spread rather than membership of the training set.
    let report = screen_json(
        &study_001_portable_model(),
        &uncertainty_probe_library(),
        &[],
    );
    let width = |ligand: &str| {
        let uncertainty = uncertainty_of(&report, ligand);
        uncertainty["upper"].as_f64().unwrap() - uncertainty["lower"].as_f64().unwrap()
    };
    let verdict = |ligand: &str| {
        report["hits"]
            .as_array()
            .unwrap()
            .iter()
            .find(|hit| hit["ligand"] == ligand)
            .unwrap()["domain_verdict"]
            .as_str()
            .unwrap()
            .to_owned()
    };

    assert_eq!(verdict("ood_above_max"), "extrapolation");
    assert_eq!(verdict("in_domain_low"), "interpolation");
    assert!(
        width("ood_above_max") < width("in_domain_low"),
        "fixture assumes the out-of-domain candidate has the narrower band: {} vs {}",
        width("ood_above_max"),
        width("in_domain_low")
    );

    // Narrower band, and still unambiguously flagged. Nothing in the report
    // presents the tighter interval as evidence of reliability.
    let flagged = uncertainty_of(&report, "ood_above_max");
    assert!(
        flagged["covers"]
            .as_str()
            .unwrap()
            .contains("extrapolation"),
        "the interval must disclaim extrapolation risk: {flagged:?}"
    );
    assert!(
        report["uncertainty_note"].as_str().is_some(),
        "the report must explain what the intervals mean"
    );
}

#[test]
fn the_bootstrap_ensemble_survives_serialization() {
    let document: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(study_001_portable_model()).unwrap())
            .unwrap();
    let ensemble = &document["uncertainty"];
    assert!(
        ensemble.is_object(),
        "the schema-2 document must carry the bootstrap ensemble"
    );

    let replicates = ensemble["replicates"].as_array().unwrap();
    let columns = ensemble["columns"].as_array().unwrap();
    let indices = ensemble["column_indices"].as_array().unwrap();
    assert_eq!(columns.len(), indices.len());
    assert_eq!(columns[0], "intercept", "the intercept column comes first");
    assert_eq!(
        ensemble["replicate_count"].as_u64().unwrap() as usize,
        replicates.len()
    );
    assert!(ensemble["seed"].as_u64().is_some(), "the seed is recorded");
    for replicate in replicates {
        assert_eq!(
            replicate.as_array().unwrap().len(),
            columns.len(),
            "every replicate must carry one coefficient per column"
        );
    }

    // Re-reading the document and screening from it reproduces the interval, so
    // the stored ensemble really is sufficient without the training data.
    let report = screen_json(&study_001_portable_model(), &study_001_library(), &[]);
    assert_eq!(
        report["uncertainty_replicates"].as_u64().unwrap() as usize,
        replicates.len()
    );
}

#[test]
fn repeated_inference_reproduces_the_same_intervals() {
    let library = uncertainty_probe_library();
    let first = screen_json(&study_001_portable_model(), &library, &[]);
    let second = screen_json(&study_001_portable_model(), &library, &[]);
    assert_eq!(first, second, "a stored ensemble must be deterministic");

    for ligand in ["centre", "in_domain_low", "ood_above_max"] {
        assert_eq!(
            uncertainty_of(&first, ligand),
            uncertainty_of(&second, ligand)
        );
    }
}

#[test]
fn a_model_without_an_ensemble_reports_no_interval_rather_than_inventing_one() {
    // The legacy schema-1 artifact never carried an ensemble.
    let report = screen_json(&study_001_model(), &study_001_library(), &["--descending"]);
    assert!(report["uncertainty_method"].is_null());
    assert!(
        report["uncertainty_unavailable"]
            .as_str()
            .unwrap()
            .contains("no bootstrap ensemble"),
        "the absence must be explained"
    );
    for hit in report["hits"].as_array().unwrap() {
        assert!(
            hit["uncertainty"].is_null(),
            "{} invented an interval",
            hit["ligand"]
        );
    }
}

fn sha256_hex_of(path: &Path) -> String {
    // Independent of StericX's own implementation: shelling out to sha256sum
    // would couple the test to the environment, so hash with a second source
    // of truth only where one is already in the repo. Here we simply assert
    // internal consistency and the documented shape.
    let bytes = std::fs::read(path).unwrap();
    let mut hasher = Sha256::new();
    hasher.update(&bytes);
    hasher.finish()
}

/// Minimal SHA-256, written independently of `src/digest.rs` so the test does
/// not confirm the implementation against itself.
struct Sha256 {
    state: [u32; 8],
    buffer: Vec<u8>,
}

impl Sha256 {
    fn new() -> Self {
        Self {
            state: [
                0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab,
                0x5be0cd19,
            ],
            buffer: Vec::new(),
        }
    }

    fn update(&mut self, bytes: &[u8]) {
        self.buffer.extend_from_slice(bytes);
    }

    fn finish(mut self) -> String {
        const K: [u32; 64] = [
            0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4,
            0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe,
            0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f,
            0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
            0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc,
            0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
            0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116,
            0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
            0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7,
            0xc67178f2,
        ];
        let bit_length = (self.buffer.len() as u64) * 8;
        self.buffer.push(0x80);
        while self.buffer.len() % 64 != 56 {
            self.buffer.push(0);
        }
        self.buffer.extend_from_slice(&bit_length.to_be_bytes());

        for chunk in self.buffer.chunks(64) {
            let mut w = [0_u32; 64];
            for (index, word) in chunk.chunks(4).enumerate() {
                w[index] = u32::from_be_bytes([word[0], word[1], word[2], word[3]]);
            }
            for index in 16..64 {
                let s0 = w[index - 15].rotate_right(7)
                    ^ w[index - 15].rotate_right(18)
                    ^ (w[index - 15] >> 3);
                let s1 = w[index - 2].rotate_right(17)
                    ^ w[index - 2].rotate_right(19)
                    ^ (w[index - 2] >> 10);
                w[index] = w[index - 16]
                    .wrapping_add(s0)
                    .wrapping_add(w[index - 7])
                    .wrapping_add(s1);
            }
            let [mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut h] = self.state;
            for index in 0..64 {
                let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
                let ch = (e & f) ^ ((!e) & g);
                let temp1 = h
                    .wrapping_add(s1)
                    .wrapping_add(ch)
                    .wrapping_add(K[index])
                    .wrapping_add(w[index]);
                let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
                let maj = (a & b) ^ (a & c) ^ (b & c);
                let temp2 = s0.wrapping_add(maj);
                h = g;
                g = f;
                f = e;
                e = d.wrapping_add(temp1);
                d = c;
                c = b;
                b = a;
                a = temp1.wrapping_add(temp2);
            }
            for (slot, value) in self.state.iter_mut().zip([a, b, c, d, e, f, g, h]) {
                *slot = slot.wrapping_add(value);
            }
        }
        self.state
            .iter()
            .map(|word| format!("{word:08x}"))
            .collect()
    }
}

#[test]
fn a_screening_run_reports_the_model_context_before_any_candidate() {
    let output = run(&[
        "screen",
        study_001_portable_model().to_str().unwrap(),
        "--library",
        study_001_library().to_str().unwrap(),
        "--top",
        "1",
    ]);
    assert_eq!(output.status, 0, "{}", output.stderr);
    let stdout = &output.stdout;
    let table_start = stdout
        .find("rank ")
        .expect("the candidate table is printed");
    let header = &stdout[..table_start];

    // Everything a reader needs to know what produced these numbers, stated
    // before the numbers themselves.
    for required in [
        "model id",
        "reaction",
        "target",
        "selected",
        "trained on",
        "validation",
        "LOO Q²",
        "stericx",
        "model sha256",
        "training data",
        "library sha256",
        "ranking",
        "uncertainty",
        "domain rule",
    ] {
        assert!(
            header.contains(required),
            "context block is missing `{required}`:\n{header}"
        );
    }
}

#[test]
fn the_report_carries_cryptographic_identity_for_every_scientific_input() {
    let report = screen_json(&study_001_portable_model(), &study_001_library(), &[]);
    let provenance = &report["provenance"];

    // Model and library hashes must be real SHA-256 of the exact bytes,
    // recomputed here by an independent implementation.
    assert_eq!(
        provenance["model_sha256"].as_str().unwrap(),
        sha256_hex_of(&study_001_portable_model()),
        "model digest does not match the file"
    );
    assert_eq!(
        provenance["library_sha256"].as_str().unwrap(),
        sha256_hex_of(&study_001_library()),
        "library digest does not match the file"
    );

    assert_eq!(
        provenance["stericx_version"],
        env!("CARGO_PKG_VERSION"),
        "the running binary must identify itself"
    );
    assert_eq!(provenance["model_schema_version"], 2);
    assert_eq!(provenance["model_id"], "mechanistically_constrained_ols");
    assert!(provenance["model_fitted_by"].as_str().is_some());

    // Training inputs, each tagged with the algorithm actually used.
    let datasets = provenance["training_datasets"].as_array().unwrap();
    assert!(!datasets.is_empty());
    for dataset in datasets {
        assert_eq!(
            dataset["algorithm"], "sha256",
            "a freshly fitted model must record cryptographic digests"
        );
        assert_eq!(dataset["digest"].as_str().unwrap().len(), 64);
    }
}

#[test]
fn two_runs_over_identical_inputs_produce_identical_provenance() {
    let first = screen_json(&study_001_portable_model(), &study_001_library(), &[]);
    let second = screen_json(&study_001_portable_model(), &study_001_library(), &[]);
    assert_eq!(
        first["provenance"], second["provenance"],
        "provenance must be a function of the inputs alone"
    );
    // No wall-clock field: a timestamp would make two identical runs differ.
    assert!(
        first["provenance"].get("timestamp").is_none()
            && first["provenance"].get("created_utc").is_none(),
        "provenance must not carry a clock: {:?}",
        first["provenance"]
    );
}

#[test]
fn a_changed_library_changes_the_recorded_digest() {
    let original = descriptor_library(&[("a", 8.0, 1.0), ("b", 9.0, 1.0)]);
    let altered = descriptor_library(&[("a", 8.0, 1.0), ("b", 9.000001, 1.0)]);

    let first = screen_json(&study_001_portable_model(), &original, &[]);
    let second = screen_json(&study_001_portable_model(), &altered, &[]);
    assert_ne!(
        first["provenance"]["library_sha256"], second["provenance"]["library_sha256"],
        "a different library must not hash the same"
    );
    // The model half is unchanged, so the difference is attributable.
    assert_eq!(
        first["provenance"]["model_sha256"],
        second["provenance"]["model_sha256"]
    );
}

#[test]
fn the_model_context_reaches_the_csv_export() {
    let output = run(&[
        "screen",
        study_001_portable_model().to_str().unwrap(),
        "--library",
        study_001_library().to_str().unwrap(),
        "--format",
        "csv",
    ]);
    assert_eq!(output.status, 0, "{}", output.stderr);
    let mut lines = output.stdout.lines();
    let header = lines.next().unwrap();
    for column in [
        "stericx_version",
        "model_id",
        "model_sha256",
        "library_sha256",
    ] {
        assert!(header.contains(column), "CSV lacks `{column}`: {header}");
    }
    // Every row is self-describing, so a filtered export stays traceable.
    let expected = sha256_hex_of(&study_001_portable_model());
    for line in lines {
        assert!(
            line.contains(&expected),
            "row lost its model digest: {line}"
        );
    }
}

#[test]
fn a_legacy_model_reports_the_context_it_has_without_inventing_the_rest() {
    let report = screen_json(&study_001_model(), &study_001_library(), &["--descending"]);
    let provenance = &report["provenance"];

    // The file itself can always be hashed, even when it records nothing.
    assert_eq!(
        provenance["model_sha256"].as_str().unwrap(),
        sha256_hex_of(&study_001_model())
    );
    assert_eq!(provenance["model_schema_version"], 1);
    // A schema-1 artifact records no identity or training digests; those must
    // read as absent rather than being filled in with a plausible value.
    assert!(provenance["model_id"].is_null());
    assert!(provenance["model_fitted_by"].is_null());
    assert!(
        provenance["training_datasets"]
            .as_array()
            .unwrap()
            .is_empty()
    );
    assert!(report["target"].is_null());
}
