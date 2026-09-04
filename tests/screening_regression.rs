//! Regression suite pinning `stericx screen` against the published studies.
//!
//! Every assertion here anchors to an artifact that is committed to this
//! repository, so the suite runs in CI with no network access and no
//! copyrighted supporting information. Where a study's inputs are *not*
//! redistributable the limitation is encoded as an explicit test rather than
//! quietly skipped — see `docs/validation/SCREENING_REGRESSION.md` for exactly
//! what is and is not covered.
//!
//! Study 001 (Ni-catalysed homo-Diels-Alder) is the only published study whose
//! fitted model can be driven through `screen`: it is a linear model over the
//! physical-organic feature space, and its model, library, frozen predictions
//! and evaluation are all committed. It therefore carries the end-to-end cases.
//!
//! Studies 007 and 009 are %Vbur(min) threshold classifiers, not fitted linear
//! models, and their per-ligand experimental data lives in copyrighted AAAS
//! supporting information the repository deliberately does not redistribute.
//! What is committed — StericX's own aggregate metrics and its own descriptor
//! database — is pinned here at the level it genuinely supports.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};
use steric_x::{PackedReactionRecord, RegressXPredictor};

const EXE: &str = env!("CARGO_BIN_EXE_stericx");

/// `stericx evaluate` accepts a frozen prediction within 1e-4 of the model's
/// own output. Screening is held to the same published tolerance: the screen
/// loop accumulates in f64 while the engine kernel is f32, so the two agree to
/// about 1e-6 in practice and 1e-4 is the documented contract.
const INFERENCE_TOLERANCE: f64 = 1.0e-4;

fn repo(path: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(path)
}

fn temp_path(extension: &str) -> PathBuf {
    static COUNTER: AtomicUsize = AtomicUsize::new(0);
    std::env::temp_dir().join(format!(
        "stericx_screen_regression_{}_{}.{extension}",
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

fn json_file(path: &Path) -> serde_json::Value {
    serde_json::from_str(&std::fs::read_to_string(path).expect("artifact is committed"))
        .expect("artifact is valid JSON")
}

// --- Study 001 committed artifacts ------------------------------------------

fn study_001_portable() -> PathBuf {
    repo("docs/study_001/stericx_portable_model.json")
}

fn study_001_legacy() -> PathBuf {
    repo("docs/study_001/stericx_model.json")
}

fn study_001_library() -> PathBuf {
    repo("data/reactions_raw.csv")
}

fn screen(model: &Path, extra: &[&str]) -> serde_json::Value {
    let library = study_001_library();
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
    serde_json::from_str(&output.stdout).expect("screen emits JSON")
}

fn hit<'a>(report: &'a serde_json::Value, ligand: &str) -> &'a serde_json::Value {
    report["hits"]
        .as_array()
        .expect("hits is an array")
        .iter()
        .find(|hit| hit["ligand"] == ligand)
        .unwrap_or_else(|| panic!("{ligand} was screened"))
}

/// The one frozen, published prediction: `stericx evaluate`'s blind holdout.
fn frozen_blind_prediction() -> (String, f64, String) {
    let text = std::fs::read_to_string(repo("docs/study_001/stericx_frozen_predictions.csv"))
        .expect("frozen predictions are committed");
    let row = text.lines().nth(1).expect("one frozen row");
    let fields = row.split(',').collect::<Vec<_>>();
    (
        fields[0].to_owned(),
        fields[3].parse().expect("predicted ddG parses"),
        fields[4].to_owned(),
    )
}

// ============================================================================
// 1. Model loading
// ============================================================================

#[test]
fn study_001_models_load_in_both_published_schema_versions() {
    // The repository publishes the same fit as a schema-1 artifact and a
    // schema-2 portable document. Both must remain loadable: the legacy one is
    // what the released studies and the Python drivers read.
    for (path, expected_schema) in [(study_001_legacy(), 1), (study_001_portable(), 2)] {
        let output = run(&["model", "inspect", path.to_str().unwrap()]);
        assert_eq!(output.status, 0, "{}: {}", path.display(), output.stderr);
        assert!(
            output
                .stdout
                .contains(&format!("schema_version={expected_schema}")),
            "{} is no longer schema {expected_schema}:\n{}",
            path.display(),
            output.stdout
        );
    }

    // The portable document must also validate clean, since it is the one a
    // third party would be handed.
    let output = run(&["model", "validate", study_001_portable().to_str().unwrap()]);
    assert_eq!(output.status, 0, "{}", output.stderr);
    assert!(output.stdout.contains("errors=0"), "{}", output.stdout);
    assert!(output.stdout.contains("warnings=0"), "{}", output.stdout);
}

#[test]
fn the_published_model_still_selects_the_descriptor_the_study_reports() {
    // Study 001's whole narrative rests on one selected term. If selection ever
    // moves, every downstream number in the study is a different claim.
    let model = json_file(&study_001_portable());
    assert_eq!(
        model["selected_features"],
        serde_json::json!(["B5_x_nbo_charge"]),
        "the published model no longer selects the descriptor STUDY_001 reports"
    );
    assert_eq!(model["training_count"], 10);
    assert_eq!(model["training_group_count"], 9);
}

// ============================================================================
// 2. Descriptor mapping
// ============================================================================

#[test]
fn the_model_drives_which_library_columns_are_required() {
    let report = screen(&study_001_portable(), &[]);
    let required = report["required_inputs"]
        .as_array()
        .unwrap()
        .iter()
        .map(|value| value.as_str().unwrap().to_owned())
        .collect::<Vec<_>>();
    // B5_x_nbo_charge is a product term, so the mapping must demand both
    // factors and neither more nor less.
    assert_eq!(required, ["sterimol_b5", "nbo_charge"]);

    // Each hit reports the descriptors its prediction actually consumed.
    let blind = hit(&report, "SIG-NIHDA-723");
    let names = blind["descriptors"]
        .as_array()
        .unwrap()
        .iter()
        .map(|d| d["name"].as_str().unwrap().to_owned())
        .collect::<Vec<_>>();
    assert_eq!(names, ["sterimol_b5", "nbo_charge"]);
}

#[test]
fn a_library_missing_a_required_descriptor_is_refused_not_guessed() {
    // The shipped Kraken database is geometry-derived and carries no donor
    // electronics, so it cannot satisfy this model. Refusing is the documented
    // behaviour; silently substituting a default would fabricate chemistry.
    let output = run(&[
        "screen",
        study_001_portable().to_str().unwrap(),
        "--library",
        repo("data/ligand_db/kraken_phosphines.csv")
            .to_str()
            .unwrap(),
    ]);
    assert_ne!(
        output.status, 0,
        "a library without nbo_charge must be refused"
    );
    assert!(
        output.stderr.contains("nbo_charge"),
        "the refusal must name the missing input: {}",
        output.stderr
    );
}

// ============================================================================
// 3. Feature scaling
// ============================================================================

#[test]
fn standardization_matches_the_constants_the_published_model_records() {
    // Recompute the standardized coordinate independently from the model's own
    // recorded mean and scale, and require screening to agree. This is what
    // makes the applicability distances comparable across ligands.
    let model = json_file(&study_001_portable());
    let geometry = &model["training_geometry"];
    let mean = geometry["means"][0].as_f64().unwrap();
    let scale = geometry["scales"][0].as_f64().unwrap();

    let report = screen(&study_001_portable(), &[]);
    let blind = hit(&report, "SIG-NIHDA-723");
    let descriptors = blind["descriptors"].as_array().unwrap();
    let b5 = descriptors[0]["value"].as_f64().unwrap();
    let nbo = descriptors[1]["value"].as_f64().unwrap();
    let standardized = (b5 * nbo - mean) / scale;

    // The nearest-training distance is measured in exactly that coordinate, so
    // it bounds how far the standardized point can be from a training point.
    let nearest = blind["nearest_training_distance"].as_f64().unwrap();
    let points = geometry["standardized_training_points"].as_array().unwrap();
    let recomputed = points
        .iter()
        .map(|row| (row[0].as_f64().unwrap() - standardized).abs())
        .fold(f64::INFINITY, f64::min);
    // The reported descriptor values are f32, so recomputing from them carries
    // single-precision rounding; 1e-6 is far tighter than any real drift.
    assert!(
        (nearest - recomputed).abs() < 1.0e-6,
        "standardization drifted: reported {nearest}, recomputed {recomputed}"
    );

    // And the recorded scale is the training standard deviation, not 1.
    assert!(scale > 0.0 && (scale - 1.0).abs() > 1.0e-6);
}

// ============================================================================
// 4. Prediction, and cross-pathway consistency
// ============================================================================

#[test]
fn screening_reproduces_the_engine_kernel_for_every_published_ligand() {
    // The generic screen pathway against the study-specific inference path.
    // `RegressXPredictor` is the kernel `fit`, `evaluate` and the study driver
    // all go through, so agreement here is agreement with the published route.
    let model = json_file(&study_001_portable());
    let mut weights = [0.0_f32; 8];
    for (index, value) in model["weights"].as_array().unwrap().iter().enumerate() {
        weights[index] = value.as_f64().unwrap() as f32;
    }
    let predictor = RegressXPredictor::new(weights);

    let report = screen(&study_001_portable(), &[]);
    let hits = report["hits"].as_array().unwrap();
    assert_eq!(hits.len(), 11, "the published library holds 11 ligands");

    for entry in hits {
        let descriptors = entry["descriptors"].as_array().unwrap();
        let b5 = descriptors[0]["value"].as_f64().unwrap() as f32;
        let nbo = descriptors[1]["value"].as_f64().unwrap() as f32;
        let record = PackedReactionRecord {
            b5,
            nbo_charge: nbo,
            ..PackedReactionRecord::default()
        };
        let expected = f64::from(predictor.predict(&record));
        let actual = entry["predicted_ddg_kcal_mol"].as_f64().unwrap();
        assert!(
            (actual - expected).abs() < INFERENCE_TOLERANCE,
            "{}: screen {actual} vs engine {expected}",
            entry["ligand"]
        );
    }
}

#[test]
fn the_frozen_blind_prediction_is_reproduced_within_the_published_tolerance() {
    // The single pre-registered, frozen prediction in the repository. If this
    // moves, a published scientific claim has moved.
    let (ligand, frozen, domain) = frozen_blind_prediction();
    let report = screen(&study_001_portable(), &[]);
    let screened = hit(&report, &ligand)["predicted_ddg_kcal_mol"]
        .as_f64()
        .unwrap();
    assert!(
        (screened - frozen).abs() < INFERENCE_TOLERANCE,
        "frozen prediction drifted: screen {screened} vs frozen {frozen}"
    );

    // The frozen artifact also records an applicability verdict, which the
    // screen pathway must not contradict.
    assert_eq!(domain, "inside_training_range");
    assert_eq!(
        hit(&report, &ligand)["applicability"],
        "inside_training_range",
        "screen disagrees with the frozen applicability verdict"
    );
}

#[test]
fn the_published_evaluation_error_still_follows_from_the_frozen_prediction() {
    // `stericx evaluate` published an MAE for the blind holdout. It is the gap
    // between the frozen prediction and the experimental value, so recomputing
    // it from committed artifacts guards the whole evaluate pathway.
    let evaluation = json_file(&repo("docs/study_001/stericx_evaluation.json"));
    let scored = &evaluation["scored_predictions"][0];
    let predicted = scored["predicted_ddg_kcal_mol"].as_f64().unwrap();
    let experimental = scored["experimental_ddg_kcal_mol"].as_f64().unwrap();
    let mae = evaluation["mae_kcal_mol"].as_f64().unwrap();
    assert!(
        ((predicted - experimental).abs() - mae).abs() < 1.0e-6,
        "the published MAE no longer follows from the published prediction"
    );

    let (_, frozen, _) = frozen_blind_prediction();
    assert!(
        (predicted - frozen).abs() < 1.0e-6,
        "evaluate and the frozen CSV disagree"
    );
}

// ============================================================================
// 5. Ranking
// ============================================================================

#[test]
fn ranking_follows_the_direction_the_published_model_records() {
    let report = screen(&study_001_portable(), &[]);
    assert_eq!(report["model_optimization"], "maximize");
    assert_eq!(report["ranking_order"], "descending");
    assert_eq!(report["ranking_overridden"], false);

    let values = report["hits"]
        .as_array()
        .unwrap()
        .iter()
        .map(|hit| hit["predicted_ddg_kcal_mol"].as_f64().unwrap())
        .collect::<Vec<_>>();
    let mut sorted = values.clone();
    sorted.sort_by(|a, b| b.total_cmp(a));
    assert_eq!(values, sorted, "hits are not in the model's stated order");

    // The best-predicted ligand in the published library, pinned.
    assert_eq!(report["hits"][0]["ligand"], "SIG-NIHDA-723");
}

#[test]
fn the_legacy_artifact_still_refuses_to_rank_without_a_direction() {
    // Deliberate behaviour introduced in v0.3: a schema-1 model records no
    // optimization direction and must not be ranked on a guess.
    let output = run(&[
        "screen",
        study_001_legacy().to_str().unwrap(),
        "--library",
        study_001_library().to_str().unwrap(),
    ]);
    assert_ne!(output.status, 0, "a directionless model must not be ranked");

    // With an explicit direction it ranks, and agrees with the portable model.
    let legacy = screen(&study_001_legacy(), &["--descending"]);
    let portable = screen(&study_001_portable(), &[]);
    let ligands = |report: &serde_json::Value| {
        report["hits"]
            .as_array()
            .unwrap()
            .iter()
            .map(|hit| hit["ligand"].as_str().unwrap().to_owned())
            .collect::<Vec<_>>()
    };
    assert_eq!(
        ligands(&legacy),
        ligands(&portable),
        "the two published artifacts rank the same library differently"
    );
}

// ============================================================================
// 6. Applicability scoring
// ============================================================================

#[test]
fn applicability_agrees_with_the_published_training_domain() {
    let model = json_file(&study_001_portable());
    let domain = &model["applicability_domain"][0];
    let minimum = domain["minimum"].as_f64().unwrap();
    let maximum = domain["maximum"].as_f64().unwrap();

    let report = screen(&study_001_portable(), &[]);
    for entry in report["hits"].as_array().unwrap() {
        let descriptors = entry["descriptors"].as_array().unwrap();
        let product =
            descriptors[0]["value"].as_f64().unwrap() * descriptors[1]["value"].as_f64().unwrap();
        let inside = product >= minimum && product <= maximum;
        let verdict = entry["domain_verdict"].as_str().unwrap();
        assert_eq!(
            verdict != "extrapolation",
            inside,
            "{}: descriptor {product} vs published range [{minimum}, {maximum}] but verdict {verdict}",
            entry["ligand"]
        );
    }

    // The published library is entirely inside the training range, which is
    // what the frozen artifact's `inside_training_range` verdict asserts.
    assert_eq!(report["domain_summary"].as_array().unwrap().len(), 1);
    assert_eq!(report["domain_summary"][0][0], "interpolation");
}

#[test]
fn the_neighbour_calibration_is_derived_from_the_published_training_points() {
    // The one applicability threshold is measured from the training set, so it
    // must remain reproducible from the training points the model carries.
    let model = json_file(&study_001_portable());
    let geometry = &model["training_geometry"];
    let points = geometry["standardized_training_points"]
        .as_array()
        .unwrap()
        .iter()
        .map(|row| row[0].as_f64().unwrap())
        .collect::<Vec<_>>();

    let mut spacing = Vec::new();
    for (index, point) in points.iter().enumerate() {
        let nearest = points
            .iter()
            .enumerate()
            .filter(|(other, _)| *other != index)
            .map(|(_, other)| (point - other).abs())
            .fold(f64::INFINITY, f64::min);
        spacing.push(nearest);
    }
    let recomputed = spacing.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let recorded = geometry["neighbor_calibration"]["threshold"]
        .as_f64()
        .unwrap();
    assert!(
        (recomputed - recorded).abs() < 1.0e-9,
        "the recorded neighbour threshold is not the training spacing: \
         recorded {recorded}, recomputed {recomputed}"
    );
}

// ============================================================================
// 7. Uncertainty serialization
// ============================================================================

#[test]
fn the_published_model_carries_a_usable_bootstrap_ensemble() {
    let model = json_file(&study_001_portable());
    let ensemble = &model["uncertainty"];
    assert!(
        ensemble.is_object(),
        "the portable model carries no ensemble"
    );

    let replicates = ensemble["replicates"].as_array().unwrap();
    assert_eq!(
        ensemble["replicate_count"].as_u64().unwrap() as usize,
        replicates.len()
    );
    assert_eq!(ensemble["columns"][0], "intercept");
    let columns = ensemble["columns"].as_array().unwrap().len();
    for replicate in replicates {
        assert_eq!(replicate.as_array().unwrap().len(), columns);
    }

    // Screening consumes it without refitting, and every interval brackets its
    // own point estimate.
    let report = screen(&study_001_portable(), &[]);
    assert_eq!(
        report["uncertainty_replicates"].as_u64().unwrap() as usize,
        replicates.len()
    );
    for entry in report["hits"].as_array().unwrap() {
        let uncertainty = &entry["uncertainty"];
        let low = uncertainty["lower"].as_f64().unwrap();
        let high = uncertainty["upper"].as_f64().unwrap();
        let central = entry["predicted_ddg_kcal_mol"].as_f64().unwrap();
        assert!(low <= central && central <= high, "{}", entry["ligand"]);
        // The name must keep disclaiming what it does not cover.
        assert_eq!(uncertainty["method"], "percentile_bootstrap_mean_response");
    }
}

#[test]
fn the_legacy_artifact_reports_no_interval_rather_than_inventing_one() {
    let report = screen(&study_001_legacy(), &["--descending"]);
    assert!(report["uncertainty_method"].is_null());
    for entry in report["hits"].as_array().unwrap() {
        assert!(entry["uncertainty"].is_null(), "{}", entry["ligand"]);
    }
}

// ============================================================================
// 8. Candidate filtering
// ============================================================================

#[test]
fn filtering_accounts_for_every_published_candidate() {
    let baseline = screen(&study_001_portable(), &[]);
    let library_size = baseline["library_size"].as_u64().unwrap();
    assert_eq!(library_size, 11);

    // --top selects after inference: everything is still screened.
    let topped = screen(&study_001_portable(), &["--top", "3"]);
    assert_eq!(topped["screened"], baseline["screened"]);
    assert_eq!(topped["returned"], 3);

    // --exclude-tested removes candidates before inference and reconciles.
    let tested = temp_path("csv");
    std::fs::write(
        &tested,
        "ligand\nSIG-NIHDA-723\nSIG-NIHDA-1058\nNOT-IN-LIBRARY\n",
    )
    .unwrap();
    let filtered = screen(
        &study_001_portable(),
        &["--exclude-tested", tested.to_str().unwrap()],
    );
    let exclusion = &filtered["exclusion"];
    assert_eq!(exclusion["excluded"], 2);
    assert_eq!(exclusion["remaining"].as_u64().unwrap() + 2, library_size);
    assert_eq!(
        exclusion["unresolved"].as_array().unwrap(),
        &vec![serde_json::json!("NOT-IN-LIBRARY")],
        "an identifier matching nothing must be reported, not dropped"
    );
    assert_eq!(filtered["library_size"], baseline["library_size"]);
}

// ============================================================================
// 9. Deterministic output
// ============================================================================

#[test]
fn every_output_format_is_byte_identical_across_runs() {
    let model = study_001_portable();
    let library = study_001_library();
    for extra in [
        vec!["--format", "json"],
        vec!["--format", "csv"],
        vec!["--format", "text"],
        vec!["--format", "json", "--diverse", "--top", "5"],
    ] {
        let mut args = vec![
            "screen",
            model.to_str().unwrap(),
            "--library",
            library.to_str().unwrap(),
        ];
        args.extend_from_slice(&extra);
        let first = run(&args);
        let second = run(&args);
        assert_eq!(first.status, 0, "{}", first.stderr);
        assert_eq!(
            first.stdout, second.stdout,
            "repeated runs differ for {extra:?}"
        );
    }
}

#[test]
fn an_exported_deck_is_reproducible_from_the_published_model() {
    let deck = temp_path("csv");
    let model = study_001_portable();
    let library = study_001_library();
    let args = [
        "screen",
        model.to_str().unwrap(),
        "--library",
        library.to_str().unwrap(),
        "--top",
        "5",
        "--export-deck",
        deck.to_str().unwrap(),
    ];
    let first = run(&args);
    assert_eq!(first.status, 0, "{}", first.stderr);
    let deck_a = std::fs::read_to_string(&deck).unwrap();
    let sidecar = deck.with_file_name(format!(
        "{}.meta.json",
        deck.file_stem().unwrap().to_string_lossy()
    ));
    let meta_a: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&sidecar).unwrap()).unwrap();

    let second = run(&args);
    assert_eq!(second.status, 0, "{}", second.stderr);
    let deck_b = std::fs::read_to_string(&deck).unwrap();
    let meta_b: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&sidecar).unwrap()).unwrap();

    assert_eq!(deck_a, deck_b, "the deck is not reproducible");
    assert_eq!(meta_a, meta_b, "the deck sidecar is not reproducible");

    // A deck must never carry a blinded experimental response, even though the
    // source library file holds one.
    let lowered = deck_a.to_lowercase();
    assert!(!lowered.contains("exp_ddg") && !lowered.contains("experimental"));
    let library = std::fs::read_to_string(study_001_library()).unwrap();
    assert!(
        library.contains("Exp_ddG_kcal_mol"),
        "fixture assumes the library carries the response column"
    );
}

// ============================================================================
// Studies 007 and 009 — what the committed artifacts genuinely support
// ============================================================================

/// Study 007 and 009 are `%Vbur(min)` threshold classifiers. The screening
/// feature space is the physical-organic set — Sterimol terms, donor
/// electronics and their products — and carries no buried-volume descriptor,
/// so these classifiers cannot be expressed as `screen` models at all.
///
/// This is asserted rather than assumed. If a buried-volume term is ever added
/// to the feature space this test fails, which is the moment to wire Studies
/// 007 and 009 into the end-to-end cases above instead of leaving them here.
#[test]
fn the_screening_feature_space_still_cannot_express_a_buried_volume_classifier() {
    let names = steric_x::model::MODEL_FEATURE_NAMES;
    for name in names {
        let lowered = name.to_lowercase();
        assert!(
            !lowered.contains("vbur") && !lowered.contains("buried"),
            "the feature space gained `{name}`: Studies 007 and 009 can now be \
             driven through `screen` and should be promoted out of this section"
        );
    }
    assert_eq!(names.len(), 8, "the published feature space changed size");
}

#[test]
fn study_007_published_metrics_are_unchanged() {
    // StericX's own committed outputs. Pinned so a descriptor-kernel change
    // that would move them cannot land without the study being re-run and the
    // numbers deliberately updated.
    let metrics = json_file(&repo("docs/study_007/crosscoupling_metrics.json"));
    assert_eq!(metrics["ligands"], 103);
    assert_eq!(metrics["descriptor_fidelity"]["n"], 479);
    let r2 = metrics["descriptor_fidelity"]["r2"].as_f64().unwrap();
    assert!(
        (r2 - 0.999_174_920_361_292_3).abs() < 1.0e-12,
        "descriptor fidelity R² moved: {r2}"
    );
    let pooled = &metrics["transferability"]["pooled"];
    assert_eq!(pooled["n"], 479);
    assert!((pooled["mcc"].as_f64().unwrap() - 0.491_774_054_911_526_15).abs() < 1.0e-12);

    // Every reaction reports its own threshold and a paper comparison.
    for key in ["I", "II", "III", "IV", "V", "RS1"] {
        let reaction = &metrics["reactions"][key];
        assert!(reaction["stericx"]["threshold"].as_f64().is_some(), "{key}");
        assert!(reaction["paper"]["thr"].as_f64().is_some(), "{key}");
    }
}

#[test]
fn study_009_published_metrics_are_unchanged() {
    let metrics = json_file(&repo("docs/study_009/pd_crosscoupling_metrics.json"));
    assert_eq!(metrics["ligands"], 132);
    assert_eq!(metrics["descriptor_fidelity"]["n"], 267);
    let r2 = metrics["descriptor_fidelity"]["r2"].as_f64().unwrap();
    assert!(
        (r2 - 0.999_360_896_081_756_4).abs() < 1.0e-12,
        "descriptor fidelity R² moved: {r2}"
    );
    assert_eq!(metrics["n_right_direction"], 6);
    for key in ["VII", "VIII", "IX", "X", "XI", "XII"] {
        assert!(
            metrics["reactions"][key]["stericx"]["threshold"]
                .as_f64()
                .is_some(),
            "{key}"
        );
    }
}

#[test]
fn the_study_thresholds_remain_inside_the_committed_descriptor_distribution() {
    // The classifier thresholds are `%Vbur(min)` values — the per-ligand minimum
    // of `percent_buried_volume` across conformers, which is what
    // studies/study_007_crosscoupling.py thresholds on. The descriptor database
    // is StericX's own committed output for all 1,541 Kraken phosphines, so a
    // kernel change that shifted buried volume would move the distribution out
    // from under every published threshold.
    //
    // This is a distribution-level guard, not a reproduction of the classifier:
    // reproducing it needs the per-ligand experimental yields, which live in
    // copyrighted supporting information this repository does not redistribute.
    let table = std::fs::read_to_string(repo("data/ligand_db/kraken_phosphines.csv"))
        .expect("the descriptor database is committed");
    let mut lines = table.lines();
    let header = lines.next().unwrap().split(',').collect::<Vec<_>>();
    let column = header
        .iter()
        .position(|name| *name == "percent_buried_volume")
        .expect("the database carries percent_buried_volume");

    let values = lines
        .filter_map(|line| line.split(',').nth(column)?.parse::<f64>().ok())
        .collect::<Vec<_>>();
    assert_eq!(values.len(), 1541, "the published library size changed");
    let low = values.iter().copied().fold(f64::INFINITY, f64::min);
    let high = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);

    for (study, path) in [
        ("007", "docs/study_007/crosscoupling_metrics.json"),
        ("009", "docs/study_009/pd_crosscoupling_metrics.json"),
    ] {
        let metrics = json_file(&repo(path));
        for (name, reaction) in metrics["reactions"].as_object().unwrap() {
            let threshold = reaction["stericx"]["threshold"].as_f64().unwrap();
            assert!(
                threshold > low && threshold < high,
                "study {study} reaction {name}: threshold {threshold} is outside the \
                 committed percent_buried_volume range [{low}, {high}]"
            );
        }
    }
}
