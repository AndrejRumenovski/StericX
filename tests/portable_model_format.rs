//! Contract tests for the versioned portable model format.
//!
//! Three properties are pinned here: a portable document round-trips without
//! drift, an existing `schema_version: 1` artifact still loads, and documents
//! that are unsupported, incomplete, or numerically inconsistent are rejected
//! rather than scored.

use steric_x::model::{
    CreationMetadata, DatasetDigest, FeatureTransform, ModelFormatError, ModelProvenance,
    PORTABLE_SCHEMA_VERSION, PortableModel, ReactionProvenance, ResponseSpec, TrainingProvenance,
};
use steric_x::{
    FitOptions, PackedReactionRecord, ReactionLabel, ScientificFitReport, train_scientific_model,
};

fn fit_report() -> ScientificFitReport {
    let records = (0..14_usize)
        .map(|index| {
            let l = 1.0 + index as f32 * 0.2;
            PackedReactionRecord {
                l,
                b1: 1.7 + 0.03 * ((index * 5) % 4) as f32,
                b5: 3.2,
                nbo_charge: -0.35,
                ir_freq: 1_650.0,
                temp_k: 298.15,
                exp_ddg: 0.6 + 1.1 * l + 0.01 * ((index * 7) % 5) as f32,
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
                format!("P{index:02}"),
                split,
                format!("group_{}", index % 4),
            )
        })
        .collect::<Vec<_>>();
    let options = FitOptions {
        bootstrap_samples: 40,
        permutation_samples: 40,
        ..FitOptions::default()
    };
    train_scientific_model(&records, &labels, options)
        .expect("fixture trains")
        .report
}

fn provenance() -> ModelProvenance {
    ModelProvenance {
        model_id: "stericx-test-0001".into(),
        stericx_version: env!("CARGO_PKG_VERSION").into(),
        record_format: "sigpack_v1".into(),
        training: TrainingProvenance {
            record_count: 12,
            group_count: 4,
            dataset_digests: vec![DatasetDigest {
                artifact: "reactions.sigpack".into(),
                algorithm: "sha256".into(),
                digest: "4af1a776378a1f1aa369e076c9b35de2afa8c60e30841ec63cd13d95ccb8f00d".into(),
                byte_count: 896,
            }],
            fit_options: FitOptions {
                bootstrap_samples: 40,
                permutation_samples: 40,
                ..FitOptions::default()
            },
        },
        reaction: ReactionProvenance {
            reaction_family: Some("Ni-catalyzed homo-Diels-Alder".into()),
            catalyst_metal: Some("Ni".into()),
            ligand_class: Some("monodentate phosphine".into()),
            source_url: Some("https://github.com/SigmanGroup/Ni-Catalyzed-hDA".into()),
            notes: vec!["Synthetic fixture; not a scientific result.".into()],
        },
    }
}

fn portable_model() -> PortableModel {
    let mut model_provenance = provenance();
    let report = fit_report();
    model_provenance.training.record_count = report.training_count;
    model_provenance.training.group_count = report.training_group_count;
    PortableModel::from_fit_report(
        report,
        ResponseSpec::transition_state_energy_difference(Some(298.15)),
        model_provenance,
        CreationMetadata {
            created_utc: "2026-08-18T12:00:00Z".into(),
            produced_by: "portable_model_format test".into(),
        },
    )
    .expect("fixture builds a portable model")
}

/// Reparses a document as free-form JSON so tests can corrupt one field.
fn as_value(model: &PortableModel) -> serde_json::Value {
    serde_json::from_str(&model.to_json().unwrap()).unwrap()
}

fn expect_malformed(value: &serde_json::Value, expected_fragment: &str) {
    let text = serde_json::to_string_pretty(value).unwrap();
    match PortableModel::from_json(&text) {
        Err(ModelFormatError::Malformed(message)) => assert!(
            message.contains(expected_fragment),
            "expected a malformed-model error mentioning {expected_fragment:?}, got: {message}"
        ),
        Err(other) => panic!("expected Malformed, got {other}"),
        Ok(_) => panic!("expected rejection for a document missing {expected_fragment:?}"),
    }
}

#[test]
fn portable_document_round_trips_without_drift() {
    let model = portable_model();

    let first = model.to_json().unwrap();
    let reloaded = PortableModel::from_json(&first).unwrap();
    let second = reloaded.to_json().unwrap();

    assert_eq!(first, second, "round trip must be byte-stable");
    assert_eq!(reloaded.schema_version(), PORTABLE_SCHEMA_VERSION);
    assert!(reloaded.is_portable());
    assert_eq!(reloaded.fit.weights, model.fit.weights);
    assert_eq!(reloaded.inference().unwrap(), model.inference().unwrap());
    assert_eq!(reloaded.provenance().unwrap(), model.provenance().unwrap());
    assert_eq!(reloaded.created, model.created);
}

#[test]
fn serialization_is_deterministic_across_repeated_writes() {
    let model = portable_model();

    let renders = (0..5).map(|_| model.to_json().unwrap()).collect::<Vec<_>>();

    assert!(
        renders.windows(2).all(|pair| pair[0] == pair[1]),
        "repeated serialization must be byte-identical"
    );
}

#[test]
fn portable_document_is_a_superset_of_the_version_1_artifact() {
    let model = portable_model();
    let value = as_value(&model);
    let object = value.as_object().unwrap();

    // Every key of the legacy artifact is still present at the document root.
    for key in [
        "schema_version",
        "model",
        "training_count",
        "training_group_count",
        "feature_names",
        "selected_feature_indices",
        "selected_features",
        "weights",
        "standardized_means",
        "standardized_scales",
        "training",
        "fixed_feature_loo",
        "fixed_feature_group_loo",
        "ridge_baseline",
        "lasso_baseline",
        "coefficient_intervals",
        "response_permutation_p_value",
        "correlation_matrix",
        "variance_inflation_factors",
        "applicability_domain",
        "notes",
    ] {
        assert!(object.contains_key(key), "missing legacy key `{key}`");
    }
    for key in ["inference", "provenance", "created"] {
        assert!(object.contains_key(key), "missing portable section `{key}`");
    }
    assert_eq!(object["schema_version"], serde_json::json!(2));

    // A legacy reader deserializes the same document as a plain fit report.
    let text = serde_json::to_string(&value).unwrap();
    let legacy: ScientificFitReport = serde_json::from_str(&text).unwrap();
    assert_eq!(legacy.weights, model.fit.weights);
    assert_eq!(legacy.selected_features, model.fit.selected_features);
}

#[test]
fn checked_in_study_001_model_still_loads_as_legacy() {
    let text = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/docs/study_001/stericx_model.json"
    ))
    .expect("Study 001 model artifact is checked in");

    let model = PortableModel::from_json(&text).expect("version 1 artifact must still load");

    assert_eq!(model.schema_version(), 1);
    assert!(!model.is_portable());
    assert_eq!(model.fit.selected_features, ["B5_x_nbo_charge"]);
    assert_eq!(model.fit.training_count, 10);

    // Legacy metadata is reported absent, never invented.
    assert!(matches!(
        model.inference(),
        Err(ModelFormatError::MissingSection {
            section: "inference"
        })
    ));
    assert!(matches!(
        model.provenance(),
        Err(ModelFormatError::MissingSection {
            section: "provenance"
        })
    ));
    assert_eq!(
        model.missing_provenance(),
        [
            "reaction_family",
            "catalyst_metal",
            "ligand_class",
            "source_url"
        ]
    );

    // The weights still drive inference even without the portable sections.
    let predictor = model.predictor();
    let record = PackedReactionRecord {
        b5: 3.4,
        nbo_charge: -0.35,
        ..PackedReactionRecord::default()
    };
    assert!(predictor.predict(&record).is_finite());
}

#[test]
fn unknown_chemistry_context_is_reported_not_filled_in() {
    let mut model = portable_model();
    let provenance = model.provenance.as_mut().unwrap();
    provenance.reaction.catalyst_metal = None;
    provenance.reaction.source_url = None;

    let reloaded = PortableModel::from_json(&model.to_json().unwrap()).unwrap();

    assert_eq!(
        reloaded.missing_provenance(),
        ["catalyst_metal", "source_url"]
    );
    let value = as_value(&reloaded);
    assert!(
        value["provenance"]["reaction"]["catalyst_metal"].is_null(),
        "an unknown field must serialize as explicit null"
    );
}

#[test]
fn future_schema_versions_are_rejected() {
    let model = portable_model();
    let mut value = as_value(&model);
    value["schema_version"] = serde_json::json!(PORTABLE_SCHEMA_VERSION + 1);
    let text = serde_json::to_string(&value).unwrap();

    match PortableModel::from_json(&text) {
        Err(ModelFormatError::UnsupportedSchemaVersion { found, maximum }) => {
            assert_eq!(found, PORTABLE_SCHEMA_VERSION + 1);
            assert_eq!(maximum, PORTABLE_SCHEMA_VERSION);
        }
        other => panic!("expected UnsupportedSchemaVersion, got {other:?}"),
    }
}

#[test]
fn version_2_documents_must_carry_every_section() {
    let model = portable_model();
    for section in ["inference", "provenance", "created"] {
        let mut value = as_value(&model);
        value.as_object_mut().unwrap().remove(section);
        let text = serde_json::to_string(&value).unwrap();

        match PortableModel::from_json(&text) {
            Err(ModelFormatError::MissingSection { section: missing }) => {
                assert_eq!(missing, section)
            }
            other => panic!("removing `{section}` should be rejected, got {other:?}"),
        }
    }
}

#[test]
fn missing_required_fit_fields_are_rejected() {
    let model = portable_model();
    for field in [
        "weights",
        "selected_feature_indices",
        "applicability_domain",
    ] {
        let mut value = as_value(&model);
        value.as_object_mut().unwrap().remove(field);
        let text = serde_json::to_string(&value).unwrap();

        match PortableModel::from_json(&text) {
            Err(ModelFormatError::Json(error)) => {
                assert!(
                    error.to_string().contains(field),
                    "error for missing `{field}` should name it: {error}"
                );
            }
            other => panic!("removing `{field}` should be rejected, got {other:?}"),
        }
    }
}

#[test]
fn malformed_coefficients_are_rejected() {
    let model = portable_model();

    // A coefficient that disagrees with the flattened weight vector.
    let mut edited = as_value(&model);
    edited["inference"]["terms"][0]["coefficient"] = serde_json::json!(99.0);
    expect_malformed(&edited, "coefficient for");

    // An intercept that disagrees with weights[0].
    let mut edited = as_value(&model);
    edited["inference"]["intercept"] = serde_json::json!(-42.0);
    expect_malformed(&edited, "inference.intercept");

    // A non-finite weight cannot be scored with.
    let mut edited = as_value(&model);
    edited["weights"][0] = serde_json::json!(f64::INFINITY);
    let text = serde_json::to_string(&edited).unwrap();
    // JSON has no infinity literal, so serde_json renders it as null and the
    // typed field rejects it before validation is even reached.
    assert!(PortableModel::from_json(&text).is_err());

    // A term pointing at a column the model never selected.
    let mut edited = as_value(&model);
    edited["inference"]["terms"][0]["feature_index"] = serde_json::json!(3);
    expect_malformed(&edited, "refers to column 3");

    // A term whose name contradicts its column.
    let mut edited = as_value(&model);
    edited["inference"]["terms"][0]["feature_name"] = serde_json::json!("B1_boltz");
    expect_malformed(&edited, "is named B1_boltz");
}

#[test]
fn malformed_scaling_information_is_rejected() {
    let model = portable_model();

    // A zero standard deviation would divide by zero during inference.
    let mut edited = as_value(&model);
    let column = model.fit.selected_feature_indices[0];
    edited["standardized_scales"][column] = serde_json::json!(0.0);
    edited["inference"]["terms"][0]["training_standard_deviation"] = serde_json::json!(0.0);
    expect_malformed(&edited, "not a positive finite number");

    // A negative standard deviation is equally unusable.
    let mut edited = as_value(&model);
    edited["standardized_scales"][column] = serde_json::json!(-2.0);
    edited["inference"]["terms"][0]["training_standard_deviation"] = serde_json::json!(-2.0);
    expect_malformed(&edited, "not a positive finite number");

    // Scaling in the inference block that drifts from the fit report.
    let mut edited = as_value(&model);
    edited["inference"]["terms"][0]["training_mean"] = serde_json::json!(1_234.5);
    expect_malformed(&edited, "training mean for");

    // An inverted applicability range.
    let mut edited = as_value(&model);
    let minimum = edited["applicability_domain"][0]["minimum"].clone();
    edited["applicability_domain"][0]["minimum"] =
        edited["applicability_domain"][0]["maximum"].clone();
    edited["applicability_domain"][0]["maximum"] = minimum;
    expect_malformed(&edited, "inverted");
}

#[test]
fn structurally_inconsistent_documents_are_rejected() {
    let model = portable_model();

    // Provenance that disagrees with the fit about the training set size.
    let mut edited = as_value(&model);
    edited["provenance"]["training"]["record_count"] = serde_json::json!(999);
    expect_malformed(&edited, "999 training rows");

    // A portable model must identify its training data.
    let mut edited = as_value(&model);
    edited["provenance"]["training"]["dataset_digests"] = serde_json::json!([]);
    expect_malformed(&edited, "must identify");

    // A digest that is not hexadecimal is not a digest.
    let mut edited = as_value(&model);
    edited["provenance"]["training"]["dataset_digests"][0]["digest"] =
        serde_json::json!("not-a-digest");
    expect_malformed(&edited, "not hexadecimal");

    // An empty model identifier makes the artifact untraceable.
    let mut edited = as_value(&model);
    edited["provenance"]["model_id"] = serde_json::json!("   ");
    expect_malformed(&edited, "model_id is empty");

    // A selected descriptor with no matching name.
    let mut edited = as_value(&model);
    edited["selected_features"][0] = serde_json::json!("ir_frequency");
    expect_malformed(&edited, "does not match column");
}

#[test]
fn feature_space_describes_every_interaction_term() {
    let model = portable_model();
    let space = &model.inference().unwrap().feature_space;

    assert_eq!(space.feature_names.len(), 8);
    assert_eq!(space.transformations.len(), 8);
    assert!(matches!(
        space.transformations[0],
        FeatureTransform::Constant
    ));
    assert_eq!(
        space.transformations[6],
        FeatureTransform::Interaction {
            factors: vec!["sterimol_b5".into(), "nbo_charge".into()],
        },
        "the B5 x NBO interaction must be reconstructable without the fitting code"
    );

    // A feature space that contradicts the fit report is rejected.
    let mut edited = as_value(&model);
    edited["inference"]["feature_space"]["feature_names"][1] = serde_json::json!("renamed");
    expect_malformed(&edited, "disagrees with feature_names");
}

#[test]
fn inference_reproduces_predictions_without_the_training_process() {
    let model = portable_model();
    let reloaded = PortableModel::from_json(&model.to_json().unwrap()).unwrap();
    let record = PackedReactionRecord {
        l: 2.4,
        b1: 1.7,
        b5: 3.2,
        nbo_charge: -0.35,
        ir_freq: 1_650.0,
        ..PackedReactionRecord::default()
    };

    // Score straight from the serialized document, using only the published
    // intercept, coefficients, and the declared feature construction.
    let inference = reloaded.inference().unwrap();
    let descriptors = |name: &str| -> f64 {
        match name {
            "sterimol_l" => f64::from(record.l),
            "sterimol_b1" => f64::from(record.b1),
            "sterimol_b5" => f64::from(record.b5),
            "nbo_charge" => f64::from(record.nbo_charge),
            "ir_frequency" => f64::from(record.ir_freq),
            other => panic!("unknown descriptor {other}"),
        }
    };
    let mut manual = inference.intercept;
    for term in &inference.terms {
        let value = match &inference.feature_space.transformations[term.feature_index] {
            FeatureTransform::Constant => 1.0,
            FeatureTransform::Descriptor { descriptor } => descriptors(descriptor),
            FeatureTransform::Interaction { factors } => {
                factors.iter().map(|factor| descriptors(factor)).product()
            }
        };
        manual += term.coefficient * value;
    }

    let engine = f64::from(reloaded.predictor().predict(&record));
    assert!(
        (manual - engine).abs() < 1.0e-5,
        "document-only inference {manual} disagreed with the engine {engine}"
    );
}

#[test]
fn json_floats_survive_a_round_trip_bit_for_bit() {
    // Guards `serde_json`'s `float_roundtrip` feature. Without it the default
    // parser shifts some seventeen-digit values by one unit in the last place,
    // which would silently rewrite published validation statistics every time a
    // model was read and written again.
    let original = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/docs/study_001/stericx_model.json"
    ))
    .expect("Study 001 model artifact is checked in");

    let parsed: ScientificFitReport = serde_json::from_str(&original).unwrap();
    let rewritten = serde_json::to_string_pretty(&parsed).unwrap();

    // Every value the checked-in document carries must come back unchanged.
    // The rewritten document may gain keys the struct has since grown, so this
    // compares what the original states rather than requiring byte equality.
    let before: serde_json::Value = serde_json::from_str(&original).unwrap();
    let after: serde_json::Value = serde_json::from_str(&rewritten).unwrap();
    assert_values_preserved(&before, &after, "$");

    for value in [
        0.9996584928554385_f64,
        0.011934137734395259,
        1.0902477383613587,
        0.36254551871085705,
    ] {
        let text = serde_json::to_string(&value).unwrap();
        let back: f64 = serde_json::from_str(&text).unwrap();
        assert_eq!(
            back.to_bits(),
            value.to_bits(),
            "{value} did not survive a JSON round trip"
        );
    }
}

/// Asserts every value reachable in `before` is present and equal in `after`.
///
/// Keys that only `after` carries are ignored, so adding an optional field to
/// the report does not fail the check, but changing a recorded number does.
fn assert_values_preserved(before: &serde_json::Value, after: &serde_json::Value, path: &str) {
    match (before, after) {
        (serde_json::Value::Object(before), serde_json::Value::Object(after)) => {
            for (key, value) in before {
                let found = after
                    .get(key)
                    .unwrap_or_else(|| panic!("{path}.{key} disappeared on re-serialization"));
                assert_values_preserved(value, found, &format!("{path}.{key}"));
            }
        }
        (serde_json::Value::Array(before), serde_json::Value::Array(after)) => {
            assert_eq!(before.len(), after.len(), "{path} changed length");
            for (index, (before, after)) in before.iter().zip(after).enumerate() {
                assert_values_preserved(before, after, &format!("{path}[{index}]"));
            }
        }
        (before, after) => assert_eq!(before, after, "{path} changed"),
    }
}
