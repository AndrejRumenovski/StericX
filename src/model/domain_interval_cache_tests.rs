//! Exact compatibility checks against the unchanged method and pre-edit C1b bits.

use super::*;

#[path = "../../tests/fixtures/c2_interval_cases.rs"]
mod frozen_cases;

fn bits(interval: Option<(f64, f64)>) -> Option<[u64; 2]> {
    interval.map(|(low, high)| [low.to_bits(), high.to_bits()])
}

#[test]
fn cached_intervals_match_frozen_c1b_endpoints_and_uncached_method() {
    let expected = include_str!("../../tests/fixtures/c2_interval_bits.jsonl")
        .lines()
        .map(|line| serde_json::from_str::<serde_json::Value>(line).unwrap())
        .collect::<Vec<_>>();
    let cases = frozen_cases::cases();
    assert_eq!(cases.len(), expected.len());
    assert_eq!(cases.len(), 46);
    for (case, expected) in cases.into_iter().zip(expected) {
        assert_eq!(case.name, expected["case"].as_str().unwrap());
        let captured: Option<[u64; 2]> =
            serde_json::from_value(expected["interval_bits"].clone()).unwrap();
        let evaluator = PredictionIntervalEvaluator::new(&case.geometry);
        assert!(evaluator.multiplier.get().is_none());
        for _ in 0..3 {
            let uncached = case
                .geometry
                .prediction_interval(case.prediction, &case.expanded);
            let cached = evaluator.prediction_interval(case.prediction, &case.expanded);
            assert_eq!(bits(uncached), captured, "uncached {}", case.name);
            assert_eq!(bits(cached), captured, "cached {}", case.name);
        }
        assert_eq!(
            evaluator.multiplier.get().is_some(),
            case.geometry.leverage(&case.expanded).is_some(),
            "initialization boundary: {}",
            case.name
        );
    }
}

#[test]
fn invalid_leverage_is_lazy_and_later_valid_candidate_initializes() {
    let geometry = frozen_cases::geometry();
    let evaluator = PredictionIntervalEvaluator::new(&geometry);
    let mut features = [0.0; MODEL_FEATURE_COUNT];
    for invalid in [f32::NAN, f32::INFINITY] {
        features[3] = invalid;
        assert!(evaluator.prediction_interval(1.0, &features).is_none());
        assert!(evaluator.multiplier.get().is_none());
    }
    features[3] = 2.0;
    assert_eq!(
        bits(evaluator.prediction_interval(1.0, &features)),
        bits(geometry.prediction_interval(1.0, &features))
    );
    assert_eq!(
        evaluator.multiplier.get().unwrap().map(f64::to_bits),
        geometry.t_multiplier().map(f64::to_bits)
    );
}

#[test]
fn unavailable_multiplier_is_cached_but_leverage_is_still_checked() {
    for observations in [1, 2] {
        let mut geometry = frozen_cases::geometry();
        geometry.observations = observations;
        assert!(geometry.validate().is_ok());
        let evaluator = PredictionIntervalEvaluator::new(&geometry);
        let mut features = [0.0; MODEL_FEATURE_COUNT];
        assert!(evaluator.multiplier.get().is_none());
        for value in [0.0, 3.0, f32::NAN, 1.0] {
            features[3] = value;
            assert_eq!(
                bits(evaluator.prediction_interval(1.0, &features)),
                bits(geometry.prediction_interval(1.0, &features))
            );
            assert_eq!(evaluator.multiplier.get(), Some(&None));
        }
    }
}

#[test]
fn nonfinite_prediction_and_endpoint_failure_do_not_cache_interval_failure() {
    let mut geometry = frozen_cases::geometry();
    geometry.residual_standard_error = f64::MAX / 4.0;
    let evaluator = PredictionIntervalEvaluator::new(&geometry);
    let features = [0.0; MODEL_FEATURE_COUNT];
    for prediction in [f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(
            evaluator
                .prediction_interval(prediction, &features)
                .is_none()
        );
        assert!(evaluator.multiplier.get().unwrap().is_some());
    }
    let cached = evaluator.prediction_interval(0.0, &features);
    assert!(cached.is_some());
    assert_eq!(
        bits(cached),
        bits(geometry.prediction_interval(0.0, &features))
    );
}

#[test]
fn candidate_arithmetic_failures_recover_with_the_same_evaluator() {
    let mut tiny_scale = frozen_cases::geometry();
    tiny_scale.scales[0] = f64::MIN_POSITIVE;
    let mut indefinite_inverse = frozen_cases::geometry();
    indefinite_inverse.xtx_inverse = vec![vec![0.1, 1.0], vec![1.0, 0.1]];
    let mut huge_residual = frozen_cases::geometry();
    huge_residual.residual_standard_error = f64::MAX / 4.0;
    for (geometry, first_feature, initialized) in [
        (tiny_scale, f32::MAX, false),
        (indefinite_inverse, -1.0, false),
        (huge_residual, 10.0, true),
    ] {
        let evaluator = PredictionIntervalEvaluator::new(&geometry);
        let mut features = [0.0; MODEL_FEATURE_COUNT];
        features[3] = first_feature;
        assert!(evaluator.prediction_interval(0.0, &features).is_none());
        assert!(geometry.prediction_interval(0.0, &features).is_none());
        assert_eq!(evaluator.multiplier.get().is_some(), initialized);
        features[3] = 0.0;
        let recovered = evaluator.prediction_interval(0.0, &features);
        assert!(recovered.is_some());
        assert_eq!(
            bits(recovered),
            bits(geometry.prediction_interval(0.0, &features))
        );
    }
    println!(
        "PredictionIntervalEvaluator stack bytes: {}",
        std::mem::size_of::<PredictionIntervalEvaluator<'_>>()
    );
}

#[test]
fn interleaved_models_and_new_contexts_keep_their_own_multiplier() {
    let mut first = frozen_cases::geometry();
    let mut second = first.clone();
    second.observations = 3;
    second.feature_indices = vec![1];
    second.scales = vec![2.0];
    let mut features = [0.0; MODEL_FEATURE_COUNT];
    features[1] = -4.0;
    features[3] = 2.0;
    {
        let first_cache = PredictionIntervalEvaluator::new(&first);
        let second_cache = PredictionIntervalEvaluator::new(&second);
        for prediction in [0.0, -2.0, 5.0] {
            for (geometry, evaluator) in [(&first, &first_cache), (&second, &second_cache)] {
                assert_eq!(
                    bits(evaluator.prediction_interval(prediction, &features)),
                    bits(geometry.prediction_interval(prediction, &features))
                );
            }
        }
        assert_ne!(
            first_cache.multiplier.get().unwrap().map(f64::to_bits),
            second_cache.multiplier.get().unwrap().map(f64::to_bits)
        );
    }
    first.observations = 122;
    let fresh = PredictionIntervalEvaluator::new(&first);
    assert!(fresh.multiplier.get().is_none());
    assert_eq!(
        bits(fresh.prediction_interval(1.0, &features)),
        bits(first.prediction_interval(1.0, &features))
    );
}

#[test]
fn unchecked_zero_parameter_edge_matches_public_method_without_endorsement() {
    let case = frozen_cases::cases()
        .into_iter()
        .find(|case| case.name == "unchecked_zero_parameters")
        .unwrap();
    assert!(case.geometry.validate().is_err());
    assert_eq!(case.geometry.leverage(&case.expanded), Some(0.0));
    let evaluator = PredictionIntervalEvaluator::new(&case.geometry);
    let result = evaluator.prediction_interval(case.prediction, &case.expanded);
    assert!(result.is_some());
    assert_eq!(
        bits(result),
        bits(
            case.geometry
                .prediction_interval(case.prediction, &case.expanded)
        )
    );
    assert!(evaluator.multiplier.get().unwrap().is_some());
}
