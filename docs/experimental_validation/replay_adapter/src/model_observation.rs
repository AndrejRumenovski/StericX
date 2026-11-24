//! Audit observation adapter only: delegates to the frozen public API.
//! These are system-under-test observations, never independent references.
use serde_json::{Value, json};
use steric_x::model::{TrainingGeometry, FeatureDomain, DomainRule, NeighborCalibration,
    assess_applicability, student_t_two_sided_quantile};

pub fn observe(v: &Value) -> Value {
    if let Some(df) = v.get("df").and_then(Value::as_f64) {
        let alpha = v.get("alpha").and_then(Value::as_f64).unwrap_or(0.05);
        return json!({"t_quantile":student_t_two_sided_quantile(alpha, df)});
    }
    if let Some(points) = v.get("points") {
        let points: Vec<Vec<f64>> = serde_json::from_value(points.clone()).expect("valid audit points");
        return json!({"calibration":NeighborCalibration::from_points(&points)});
    }
    let geometry: TrainingGeometry = serde_json::from_value(v["geometry"].clone()).expect("valid audit geometry");
    let features: [f32; 8] = serde_json::from_value(v["features"].clone()).expect("eight finite audit features");
    let ranges: Vec<FeatureDomain> = serde_json::from_value(v.get("ranges").cloned().unwrap_or(json!([]))).expect("valid ranges");
    let rule: DomainRule = serde_json::from_value(v.get("rule").cloned().unwrap_or(json!("max_neighbor"))).expect("valid domain rule");
    let prediction=v.get("prediction").and_then(Value::as_f64).unwrap_or(0.0);
    json!({"design":geometry.design_vector(&features),"leverage":geometry.leverage(&features),
        "mahalanobis":geometry.mahalanobis_distance(&features),"prediction_interval":geometry.prediction_interval(prediction,&features),
        "confidence_interval":geometry.confidence_interval(prediction,&features),"t_multiplier":geometry.t_multiplier(),
        "assessment":assess_applicability(Some(&geometry),&ranges,&geometry.feature_indices,&features,rule)})
}
