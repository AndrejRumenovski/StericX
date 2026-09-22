use crate::model::{MODEL_FEATURE_COUNT, TrainingGeometry};

pub struct Case {
    pub name: String,
    pub geometry: TrainingGeometry,
    pub expanded: [f32; MODEL_FEATURE_COUNT],
    pub prediction: f64,
}

pub fn geometry() -> TrainingGeometry {
    TrainingGeometry {
        feature_indices: vec![3],
        means: vec![0.0],
        scales: vec![1.0],
        xtx_inverse: vec![vec![0.1, 0.0], vec![0.0, 0.1]],
        observations: 10,
        parameters: 2,
        residual_standard_error: 0.5,
        warning_leverage: 0.6,
        standardized_training_points: vec![],
        training_labels: vec![],
        neighbor_calibration: None,
    }
}

pub fn cases() -> Vec<Case> {
    let mut out = Vec::new();
    let mut add = |name: &str, geometry, feature, prediction| {
        let mut expanded = [0.0; MODEL_FEATURE_COUNT];
        expanded[3] = feature;
        out.push(Case {
            name: name.to_owned(),
            geometry,
            expanded,
            prediction,
        });
    };
    for n in [3, 4, 7, 10, 32, 122, 1002] {
        let mut g = geometry();
        g.observations = n;
        add(&format!("df_{}_center", n - 2), g.clone(), 0.0, 1.0);
        add(&format!("df_{}_far", n - 2), g, 3.0, -2.0);
    }
    for (name, s) in [("positive_zero", 0.0), ("negative_zero", -0.0)] {
        let mut g = geometry();
        g.residual_standard_error = s;
        add(&format!("{name}_positive_prediction"), g.clone(), 0.0, 0.0);
        add(&format!("{name}_negative_prediction"), g, 3.0, -0.0);
    }
    for n in [0, 1, 2] {
        let mut g = geometry();
        g.observations = n;
        add(&format!("observations_{n}"), g, 0.0, 1.0);
    }
    for (name, s) in [
        ("scale_zero", 0.0),
        ("scale_negative", -1.0),
        ("scale_nan", f64::NAN),
        ("scale_infinity", f64::INFINITY),
    ] {
        let mut g = geometry();
        g.scales[0] = s;
        add(name, g, 1.0, 1.0);
    }
    for (name, s) in [
        ("residual_negative", -1.0),
        ("residual_nan", f64::NAN),
        ("residual_infinity", f64::INFINITY),
    ] {
        let mut g = geometry();
        g.residual_standard_error = s;
        add(name, g, 1.0, 1.0);
    }
    let mut g = geometry();
    g.feature_indices[0] = MODEL_FEATURE_COUNT;
    add("feature_out_of_bounds", g, 1.0, 1.0);
    let mut g = geometry();
    g.feature_indices.push(3);
    g.means.push(0.0);
    g.scales.push(1.0);
    g.parameters = 3;
    g.xtx_inverse = vec![
        vec![0.1, 0.0, 0.0],
        vec![0.0, 0.1, 0.0],
        vec![0.0, 0.0, 0.1],
    ];
    add("duplicate_feature", g, 1.0, 1.0);
    let mut g = geometry();
    g.xtx_inverse[1].pop();
    add("malformed_inverse", g, 1.0, 1.0);
    let mut g = geometry();
    g.xtx_inverse[1][1] = f64::NAN;
    add("nonfinite_inverse", g, 1.0, 1.0);
    let mut g = geometry();
    g.standardized_training_points = vec![vec![0.0]];
    add("inconsistent_training_points", g, 1.0, 1.0);
    let mut g = geometry();
    g.feature_indices.clear();
    g.means.clear();
    g.scales.clear();
    g.xtx_inverse.clear();
    g.parameters = 0;
    add("unchecked_zero_parameters", g, 1.0, 1.0);
    add("selected_nan", geometry(), f32::NAN, 1.0);
    add("selected_infinity", geometry(), f32::INFINITY, 1.0);
    let mut g = geometry();
    g.scales[0] = f64::MIN_POSITIVE;
    add("leverage_overflow", g.clone(), f32::MAX, 1.0);
    add("after_leverage_overflow", g, 0.0, 1.0);
    let mut g = geometry();
    g.xtx_inverse = vec![vec![0.1, 1.0], vec![1.0, 0.1]];
    add("finite_negative_leverage", g, -1.0, 1.0);
    add("prediction_nan", geometry(), 0.0, f64::NAN);
    add("prediction_infinity", geometry(), 0.0, f64::INFINITY);
    let mut g = geometry();
    g.residual_standard_error = f64::MAX;
    add("half_width_overflow", g, 0.0, 0.0);
    let mut g = geometry();
    g.residual_standard_error = f64::MAX / 4.0;
    add("endpoint_overflow", g.clone(), 0.0, f64::MAX);
    add("after_endpoint_overflow", g, 0.0, 0.0);
    let mut unused = [0.0; MODEL_FEATURE_COUNT];
    unused[1] = f32::NAN;
    out.push(Case {
        name: "unused_nan".to_owned(),
        geometry: geometry(),
        expanded: unused,
        prediction: 1.0,
    });
    unused[1] = f32::INFINITY;
    out.push(Case {
        name: "unused_infinity".to_owned(),
        geometry: geometry(),
        expanded: unused,
        prediction: 1.0,
    });
    out
}
