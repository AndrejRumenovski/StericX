use super::*;

fn pair() -> [BuriedVolumeParams; 2] {
    [
        BuriedVolumeParams {
            buried_volume: 10.0,
            qvbur_min: 2.0,
            qvbur_max: 8.0,
            max_delta_qvbur: 6.0,
            near_vbur: 4.0,
            far_vbur: 6.0,
            ..BuriedVolumeParams::default()
        },
        BuriedVolumeParams {
            buried_volume: 30.0,
            qvbur_min: 6.0,
            qvbur_max: 24.0,
            max_delta_qvbur: 18.0,
            near_vbur: 12.0,
            far_vbur: 18.0,
            ..BuriedVolumeParams::default()
        },
    ]
}

#[test]
fn common_weight_scale_preserves_analytic_means_and_extrema() {
    for weights in [
        [1.0, 3.0],
        [1e-10, 3e-10],
        [1e38, 3e38],
        [f32::from_bits(1), f32::from_bits(3)],
    ] {
        let result = BuriedVolumeCalculator::aggregate(&pair(), &weights).unwrap();
        assert!((result.vbur_boltz - 25.0).abs() < 4e-6);
        assert!((result.qvbur_min_boltz - 5.0).abs() < 1e-6);
        assert_eq!(result.vbur_min, 10.0);
        assert_eq!(result.vbur_max, 30.0);
        assert_eq!(result.vbur_delta, 20.0);
        assert_eq!(result.max_delta_qvbur_vburminconf, 6.0);
    }
}

#[test]
fn invalid_probabilities_and_nonfinite_descriptor_fields_are_rejected() {
    for weights in [
        [0.0, 0.0],
        [-0.5, 1.5],
        [f32::NAN, 1.0],
        [f32::INFINITY, 1.0],
    ] {
        assert!(BuriedVolumeCalculator::aggregate(&pair(), &weights).is_err());
    }
    let setters: [fn(&mut BuriedVolumeParams, f32); 9] = [
        |p, v| p.buried_volume = v,
        |p, v| p.percent_buried_volume = v,
        |p, v| p.qvbur_min = v,
        |p, v| p.qvbur_max = v,
        |p, v| p.max_delta_qvbur = v,
        |p, v| p.ovbur_min = v,
        |p, v| p.ovbur_max = v,
        |p, v| p.near_vbur = v,
        |p, v| p.far_vbur = v,
    ];
    for set in setters {
        for value in [f32::NAN, f32::INFINITY, f32::NEG_INFINITY] {
            let mut conformers = pair();
            set(&mut conformers[1], value);
            // Invalid zero-weight states must not be silently accepted either.
            assert!(BuriedVolumeCalculator::aggregate(&conformers, &[1.0, 0.0]).is_err());
        }
    }
}

#[test]
fn tiny_weight_is_not_rounded_before_multiplication() {
    let mut conformers = pair();
    conformers[0].buried_volume = f32::MAX;
    conformers[1].buried_volume = 0.0;
    let result = BuriedVolumeCalculator::aggregate(&conformers, &[f32::from_bits(1), 2.0]).unwrap();
    let reference = (f64::from(f32::MAX) * f64::from(f32::from_bits(1)) / 2.0) as f32;
    assert_eq!(result.vbur_boltz, reference);
    assert!(result.vbur_boltz > 0.0);
}
