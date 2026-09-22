use steric_x::EyringKineticLink;
use steric_x::kinetics::RateConstantError;

// Decimal references use exact SI kB, NA, h and 4184 J/kcal at the actual
// represented f32 temperature 298.149993896484375 K, precision 80 digits.
#[test]
fn absolute_barriers_have_independent_rates_even_at_identical_selectivity() {
    for (barrier, reference) in [
        (0.0, 6_212_437_864_443.485_f64),
        (10.0, 290_543.418_103_518_5),
        (11.0, 53_728.651_401_089_606),
        (20.0, 0.013_588_140_379_869_692),
        (21.0, 0.002_512_782_641_660_009_7),
        (50.0, 1.389_974_955_735_155_5e-24),
    ] {
        let actual = EyringKineticLink::calculate_rate_constant_checked(barrier, 298.15).unwrap();
        assert!((f64::from(actual) / reference - 1.0).abs() < 6e-8);
    }
    let first = EyringKineticLink::calculate_rate_constant(10.0, 298.15)
        / EyringKineticLink::calculate_rate_constant(11.0, 298.15);
    let second = EyringKineticLink::calculate_rate_constant(20.0, 298.15)
        / EyringKineticLink::calculate_rate_constant(21.0, 298.15);
    assert!((first / second - 1.0).abs() < 2e-7);
    assert!(
        EyringKineticLink::calculate_rate_constant(10.0, 298.15)
            > 1e7 * EyringKineticLink::calculate_rate_constant(20.0, 298.15)
    );
}

#[test]
fn tiny_selectivity_and_minority_are_not_erased_by_subtraction() {
    let ee = EyringKineticLink::calculate_enantiomeric_excess(1e-12, 298.15);
    assert!((f64::from(ee) / 8.439_033_067_072_845e-11 - 1.0).abs() < 6e-8);
    for (difference, expected) in [
        (20.0, 2.187_247_691_866_757_4e-13),
        (50.0, 2.237_406_612_451_761_4e-35),
    ] {
        let (major, minor) = EyringKineticLink::calculate_enantiomeric_ratio(difference, 298.15);
        assert_eq!(major, 100.0);
        assert!((f64::from(minor) / expected - 1.0).abs() < 6e-8);
        let positive = EyringKineticLink::product_ratio(difference, 298.15);
        let negative = EyringKineticLink::product_ratio(-difference, 298.15);
        assert_eq!(positive.percent_s, negative.percent_r);
        assert_eq!(positive.percent_r, negative.percent_s);
        assert_eq!(positive.ee_percent, negative.ee_percent);
    }
}

#[test]
fn extreme_rates_have_explicit_range_errors_and_finite_log_rates() {
    assert_eq!(
        EyringKineticLink::calculate_rate_constant_checked(100.0, 298.15),
        Err(RateConstantError::Underflow)
    );
    assert_eq!(
        EyringKineticLink::calculate_rate_constant_checked(-100.0, 298.15),
        Err(RateConstantError::Overflow)
    );
    for barrier in [-f32::MAX, 0.0, f32::MAX] {
        for temperature in [f32::from_bits(1), 298.15, f32::MAX] {
            assert!(
                EyringKineticLink::calculate_log_rate_constant(barrier, temperature)
                    .unwrap()
                    .is_finite()
            );
        }
    }
    for (barrier, temperature) in [
        (f32::NAN, 298.15),
        (f32::INFINITY, 298.15),
        (1.0, 0.0),
        (1.0, -1.0),
        (1.0, f32::INFINITY),
    ] {
        assert_eq!(
            EyringKineticLink::calculate_rate_constant_checked(barrier, temperature),
            Err(RateConstantError::InvalidInput)
        );
        assert!(EyringKineticLink::calculate_enantiomeric_excess(barrier, temperature).is_nan());
    }
}

#[test]
fn difference_only_cli_does_not_invent_an_absolute_rate() {
    let result = std::process::Command::new(env!("CARGO_BIN_EXE_stericx"))
        .args(["simulate", "--ddg", "1", "--temp", "298.15"])
        .output()
        .unwrap();
    assert!(result.status.success());
    let output = String::from_utf8(result.stdout).unwrap();
    assert!(!output.contains("rate_constant"));
    assert!(output.contains("ddg_convention=G_dagger_S_minus_G_dagger_R"));
    assert!(
        output.contains("selectivity_assumption=competing_irreversible_pathways_equal_prefactors")
    );
    assert!(output.contains("ee_percent="));
}
