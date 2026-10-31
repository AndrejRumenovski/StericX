use steric_x::ensemble::{EnergyKind, EnergyMetadata, EnergyUnit, Population};

fn energies(values: &[f64]) -> Vec<EnergyMetadata> {
    values
        .iter()
        .enumerate()
        .map(|(i, value)| EnergyMetadata {
            value: *value,
            unit: EnergyUnit::KcalMol,
            kind: EnergyKind::Electronic,
            model: "analytical test".into(),
            provenance: format!("analytical state {i}"),
            temperature_k: None,
            degeneracy: 1,
        })
        .collect()
}

#[test]
fn two_state_analytic_ratios_and_temperature_response() {
    // Independent SI R, not imported from the implementation.
    let r = 1.380_649e-23 * 6.022_140_76e23 / 4184.0;
    let temperature = 298.15;
    let states = energies(&[0.0, r * temperature * 2.0_f64.ln()]);
    let population = Population::boltzmann(&states, temperature, "enumerated states").unwrap();
    assert!((population.weights()[0] - 2.0 / 3.0).abs() < 8.0 * f64::EPSILON);
    let hotter = Population::boltzmann(&states, 2.0 * temperature, "enumerated states").unwrap();
    assert!(
        (hotter.weights()[0] - 2.0_f64.sqrt() / (1.0 + 2.0_f64.sqrt())).abs() < 8.0 * f64::EPSILON
    );
    assert!(hotter.weights()[1] > population.weights()[1]);
}

#[test]
fn multiplicity_and_explicit_degeneracy_are_not_merged() {
    let mut grouped = energies(&[0.0, 0.0, 0.0]);
    for (state, g) in grouped.iter_mut().zip([1, 2, 4]) {
        state.degeneracy = g;
    }
    let population = Population::boltzmann(
        &grouped,
        300.0,
        "three groups representing seven microstates",
    )
    .unwrap();
    for (actual, g) in population.weights().iter().zip([1.0, 2.0, 4.0]) {
        assert!((actual - g / 7.0).abs() < 8.0 * f64::EPSILON);
    }
    let explicit = Population::boltzmann(
        &energies(&[0.0; 7]),
        300.0,
        "seven distinct equal-energy states",
    )
    .unwrap();
    assert_eq!(explicit.weights().len(), 7);
    for actual in explicit.weights() {
        assert_eq!(*actual, 1.0 / 7.0);
    }
}

#[test]
fn invalid_thermodynamic_metadata_and_nonfinite_inputs_fail() {
    let valid = energies(&[0.0, 1.0]);
    for temperature in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(Population::boltzmann(&valid, temperature, "test").is_err());
    }
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(Population::boltzmann(&energies(&[0.0, value]), 300.0, "test").is_err());
    }
    assert!(Population::boltzmann(&[], 300.0, "test").is_err());
    assert!(Population::boltzmann(&valid, 300.0, "").is_err());
    let mut free = valid.clone();
    for state in &mut free {
        state.kind = EnergyKind::GibbsFreeEnergy;
        state.temperature_k = Some(298.15);
    }
    assert!(
        Population::boltzmann(&free, 298.15, "free energies, no double-counted entropy").is_ok()
    );
    assert!(Population::boltzmann(&free, 300.0, "test").is_err());
    for field in ["unit", "kind", "model", "provenance", "degeneracy"] {
        let mut invalid = valid.clone();
        match field {
            "unit" => invalid[1].unit = EnergyUnit::Hartree,
            "kind" => invalid[1].kind = EnergyKind::ForceField,
            "model" => invalid[1].model = "different model".into(),
            "provenance" => invalid[1].provenance.clear(),
            "degeneracy" => invalid[1].degeneracy = 0,
            _ => unreachable!(),
        }
        assert!(
            Population::boltzmann(&invalid, 300.0, "test").is_err(),
            "{field}"
        );
    }
}

#[test]
fn extreme_inputs_remain_normalized_without_overflow() {
    let cold = Population::boltzmann(&energies(&[0.0, 1.0]), f64::from_bits(1), "test").unwrap();
    assert_eq!(cold.weights(), [1.0, 0.0]);
    let mut extremes = energies(&[-f64::MAX, f64::MAX]);
    for state in &mut extremes {
        state.unit = EnergyUnit::KjMol;
    }
    let population = Population::boltzmann(
        &extremes,
        f64::MAX,
        "arithmetic stress, not physical temperatures",
    )
    .unwrap();
    // Delta-E over T is exactly 2 in this analytical extreme case.
    let expected = (-2.0_f64 / (1.380_649e-23 * 6.022_140_76e23 / 1000.0)).exp();
    assert!(population.weights()[1] > 0.0);
    assert!((population.weights()[1] / expected - 1.0).abs() < 2e-12);
}
