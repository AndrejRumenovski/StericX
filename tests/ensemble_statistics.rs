use steric_x::ensemble::{
    DescriptorOptions, DescriptorOutcome, EnergyKind, EnergyMetadata, EnergyUnit, Ensemble,
    EnsembleDescriptors, MissingConformerPolicy, Population, StatisticsOptions, SterimolAxis,
};
use steric_x::geometry::SterimolParams;

const SDF: &str = "state\n  StericX\n\n  4  3  0  0  0  0            999 V2000\n0 0 0 P\n1.4 0 0.45 C\n-0.7 1.212 0.45 C\n-0.7 -1.212 0.45 C\n  1  2  1  0\n  1  3  1  0\n  1  4  1  0\nM  END\n$$$$\n";

fn descriptors(values: &[Option<f32>]) -> EnsembleDescriptors {
    let ensemble = Ensemble::from_sdf_str(&SDF.repeat(values.len()), "states.sdf").unwrap();
    let mut options = DescriptorOptions {
        axis: SterimolAxis::Coordination,
        ..Default::default()
    };
    options.volume.density = 0.2;
    let mut result = ensemble.descriptors(&options);
    for (row, value) in result.conformers.iter_mut().zip(values) {
        row.sterimol = match value {
            Some(v) => DescriptorOutcome::Success {
                value: SterimolParams {
                    l: *v,
                    b1: *v,
                    b5: *v,
                },
            },
            None => DescriptorOutcome::Error {
                message: "deliberate missing value".into(),
            },
        };
    }
    result
}

fn close(actual: Option<f64>, expected: f64) {
    assert!((actual.unwrap() - expected).abs() <= 8.0 * f64::EPSILON * expected.abs().max(1.0));
}

#[test]
fn independently_known_statistics_and_supplied_weights() {
    let result = descriptors(&[Some(1.0), Some(2.0), Some(4.0), Some(7.0)]);
    let weights = Population::supplied(&[1.0, 2.0, 3.0, 4.0], "declared test counts").unwrap();
    let stats = result
        .statistics(&StatisticsOptions::default(), Some(&weights))
        .unwrap();
    assert_eq!(stats.descriptors.len(), 14);
    let b5 = &stats.descriptors["sterimol_b5"];
    close(b5.mean, 3.5);
    close(b5.standard_deviation, 5.25_f64.sqrt());
    close(b5.median, 3.0);
    close(b5.minimum, 1.0);
    close(b5.maximum, 7.0);
    close(b5.weighted_mean, 4.5);
    assert_eq!(b5.boltzmann_mean, None);
    assert_eq!(b5.minimum_conformers[0].index, 0);
    assert_eq!(b5.maximum_conformers[0].index, 3);
    assert!(b5.complete);
    assert_eq!(stats.regions.len(), 36);
    for region in stats.regions.values() {
        assert!(region.complete);
        close(region.standard_deviation, 0.0);
    }
}

#[test]
fn failures_require_explicit_conditioning_and_expose_lost_mass() {
    let result = descriptors(&[Some(1.0), None, Some(4.0), Some(7.0)]);
    let weights = Population::supplied(&[1.0, 2.0, 3.0, 4.0], "test counts").unwrap();
    let strict = result
        .statistics(&StatisticsOptions::default(), Some(&weights))
        .unwrap();
    let b5 = &strict.descriptors["sterimol_b5"];
    assert_eq!((b5.total_count, b5.successful_count), (4, 3));
    assert_eq!(b5.missing_conformers[0].index, 1);
    assert!(b5.mean.is_none() && b5.minimum.is_none() && b5.weighted_mean.is_none());
    assert!(b5.unavailable_reason.is_some());
    close(b5.retained_population_mass, 0.8);
    // An independent family remains complete.
    assert!(strict.descriptors["buried_volume"].mean.is_some());
    let conditional = StatisticsOptions {
        missing: MissingConformerPolicy::AvailableCases,
    };
    let stats = result.statistics(&conditional, Some(&weights)).unwrap();
    let b5 = &stats.descriptors["sterimol_b5"];
    close(b5.mean, 4.0);
    close(b5.standard_deviation, 6.0_f64.sqrt());
    close(b5.median, 4.0);
    close(b5.weighted_mean, 41.0 / 8.0);
    assert!(!b5.complete);
    let zero_mass = Population::supplied(&[0.0, 1.0, 0.0, 0.0], "test zero retained mass").unwrap();
    let stats = result.statistics(&conditional, Some(&zero_mass)).unwrap();
    assert!(stats.descriptors["sterimol_b5"].weighted_mean.is_none());
    close(
        stats.descriptors["sterimol_b5"].retained_population_mass,
        0.0,
    );
}

#[test]
fn ties_singletons_nonfinite_values_and_empty_ensembles_are_explicit() {
    let result = descriptors(&[Some(1.0), Some(1.0), Some(2.0), Some(2.0)]);
    let stats = result
        .statistics(&StatisticsOptions::default(), None)
        .unwrap();
    let b5 = &stats.descriptors["sterimol_b5"];
    assert_eq!(
        b5.minimum_conformers
            .iter()
            .map(|c| c.index)
            .collect::<Vec<_>>(),
        [0, 1]
    );
    assert_eq!(
        b5.maximum_conformers
            .iter()
            .map(|c| c.index)
            .collect::<Vec<_>>(),
        [2, 3]
    );
    assert!(b5.weighted_mean.is_none() && b5.boltzmann_mean.is_none());
    let single = descriptors(&[Some(f32::MAX)])
        .statistics(&StatisticsOptions::default(), None)
        .unwrap();
    close(single.descriptors["sterimol_b5"].mean, f64::from(f32::MAX));
    close(single.descriptors["sterimol_b5"].standard_deviation, 0.0);
    let invalid = descriptors(&[Some(f32::NAN)])
        .statistics(&StatisticsOptions::default(), None)
        .unwrap();
    assert_eq!(invalid.descriptors["sterimol_b5"].successful_count, 0);
    let mut empty = result;
    empty.conformers.clear();
    assert!(
        empty
            .statistics(&StatisticsOptions::default(), None)
            .is_err()
    );
}

#[test]
fn population_validation_scale_invariance_and_input_length() {
    for values in [
        &[][..],
        &[0.0, 0.0],
        &[-1.0, 1.0],
        &[f64::NAN],
        &[f64::INFINITY],
    ] {
        assert!(Population::supplied(values, "test").is_err());
    }
    assert!(Population::supplied(&[1.0], " ").is_err());
    for scale in [1e-300, 1.0, 1e300] {
        let weights = Population::supplied(&[scale, 2.0 * scale], "test").unwrap();
        close(Some(weights.weights()[0]), 1.0 / 3.0);
        close(Some(weights.weights()[1]), 2.0 / 3.0);
    }
    let result = descriptors(&[Some(1.0), Some(2.0)]);
    let wrong_length = Population::supplied(&[1.0], "test").unwrap();
    assert!(
        result
            .statistics(&StatisticsOptions::default(), Some(&wrong_length))
            .is_err()
    );
}

#[test]
fn correspondence_checks_include_later_supplied_graphs() {
    let mut result = descriptors(&[Some(1.0), Some(2.0), Some(3.0)]);
    if let DescriptorOutcome::Success { value } = &mut result.conformers[0].geometry {
        value.bonds = None;
    }
    let stats = result
        .statistics(&StatisticsOptions::default(), None)
        .unwrap();
    assert!(stats.regions.is_empty());
    assert!(stats.regional_unavailable_reason.is_some());
    if let DescriptorOutcome::Success { value } = &mut result.conformers[2].geometry {
        value.bonds = Some(vec![[0, 1], [1, 2], [2, 3]]);
    }
    assert!(
        result
            .statistics(&StatisticsOptions::default(), None)
            .unwrap_err()
            .contains("connectivity")
    );
    if let DescriptorOutcome::Success { value } = &mut result.conformers[2].geometry {
        value.elements[0] = "N".into();
    }
    assert!(
        result
            .statistics(&StatisticsOptions::default(), None)
            .unwrap_err()
            .contains("identities")
    );
}

#[test]
fn energy_weighting_requires_opt_in_and_matches_ordered_provenance() {
    let mut ensemble = Ensemble::from_sdf_str(&SDF.repeat(2), "energy.sdf").unwrap();
    assert!(ensemble.boltzmann_population(298.15, "test").is_err());
    for i in 0..2 {
        ensemble
            .set_energy(
                i,
                EnergyMetadata {
                    value: 0.0,
                    unit: EnergyUnit::KcalMol,
                    kind: EnergyKind::Electronic,
                    model: "test method".into(),
                    provenance: format!("state_{i}.log"),
                    temperature_k: None,
                    degeneracy: 1,
                },
            )
            .unwrap();
    }
    let mut options = DescriptorOptions {
        axis: SterimolAxis::Coordination,
        ..Default::default()
    };
    options.volume.density = 0.2;
    let mut result = ensemble.descriptors(&options);
    let stats = result
        .statistics(&StatisticsOptions::default(), None)
        .unwrap();
    assert!(stats.descriptors["sterimol_b5"].boltzmann_mean.is_none());
    let population = ensemble
        .boltzmann_population(
            298.15,
            "distinct enumerated states, electronic-energy approximation",
        )
        .unwrap();
    let stats = result
        .statistics(&StatisticsOptions::default(), Some(&population))
        .unwrap();
    close(
        stats.descriptors["sterimol_b5"].boltzmann_mean,
        stats.descriptors["sterimol_b5"].mean.unwrap(),
    );
    result.conformers.swap(0, 1);
    assert!(
        result
            .statistics(&StatisticsOptions::default(), Some(&population))
            .unwrap_err()
            .contains("provenance")
    );
}

#[test]
fn sensitivity_reports_preserve_ranges_extremizers_and_csv_identity() {
    let mut result = descriptors(&[Some(1.0), Some(1.0), Some(7.0)]);
    result.conformers[0].identifier = "comma, quote\" and\nnewline".into();
    let statistics = result
        .statistics(&StatisticsOptions::default(), None)
        .unwrap();
    let report = statistics.sensitivity("ligand,\"A\"");
    assert_eq!(report.ranges["sterimol_b5"], Some(6.0));
    assert_eq!(
        report.statistics.descriptors["sterimol_b5"]
            .minimum_conformers
            .len(),
        2
    );
    for range in report.regional_ranges.values() {
        assert_eq!(*range, Some(0.0));
    }
    let text = report.to_csv().unwrap();
    let mut reader = csv::Reader::from_reader(text.as_bytes());
    let headers = reader.headers().unwrap().clone();
    let column = |name: &str| headers.iter().position(|h| h == name).unwrap();
    let rows = reader.records().collect::<Result<Vec<_>, _>>().unwrap();
    assert_eq!(rows.len(), 50);
    let b5 = rows
        .iter()
        .find(|r| &r[column("descriptor")] == "sterimol_b5")
        .unwrap();
    assert_eq!(&b5[column("ligand_id")], "ligand,\"A\"");
    assert_eq!(&b5[column("range")], "6");
    assert_eq!(&b5[column("boltzmann_mean")], "");
    let ids: serde_json::Value = serde_json::from_str(&b5[column("minimum_conformers")]).unwrap();
    assert_eq!(ids[0]["identifier"], result.conformers[0].identifier);
    assert_eq!(ids[1]["index"], 1);
    result.conformers[1].sterimol = DescriptorOutcome::Error {
        message: "failed state".into(),
    };
    let incomplete = result
        .statistics(&StatisticsOptions::default(), None)
        .unwrap()
        .sensitivity("incomplete");
    assert_eq!(incomplete.ranges["sterimol_b5"], None);
    assert_eq!(
        incomplete.statistics.descriptors["sterimol_b5"].missing_conformers[0].index,
        1
    );
}
