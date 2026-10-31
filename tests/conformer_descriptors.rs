use steric_x::ensemble::{
    DescriptorOptions, EnergyKind, EnergyMetadata, EnergyUnit, Ensemble, SterimolAxis,
};
use steric_x::geometry::{BuriedVolumeCalculator, coordination_center_with_neighbors};

const SDF: &str = "same title\n  StericX\n\n  4  3  0  0  0  0            999 V2000\n0 0 0 P\n1.4 0 0.45 C\n-0.7 1.212 0.45 C\n-0.7 -1.212 0.45 C\n  1  2  1  0\n  1  3  1  0\n  1  4  1  0\nM  END\n>  <ENERGY>\n-123.45\n\n>  <ENERGY>\nunknown units\n\n$$$$\n";

fn coordination() -> DescriptorOptions {
    DescriptorOptions {
        axis: SterimolAxis::Coordination,
        ..Default::default()
    }
}

#[test]
fn malformed_middle_records_and_duplicate_identifiers_are_retained() {
    let bad = "bad record\n$$$$\n";
    let text = format!("{SDF}{bad}{SDF}");
    let ensemble = Ensemble::from_sdf_str(&text, "three.sdf").unwrap();
    let results = ensemble.descriptors(&coordination());
    assert_eq!(results.conformers.len(), 3);
    for (i, r) in results.conformers.iter().enumerate() {
        assert_eq!(r.index, i);
    }
    assert_eq!(
        results.conformers[0].identifier,
        results.conformers[2].identifier
    );
    assert_eq!(results.conformers[0].source_geometry, SDF);
    assert_eq!(results.conformers[1].source_geometry, bad);
    assert!(results.conformers[0].sterimol.value().is_some());
    assert!(results.conformers[1].geometry.error().is_some());
    assert!(results.conformers[1].sterimol.error().is_some());
    assert!(results.conformers[1].buried_volume.error().is_some());
    assert!(results.conformers[2].buried_volume.value().is_some());
    let properties = &results.conformers[0].properties;
    assert_eq!(properties.len(), 2);
    assert_eq!(properties[0].name, "ENERGY");
    assert_eq!(properties[1].name, "ENERGY");
    assert_eq!(properties[0].value, "-123.45");
    assert!(results.conformers[0].energy.is_none()); // No invented units.
    assert_eq!(
        ensemble.conformers[0]
            .frame
            .as_ref()
            .unwrap()
            .bonds
            .as_ref()
            .unwrap(),
        &[[0, 1], [0, 2], [0, 3]]
    );
    let json = serde_json::to_value(&results).unwrap();
    assert_eq!(json["conformers"][1]["sterimol"]["status"], "error");
}

#[test]
fn source_blank_titles_crlf_and_explicit_empty_records_are_not_dropped() {
    let blank = SDF.replacen("same title", "", 1).replace('\n', "\r\n");
    let ensemble = Ensemble::from_sdf_str(&format!("{blank}$$$$\n{blank}"), "blank.sdf").unwrap();
    assert_eq!(ensemble.conformers.len(), 3);
    assert_eq!(ensemble.conformers[0].identifier, "");
    assert!(ensemble.conformers[0].frame.is_ok());
    assert!(ensemble.conformers[1].frame.is_err());
    assert!(ensemble.conformers[2].frame.is_ok());
    assert_eq!(ensemble.conformers[0].source_geometry, blank);
}

#[test]
fn explicit_bond_support_survives_unavailable_trivalent_descriptors() {
    let ensemble =
        Ensemble::from_xyz_str("2\nethane axis\nC 0 0 0\nC 0 0 1.5\n", "axis.xyz").unwrap();
    let options = DescriptorOptions {
        donor_index: Some(0),
        reference_index: Some(1),
        ..Default::default()
    };
    let result = ensemble.descriptors(&options);
    let row = &result.conformers[0];
    // Independent support of a single radius-1.7 sphere 1.5 Å along the axis.
    let sterimol = row.sterimol.value().unwrap();
    assert!((sterimol.l - 3.2).abs() < 3e-7);
    assert_eq!(sterimol.b1, 1.7);
    assert_eq!(sterimol.b5, 1.7);
    assert!(row.buried_volume.error().is_some());
    assert!(row.pyramidalization.error().is_some());
}

#[test]
fn failed_volume_settings_do_not_erase_other_descriptor_results() {
    let ensemble = Ensemble::from_sdf_str(SDF, "valid.sdf").unwrap();
    let mut options = DescriptorOptions {
        donor_index: Some(0),
        reference_index: Some(1),
        ..Default::default()
    };
    options.volume.sphere_radius = -1.0;
    let results = ensemble.descriptors(&options);
    let row = &results.conformers[0];
    assert!(row.sterimol.value().is_some());
    assert!(row.pyramidalization.value().is_some());
    assert!(row.buried_volume.error().unwrap().contains("sphere radius"));
}

#[test]
fn regional_output_preserves_native_summary_and_defines_every_frame() {
    let ensemble = Ensemble::from_sdf_str(SDF, "valid.sdf").unwrap();
    let options = coordination();
    let results = ensemble.descriptors(&options);
    let details = results.conformers[0].buried_volume.value().unwrap();
    let molecule = &ensemble.conformers[0].frame.as_ref().unwrap().molecule;
    let legacy =
        BuriedVolumeCalculator::compute_with_neighbors(molecule, 0, [1, 2, 3], options.volume)
            .unwrap();
    assert_eq!(details.summary, legacy);
    assert_eq!(
        details.center,
        coordination_center_with_neighbors(molecule, 0, [1, 2, 3], options.volume)
            .unwrap()
            .to_array()
    );
    let mut planes = details
        .orientations
        .iter()
        .map(|r| r.plane_atom_index)
        .collect::<Vec<_>>();
    planes.sort();
    assert_eq!(planes, [1, 2, 3]);
    for region in &details.orientations {
        assert_eq!(region.quadrants.len(), 4);
        assert_eq!(region.octants.len(), 8);
        for axis in region.axes {
            assert!((axis.into_iter().map(|x| x * x).sum::<f32>() - 1.0).abs() < 5e-7);
        }
        for q in region.quadrants {
            assert!(q >= legacy.qvbur_min && q <= legacy.qvbur_max);
        }
        for o in region.octants {
            assert!(o >= legacy.ovbur_min && o <= legacy.ovbur_max);
        }
    }
}

#[test]
fn typed_energy_is_preserved_and_invalid_metadata_is_explicit() {
    let mut ensemble = Ensemble::from_sdf_str(SDF, "energy.sdf").unwrap();
    let energy = EnergyMetadata {
        value: -123.45,
        unit: EnergyUnit::Hartree,
        kind: EnergyKind::Electronic,
        model: "declared electronic method / gas".into(),
        provenance: "calculation.log".into(),
        temperature_k: None,
        degeneracy: 1,
    };
    ensemble.set_energy(0, energy.clone()).unwrap();
    assert_eq!(
        ensemble.descriptors(&coordination()).conformers[0].energy,
        Some(energy.clone())
    );
    let mut invalid = energy;
    invalid.value = f64::NAN;
    assert!(ensemble.set_energy(0, invalid.clone()).is_err());
    // Even a direct unchecked caller mutation cannot serialize NaN as energy null
    // without the accompanying explicit error.
    ensemble.conformers[0].energy = Some(invalid);
    let result = ensemble.descriptors(&coordination());
    assert!(result.conformers[0].energy_error.is_some());
    assert!(result.conformers[0].energy.is_none());
    assert!(result.conformers[0].sterimol.value().is_some());
}

#[test]
fn multi_xyz_preserves_bad_frames_and_unknown_boundaries() {
    let good = "1\none\nH 0 0 0\n";
    let bad = "1\ninvalid coordinate\nH NaN 0 0\n";
    let tail = format!("unknown count\nremaining text\n{good}");
    let ensemble =
        Ensemble::from_xyz_str(&format!("{good}{bad}{good}{tail}"), "multi.xyz").unwrap();
    assert_eq!(ensemble.conformers.len(), 4);
    assert!(ensemble.conformers[0].frame.is_ok());
    assert!(ensemble.conformers[1].frame.is_err());
    assert!(ensemble.conformers[2].frame.is_ok());
    assert!(ensemble.conformers[3].frame.is_err());
    assert_eq!(ensemble.conformers[3].source_geometry, tail);
}
