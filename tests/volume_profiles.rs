use steric_x::density::*;
use steric_x::ensemble::{DescriptorOptions, Ensemble};
#[test]
fn hard_sphere_scan_keeps_center_order_duplicates_and_failures() {
    let ensemble = Ensemble::from_xyz_str(
        "4\nprofile\nP 0 0 0\nC 1.4 0 .45\nC -.7 1.212 .45\nC -.7 -1.212 .45\n",
        "profile.xyz",
    )
    .unwrap();
    let options = DescriptorOptions::default();
    let conformer = &ensemble.conformers[0];
    let profile = conformer
        .buried_volume_scan(&[3.5, -1., 2.5, 3.5], &options)
        .unwrap();
    assert_eq!(profile.samples.len(), 4);
    assert!(profile.samples[1].outcome.error().is_some());
    assert_eq!(profile.samples[0].outcome, profile.samples[3].outcome);
    for sample in &profile.samples {
        if let Some(value) = sample.outcome.value() {
            assert_eq!(Some(value.center), profile.center_angstrom);
        }
    }
    let mut bad = conformer.clone();
    bad.frame = Err("bad geometry".into());
    let p = bad.buried_volume_scan(&[2., 3.], &options).unwrap();
    assert!(
        p.samples
            .iter()
            .all(|s| s.outcome.error() == Some("bad geometry"))
    );
    assert!(conformer.buried_volume_scan(&[f64::NAN], &options).is_err());
    assert!(conformer.buried_volume_scan(&[], &options).is_err());
}
#[test]
fn density_scan_retains_out_of_grid_positions() {
    let g = DensityGrid::from_cube_str(
        "grid",
        "uniform\ndensity\n1 -2 -2 -2\n2 4 0 0\n2 0 4 0\n2 0 0 4\n1 1 0 0 0\n2 2 2 2 2 2 2 2",
        DensityDeclaration {
            field: DensityField::TotalElectronDensity,
            unit: DensityUnit::ElectronsPerBohrCubed,
            provenance: "enumerable".into(),
        },
    )
    .unwrap();
    // All these radii have a defined failure: no nodes, incomplete domain, or invalid radius.
    let p = g
        .buried_volume_scan([0.; 3], &[0.1, 3., -1., 0.1], 1.)
        .unwrap();
    assert_eq!(p.samples.len(), 4);
    assert!(p.samples.iter().all(|s| s.outcome.error().is_some()));
    assert_eq!(
        p.samples[0].requested_radius_angstrom,
        p.samples[3].requested_radius_angstrom
    );
}
