use steric_x::density::*;
fn grid(shift: f64, mode: &str) -> DensityGrid {
    let mut s = format!(
        "analytic\ntotal density\n1 {} -2 -2\n5 1 0 0\n5 0 1 0\n5 0 0 1\n1 1 {shift} 0 0\n",
        shift - 2.
    );
    for i in 0..5 {
        for _j in 0..5 {
            for _k in 0..5 {
                let rho = match mode {
                    "full" => 2.,
                    "empty" => 0.,
                    _ => {
                        if i > 2 {
                            2.
                        } else {
                            1.
                        }
                    }
                };
                s.push_str(&format!("{rho} "));
            }
        }
    }
    DensityGrid::from_cube_str(
        "analytic",
        &s,
        DensityDeclaration {
            field: DensityField::TotalElectronDensity,
            unit: DensityUnit::ElectronsPerBohrCubed,
            provenance: "enumerable field".into(),
        },
    )
    .unwrap()
}
#[test]
fn independently_enumerable_occupation_and_strict_threshold() {
    for (mode, count) in [("full", 7), ("empty", 0), ("half", 1)] {
        let r = grid(0., mode)
            .buried_volume([0.; 3], 1.01 * BOHR_ANGSTROM, 1.)
            .unwrap();
        assert_eq!(r.sphere_nodes, 7);
        assert_eq!(r.occupied_nodes, count);
        assert!((r.buried_volume_angstrom3 - count as f64 * BOHR_ANGSTROM.powi(3)).abs() < 1e-14);
        assert_eq!(r.percent_buried_volume, 100. * count as f64 / 7.);
    }
}
#[test]
fn translation_and_repeatability_preserve_counts() {
    let a = grid(0., "half")
        .buried_volume([0.; 3], 1.01 * BOHR_ANGSTROM, 1.)
        .unwrap();
    let g = grid(10., "half");
    let b = g
        .buried_volume([10. * BOHR_ANGSTROM, 0., 0.], 1.01 * BOHR_ANGSTROM, 1.)
        .unwrap();
    assert_eq!(a.occupied_nodes, b.occupied_nodes);
    assert_eq!(a.sphere_nodes, b.sphere_nodes);
    let again = g
        .buried_volume([10. * BOHR_ANGSTROM, 0., 0.], 1.01 * BOHR_ANGSTROM, 1.)
        .unwrap();
    assert_eq!(
        serde_json::to_string(&b).unwrap(),
        serde_json::to_string(&again).unwrap()
    );
}
#[test]
fn unsupported_domains_and_parameters_fail_explicitly() {
    let g = grid(0., "full");
    for (c, r, t) in [
        ([0.; 3], 3. * BOHR_ANGSTROM, 1.),
        ([0.; 3], 0., 1.),
        ([0.; 3], 1., 0.),
        ([0.; 3], 1., f64::NAN),
        ([f64::INFINITY, 0., 0.], 1., 1.),
        ([0.1; 3], 1e-10, 1.),
    ] {
        assert!(g.buried_volume(c, r, t).is_err());
    }
    let text =
        "anisotropic\ndensity\n1 -2 -2 -2\n2 4 0 0\n2 0 5 0\n2 0 0 4\n1 1 0 0 0\n2 2 2 2 2 2 2 2";
    let a = DensityGrid::from_cube_str("fixture", text, g.metadata().declaration.clone()).unwrap();
    assert!(
        a.buried_volume([0.; 3], 0.5, 1.)
            .unwrap_err()
            .to_string()
            .contains("isotropic")
    );
}
