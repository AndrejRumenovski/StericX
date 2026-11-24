use steric_x::geometry::{
    AtomicSphere, CertifiedVolumeOptions, sphere_union_volume_certified_radius_profile,
    sphere_union_volume_certified_regions,
};
fn options() -> CertifiedVolumeOptions {
    CertifiedVolumeOptions {
        absolute_tolerance: 0.02,
        max_cells: 262144,
    }
}
#[test]
fn profile_preserves_order_errors_duplicates_and_standalone_certificates() {
    let atoms = [AtomicSphere {
        center: [0.25, 0., 0.5],
        radius: 1.,
    }];
    let radii = [2., 1., 2., 0., f64::NAN, 1.5];
    let p =
        sphere_union_volume_certified_radius_profile([0.; 3], &atoms, &radii, options()).unwrap();
    assert_eq!(p.entries.len(), radii.len());
    assert_eq!(p.distinct_integrations, 3);
    assert!(!p.all_targets_met);
    assert_eq!(p.entries[2].reused_from_index, Some(0));
    assert_eq!(p.entries[2].result.as_ref().unwrap().elapsed_seconds, 0.);
    assert!(p.entries[3].result.is_none() && p.entries[3].error.is_some());
    assert!(p.entries[4].error.is_some());
    let mut cells = 0;
    for (i, e) in p.entries.iter().enumerate() {
        assert_eq!(e.radius_bits, radii[i].to_bits());
        if let Some(r) = &e.result {
            let reference = sphere_union_volume_certified_regions(
                AtomicSphere {
                    center: [0.; 3],
                    radius: radii[i],
                },
                &atoms,
                options(),
            )
            .unwrap();
            for (a, b) in r.octants.iter().zip(reference.octants.iter()) {
                assert_eq!(a.lower_bound_angstrom3, b.lower_bound_angstrom3);
                assert_eq!(a.upper_bound_angstrom3, b.upper_bound_angstrom3);
            }
            if e.reused_from_index.is_none() {
                cells += r.cells_evaluated;
            }
        }
    }
    assert_eq!(p.cells_evaluated, cells);
}
#[test]
fn occupied_volume_is_monotone_but_percentage_can_decrease() {
    let atoms = [AtomicSphere {
        center: [0.; 3],
        radius: 1.,
    }];
    let p = sphere_union_volume_certified_radius_profile(
        [0.; 3],
        &atoms,
        &[0.5, 1., 2.],
        CertifiedVolumeOptions {
            absolute_tolerance: 0.001,
            ..options()
        },
    )
    .unwrap();
    assert!(p.all_targets_met);
    let mut previous = 0.;
    for e in &p.entries {
        let r = e.result.as_ref().unwrap();
        let exact = 4. * std::f64::consts::PI * e.radius_angstrom.min(1.).powi(3) / 3.;
        assert!(r.total.lower_bound_angstrom3 <= exact && exact <= r.total.upper_bound_angstrom3);
        assert!(previous <= r.total.upper_bound_angstrom3);
        previous = r.total.lower_bound_angstrom3;
    }
    assert!(
        p.entries[2]
            .result
            .as_ref()
            .unwrap()
            .total
            .percent_buried_volume
            .unwrap()
            .upper_bound
            < 13.
    );
    assert!(
        p.entries[0]
            .result
            .as_ref()
            .unwrap()
            .total
            .percent_buried_volume
            .unwrap()
            .lower_bound
            > 98.
    );
}
#[test]
fn profile_validates_common_inputs_even_with_no_radii() {
    assert!(
        sphere_union_volume_certified_radius_profile([f64::NAN, 0., 0.], &[], &[], options())
            .is_err()
    );
    let p = sphere_union_volume_certified_radius_profile([0.; 3], &[], &[], options()).unwrap();
    assert!(p.entries.is_empty());
    assert_eq!(p.distinct_integrations, 0);
}
