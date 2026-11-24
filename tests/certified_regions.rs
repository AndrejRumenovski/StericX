use steric_x::geometry::{
    AtomicSphere, CertifiedVolumeOptions, CertifiedVolumeStatus,
    sphere_union_volume_certified_regions,
};
fn s(c: [f64; 3], r: f64) -> AtomicSphere {
    AtomicSphere {
        center: c,
        radius: r,
    }
}
fn opts(t: f64, b: usize) -> CertifiedVolumeOptions {
    CertifiedVolumeOptions {
        absolute_tolerance: t,
        max_cells: b,
    }
}
fn enclosed(m: steric_x::geometry::CertifiedRegionalMeasure, x: f64) {
    assert!(m.lower_bound_angstrom3 >= 0.0);
    assert!(
        m.lower_bound_angstrom3 <= x && x <= m.upper_bound_angstrom3,
        "{m:?}, exact {x}"
    );
    assert!((m.volume_angstrom3 - x).abs() <= m.absolute_error_bound_angstrom3);
}
#[test]
fn central_ball_has_analytic_eighths_and_quarters() {
    let r = sphere_union_volume_certified_regions(
        s([0.; 3], 2.),
        &[s([0.; 3], 1.)],
        opts(0.02, 262144),
    )
    .unwrap();
    let v = 4. * std::f64::consts::PI / 3.;
    assert!(r.target_met, "{r:?}");
    enclosed(r.total, v);
    for q in r.quadrants {
        enclosed(q, v / 4.);
    }
    for o in r.octants {
        enclosed(o, v / 8.);
        let p = o.percent_buried_volume.unwrap();
        assert!(p.lower_bound <= 12.5 && 12.5 <= p.upper_bound);
    }
    assert_eq!(r.cells_evaluated, 4 + 2 * r.subdivisions);
}
#[test]
fn sphere_plane_cap_and_region_labels_are_correct() {
    let r = sphere_union_volume_certified_regions(
        s([0.; 3], 3.),
        &[s([0., 0., 0.5], 1.)],
        opts(0.02, 262144),
    )
    .unwrap();
    let pi = std::f64::consts::PI;
    // Negative cap h=1/2: pi*h²*(r-h/3)=5*pi/24. Positive remainder=9*pi/8.
    for i in 0..4 {
        enclosed(r.octants[i], 9. * pi / 32.);
        enclosed(r.octants[i + 4], 5. * pi / 96.);
    }
    assert!(r.target_met);
    let r = sphere_union_volume_certified_regions(
        s([0.; 3], 3.),
        &[s([1., 1., 1.], 0.5)],
        opts(0.01, 262144),
    )
    .unwrap();
    enclosed(r.octants[0], pi / 6.);
    for &octant in r.octants.iter().skip(1) {
        enclosed(octant, 0.);
    }
}
#[test]
fn nested_budgets_and_exact_rigid_transform_preserve_enclosures() {
    let region = s([0.; 3], 2.);
    let atoms = [s([0.25, 0.5, -0.25], 1.), s([-0.5, 0.25, 0.5], 1.)];
    let mut prev = [(0., f64::INFINITY); 8];
    for b in [4, 6, 32, 1024] {
        let r = sphere_union_volume_certified_regions(region, &atoms, opts(1e-12, b)).unwrap();
        assert!(!r.target_met);
        assert_eq!(r.status, CertifiedVolumeStatus::CellBudgetExhausted);
        assert!(r.cells_evaluated <= b);
        for (i, m) in r.octants.iter().enumerate() {
            assert!(prev[i].0 <= m.lower_bound_angstrom3 && m.upper_bound_angstrom3 <= prev[i].1);
            prev[i] = (m.lower_bound_angstrom3, m.upper_bound_angstrom3);
        }
    }
    let a = sphere_union_volume_certified_regions(region, &atoms, opts(0.05, 65536)).unwrap();
    let shift = |mut s: AtomicSphere| {
        s.center = s.center.map(|x| x + 1024.);
        s
    };
    let b = sphere_union_volume_certified_regions(
        shift(region),
        &[shift(atoms[1]), shift(atoms[0]), shift(atoms[0])],
        opts(0.05, 65536),
    )
    .unwrap();
    for (a, b) in a.octants.iter().zip(b.octants.iter()) {
        assert_eq!(a.lower_bound_angstrom3, b.lower_bound_angstrom3);
        assert_eq!(a.upper_bound_angstrom3, b.upper_bound_angstrom3);
    }
}
#[test]
fn empty_full_and_invalid_inputs_have_explicit_status() {
    let empty = sphere_union_volume_certified_regions(s([0.; 3], 1.), &[], opts(0.1, 4)).unwrap();
    enclosed(empty.total, 0.);
    for q in empty.quadrants {
        enclosed(q, 0.);
    }
    assert!(empty.target_met);
    let full = sphere_union_volume_certified_regions(
        s([0.; 3], 1.),
        &[s([0.; 3], 2.)],
        opts(0.02, 262144),
    )
    .unwrap();
    for o in full.octants {
        enclosed(o, std::f64::consts::PI / 6.);
        let p = o.percent_buried_volume.unwrap();
        assert_eq!(p.upper_bound, 100.);
    }
    assert!(sphere_union_volume_certified_regions(s([0.; 3], 1.), &[], opts(0.1, 3)).is_err());
    assert!(
        sphere_union_volume_certified_regions(s([0.; 3], f64::NAN), &[], opts(0.1, 4)).is_err()
    );
}

#[test]
fn molecular_frame_matches_existing_axes_and_declares_local_inputs() {
    use glam::Vec3;
    use steric_x::geometry::{Atom, BuriedVolumeCalculator, BuriedVolumeConfig, Molecule};
    let molecule = Molecule {
        atoms: vec![
            Atom::new("P", Vec3::ZERO),
            Atom::new("C", Vec3::new(1., 0., 1.)),
            Atom::new("C", Vec3::new(-0.5, 1., 1.)),
            Atom::new("H", Vec3::new(-0.5, -1., 1.)),
        ],
    };
    let center = Vec3::new(0., 0., -2.1);
    let config = BuriedVolumeConfig::default();
    let detail = BuriedVolumeCalculator::compute_detailed_with_center_and_neighbors(
        &molecule,
        0,
        [1, 2, 3],
        center,
        config,
    )
    .unwrap();
    for orientation in detail.orientations {
        let frame = BuriedVolumeCalculator::prepare_certified_frame(
            &molecule,
            0,
            orientation.plane_atom_index,
            center,
            f64::from(config.radii_scale),
            false,
        )
        .unwrap();
        assert_eq!(frame.axes, orientation.axes);
        assert_eq!(frame.atoms.len(), 3);
        assert!(frame.atoms[0].center[2] < 0.);
        assert_eq!(
            frame.atoms[0].radius,
            f64::from(molecule.atoms[0].vdw_radius) * f64::from(config.radii_scale)
        );
        let r = frame.compute_regions(3.5, opts(0.2, 65536)).unwrap();
        assert!(r.target_met);
    }
    assert!(
        BuriedVolumeCalculator::prepare_certified_frame(&molecule, 10, 1, center, 1.17, false)
            .is_err()
    );
}
