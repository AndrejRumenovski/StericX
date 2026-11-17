use steric_x::geometry::{
    AtomicSphere, CertifiedVolumeOptions, CertifiedVolumeStatus, sphere_union_volume_certified,
};
fn sphere(center: [f64; 3], radius: f64) -> AtomicSphere {
    AtomicSphere { center, radius }
}
fn options(tolerance: f64, max_cells: usize) -> CertifiedVolumeOptions {
    CertifiedVolumeOptions {
        absolute_tolerance: tolerance,
        max_cells,
    }
}

#[test]
fn independent_analytic_volumes_and_zero_measure_contacts_are_enclosed() {
    let region = sphere([0.0; 3], 3.0);
    let cases = [
        (vec![], 0.0),
        (
            vec![sphere([0.0; 3], 1.0)],
            4.0 * std::f64::consts::PI / 3.0,
        ),
        (vec![sphere([4.0, 0.0, 0.0], 1.0)], 0.0),
        (vec![sphere([0.0; 3], 4.0)], 36.0 * std::f64::consts::PI),
        (
            vec![sphere([-0.5, 0.0, 0.0], 1.0), sphere([0.5, 0.0, 0.0], 1.0)],
            9.0 * std::f64::consts::PI / 4.0,
        ),
    ];
    for (atoms, reference) in cases {
        let r = sphere_union_volume_certified(region, &atoms, options(1e-10, 32767)).unwrap();
        assert!(
            r.lower_bound_angstrom3 <= reference && reference <= r.upper_bound_angstrom3,
            "{r:?}"
        );
        assert!(r.target_met);
    }
}

#[test]
fn many_sphere_enclosures_nest_and_budget_failure_is_explicit() {
    // Triple intersections are empty: pair-only inclusion-exclusion is exact.
    let region = sphere([0.0; 3], 4.0);
    let atoms = [-1.5, 0.0, 1.5].map(|x| sphere([x, 0.0, 0.0], 1.0));
    let reference = (4.0 - 2.0 * 11.0 / 96.0) * std::f64::consts::PI;
    let mut previous = (0.0, f64::INFINITY);
    for budget in [1, 2, 31, 1023, 16383] {
        let r = sphere_union_volume_certified(region, &atoms, options(1e-12, budget)).unwrap();
        assert!(!r.target_met);
        assert_eq!(r.status, CertifiedVolumeStatus::CellBudgetExhausted);
        assert!(r.cells_evaluated <= budget);
        assert_eq!(r.cells_evaluated, 1 + 2 * r.subdivisions);
        assert!(previous.0 <= r.lower_bound_angstrom3 && r.upper_bound_angstrom3 <= previous.1);
        assert!(r.lower_bound_angstrom3 <= reference && reference <= r.upper_bound_angstrom3);
        previous = (r.lower_bound_angstrom3, r.upper_bound_angstrom3);
    }
    let r = sphere_union_volume_certified(region, &atoms, options(0.05, 262143)).unwrap();
    assert!(r.target_met && r.absolute_error_bound_angstrom3 <= 0.05);
    assert_eq!(r.status, CertifiedVolumeStatus::ToleranceAchieved);
    assert!((r.volume_angstrom3 - reference).abs() <= r.absolute_error_bound_angstrom3);
    assert!(
        r.relative_error_bound.unwrap()
            >= r.absolute_error_bound_angstrom3 / r.lower_bound_angstrom3
    );
}

#[test]
fn permutation_duplicates_and_exact_translation_preserve_bounds() {
    let atoms = [
        sphere([0.0, 0.25, -0.5], 1.0),
        sphere([1.5, -0.25, 0.0], 1.0),
        sphere([-1.5, 0.25, 0.0], 1.0),
    ];
    let opts = options(0.05, 262143);
    let a = sphere_union_volume_certified(sphere([0.0; 3], 4.0), &atoms, opts).unwrap();
    let b = sphere_union_volume_certified(
        sphere([0.0; 3], 4.0),
        &[atoms[2], atoms[0], atoms[1], atoms[0]],
        opts,
    )
    .unwrap();
    let shift = |s: AtomicSphere| sphere(s.center.map(|x| x + 1024.0), s.radius);
    let c = sphere_union_volume_certified(shift(sphere([0.0; 3], 4.0)), &atoms.map(shift), opts)
        .unwrap();
    for r in [b, c] {
        assert_eq!(
            a.lower_bound_angstrom3.to_bits(),
            r.lower_bound_angstrom3.to_bits()
        );
        assert_eq!(
            a.upper_bound_angstrom3.to_bits(),
            r.upper_bound_angstrom3.to_bits()
        );
        assert_eq!(a.cells_evaluated, r.cells_evaluated);
    }
}

#[test]
fn numerical_limits_and_invalid_inputs_do_not_claim_false_success() {
    for scale in [f64::from_bits(1), 1e-150, 1e-100, 1e100] {
        let r = sphere_union_volume_certified(
            sphere([0.0; 3], scale),
            &[sphere([0.0; 3], scale)],
            options(f64::from_bits(1), 1),
        )
        .unwrap();
        assert!(r.lower_bound_angstrom3.is_finite() && r.upper_bound_angstrom3.is_finite());
        assert!(
            r.lower_bound_angstrom3 >= 0.0 && r.lower_bound_angstrom3 <= r.upper_bound_angstrom3
        );
        assert_eq!(
            r.target_met,
            r.absolute_error_bound_angstrom3 <= f64::from_bits(1)
        );
    }
    let region = sphere([0.0; 3], 1.0);
    for opts in [
        options(0.0, 1),
        options(f64::NAN, 1),
        options(f64::INFINITY, 1),
        options(0.05, 0),
    ] {
        assert!(sphere_union_volume_certified(region, &[], opts).is_err());
    }
    for invalid in [
        sphere([f64::NAN, 0.0, 0.0], 1.0),
        sphere([0.0; 3], -1.0),
        sphere([0.0; 3], f64::INFINITY),
    ] {
        assert!(sphere_union_volume_certified(region, &[invalid], options(0.05, 1)).is_err());
    }
}

#[test]
fn molecular_wrapper_preserves_configured_radii_center_and_hydrogen_policy() {
    use glam::Vec3;
    use steric_x::geometry::{Atom, BuriedVolumeCalculator, Molecule};
    let molecule = Molecule {
        atoms: vec![
            Atom::new("C", Vec3::new(0.25, -0.5, 0.0)),
            Atom::new("H", Vec3::new(2.25, 0.0, 0.0)),
        ],
    };
    let center = [0.125, -0.25, 0.5];
    let scale = f64::from(1.17_f32);
    for include_hydrogens in [false, true] {
        let expected: Vec<_> = molecule
            .atoms
            .iter()
            .filter(|a| include_hydrogens || a.element != "H")
            .map(|a| {
                sphere(
                    a.position.as_dvec3().to_array(),
                    f64::from(a.vdw_radius) * scale,
                )
            })
            .collect();
        let explicit =
            sphere_union_volume_certified(sphere(center, 3.5), &expected, options(0.05, 262143))
                .unwrap();
        let wrapper = BuriedVolumeCalculator::compute_certified_with_center(
            &molecule,
            center,
            3.5,
            scale,
            include_hydrogens,
            options(0.05, 262143),
        )
        .unwrap();
        assert_eq!(
            explicit.lower_bound_angstrom3.to_bits(),
            wrapper.lower_bound_angstrom3.to_bits()
        );
        assert_eq!(
            explicit.upper_bound_angstrom3.to_bits(),
            wrapper.upper_bound_angstrom3.to_bits()
        );
        assert_eq!(explicit.cells_evaluated, wrapper.cells_evaluated);
        assert!(wrapper.target_met);
    }
    assert!(
        BuriedVolumeCalculator::compute_certified_with_center(
            &molecule,
            center,
            3.5,
            0.0,
            false,
            options(0.05, 1),
        )
        .is_err()
    );
}
