use steric_x::geometry::{
    AtomicSphere, CertifiedVolumeOptions, sphere_union_volume_certified,
    sphere_union_volume_certified_geometric,
};
fn s(x: f64, y: f64, z: f64, radius: f64) -> AtomicSphere {
    AtomicSphere {
        center: [x, y, z],
        radius,
    }
}
fn o(t: f64, budget: usize) -> CertifiedVolumeOptions {
    CertifiedVolumeOptions {
        absolute_tolerance: t,
        max_cells: budget,
    }
}

#[test]
fn contained_triangle_free_unions_include_tangent_pairs_and_exact_references() {
    let region = s(0., 0., 0., 4.);
    for (atoms, exact) in [
        (
            vec![s(-1., 0., 0., 1.), s(0., 0., 0., 1.), s(1., 0., 0., 1.)],
            19. * std::f64::consts::PI / 6.,
        ),
        (
            vec![s(-2., 0., 0., 0.5), s(0., 0., 0., 0.5), s(2., 0., 0., 0.5)],
            std::f64::consts::PI / 2.,
        ),
    ] {
        let r = sphere_union_volume_certified_geometric(region, &atoms, o(1e-5, 1)).unwrap();
        assert!(r.target_met && r.analytical && r.cells_evaluated == 0);
        assert!(r.lower_bound_angstrom3 <= exact && exact <= r.upper_bound_angstrom3);
        let mut reversed = atoms.clone();
        reversed.reverse();
        reversed.push(atoms[0]);
        let b = sphere_union_volume_certified_geometric(region, &reversed, o(1e-5, 1)).unwrap();
        assert_eq!(r.lower_bound_angstrom3, b.lower_bound_angstrom3);
        assert_eq!(r.upper_bound_angstrom3, b.upper_bound_angstrom3);
    }
}

#[test]
fn possible_triples_and_clipped_balls_use_original_certified_fallback() {
    for (region, atoms) in [
        (
            s(0., 0., 0., 4.),
            vec![
                s(-0.25, 0., 0., 1.),
                s(0., 0.25, 0., 1.),
                s(0.25, 0., 0., 1.),
            ],
        ),
        (
            s(0., 0., 0., 1.),
            vec![s(-0.5, 0., 0., 1.), s(0., 0.5, 0., 1.), s(0.5, 0., 0., 1.)],
        ),
    ] {
        let opts = o(1e-6, 1023);
        let a = sphere_union_volume_certified(region, &atoms, opts).unwrap();
        let b = sphere_union_volume_certified_geometric(region, &atoms, opts).unwrap();
        assert!(!b.analytical);
        assert_eq!(a.lower_bound_angstrom3, b.lower_bound_angstrom3);
        assert_eq!(a.upper_bound_angstrom3, b.upper_bound_angstrom3);
        assert_eq!(a.cells_evaluated, b.cells_evaluated);
    }
}

#[test]
fn subnormal_absolute_target_does_not_require_a_false_zero_volume() {
    let u = f64::from_bits(1);
    for r in [u, 1e-200, 1e-150, f64::from_bits((1023_u64 - 359) << 52)] {
        let result =
            sphere_union_volume_certified_geometric(s(0., 0., 0., r), &[s(r, 0., 0., r)], o(u, 1))
                .unwrap();
        assert!(result.target_met);
        assert_eq!(result.lower_bound_angstrom3, 0.);
        assert_eq!(result.upper_bound_angstrom3, u);
        assert_eq!(result.absolute_error_bound_angstrom3, u);
        assert_eq!(result.relative_error_bound, None);
        let tangent = sphere_union_volume_certified_geometric(
            s(0., 0., 0., r),
            &[s(2. * r, 0., 0., r)],
            o(u, 1),
        )
        .unwrap();
        assert_eq!(tangent.upper_bound_angstrom3, 0.);
    }
}

#[test]
fn invalid_inputs_are_rejected_before_any_analytical_shortcut() {
    for options in [o(0., 1), o(f64::NAN, 1), o(0.05, 0)] {
        assert!(sphere_union_volume_certified_geometric(s(0., 0., 0., 1.), &[], options).is_err());
    }
    assert!(sphere_union_volume_certified_geometric(s(0., 0., 0., 0.), &[], o(0.05, 1)).is_err());
    assert!(
        sphere_union_volume_certified_geometric(
            s(0., 0., 0., 1.),
            &[s(f64::NAN, 0., 0., 1.)],
            o(0.05, 1)
        )
        .is_err()
    );
}
