use glam::Vec3;
use steric_x::{Atom, Molecule, SterimolCalculator};

#[test]
fn explicit_sampling_bounds_the_independent_three_sphere_envelope() {
    // Equal spheres at 0 and +/-v have exact transverse minimum equal to
    // their radius: choose a direction perpendicular to v. No Kraken data.
    let angle = 0.5_f32.to_radians();
    let v = Vec3::new(3.0 * angle.cos(), 3.0 * angle.sin(), 0.0);
    let molecule = Molecule {
        atoms: [Vec3::ZERO, v, -v]
            .into_iter()
            .map(|position| Atom::new("C", position))
            .collect(),
    };
    let dummy = Vec3::new(0.0, 0.0, -2.28);
    let baseline = SterimolCalculator::compute_with_dummy(&molecule, 0, dummy).unwrap();
    let sampled =
        SterimolCalculator::compute_with_dummy_with_sampling(&molecule, 0, dummy, 3600).unwrap();
    let exact = molecule.atoms[0].vdw_radius;
    let roundoff = 64.0 * f32::EPSILON * 3.0;
    let bound = 3.0 * std::f32::consts::PI / 3599.0 + roundoff;
    assert!(sampled.b1 >= exact - roundoff);
    assert!(sampled.b1 - exact <= bound);
    assert!(baseline.b1 - exact > 0.02);
    assert_eq!(sampled.l.to_bits(), baseline.l.to_bits());
    assert_eq!(sampled.b5.to_bits(), baseline.b5.to_bits());
}

#[test]
fn optional_sampling_validates_count_and_geometry() {
    let molecule = Molecule {
        atoms: vec![Atom::new("N", Vec3::ZERO)],
    };
    let dummy = Vec3::new(0.0, 0.0, -2.28);
    for count in [0, 1] {
        assert!(
            SterimolCalculator::compute_with_dummy_with_sampling(&molecule, 0, dummy, count)
                .is_err()
        );
    }
    for count in [2, 360, 3600] {
        let result =
            SterimolCalculator::compute_with_dummy_with_sampling(&molecule, 0, dummy, count)
                .unwrap();
        assert_eq!(result.b1, molecule.atoms[0].vdw_radius);
        assert_eq!(result.b5, molecule.atoms[0].vdw_radius);
    }
    assert!(
        SterimolCalculator::compute_with_dummy_with_sampling(&molecule, 1, dummy, 3600).is_err()
    );
    assert!(
        SterimolCalculator::compute_with_dummy_with_sampling(&molecule, 0, Vec3::ZERO, 3600)
            .is_err()
    );
}

#[test]
fn sampled_envelope_preserves_atom_order_and_bounds_rigid_transform_error() {
    let angle = 0.5_f32.to_radians();
    let v = Vec3::new(3.0 * angle.cos(), 3.0 * angle.sin(), 0.0);
    for i in 0..12 {
        let rotation =
            glam::Quat::from_axis_angle(Vec3::new(1.0, 2.0, 3.0).normalize(), i as f32 * 0.431);
        let shift = Vec3::new(17.0, -23.0, 11.0) * i as f32;
        let positions = [Vec3::ZERO, v, -v].map(|p| rotation * p + shift);
        let dummy = rotation * Vec3::new(0.0, 0.0, -2.28) + shift;
        let mut reference = None;
        for order in [[0, 1, 2], [2, 0, 1], [1, 2, 0]] {
            let molecule = Molecule {
                atoms: order
                    .iter()
                    .map(|&j| Atom::new("C", positions[j]))
                    .collect(),
            };
            let donor = order.iter().position(|&j| j == 0).unwrap();
            let p =
                SterimolCalculator::compute_with_dummy_with_sampling(&molecule, donor, dummy, 3600)
                    .unwrap();
            let bits = [p.l.to_bits(), p.b1.to_bits(), p.b5.to_bits()];
            if let Some(previous) = reference {
                assert_eq!(bits, previous);
            }
            reference = Some(bits);
            let roundoff = 64.0 * f32::EPSILON * (shift.length() + 5.0);
            assert!(p.b1 >= 1.7 - roundoff);
            assert!(p.b1 <= 1.7 + 3.0 * std::f32::consts::PI / 3599.0 + roundoff);
            assert!((p.l - 3.98).abs() <= roundoff);
            assert!((p.b5 - 4.7).abs() <= roundoff);
        }
    }
}
