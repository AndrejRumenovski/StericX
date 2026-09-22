//! Independent audit witnesses for the geometry corrections, not golden SUT outputs.
use glam::Vec3;
use serde_json::Value;
use steric_x::geometry::{coordination_center_with_neighbors, parse_sdf_with_topology};
use steric_x::{
    Atom, BuriedVolumeCalculator as BV, BuriedVolumeConfig, Molecule,
    PyramidalizationCalculator as Pyr, SterimolCalculator as Sterimol,
};

fn molecule(request: &Value) -> Molecule {
    Molecule {
        atoms: request["atoms"]
            .as_array()
            .unwrap()
            .iter()
            .map(|a| {
                let p = &a["position"];
                let mut atom = Atom::new(
                    a["element"].as_str().unwrap(),
                    Vec3::new(
                        p[0].as_f64().unwrap() as f32,
                        p[1].as_f64().unwrap() as f32,
                        p[2].as_f64().unwrap() as f32,
                    ),
                );
                if let Some(r) = a["radius"].as_f64() {
                    atom.vdw_radius = r as f32;
                }
                atom
            })
            .collect(),
    }
}

fn request(id: &str) -> Value {
    include_str!("../docs/scientific_remediation/geometry/focused_v1/requests.jsonl")
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .find(|r| r["id"] == id)
        .unwrap()
}

#[test]
fn audited_near_axis_supports_match_direct_geometry_and_exact_rotation() {
    for id in [
        "KRAKEN:963:48394__controlled",
        "KRAKEN:963:48394__minimal",
        "KRAKEN:1075:49864__controlled",
        "KRAKEN:1075:49864__minimal",
    ] {
        let r = request(id);
        let m = molecule(&r);
        let d = r["donor"].as_u64().unwrap() as usize;
        let c = &r["center"];
        let center = Vec3::new(
            c[0].as_f64().unwrap() as f32,
            c[1].as_f64().unwrap() as f32,
            c[2].as_f64().unwrap() as f32,
        );
        let got = Sterimol::compute_with_dummy(&m, d, center).unwrap();
        // Independent scalar projection/perpendicular distance; no rotation arc.
        let axis = (m.atoms[d].position.as_dvec3() - center.as_dvec3()).normalize();
        let mut expected_l = f64::NEG_INFINITY;
        let mut expected_b5 = 0.0_f64;
        for atom in &m.atoms {
            let displacement = atom.position.as_dvec3() - center.as_dvec3();
            let height = displacement.dot(axis);
            expected_l = expected_l.max(height + f64::from(atom.vdw_radius));
            expected_b5 = expected_b5
                .max((displacement - height * axis).length() + f64::from(atom.vdw_radius));
        }
        // Fixed float-rounding budget, far below the frozen 0.001–0.004 Å defects.
        assert!((f64::from(got.l) - expected_l).abs() < 2.0e-6, "{id}");
        assert!((f64::from(got.b5) - expected_b5).abs() < 2.0e-6, "{id}");
        let rotate = |p: Vec3| Vec3::new(p.x, -p.z, p.y);
        let rotated = Molecule {
            atoms: m
                .atoms
                .iter()
                .cloned()
                .map(|mut a| {
                    a.position = rotate(a.position);
                    a
                })
                .collect(),
        };
        let other = Sterimol::compute_with_dummy(&rotated, d, rotate(center)).unwrap();
        assert_eq!(got.l, other.l);
        assert_eq!(got.b5, other.b5);
        // B1's retained one-degree azimuthal scan has an explicit discretization
        // limit; this test does not pretend its phase is continuously invariant.
    }
}

#[test]
fn nearly_collinear_pyramidalization_matches_frozen_independent_equations() {
    let r = request("nearly_collinear");
    let got = Pyr::compute(&molecule(&r), 0, [1, 2, 3]).unwrap();
    // Values predate this correction: sealed geometry/reference_results.jsonl.
    assert!((f64::from(got.pyr_alpha) - 29.999_924_973_640_23).abs() < 2.0e-6);
    assert!((f64::from(got.pyr_p) - 3.086_419_760_576_517_5e-11).abs() < 2.0e-18);
}

#[test]
fn short_nonzero_axes_and_finite_extreme_coordinates_do_not_become_zero_or_inf() {
    for id in [
        "short_axis_0.0003452669479884207",
        "short_axis_0.0003452669770922512",
        "short_axis_0.00034526700619608164",
    ] {
        let r = request(id);
        let got = Sterimol::compute(&molecule(&r), 0, 1).unwrap();
        assert!(got.l > 1.7 && got.b5 > 3.0);
    }
    let r = request("huge_finite_coordinate");
    let got = Sterimol::compute(&molecule(&r), 0, 1).unwrap();
    assert!(got.b5.is_finite() && got.b5 > 1.0e37);
    let mut m = Molecule {
        atoms: vec![
            Atom::new("H", Vec3::ZERO),
            Atom::new("C", Vec3::Z * f32::from_bits(1)),
        ],
    };
    assert!(Sterimol::compute(&m, 0, 1).unwrap().l > 0.0);
    m.atoms[1].position = Vec3::ZERO;
    assert!(Sterimol::compute(&m, 0, 1).is_err());
    m.atoms[1].position = Vec3::splat(f32::MAX);
    assert!(Sterimol::compute(&m, 0, 1).is_err());
}

#[test]
fn valid_symmetric_volume_and_permuted_triad_are_accepted() {
    let r = request("phosphine_PH3");
    let m = molecule(&r);
    let config = BuriedVolumeConfig {
        center_distance: 2.28,
        ..Default::default()
    };
    let base = BV::compute_with_neighbors(&m, 0, [1, 2, 3], config).unwrap();
    assert!((f64::from(base.buried_volume) - 31.937_214_517_315_244).abs() < 1.0e-5);
    assert_eq!(base.max_delta_qvbur, 0.0);
    for order in [[1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]] {
        assert_eq!(
            base,
            BV::compute_with_neighbors(&m, 0, order, config).unwrap()
        );
    }
}

#[test]
fn planar_sign_ties_require_an_explicit_center() {
    let m = Molecule {
        atoms: vec![
            Atom::new("P", Vec3::ZERO),
            Atom::new("C", Vec3::X),
            Atom::new("C", Vec3::new(-0.5, 0.866_025_4, 0.0)),
            Atom::new("C", Vec3::new(-0.5, -0.866_025_4, 0.0)),
        ],
    };
    let config = BuriedVolumeConfig::default();
    for n in [[1, 2, 3], [3, 2, 1], [2, 1, 3]] {
        assert!(
            coordination_center_with_neighbors(&m, 0, n, config)
                .unwrap_err()
                .to_string()
                .contains("ambiguous")
        );
        assert!(
            BV::compute_with_center_and_neighbors(&m, 0, n, Vec3::new(0.0, 0.0, -2.1), config)
                .is_ok()
        );
    }
}

#[test]
fn rounded_planar_rotations_remain_ambiguous_at_input_precision() {
    let mut checked = 0;
    for line in
        include_str!("../docs/scientific_remediation/geometry/planar_precision/requests.jsonl")
            .lines()
    {
        let r: Value = serde_json::from_str(line).unwrap();
        let m = molecule(&r);
        let donor = r["donor"].as_u64().unwrap() as usize;
        let neighbors: [usize; 3] = r["neighbors"]
            .as_array()
            .unwrap()
            .iter()
            .map(|n| n.as_u64().unwrap() as usize)
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        let error =
            coordination_center_with_neighbors(&m, donor, neighbors, BuriedVolumeConfig::default());
        assert!(
            error.is_err(),
            "{} falsely selected a planar sign: {error:?}",
            r["id"]
        );
        checked += 1;
    }
    assert_eq!(checked, 124);
    // A clearly resolved nearby obstruction remains usable. This tests the
    // clearance decision, not just rejecting every planar configuration.
    let mut m = Molecule {
        atoms: vec![
            Atom::new("P", Vec3::ZERO),
            Atom::new("C", Vec3::new(1.8, 0.0, 0.0)),
            Atom::new("C", Vec3::new(-0.9, 1.558_845_8, 0.0)),
            Atom::new("C", Vec3::new(-0.9, -1.558_845_8, 0.0)),
        ],
    };
    m.atoms.push(Atom::new("C", Vec3::new(0.0, 0.0, 2.5)));
    let center =
        coordination_center_with_neighbors(&m, 0, [1, 2, 3], BuriedVolumeConfig::default())
            .unwrap();
    assert!(center.z < 0.0);
}

#[test]
fn explicit_sdf_bonds_do_not_include_a_close_nonbonded_contact() {
    let sdf = "contact\nstericx\n\n  5  3  0  0  0  0            999 V2000\n0 0 0 P\n1.5 0 0.6 C\n-0.7 1.2 0.6 C\n-0.7 -1.2 0.6 C\n0 0 2.5 P\n  1  2  1  0\n  1  3  1  0\n  1  4  1  0\nM  END\n$$$$\n";
    let frames = parse_sdf_with_topology(sdf).unwrap();
    let frame = &frames[0];
    assert_eq!(frame.bonded_neighbors(0).unwrap(), vec![1, 2, 3]);
    assert_eq!(steric_x::bonded_neighbors(&frame.molecule, 0).len(), 4);
    assert!(BV::compute(&frame.molecule, 0, 1, Default::default()).is_err());
    assert!(BV::compute_with_neighbors(&frame.molecule, 0, [1, 2, 3], Default::default()).is_ok());
    assert!(parse_sdf_with_topology(&sdf.replace("  1  4  1  0", "  1  6  1  0")).is_err());
    assert!(parse_sdf_with_topology(&sdf.replace("  1  4  1  0", "  1  2  1  0")).is_err());
}

#[test]
fn all_eight_audited_topology_failures_use_the_independently_verified_graph() {
    let mut topology_failures = 0;
    for line in
        include_str!("../docs/scientific_remediation/geometry/topology_witnesses.jsonl").lines()
    {
        let r: Value = serde_json::from_str(line).unwrap();
        let id = r["id"].as_str().unwrap();
        let conformer = id.rsplit(':').next().unwrap();
        let file = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(format!(
            "docs/scientific_accuracy_audit/kraken/primary/conformers/{conformer}.body"
        ));
        let frames = parse_sdf_with_topology(&std::fs::read_to_string(file).unwrap()).unwrap();
        let frame = &frames[0];
        let donor = r["donor"].as_u64().unwrap() as usize;
        let expected: Vec<usize> = r["neighbors"]
            .as_array()
            .unwrap()
            .iter()
            .map(|n| n.as_u64().unwrap() as usize)
            .collect();
        assert_eq!(frame.bonded_neighbors(donor).unwrap(), expected, "{id}");
        let config = BuriedVolumeConfig {
            center_distance: 2.28,
            ..Default::default()
        };
        let neighbors: [usize; 3] = expected.try_into().unwrap();
        let got = BV::compute_with_neighbors(&frame.molecule, donor, neighbors, config).unwrap();
        assert!(
            got.buried_volume > 0.0 && got.buried_volume.is_finite(),
            "{id}"
        );
        if id != "KRAKEN:1299:54318" {
            assert_eq!(
                steric_x::bonded_neighbors(&frame.molecule, donor).len(),
                4,
                "{id}"
            );
            assert!(BV::compute(&frame.molecule, donor, neighbors[0], config).is_err());
            topology_failures += 1;
        } else {
            assert_eq!(got.max_delta_qvbur, 0.0);
            assert!(BV::compute(&frame.molecule, donor, neighbors[0], config).is_ok());
        }
    }
    assert_eq!(topology_failures, 8);
}
