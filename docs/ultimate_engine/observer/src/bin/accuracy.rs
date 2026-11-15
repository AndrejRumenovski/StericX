//! Actual public finite-grid API on declared analytical geometries.
use glam::Vec3;
use serde::Deserialize;
use serde_json::json;
use std::{error::Error, fs, hint::black_box, time::Instant};
use steric_x::geometry::{Atom, BuriedVolumeCalculator, BuriedVolumeConfig, Molecule};

#[derive(Deserialize)]
struct Case {
    id: String,
    radius: f32,
    centers: Vec<[f32; 3]>,
    radii: Vec<f32>,
}

fn main() -> Result<(), Box<dyn Error>> {
    let path = std::env::args().nth(1).ok_or("input cases required")?;
    let cases: Vec<Case> = serde_json::from_str(&fs::read_to_string(path)?)?;
    let mut rows = Vec::new();
    for case in cases {
        assert_eq!(case.centers.len(), case.radii.len());
        // Four excluded H atoms supply a valid frame. Only explicit carbon
        // spheres contribute occupancy. All center and radius choices are input.
        let mut atoms = vec![
            Atom::new("H", Vec3::new(0., 0., -2.)),
            Atom::new("H", Vec3::new(1., 0., -2.)),
            Atom::new("H", Vec3::new(0., 1., -2.)),
            Atom::new("H", Vec3::new(-1., -1., -2.)),
        ];
        for (center, radius) in case.centers.iter().zip(&case.radii) {
            atoms.push(Atom {
                element: "C".into(),
                position: Vec3::from_array(*center),
                vdw_radius: *radius,
            });
        }
        let molecule = Molecule { atoms };
        for density in [0.08, 0.01, 0.00125, 0.00015625] {
            let config = BuriedVolumeConfig {
                sphere_radius: case.radius,
                density,
                radii_scale: 1.,
                include_hydrogens: false,
                ..Default::default()
            };
            let run = || {
                BuriedVolumeCalculator::compute_detailed_with_center_and_neighbors(
                    &molecule,
                    0,
                    [1, 2, 3],
                    Vec3::ZERO,
                    config,
                )
                .unwrap()
            };
            black_box(run());
            let mut times = Vec::new();
            for _ in 0..5 {
                let start = Instant::now();
                black_box(run());
                times.push(start.elapsed().as_nanos());
            }
            rows.push(
                json!({"id":case.id,"radius":case.radius,"centers":case.centers,
                "radii":case.radii,"density":density,"times_ns":times,"result":run()}),
            );
        }
    }
    println!("{}", json!(rows));
    Ok(())
}
