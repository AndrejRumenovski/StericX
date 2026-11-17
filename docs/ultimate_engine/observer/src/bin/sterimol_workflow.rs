//! Matched in-memory Sterimol workflow; reference mathematics lives elsewhere.
use glam::Vec3;
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::{fs, hint::black_box, time::Instant};
use steric_x::geometry::{Atom, Molecule, SterimolCalculator};

#[derive(Deserialize)]
struct InputAtom {
    position: [f32; 3],
    radius: f32,
}
#[derive(Deserialize)]
struct Input {
    id: String,
    atoms: Vec<InputAtom>,
    base: usize,
    dummy: [f32; 3],
}
struct Case {
    id: String,
    molecule: Molecule,
    base: usize,
    dummy: Vec3,
}
#[derive(Serialize, PartialEq)]
struct Row<'a> {
    id: &'a str,
    l: f64,
    b1: f64,
    b5: f64,
}
fn batch<'a>(cases: &'a [Case], method: &str) -> Vec<Row<'a>> {
    cases
        .iter()
        .map(|case| {
            let (l, b1, b5) = if method == "continuous" {
                let p = SterimolCalculator::compute_continuous_with_dummy(
                    &case.molecule,
                    case.base,
                    case.dummy,
                )
                .unwrap();
                (p.l, p.b1, p.b5)
            } else {
                let p = SterimolCalculator::compute_with_dummy_with_sampling(
                    &case.molecule,
                    case.base,
                    case.dummy,
                    361,
                )
                .unwrap();
                (f64::from(p.l), f64::from(p.b1), f64::from(p.b5))
            };
            Row {
                id: &case.id,
                l,
                b1,
                b5,
            }
        })
        .collect()
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input: Vec<Input> = serde_json::from_str(&fs::read_to_string(&args[1])?)?;
    let method = &args[2];
    assert!(method == "continuous" || method == "sampled361");
    let iterations: usize = args[3].parse()?;
    assert!(iterations > 0);
    let cases: Vec<Case> = input
        .into_iter()
        .map(|c| Case {
            id: c.id,
            base: c.base,
            dummy: Vec3::from_array(c.dummy),
            molecule: Molecule {
                atoms: c
                    .atoms
                    .into_iter()
                    .map(|a| Atom {
                        element: "C".into(),
                        position: Vec3::from_array(a.position),
                        vdw_radius: a.radius,
                    })
                    .collect(),
            },
        })
        .collect();
    let expected = batch(&cases, method);
    let started = Instant::now();
    for _ in 0..iterations {
        black_box(batch(&cases, method));
    }
    let elapsed_ns = started.elapsed().as_nanos();
    let rows = batch(&cases, method);
    assert!(rows == expected);
    println!(
        "{}",
        json!({"rows": rows, "elapsed_ns": elapsed_ns, "iterations": iterations})
    );
    Ok(())
}
