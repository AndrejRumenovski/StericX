//! Measurement adapter calling public production APIs, not a reference solver.
use glam::Vec3;
use serde_json::{Value, json};
use std::io::{self, BufRead, Write};
use steric_x::{Atom, Molecule, SterimolCalculator};

fn vector(v: &Value) -> Vec3 {
    Vec3::new(
        v[0].as_f64().unwrap() as f32,
        v[1].as_f64().unwrap() as f32,
        v[2].as_f64().unwrap() as f32,
    )
}

fn main() {
    let mut out = io::BufWriter::new(io::stdout().lock());
    for line in io::stdin().lock().lines() {
        let r: Value = serde_json::from_str(&line.unwrap()).unwrap();
        let molecule = Molecule {
            atoms: r["atoms"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| {
                    let mut a = Atom::new(v["element"].as_str().unwrap(), vector(&v["position"]));
                    if let Some(radius) = v["radius"].as_f64() {
                        a.vdw_radius = radius as f32;
                    }
                    a
                })
                .collect(),
        };
        let donor = r["donor"].as_u64().unwrap() as usize;
        let center = vector(&r["center"]);
        let result = if let Some(count) = r["n_rot_vectors"].as_u64() {
            SterimolCalculator::compute_with_dummy_with_sampling(
                &molecule,
                donor,
                center,
                count as usize,
            )
        } else {
            SterimolCalculator::compute_with_dummy(&molecule, donor, center)
        };
        let value = match result {
            Ok(p) => json!({"id": r["id"], "l": p.l + 0.4_f32, "b1": p.b1, "b5": p.b5,
                "bits": [(p.l + 0.4_f32).to_bits(), p.b1.to_bits(), p.b5.to_bits()]}),
            Err(e) => json!({"id": r["id"], "error": e.to_string()}),
        };
        writeln!(out, "{value}").unwrap();
        out.flush().unwrap();
    }
}
