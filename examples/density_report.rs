//! Thin JSON-lines driver for density validation and reproducible research.
use serde::Deserialize;
use std::io::{self, BufRead};
use steric_x::density::{DensityDeclaration, DensityGrid};
#[derive(Deserialize)]
struct Request {
    path: String,
    declaration: DensityDeclaration,
    center_angstrom: [f64; 3],
    radius_angstrom: f64,
    threshold_electrons_per_bohr3: f64,
}
fn main() {
    for line in io::stdin().lock().lines() {
        let result = line.map_err(|e| e.to_string()).and_then(|line| {
            let r: Request = serde_json::from_str(&line).map_err(|e| e.to_string())?;
            let grid =
                DensityGrid::from_cube_file(&r.path, r.declaration).map_err(|e| e.to_string())?;
            grid.buried_volume(
                r.center_angstrom,
                r.radius_angstrom,
                r.threshold_electrons_per_bohr3,
            )
            .map_err(|e| e.to_string())
        });
        let output = match result {
            Ok(value) => serde_json::json!({"status":"success","value":value}),
            Err(error) => serde_json::json!({"status":"error","message":error}),
        };
        println!(
            "{}",
            serde_json::to_string(&output).expect("finite validated result")
        );
    }
}
