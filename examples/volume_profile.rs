//! Reproducible fixed-center profiles; all scientific work stays in the library.
use serde::Deserialize;
use std::io::{self, BufRead};
use steric_x::{
    density::{DensityDeclaration, DensityGrid},
    ensemble::{DescriptorOptions, Ensemble},
};
#[derive(Deserialize)]
#[serde(tag = "mode", rename_all = "snake_case")]
enum Request {
    HardSphere {
        path: String,
        options: DescriptorOptions,
        radii_angstrom: Vec<f64>,
    },
    Density {
        path: String,
        declaration: DensityDeclaration,
        center_angstrom: [f64; 3],
        radii_angstrom: Vec<f64>,
        threshold_electrons_per_bohr3: f64,
    },
}
fn run(request: Request) -> Result<serde_json::Value, String> {
    match request {
        Request::HardSphere {
            path,
            options,
            radii_angstrom,
        } => {
            let ensemble = Ensemble::from_file(path)?;
            let profiles = ensemble
                .conformers
                .iter()
                .map(|c| c.buried_volume_scan(&radii_angstrom, &options))
                .collect::<Result<Vec<_>, _>>()?;
            serde_json::to_value(profiles).map_err(|e| e.to_string())
        }
        Request::Density {
            path,
            declaration,
            center_angstrom,
            radii_angstrom,
            threshold_electrons_per_bohr3,
        } => {
            let grid = DensityGrid::from_cube_file(path, declaration).map_err(|e| e.to_string())?;
            let profile = grid
                .buried_volume_scan(
                    center_angstrom,
                    &radii_angstrom,
                    threshold_electrons_per_bohr3,
                )
                .map_err(|e| e.to_string())?;
            serde_json::to_value(profile).map_err(|e| e.to_string())
        }
    }
}
fn main() {
    for line in io::stdin().lock().lines() {
        let result = line
            .map_err(|e| e.to_string())
            .and_then(|s| serde_json::from_str(&s).map_err(|e| e.to_string()))
            .and_then(run);
        println!(
            "{}",
            match result {
                Ok(v) => serde_json::json!({"status":"success","value":v}),
                Err(e) => serde_json::json!({"status":"error","message":e}),
            }
        );
    }
}
