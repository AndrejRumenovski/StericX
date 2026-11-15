//! Source-bound parser observation, not a numerical reference.
use std::{error::Error, fs};
use steric_x::density::{DensityDeclaration, DensityField, DensityGrid, DensityUnit};

fn main() -> Result<(), Box<dyn Error>> {
    let path = std::env::args().nth(1).ok_or("cube path required")?;
    let declaration = DensityDeclaration {
        field: DensityField::TotalElectronDensity,
        unit: DensityUnit::ElectronsPerBohrCubed,
        provenance: "parser equivalence fixture".into(),
    };
    let result = match DensityGrid::from_cube_str("frozen.cube", &fs::read_to_string(path)?, declaration) {
        Ok(grid) => serde_json::json!({"metadata":grid.metadata(),"values":grid.values()}),
        Err(error) => serde_json::json!({"error":error.to_string()}),
    };
    println!("{result}");
    Ok(())
}
