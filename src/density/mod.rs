//! Scalar total-electron-density grids, separate from atomic hard spheres.
//!
//! Cube files do not identify their physical field or units reliably. Loading
//! requires an explicit declaration; no density is inferred from geometry.

use serde::{Deserialize, Serialize};
use std::{error::Error, fmt, path::Path};

/// CODATA 2022 atomic unit of length in Å. Source: NIST allascii.txt.
pub const BOHR_ANGSTROM: f64 = 0.529_177_210_544;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DensityError(pub String);
impl fmt::Display for DensityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}
impl Error for DensityError {}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DensityField {
    TotalElectronDensity,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DensityUnit {
    ElectronsPerBohrCubed,
}

/// Caller assertion about the scalar field, not a deduction from cube comments.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct DensityDeclaration {
    pub field: DensityField,
    pub unit: DensityUnit,
    /// Quantum method/basis and original calculation/file identifier.
    pub provenance: String,
}

#[derive(Clone, Debug, Serialize)]
pub struct CubeAtom {
    pub atomic_number: u8,
    pub charge: f64,
    pub position_bohr: [f64; 3],
    pub position_angstrom: [f64; 3],
}
#[derive(Clone, Debug, Serialize)]
pub struct GridMetadata {
    pub source: String,
    pub comments: [String; 2],
    pub declaration: DensityDeclaration,
    pub dimensions: [usize; 3],
    pub origin_bohr: [f64; 3],
    pub vectors_bohr: [[f64; 3]; 3],
    pub origin_angstrom: [f64; 3],
    pub vectors_angstrom: [[f64; 3]; 3],
    pub bohr_to_angstrom: f64,
    pub atoms: Vec<CubeAtom>,
    pub ordering: String,
}

/// Validated scalar cube. Internal arrays cannot be mutated after validation.
#[derive(Clone, Debug)]
pub struct DensityGrid {
    metadata: GridMetadata,
    values: Vec<f64>,
}

fn fail(message: impl Into<String>) -> DensityError {
    DensityError(message.into())
}
fn number(token: &str) -> Result<f64, DensityError> {
    let n = token
        .replace(['D', 'd'], "E")
        .parse::<f64>()
        .map_err(|_| fail(format!("invalid cube number {token:?}")))?;
    if !n.is_finite() {
        return Err(fail("cube numbers must be finite"));
    }
    Ok(n)
}
fn xyz(tokens: &[&str]) -> Result<[f64; 3], DensityError> {
    Ok([number(tokens[0])?, number(tokens[1])?, number(tokens[2])?])
}
fn angstrom(v: [f64; 3]) -> Result<[f64; 3], DensityError> {
    let result = v.map(|x| x * BOHR_ANGSTROM);
    if result.iter().any(|x| !x.is_finite()) {
        return Err(fail("coordinate conversion overflow"));
    }
    Ok(result)
}
fn determinant(v: [[f64; 3]; 3]) -> f64 {
    v[0][0] * (v[1][1] * v[2][2] - v[1][2] * v[2][1])
        - v[0][1] * (v[1][0] * v[2][2] - v[1][2] * v[2][0])
        + v[0][2] * (v[1][0] * v[2][1] - v[1][1] * v[2][0])
}

impl DensityGrid {
    pub fn from_cube_file(
        path: impl AsRef<Path>,
        declaration: DensityDeclaration,
    ) -> Result<Self, DensityError> {
        let path = path.as_ref();
        let text =
            std::fs::read_to_string(path).map_err(|e| fail(format!("{}: {e}", path.display())))?;
        Self::from_cube_str(&path.display().to_string(), &text, declaration)
    }

    /// Read positive-count Gaussian scalar cubes, with coordinates in Bohr.
    /// Negative-count unit extensions and molecular-orbital cubes are rejected.
    pub fn from_cube_str(
        source: &str,
        text: &str,
        declaration: DensityDeclaration,
    ) -> Result<Self, DensityError> {
        if source.trim().is_empty() || declaration.provenance.trim().is_empty() {
            return Err(fail("source and density provenance must be explicit"));
        }
        let mut lines = text.lines();
        let comments = [
            lines
                .next()
                .ok_or_else(|| fail("missing cube title"))?
                .to_owned(),
            lines
                .next()
                .ok_or_else(|| fail("missing cube description"))?
                .to_owned(),
        ];
        let header: Vec<_> = lines
            .next()
            .ok_or_else(|| fail("missing atom count/origin"))?
            .split_whitespace()
            .collect();
        if !(header.len() == 4 || header.len() == 5 && header[4] == "1") {
            return Err(fail("cube header must declare one scalar field"));
        }
        let n_atoms = header[0]
            .parse::<usize>()
            .map_err(|_| fail("positive atom count required; orbital cubes unsupported"))?;
        if n_atoms == 0 {
            return Err(fail("positive atom count required"));
        }
        let origin_bohr = xyz(&header[1..4])?;
        let mut dimensions = [0; 3];
        let mut vectors_bohr = [[0.; 3]; 3];
        for axis in 0..3 {
            let tokens: Vec<_> = lines
                .next()
                .ok_or_else(|| fail("missing grid vector"))?
                .split_whitespace()
                .collect();
            if tokens.len() != 4 {
                return Err(fail(
                    "grid record requires count and three vector components",
                ));
            }
            dimensions[axis] = tokens[0].parse::<usize>().map_err(|_| {
                fail("positive grid counts required; ambiguous unit extensions unsupported")
            })?;
            if dimensions[axis] < 2 {
                return Err(fail("each grid dimension must have at least two nodes"));
            }
            vectors_bohr[axis] = xyz(&tokens[1..4])?;
        }
        let size = dimensions
            .into_iter()
            .try_fold(1usize, |n, d| n.checked_mul(d))
            .ok_or_else(|| fail("grid size overflow"))?;
        let vectors_angstrom = [
            angstrom(vectors_bohr[0])?,
            angstrom(vectors_bohr[1])?,
            angstrom(vectors_bohr[2])?,
        ];
        let cell = determinant(vectors_angstrom).abs();
        if !cell.is_finite() || cell == 0. {
            return Err(fail("grid vectors must span a finite nonzero cell"));
        }
        let origin_angstrom = angstrom(origin_bohr)?;
        // Bounds must remain representable, even for an oblique input grid.
        for axis in 0..3 {
            let extent: f64 = (0..3)
                .map(|v| vectors_angstrom[v][axis].abs() * (dimensions[v] - 1) as f64)
                .sum();
            if !(origin_angstrom[axis].abs() + extent).is_finite() {
                return Err(fail("grid coordinate extent overflow"));
            }
        }
        let mut atoms = Vec::new();
        for _ in 0..n_atoms {
            let tokens: Vec<_> = lines
                .next()
                .ok_or_else(|| fail("truncated atom records"))?
                .split_whitespace()
                .collect();
            if tokens.len() != 5 {
                return Err(fail("atom record requires number, charge and coordinates"));
            }
            let atomic_number = tokens[0]
                .parse::<u8>()
                .map_err(|_| fail("invalid atomic number"))?;
            if atomic_number > 118 {
                return Err(fail("atomic number exceeds 118"));
            }
            let position_bohr = xyz(&tokens[2..5])?;
            atoms.push(CubeAtom {
                atomic_number,
                charge: number(tokens[1])?,
                position_bohr,
                position_angstrom: angstrom(position_bohr)?,
            });
        }
        let mut values = Vec::new();
        for token in lines.flat_map(str::split_whitespace) {
            if values.len() == size {
                return Err(fail(
                    "excess cube values; multiple scalar fields unsupported",
                ));
            }
            let value = number(token)?;
            if value < 0. {
                return Err(fail("total electron density must be nonnegative"));
            }
            values.push(value);
        }
        if values.len() != size {
            return Err(fail(format!(
                "expected {size} density values, read {}",
                values.len()
            )));
        }
        Ok(Self {
            metadata: GridMetadata {
                source: source.into(),
                comments,
                declaration,
                dimensions,
                origin_bohr,
                vectors_bohr,
                origin_angstrom,
                vectors_angstrom,
                bohr_to_angstrom: BOHR_ANGSTROM,
                atoms,
                ordering: "x_slow_y_middle_z_fast".into(),
            },
            values,
        })
    }
    pub fn metadata(&self) -> &GridMetadata {
        &self.metadata
    }
    pub fn values(&self) -> &[f64] {
        &self.values
    }
    pub fn value(&self, index: [usize; 3]) -> Option<f64> {
        let [nx, ny, nz] = self.metadata.dimensions;
        let [i, j, k] = index;
        (i < nx && j < ny && k < nz).then(|| self.values[(i * ny + j) * nz + k])
    }
}
