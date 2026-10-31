use super::{DensityError, DensityGrid, GridMetadata, determinant, fail};
use serde::Serialize;

/// Node quadrature of a density superlevel set inside an explicit sphere.
/// `buried_volume_angstrom3` is count × cell volume, not analytic normalization.
#[derive(Clone, Debug, Serialize)]
pub struct DensityVolume {
    pub method: String,
    pub center_angstrom: [f64; 3],
    pub radius_angstrom: f64,
    pub threshold_electrons_per_bohr3: f64,
    pub occupied_nodes: usize,
    pub sphere_nodes: usize,
    pub cell_volume_angstrom3: f64,
    pub buried_volume_angstrom3: f64,
    pub sampled_sphere_volume_angstrom3: f64,
    pub analytic_sphere_volume_angstrom3: f64,
    pub sphere_sampling_error_percent: f64,
    pub percent_buried_volume: f64,
    pub grid: GridMetadata,
}

impl DensityGrid {
    /// Deterministic strict-threshold occupation. No interpolation, padding,
    /// radius scaling or atomic hard spheres enter this separate calculation.
    pub fn buried_volume(
        &self,
        center_angstrom: [f64; 3],
        radius_angstrom: f64,
        threshold_electrons_per_bohr3: f64,
    ) -> Result<DensityVolume, DensityError> {
        if center_angstrom.iter().any(|x| !x.is_finite())
            || !radius_angstrom.is_finite()
            || radius_angstrom <= 0.
        {
            return Err(fail(
                "center must be finite and sphere radius finite and positive",
            ));
        }
        if !threshold_electrons_per_bohr3.is_finite() || threshold_electrons_per_bohr3 <= 0. {
            return Err(fail(
                "density threshold must be explicit, finite and positive in electrons/Bohr^3",
            ));
        }
        let m = &self.metadata;
        let v = m.vectors_bohr;
        let spacing_bohr = v[0][0];
        if spacing_bohr <= 0.
            || (0..3).any(|i| (0..3).any(|j| v[i][j] != if i == j { spacing_bohr } else { 0. }))
        {
            return Err(fail(
                "validated density analysis requires positive axis-aligned isotropic grid vectors; original vectors retained",
            ));
        }
        let spacing = m.vectors_angstrom[0][0];
        for (axis, &center) in center_angstrom.iter().enumerate() {
            let lo = m.origin_angstrom[axis];
            let hi = lo + (m.dimensions[axis] - 1) as f64 * spacing;
            if center - radius_angstrom < lo || center + radius_angstrom > hi {
                return Err(fail(format!(
                    "sphere extends outside density grid on axis {axis}; no extrapolation or padding"
                )));
            }
        }
        let cell = determinant(m.vectors_angstrom).abs();
        let analytic = 4. * std::f64::consts::PI * radius_angstrom.powi(3) / 3.;
        if !analytic.is_finite() || analytic == 0. {
            return Err(fail("sphere volume is outside finite numerical range"));
        }
        let radius_squared = radius_angstrom * radius_angstrom;
        let mut sphere_nodes = 0;
        let mut occupied_nodes = 0;
        let [nx, ny, nz] = m.dimensions;
        for i in 0..nx {
            let x = m.origin_angstrom[0] + i as f64 * spacing - center_angstrom[0];
            for j in 0..ny {
                let y = m.origin_angstrom[1] + j as f64 * spacing - center_angstrom[1];
                for k in 0..nz {
                    let z = m.origin_angstrom[2] + k as f64 * spacing - center_angstrom[2];
                    if (x * x + y * y) + z * z <= radius_squared {
                        sphere_nodes += 1;
                        if self.values[(i * ny + j) * nz + k] > threshold_electrons_per_bohr3 {
                            occupied_nodes += 1;
                        }
                    }
                }
            }
        }
        if sphere_nodes == 0 {
            return Err(fail("sphere contains no sampled density nodes"));
        }
        let sampled = sphere_nodes as f64 * cell;
        let buried = occupied_nodes as f64 * cell;
        let sampling_error = (sampled - analytic) / analytic * 100.;
        if !sampled.is_finite() || !buried.is_finite() || !sampling_error.is_finite() {
            return Err(fail("density integration result overflow"));
        }
        Ok(DensityVolume {
            method: "electron_density_node_quadrature".into(),
            center_angstrom,
            radius_angstrom,
            threshold_electrons_per_bohr3,
            occupied_nodes,
            sphere_nodes,
            cell_volume_angstrom3: cell,
            buried_volume_angstrom3: buried,
            sampled_sphere_volume_angstrom3: sampled,
            analytic_sphere_volume_angstrom3: analytic,
            sphere_sampling_error_percent: sampling_error,
            percent_buried_volume: 100. * occupied_nodes as f64 / sphere_nodes as f64,
            grid: m.clone(),
        })
    }
}
