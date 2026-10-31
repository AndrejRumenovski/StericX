//! Buried-volume radius profiles at a fixed center, following DBSTEP's scan variable.
use crate::density::{DensityError, DensityGrid, DensityVolume};
use crate::ensemble::{Conformer, DescriptorOptions, DescriptorOutcome};
use crate::geometry::{
    BuriedVolumeCalculator, BuriedVolumeDetails, coordination_center_with_neighbors,
    detect_trivalent_donor,
};
use glam::Vec3;
use serde::Serialize;

#[derive(Clone, Debug, Serialize)]
pub struct RadiusSample<T> {
    pub requested_radius_angstrom: f64,
    /// Actual engine radius (hard-sphere uses f32; density uses f64).
    pub evaluated_radius_angstrom: Option<f64>,
    pub outcome: DescriptorOutcome<T>,
}
#[derive(Clone, Debug, Serialize)]
pub struct HardSphereProfile {
    pub definition: String,
    pub source: String,
    pub source_geometry: String,
    pub conformer_index: usize,
    pub conformer_identifier: String,
    pub options: DescriptorOptions,
    pub center_angstrom: Option<[f32; 3]>,
    pub samples: Vec<RadiusSample<BuriedVolumeDetails>>,
}
#[derive(Clone, Debug, Serialize)]
pub struct DensityProfile {
    pub definition: String,
    pub center_angstrom: [f64; 3],
    pub threshold_electrons_per_bohr3: f64,
    pub samples: Vec<RadiusSample<DensityVolume>>,
}
fn valid_radii(radii: &[f64]) -> Result<(), String> {
    if radii.is_empty() || radii.iter().any(|r| !r.is_finite()) {
        return Err("radius profile requires a nonempty list of finite radii".into());
    }
    Ok(())
}
fn native_radius(r: f64) -> Result<f32, String> {
    let value = r as f32;
    if r <= 0. || !value.is_finite() || value <= 0. {
        return Err("radius must be positive and representable as a finite positive f32".into());
    }
    Ok(value)
}
impl Conformer {
    /// Vary only the sampling-sphere radius; center and molecular frame stay fixed.
    /// Order, duplicate radii and failed finite-radius calculations are retained.
    pub fn buried_volume_scan(
        &self,
        radii_angstrom: &[f64],
        options: &DescriptorOptions,
    ) -> Result<HardSphereProfile, String> {
        valid_radii(radii_angstrom)?;
        let context = (|| {
            let frame = self.frame.as_ref().map_err(Clone::clone)?;
            let t = detect_trivalent_donor(
                &frame.molecule,
                frame.bonds.as_deref(),
                &options.donor_element,
                options.donor_index,
                None,
                true,
            )?;
            let mut config = options.volume;
            if let Some(r) = radii_angstrom.iter().find_map(|r| native_radius(*r).ok()) {
                config.sphere_radius = r;
            }
            let center = match options.center {
                Some(c) => Vec3::from_array(c),
                None => coordination_center_with_neighbors(
                    &frame.molecule,
                    t.donor_idx,
                    t.substituents,
                    config,
                )
                .map_err(|e| e.to_string())?,
            };
            if !center.is_finite() {
                return Err("profile center must be finite".to_string());
            }
            Ok((frame, t, center))
        })();
        let samples = radii_angstrom
            .iter()
            .map(|&r| {
                let converted = native_radius(r);
                let evaluated_radius_angstrom = converted.as_ref().ok().map(|v| f64::from(*v));
                let result = converted.and_then(|radius| {
                    let (frame, t, center) = context.as_ref().map_err(Clone::clone)?;
                    let mut config = options.volume;
                    config.sphere_radius = radius;
                    BuriedVolumeCalculator::compute_detailed_with_center_and_neighbors(
                        &frame.molecule,
                        t.donor_idx,
                        t.substituents,
                        *center,
                        config,
                    )
                    .map_err(|e| e.to_string())
                });
                RadiusSample {
                    requested_radius_angstrom: r,
                    evaluated_radius_angstrom,
                    outcome: result.into(),
                }
            })
            .collect();
        Ok(HardSphereProfile {
            definition: "sampling_sphere_radius_at_fixed_center".into(),
            source: self.source.clone(),
            source_geometry: self.source_geometry.clone(),
            conformer_index: self.index,
            conformer_identifier: self.identifier.clone(),
            options: options.clone(),
            center_angstrom: context.as_ref().ok().map(|(_, _, c)| c.to_array()),
            samples,
        })
    }
}
impl DensityGrid {
    /// Radius profile on the unchanged input density grid. Out-of-domain radii
    /// produce individual errors; they are never padded or omitted.
    pub fn buried_volume_scan(
        &self,
        center_angstrom: [f64; 3],
        radii_angstrom: &[f64],
        threshold_electrons_per_bohr3: f64,
    ) -> Result<DensityProfile, DensityError> {
        valid_radii(radii_angstrom).map_err(DensityError)?;
        if center_angstrom.iter().any(|x| !x.is_finite())
            || !threshold_electrons_per_bohr3.is_finite()
            || threshold_electrons_per_bohr3 <= 0.
        {
            return Err(DensityError(
                "profile center and positive threshold must be finite".into(),
            ));
        }
        let samples = radii_angstrom
            .iter()
            .map(|&r| RadiusSample {
                requested_radius_angstrom: r,
                evaluated_radius_angstrom: (r > 0.).then_some(r),
                outcome: self
                    .buried_volume(center_angstrom, r, threshold_electrons_per_bohr3)
                    .map_err(|e| e.to_string())
                    .into(),
            })
            .collect();
        Ok(DensityProfile {
            definition: "sampling_sphere_radius_at_fixed_center".into(),
            center_angstrom,
            threshold_electrons_per_bohr3,
            samples,
        })
    }
}
