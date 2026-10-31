//! Ordered conformer records calling the existing native descriptor kernels.
//!
//! Parsing and descriptor failures are data: each input record keeps its index,
//! title, original source and individual calculation outcomes. Loading an SDF
//! does not infer energies, populations, atom mapping or ensemble completeness.

mod input;
mod population;
mod sensitivity;
mod statistics;

use crate::geometry::{
    BuriedVolumeCalculator, BuriedVolumeConfig, BuriedVolumeDetails, MolecularFrame,
    PyramidalizationCalculator, PyramidalizationParams, SterimolCalculator, SterimolParams,
    coordination_center_with_neighbors, detect_trivalent_donor,
};
use glam::Vec3;
use serde::{Deserialize, Serialize};

pub use input::{Conformer, Ensemble, SourceProperty};
pub use population::{Population, PopulationOrigin};
pub use sensitivity::SensitivityReport;
pub use statistics::{
    ConformerIdentity, DescriptorStatistics, EnsembleStatistics, MissingConformerPolicy,
    StatisticsOptions,
};

/// Energy units must be declared by the source or caller, never guessed.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EnergyUnit {
    KcalMol,
    KjMol,
    Hartree,
    ElectronVolt,
}

/// These energy quantities support different physical interpretations.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EnergyKind {
    Electronic,
    ForceField,
    GibbsFreeEnergy,
}

/// Explicit, unmodified per-state energy provenance. This metadata alone does
/// not establish equilibrium populations or that a conformer search is complete.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct EnergyMetadata {
    pub value: f64,
    pub unit: EnergyUnit,
    pub kind: EnergyKind,
    /// Method, solvent and standard-state convention sufficient to identify the
    /// common energy model. Different models must not be combined implicitly.
    pub model: String,
    /// Original calculation/file/property identifier.
    pub provenance: String,
    /// Temperature at which free energies were evaluated, if applicable (K).
    pub temperature_k: Option<f64>,
    /// Number of equal-energy microstates represented by this record. Explicitly
    /// enumerated states use one; equal energies alone do not imply multiplicity.
    pub degeneracy: u32,
}

impl EnergyMetadata {
    pub fn validate(&self) -> Result<(), String> {
        if !self.value.is_finite() || self.degeneracy == 0 {
            return Err("energy must be finite and degeneracy positive".into());
        }
        if self.model.trim().is_empty() || self.provenance.trim().is_empty() {
            return Err("energy model and provenance must be explicit".into());
        }
        if self
            .temperature_k
            .is_some_and(|t| !t.is_finite() || t <= 0.0)
        {
            return Err("energy temperature must be finite and positive".into());
        }
        if self.kind == EnergyKind::GibbsFreeEnergy && self.temperature_k.is_none() {
            return Err("Gibbs free energy requires its evaluation temperature".into());
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SterimolAxis {
    /// Attachment atom is excluded; axis points from donor to reference atom.
    Bond,
    /// Virtual center to donor; all real atoms contribute.
    Coordination,
}

/// Explicit settings returned with the results. Defaults match the descriptor
/// CLI: bond-axis Sterimol and a 2.28 Å virtual-center distance for buried volume.
/// The lower-level `BuriedVolumeConfig::default()` separately retains 2.1 Å.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct DescriptorOptions {
    pub donor_element: String,
    pub donor_index: Option<usize>,
    pub reference_index: Option<usize>,
    pub axis: SterimolAxis,
    pub volume: BuriedVolumeConfig,
    /// An explicit source-coordinate center overrides geometric inference.
    pub center: Option<[f32; 3]>,
    /// None means the legacy 360-direction scan. Explicit endpoint-inclusive
    /// sampling is supported by the existing coordination-axis kernel.
    pub angular_samples: Option<usize>,
    /// None uses raw L for bond axes and +0.40 Å for coordination axes.
    pub l_correction: Option<f32>,
}

impl Default for DescriptorOptions {
    fn default() -> Self {
        Self {
            donor_element: "P".into(),
            donor_index: None,
            reference_index: None,
            axis: SterimolAxis::Bond,
            volume: BuriedVolumeConfig {
                center_distance: 2.28,
                ..Default::default()
            },
            center: None,
            angular_samples: None,
            l_correction: None,
        }
    }
}

/// A calculation failure is retained alongside other successful descriptors.
#[derive(Clone, Debug, PartialEq, Serialize)]
#[serde(tag = "status", rename_all = "snake_case")]
pub enum DescriptorOutcome<T> {
    Success { value: T },
    Error { message: String },
}

impl<T> From<Result<T, String>> for DescriptorOutcome<T> {
    fn from(result: Result<T, String>) -> Self {
        match result {
            Ok(value) => Self::Success { value },
            Err(message) => Self::Error { message },
        }
    }
}

impl<T> DescriptorOutcome<T> {
    pub fn value(&self) -> Option<&T> {
        match self {
            Self::Success { value } => Some(value),
            Self::Error { .. } => None,
        }
    }
    pub fn error(&self) -> Option<&str> {
        match self {
            Self::Error { message } => Some(message),
            Self::Success { .. } => None,
        }
    }
}

/// Parsed f32 engine input in unchanged atom order, in addition to original text.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct GeometrySnapshot {
    pub elements: Vec<String>,
    pub coordinates_angstrom: Vec<[f32; 3]>,
    pub radii_angstrom: Vec<f32>,
    pub bonds: Option<Vec<[usize; 2]>>,
}

impl From<&MolecularFrame> for GeometrySnapshot {
    fn from(frame: &MolecularFrame) -> Self {
        Self {
            elements: frame
                .molecule
                .atoms
                .iter()
                .map(|a| a.element.clone())
                .collect(),
            coordinates_angstrom: frame
                .molecule
                .atoms
                .iter()
                .map(|a| a.position.to_array())
                .collect(),
            radii_angstrom: frame.molecule.atoms.iter().map(|a| a.vdw_radius).collect(),
            bonds: frame.bonds.clone(),
        }
    }
}

/// One result per source record, including malformed/unsuccessful records.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct ConformerDescriptors {
    /// Unique positional identity within this source (zero-based).
    pub index: usize,
    /// Original title; duplicate and blank titles are retained verbatim.
    pub identifier: String,
    pub source: String,
    pub source_format: String,
    pub source_geometry: String,
    pub properties: Vec<SourceProperty>,
    pub geometry: DescriptorOutcome<GeometrySnapshot>,
    pub energy: Option<EnergyMetadata>,
    pub energy_error: Option<String>,
    pub sterimol: DescriptorOutcome<SterimolParams>,
    pub buried_volume: DescriptorOutcome<BuriedVolumeDetails>,
    pub pyramidalization: DescriptorOutcome<PyramidalizationParams>,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct EnsembleDescriptors {
    pub options: DescriptorOptions,
    pub conformers: Vec<ConformerDescriptors>,
}

impl Ensemble {
    /// Calculate in source order. No record is dropped or implicitly reweighted.
    pub fn descriptors(&self, options: &DescriptorOptions) -> EnsembleDescriptors {
        EnsembleDescriptors {
            options: options.clone(),
            conformers: self
                .conformers
                .iter()
                .map(|c| c.descriptors(options))
                .collect(),
        }
    }
}

impl Conformer {
    pub fn descriptors(&self, options: &DescriptorOptions) -> ConformerDescriptors {
        let calculate = || -> Result<_, String> {
            let frame = self.frame.as_ref().map_err(Clone::clone)?;
            let topology = detect_trivalent_donor(
                &frame.molecule,
                frame.bonds.as_deref(),
                &options.donor_element,
                options.donor_index,
                None,
                true,
            );
            let center = topology.as_ref().map_err(Clone::clone).and_then(|t| {
                options.center.map(Vec3::from_array).map_or_else(
                    || {
                        coordination_center_with_neighbors(
                            &frame.molecule,
                            t.donor_idx,
                            t.substituents,
                            options.volume,
                        )
                        .map_err(|e| e.to_string())
                    },
                    Ok,
                )
            });
            let sterimol = (|| {
                let mut params = match options.axis {
                    SterimolAxis::Bond => {
                        if options.angular_samples.is_some() {
                            return Err(
                                "explicit angular sampling requires coordination axis".into()
                            );
                        }
                        // Explicit axis atoms define Sterimol even when a trivalent
                        // regional/pyramidalization frame is unavailable.
                        let (donor, reference) =
                            match (options.donor_index, options.reference_index) {
                                (Some(a), Some(b)) => (a, b),
                                _ => {
                                    let t = detect_trivalent_donor(
                                        &frame.molecule,
                                        frame.bonds.as_deref(),
                                        &options.donor_element,
                                        options.donor_index,
                                        options.reference_index,
                                        false,
                                    )?;
                                    (t.donor_idx, t.reference_idx)
                                }
                            };
                        SterimolCalculator::compute(&frame.molecule, donor, reference)
                            .map_err(|e| e.to_string())?
                    }
                    SterimolAxis::Coordination => {
                        let t = topology.as_ref().map_err(Clone::clone)?;
                        let center = *center.as_ref().map_err(Clone::clone)?;
                        match options.angular_samples {
                            Some(n) => SterimolCalculator::compute_with_dummy_with_sampling(
                                &frame.molecule,
                                t.donor_idx,
                                center,
                                n,
                            ),
                            None => SterimolCalculator::compute_with_dummy(
                                &frame.molecule,
                                t.donor_idx,
                                center,
                            ),
                        }
                        .map_err(|e| e.to_string())?
                    }
                };
                params.l += options.l_correction.unwrap_or(match options.axis {
                    SterimolAxis::Bond => 0.0,
                    SterimolAxis::Coordination => 0.40,
                });
                if !params.l.is_finite() {
                    return Err("corrected Sterimol L exceeds finite f32 range".into());
                }
                Ok(params)
            })();
            let volume = (|| {
                let t = topology.as_ref().map_err(Clone::clone)?;
                BuriedVolumeCalculator::compute_detailed_with_center_and_neighbors(
                    &frame.molecule,
                    t.donor_idx,
                    t.substituents,
                    *center.as_ref().map_err(Clone::clone)?,
                    options.volume,
                )
                .map_err(|e| e.to_string())
            })();
            let pyramid = topology.as_ref().map_err(Clone::clone).and_then(|t| {
                PyramidalizationCalculator::compute(&frame.molecule, t.donor_idx, t.substituents)
                    .map_err(|e| e.to_string())
            });
            Ok((sterimol, volume, pyramid))
        };
        let (sterimol, volume, pyramid) =
            calculate().unwrap_or_else(|e| (Err(e.clone()), Err(e.clone()), Err(e)));
        let energy_error = self.energy.as_ref().and_then(|e| e.validate().err());
        ConformerDescriptors {
            index: self.index,
            identifier: self.identifier.clone(),
            source: self.source.clone(),
            source_format: self.source_format.clone(),
            source_geometry: self.source_geometry.clone(),
            properties: self.properties.clone(),
            geometry: self
                .frame
                .as_ref()
                .map(GeometrySnapshot::from)
                .map_err(Clone::clone)
                .into(),
            energy: if energy_error.is_none() {
                self.energy.clone()
            } else {
                None
            },
            energy_error,
            sterimol: sterimol.into(),
            buried_volume: volume.into(),
            pyramidalization: pyramid.into(),
        }
    }
}
