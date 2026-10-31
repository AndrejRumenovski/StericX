use super::{ConformerDescriptors, EnsembleDescriptors, Population, PopulationOrigin};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MissingConformerPolicy {
    /// Any unsuccessful value leaves all statistics for that field unavailable.
    #[default]
    RequireComplete,
    /// Explicit conditional statistics of successful records; original failures
    /// and retained population mass remain visible in every field's result.
    AvailableCases,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct StatisticsOptions {
    pub missing: MissingConformerPolicy,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct ConformerIdentity {
    pub index: usize,
    pub identifier: String,
    pub source: String,
}

impl From<&ConformerDescriptors> for ConformerIdentity {
    fn from(row: &ConformerDescriptors) -> Self {
        Self {
            index: row.index,
            identifier: row.identifier.clone(),
            source: row.source.clone(),
        }
    }
}

/// Descriptive statistics of the supplied finite population of conformers.
/// Standard deviation uses divisor n (ddof=0), not n-1. All statistics are f64
/// reductions of native f32 descriptors. Extrema retain every exactly tied ID.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct DescriptorStatistics {
    pub unit: String,
    pub total_count: usize,
    pub successful_count: usize,
    pub missing_conformers: Vec<ConformerIdentity>,
    pub complete: bool,
    pub unavailable_reason: Option<String>,
    pub mean: Option<f64>,
    pub minimum: Option<f64>,
    pub maximum: Option<f64>,
    pub standard_deviation: Option<f64>,
    pub median: Option<f64>,
    pub weighted_mean: Option<f64>,
    /// Present only for weights derived by the energy/provenance-checked API.
    pub boltzmann_mean: Option<f64>,
    /// Weight carried by successful states before any explicit conditioning.
    pub retained_population_mass: Option<f64>,
    pub minimum_conformers: Vec<ConformerIdentity>,
    pub maximum_conformers: Vec<ConformerIdentity>,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct EnsembleStatistics {
    pub descriptor_options: super::DescriptorOptions,
    pub options: StatisticsOptions,
    pub population: Option<Population>,
    /// Scalar statistics assume supplied records represent the intended species.
    /// This check covers element order and adjacency, not bond orders, charges,
    /// stereochemistry or experimental identity.
    pub structural_correspondence: String,
    pub descriptors: BTreeMap<String, DescriptorStatistics>,
    /// Molecular-frame region keys identify the original plane atom and bin.
    /// They are unavailable without verified ordered adjacency in every record.
    pub regions: BTreeMap<String, DescriptorStatistics>,
    pub regional_unavailable_reason: Option<String>,
}

impl EnsembleDescriptors {
    pub fn statistics(
        &self,
        options: &StatisticsOptions,
        population: Option<&Population>,
    ) -> Result<EnsembleStatistics, String> {
        if self.conformers.is_empty() {
            return Err("cannot summarize an empty ensemble".into());
        }
        if let Some(pop) = population {
            if pop.weights().len() != self.conformers.len() {
                return Err("weights must match the complete ordered conformer list".into());
            }
            if let PopulationOrigin::Boltzmann { energies, .. } = pop.origin() {
                for (row, energy) in self.conformers.iter().zip(energies) {
                    if row.energy.as_ref() != Some(energy) {
                        return Err(format!(
                            "population energy provenance does not match conformer {}",
                            row.index
                        ));
                    }
                }
            }
        }
        let geometries = self
            .conformers
            .iter()
            .filter_map(|r| r.geometry.value())
            .collect::<Vec<_>>();
        if let Some(first) = geometries.first() {
            let first_bonds = geometries.iter().find_map(|g| g.bonds.as_ref());
            for frame in &geometries {
                if frame.elements != first.elements {
                    return Err("ensemble atom identities/order differ; map corresponding atoms before aggregation".into());
                }
                if let (Some(a), Some(b)) = (first_bonds, &frame.bonds) {
                    if a != b {
                        return Err("ensemble supplied connectivity differs".into());
                    }
                }
            }
        }
        let ordered_graph = geometries.len() == self.conformers.len()
            && geometries.iter().all(|g| g.bonds.is_some());
        let fields = [
            ("sterimol_l", "angstrom"),
            ("sterimol_b1", "angstrom"),
            ("sterimol_b5", "angstrom"),
            ("buried_volume", "angstrom^3"),
            ("percent_buried_volume", "percent"),
            ("qvbur_min", "angstrom^3"),
            ("qvbur_max", "angstrom^3"),
            ("max_delta_qvbur", "angstrom^3"),
            ("ovbur_min", "angstrom^3"),
            ("ovbur_max", "angstrom^3"),
            ("near_vbur", "angstrom^3"),
            ("far_vbur", "angstrom^3"),
            ("pyr_p", "dimensionless"),
            ("pyr_alpha", "degree"),
        ];
        let mut descriptors = BTreeMap::new();
        for (name, unit) in fields {
            let values: Vec<Option<f32>> =
                self.conformers.iter().map(|r| scalar(r, name)).collect();
            descriptors.insert(
                name.into(),
                summarize(&self.conformers, &values, unit, options, population),
            );
        }
        let mut regions = BTreeMap::new();
        if ordered_graph {
            let mut columns: BTreeMap<String, Vec<Option<f32>>> = BTreeMap::new();
            for (i, row) in self.conformers.iter().enumerate() {
                if let Some(volume) = row.buried_volume.value() {
                    for plane in &volume.orientations {
                        for (kind, values) in [
                            ("quadrant", plane.quadrants.as_slice()),
                            ("octant", plane.octants.as_slice()),
                        ] {
                            for (bin, value) in values.iter().enumerate() {
                                let key =
                                    format!("plane_{}.{}_{}", plane.plane_atom_index, kind, bin);
                                columns
                                    .entry(key)
                                    .or_insert_with(|| vec![None; self.conformers.len()])[i] =
                                    Some(*value);
                            }
                        }
                    }
                }
            }
            for (name, values) in columns {
                regions.insert(
                    name,
                    summarize(&self.conformers, &values, "angstrom^3", options, population),
                );
            }
        }
        let regional_unavailable_reason = if !ordered_graph {
            Some("regional statistics require complete matching element order and explicit adjacency; raw per-conformer regions remain available".into())
        } else if regions.is_empty() {
            Some("no conformer has a defined regional volume frame".into())
        } else {
            None
        };
        Ok(EnsembleStatistics {
            descriptor_options: self.options.clone(),
            options: options.clone(),
            population: population.cloned(),
            structural_correspondence: if ordered_graph {
                "matching_element_order_and_adjacency"
            } else {
                "element_order_checked_where_available_topology_unverified"
            }
            .into(),
            descriptors,
            regions,
            regional_unavailable_reason,
        })
    }
}

fn scalar(row: &ConformerDescriptors, name: &str) -> Option<f32> {
    match name {
        "sterimol_l" => row.sterimol.value().map(|v| v.l),
        "sterimol_b1" => row.sterimol.value().map(|v| v.b1),
        "sterimol_b5" => row.sterimol.value().map(|v| v.b5),
        "pyr_p" => row.pyramidalization.value().map(|v| v.pyr_p),
        "pyr_alpha" => row.pyramidalization.value().map(|v| v.pyr_alpha),
        _ => row.buried_volume.value().map(|v| match name {
            "buried_volume" => v.summary.buried_volume,
            "percent_buried_volume" => v.summary.percent_buried_volume,
            "qvbur_min" => v.summary.qvbur_min,
            "qvbur_max" => v.summary.qvbur_max,
            "max_delta_qvbur" => v.summary.max_delta_qvbur,
            "ovbur_min" => v.summary.ovbur_min,
            "ovbur_max" => v.summary.ovbur_max,
            "near_vbur" => v.summary.near_vbur,
            "far_vbur" => v.summary.far_vbur,
            _ => unreachable!("declared scalar descriptor"),
        }),
    }
}

fn summarize(
    rows: &[ConformerDescriptors],
    values: &[Option<f32>],
    unit: &str,
    options: &StatisticsOptions,
    population: Option<&Population>,
) -> DescriptorStatistics {
    let available: Vec<(usize, f64)> = values
        .iter()
        .enumerate()
        .filter_map(|(i, v)| v.filter(|v| v.is_finite()).map(|v| (i, f64::from(v))))
        .collect();
    let missing = rows
        .iter()
        .zip(values)
        .filter(|(_, v)| !v.is_some_and(f32::is_finite))
        .map(|(r, _)| r.into())
        .collect();
    let mass = population.map(|p| available.iter().map(|(i, _)| p.weights()[*i]).sum::<f64>());
    let mut result = DescriptorStatistics {
        unit: unit.into(),
        total_count: rows.len(),
        successful_count: available.len(),
        missing_conformers: missing,
        complete: available.len() == rows.len(),
        unavailable_reason: None,
        mean: None,
        minimum: None,
        maximum: None,
        standard_deviation: None,
        median: None,
        weighted_mean: None,
        boltzmann_mean: None,
        retained_population_mass: mass,
        minimum_conformers: Vec::new(),
        maximum_conformers: Vec::new(),
    };
    if available.is_empty()
        || (!result.complete && options.missing == MissingConformerPolicy::RequireComplete)
    {
        result.unavailable_reason = Some(
            "one or more conformer descriptors unavailable; no complete-ensemble statistic".into(),
        );
        return result;
    }
    // f32 descriptor inputs make both their f64 sum and squared deviations
    // representable for any realizable in-memory ensemble.
    let n = available.len() as f64;
    let mean = available.iter().map(|(_, v)| *v).sum::<f64>() / n;
    let variance = available
        .iter()
        .map(|(_, v)| (v - mean).powi(2))
        .sum::<f64>()
        / n;
    let mut sorted = available.iter().map(|(_, v)| *v).collect::<Vec<_>>();
    sorted.sort_by(f64::total_cmp);
    let minimum = sorted[0];
    let maximum = sorted[sorted.len() - 1];
    result.mean = Some(mean);
    result.minimum = Some(minimum);
    result.maximum = Some(maximum);
    result.standard_deviation = Some(variance.sqrt());
    result.median = Some((sorted[(sorted.len() - 1) / 2] + sorted[sorted.len() / 2]) / 2.0);
    result.minimum_conformers = available
        .iter()
        .filter(|(_, v)| *v == minimum)
        .map(|(i, _)| (&rows[*i]).into())
        .collect();
    result.maximum_conformers = available
        .iter()
        .filter(|(_, v)| *v == maximum)
        .map(|(i, _)| (&rows[*i]).into())
        .collect();
    if let (Some(pop), Some(mass)) = (population, mass) {
        if mass > 0.0 {
            let weighted = available
                .iter()
                .map(|(i, v)| v * (pop.weights()[*i] / mass))
                .sum();
            result.weighted_mean = Some(weighted);
            if pop.is_boltzmann() {
                result.boltzmann_mean = Some(weighted);
            }
        }
    }
    result
}
