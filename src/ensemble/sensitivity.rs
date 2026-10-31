use super::{DescriptorStatistics, EnergyKind, EnsembleStatistics, PopulationOrigin};
use serde::Serialize;

/// Physically interpretable spread of each supplied ensemble descriptor.
/// No universal flexibility score, sampling correction, or uncertainty estimate.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct SensitivityReport {
    pub ligand_id: String,
    pub statistics: EnsembleStatistics,
    /// Maximum minus minimum in the descriptor's original unit. Missing when
    /// the statistics policy does not permit an ensemble comparison.
    pub ranges: std::collections::BTreeMap<String, Option<f64>>,
    pub regional_ranges: std::collections::BTreeMap<String, Option<f64>>,
}

impl EnsembleStatistics {
    pub fn sensitivity(&self, ligand_id: &str) -> SensitivityReport {
        let range = |s: &DescriptorStatistics| s.maximum.zip(s.minimum).map(|(max, min)| max - min);
        SensitivityReport {
            ligand_id: ligand_id.into(),
            statistics: self.clone(),
            ranges: self
                .descriptors
                .iter()
                .map(|(k, s)| (k.clone(), range(s)))
                .collect(),
            regional_ranges: self
                .regions
                .iter()
                .map(|(k, s)| (k.clone(), range(s)))
                .collect(),
        }
    }
}

impl SensitivityReport {
    /// Compact table for research workflows. Full settings and energy metadata
    /// remain in the serializable report; archive its JSON alongside this CSV.
    pub fn to_csv(&self) -> Result<String, String> {
        let mut writer = csv::Writer::from_writer(Vec::new());
        writer
            .write_record([
                "ligand_id",
                "descriptor",
                "regional",
                "unit",
                "total_count",
                "successful_count",
                "complete",
                "minimum",
                "maximum",
                "range",
                "standard_deviation",
                "mean",
                "median",
                "weighted_mean",
                "boltzmann_mean",
                "retained_population_mass",
                "minimum_conformers",
                "maximum_conformers",
                "missing_conformers",
                "unavailable_reason",
                "population_origin",
                "temperature_k",
            ])
            .map_err(|e| e.to_string())?;
        let (origin, temperature): (String, Option<f64>) =
            match self.statistics.population.as_ref().map(|p| p.origin()) {
                None => ("unavailable".into(), None),
                Some(PopulationOrigin::Supplied { .. }) => ("supplied".into(), None),
                Some(PopulationOrigin::Boltzmann {
                    energy_kind,
                    temperature_k,
                    ..
                }) => {
                    let label = match energy_kind {
                        EnergyKind::Electronic => "boltzmann_electronic",
                        EnergyKind::ForceField => "boltzmann_force_field",
                        EnergyKind::GibbsFreeEnergy => "boltzmann_gibbs_free_energy",
                    };
                    (label.into(), Some(*temperature_k))
                }
            };
        let number = |value: Option<f64>| value.map(|v| v.to_string()).unwrap_or_default();
        for (regional, fields, ranges) in [
            (false, &self.statistics.descriptors, &self.ranges),
            (true, &self.statistics.regions, &self.regional_ranges),
        ] {
            for (name, s) in fields {
                writer
                    .write_record([
                        self.ligand_id.clone(),
                        name.clone(),
                        regional.to_string(),
                        s.unit.clone(),
                        s.total_count.to_string(),
                        s.successful_count.to_string(),
                        s.complete.to_string(),
                        number(s.minimum),
                        number(s.maximum),
                        number(ranges[name]),
                        number(s.standard_deviation),
                        number(s.mean),
                        number(s.median),
                        number(s.weighted_mean),
                        number(s.boltzmann_mean),
                        number(s.retained_population_mass),
                        serde_json::to_string(&s.minimum_conformers).map_err(|e| e.to_string())?,
                        serde_json::to_string(&s.maximum_conformers).map_err(|e| e.to_string())?,
                        serde_json::to_string(&s.missing_conformers).map_err(|e| e.to_string())?,
                        s.unavailable_reason.clone().unwrap_or_default(),
                        origin.clone(),
                        number(temperature),
                    ])
                    .map_err(|e| e.to_string())?;
            }
        }
        let bytes = writer.into_inner().map_err(|e| e.to_string())?;
        String::from_utf8(bytes).map_err(|e| e.to_string())
    }
}
