use super::{EnergyKind, EnergyMetadata, EnergyUnit, Ensemble};
use crate::kinetics::GAS_CONSTANT_KCAL_MOL_K;
use serde::Serialize;

/// Provenance is part of a population result, not inferred from its numbers.
#[derive(Clone, Debug, PartialEq, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum PopulationOrigin {
    Supplied {
        provenance: String,
    },
    Boltzmann {
        temperature_k: f64,
        energy_kind: EnergyKind,
        energies: Vec<EnergyMetadata>,
        assumptions: String,
    },
}

/// Validated normalized weights in input order. No implicit uniform fallback.
/// Constructors keep the raw input interpretation alongside the normalized data.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct Population {
    weights: Vec<f64>,
    origin: PopulationOrigin,
}

impl Population {
    pub fn supplied(weights: &[f64], provenance: &str) -> Result<Self, String> {
        if provenance.trim().is_empty() {
            return Err("supplied weights require provenance and interpretation".into());
        }
        Ok(Self {
            weights: normalize(weights)?,
            origin: PopulationOrigin::Supplied {
                provenance: provenance.into(),
            },
        })
    }

    /// Normalize g_i exp[-(E_i-E_min)/(RT)] via shifted log partition terms.
    /// All states must use one declared energy model, kind and unit. Mixed units
    /// require an explicit caller conversion rather than an implicit absolute-
    /// energy conversion that can destroy small differences at large offsets.
    /// Electronic/force-field weighting is an approximation, not a free-energy
    /// equilibrium population. Free energies must be evaluated at exactly T.
    pub fn boltzmann(
        energies: &[EnergyMetadata],
        temperature_k: f64,
        assumptions: &str,
    ) -> Result<Self, String> {
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err("population temperature must be finite and positive (K)".into());
        }
        if assumptions.trim().is_empty() {
            return Err("declare population/state-membership and degeneracy assumptions".into());
        }
        let first = energies
            .first()
            .ok_or("energies must be complete and nonempty")?;
        for energy in energies {
            energy.validate()?;
            if energy.unit != first.unit || energy.kind != first.kind || energy.model != first.model
            {
                return Err("conformer energies must share units, kind and energy model".into());
            }
            if energy.kind == EnergyKind::GibbsFreeEnergy
                && energy.temperature_k != Some(temperature_k)
            {
                return Err(
                    "free-energy evaluation temperature differs from population temperature".into(),
                );
            }
        }
        // SI exact N_A and electron charge; CODATA 2022 E_h in joules.
        // Sources: https://physics.nist.gov/cuu/Constants/Table/allascii.txt
        // The thermochemical kcal is exactly 4184 J. Keep the existing R.
        let kcal_per_unit = match first.unit {
            EnergyUnit::KcalMol => 1.0,
            EnergyUnit::KjMol => 1000.0 / 4184.0,
            EnergyUnit::Hartree => 4.359_744_722_206e-18 * 6.022_140_76e23 / 4184.0,
            EnergyUnit::ElectronVolt => 1.602_176_634e-19 * 6.022_140_76e23 / 4184.0,
        };
        let gas_constant = GAS_CONSTANT_KCAL_MOL_K / kcal_per_unit;
        let minimum = energies
            .iter()
            .map(|e| e.value)
            .fold(f64::INFINITY, f64::min);
        let logs: Vec<f64> = energies
            .iter()
            .map(|e| {
                // Dividing in stages avoids RT underflow at positive subnormal T.
                // Opposite-sign extreme energies can overflow their difference
                // even when division by a large T yields a finite exponent.
                let difference = e.value - minimum;
                let scaled = if difference.is_infinite() {
                    e.value / temperature_k - minimum / temperature_k
                } else {
                    difference / temperature_k
                };
                f64::from(e.degeneracy).ln() - scaled / gas_constant
            })
            .collect();
        let maximum = logs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let raw: Vec<f64> = logs.iter().map(|log| (log - maximum).exp()).collect();
        Ok(Self {
            weights: normalize(&raw)?,
            origin: PopulationOrigin::Boltzmann {
                temperature_k,
                energy_kind: first.kind,
                energies: energies.to_vec(),
                assumptions: assumptions.into(),
            },
        })
    }

    pub fn weights(&self) -> &[f64] {
        &self.weights
    }
    pub fn origin(&self) -> &PopulationOrigin {
        &self.origin
    }
    pub fn is_boltzmann(&self) -> bool {
        matches!(self.origin, PopulationOrigin::Boltzmann { .. })
    }
}

impl Ensemble {
    /// Requires declared metadata for every record, including unsuccessful
    /// geometries. Missing states/energies do not become equal populations.
    pub fn boltzmann_population(
        &self,
        temperature_k: f64,
        assumptions: &str,
    ) -> Result<Population, String> {
        let energies = self
            .conformers
            .iter()
            .map(|c| {
                c.energy
                    .clone()
                    .ok_or_else(|| format!("conformer {} has no declared energy metadata", c.index))
            })
            .collect::<Result<Vec<_>, _>>()?;
        Population::boltzmann(&energies, temperature_k, assumptions)
    }
}

fn normalize(values: &[f64]) -> Result<Vec<f64>, String> {
    if values.is_empty() || values.iter().any(|v| !v.is_finite() || *v < 0.0) {
        return Err("weights must be finite, non-negative and nonempty".into());
    }
    let scale = values.iter().copied().fold(0.0_f64, f64::max);
    if scale == 0.0 {
        return Err("weights must have a positive total".into());
    }
    let scaled: Vec<f64> = values.iter().map(|w| w / scale).collect();
    // Compensate accumulation so thousands of unequal conformer weights do not
    // progressively consume the normalization error budget.
    let mut sum = 0.0;
    let mut correction = 0.0;
    for value in &scaled {
        let adjusted = value - correction;
        let next = sum + adjusted;
        correction = (next - sum) - adjusted;
        sum = next;
    }
    Ok(scaled.iter().map(|w| w / sum).collect())
}
