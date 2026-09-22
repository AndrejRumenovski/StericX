/// Boltzmann constant in joules per kelvin.
pub const BOLTZMANN_CONSTANT_J_K: f64 = 1.380_649e-23;
/// Planck constant in joule-seconds.
pub const PLANCK_CONSTANT_J_S: f64 = 6.626_070_15e-34;
/// Molar gas constant in kcal mol⁻¹ K⁻¹.
pub const GAS_CONSTANT_KCAL_MOL_K: f64 = 0.001_987_204_258_640_831_6;

/// Failure to represent an absolute transition-state-theory rate as finite,
/// nonzero f32. Log rates remain available for valid extreme inputs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RateConstantError {
    InvalidInput,
    Underflow,
    Overflow,
}

impl std::fmt::Display for RateConstantError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::InvalidInput => {
                "activation free energy must be finite and temperature positive and finite"
            }
            Self::Underflow => "rate rounds to zero in f32; use calculate_log_rate_constant",
            Self::Overflow => "rate overflows f32; use calculate_log_rate_constant",
        })
    }
}

impl std::error::Error for RateConstantError {}

/// Predicted enantiomeric product distribution at one temperature.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ProductRatio {
    /// Temperature in kelvin.
    pub temp_k: f32,
    /// Percentage assigned to the R product.
    pub percent_r: f32,
    /// Percentage assigned to the S product.
    pub percent_s: f32,
    /// Absolute predicted enantiomeric excess.
    pub ee_percent: f32,
}

/// Transition-state-theory link from barriers to rates and selectivity.
#[derive(Clone, Copy, Debug, Default)]
pub struct EyringKineticLink;

impl EyringKineticLink {
    /// Calculates `ln(k / s^-1)` for an absolute activation free energy ΔG‡.
    /// Assumes unit transmission coefficient and a unimolecular rate prefactor.
    /// A barrier difference ΔΔG‡ alone does not determine an absolute rate.
    pub fn calculate_log_rate_constant(
        activation_free_energy_kcal: f32,
        temp_k: f32,
    ) -> Result<f64, RateConstantError> {
        if !valid_inputs(activation_free_energy_kcal, temp_k) {
            return Err(RateConstantError::InvalidInput);
        }
        let temperature = f64::from(temp_k);
        Ok(
            (BOLTZMANN_CONSTANT_J_K / PLANCK_CONSTANT_J_S).ln() + temperature.ln()
                - f64::from(activation_free_energy_kcal) / (GAS_CONSTANT_KCAL_MOL_K * temperature),
        )
    }

    /// Calculates an absolute Eyring rate, reporting f32 underflow/overflow.
    pub fn calculate_rate_constant_checked(
        activation_free_energy_kcal: f32,
        temp_k: f32,
    ) -> Result<f32, RateConstantError> {
        let rate =
            Self::calculate_log_rate_constant(activation_free_energy_kcal, temp_k)?.exp() as f32;
        if rate == 0.0 {
            Err(RateConstantError::Underflow)
        } else if !rate.is_finite() {
            Err(RateConstantError::Overflow)
        } else {
            Ok(rate)
        }
    }

    /// Calculates `k = (k_B T / h) exp(-ΔG‡ / RT)` from an absolute barrier.
    ///
    /// Units are kcal/mol, kelvin, and s^-1; transmission coefficient is one.
    /// Invalid inputs return NaN, and f32 range limits return zero or infinity.
    /// Use the checked or log-rate API when these limits must be distinguished.
    #[must_use]
    pub fn calculate_rate_constant(activation_free_energy_kcal: f32, temp_k: f32) -> f32 {
        match Self::calculate_rate_constant_checked(activation_free_energy_kcal, temp_k) {
            Ok(rate) => rate,
            Err(RateConstantError::InvalidInput) => f32::NAN,
            Err(RateConstantError::Underflow) => 0.0,
            Err(RateConstantError::Overflow) => f32::INFINITY,
        }
    }

    /// Returns normalized `(major_percent, minor_percent)` from ΔΔG‡.
    ///
    /// This evaluates `major/minor = exp(|ΔΔG‡| / RT)` with a numerically
    /// stable logistic form. The percentages always sum to 100; the sign of
    /// `ddg_kcal` selects identity rather than changing the major/minor values.
    /// Invalid inputs return `(NaN, NaN)`.
    #[must_use]
    pub fn calculate_enantiomeric_ratio(ddg_kcal: f32, temp_k: f32) -> (f32, f32) {
        if !valid_inputs(ddg_kcal, temp_k) {
            return (f32::NAN, f32::NAN);
        }
        let log_ratio = f64::from(ddg_kcal).abs() / (GAS_CONSTANT_KCAL_MOL_K * f64::from(temp_k));
        let minor_ratio = (-log_ratio).exp();
        let minor_percent = (100.0 * minor_ratio / (1.0 + minor_ratio)) as f32;
        (100.0 - minor_percent, minor_percent)
    }

    /// Calculates absolute enantiomeric excess from ΔΔG‡.
    #[must_use]
    pub fn calculate_enantiomeric_excess(ddg_kcal: f32, temp_k: f32) -> f32 {
        if !valid_inputs(ddg_kcal, temp_k) {
            return f32::NAN;
        }
        let half_log_ratio =
            f64::from(ddg_kcal).abs() / (2.0 * GAS_CONSTANT_KCAL_MOL_K * f64::from(temp_k));
        (100.0 * half_log_ratio.tanh()) as f32
    }

    /// Converts ΔΔG‡ into an R:S product distribution.
    ///
    /// The sign convention is `ΔΔG‡ = G‡(S) - G‡(R)`, so positive values favor
    /// the R product. Equal barriers produce a 50:50 mixture.
    #[must_use]
    pub fn product_ratio(ddg_kcal: f32, temp_k: f32) -> ProductRatio {
        if !valid_inputs(ddg_kcal, temp_k) {
            return ProductRatio {
                temp_k,
                percent_r: f32::NAN,
                percent_s: f32::NAN,
                ee_percent: f32::NAN,
            };
        }
        let (major_percent, minor_percent) = Self::calculate_enantiomeric_ratio(ddg_kcal, temp_k);
        let (percent_r, percent_s) = if ddg_kcal >= 0.0 {
            (major_percent, minor_percent)
        } else {
            (minor_percent, major_percent)
        };
        ProductRatio {
            temp_k,
            percent_r,
            percent_s,
            ee_percent: Self::calculate_enantiomeric_excess(ddg_kcal, temp_k),
        }
    }

    /// Simulates an inclusive, evenly spaced temperature range.
    ///
    /// Invalid ranges or a non-positive step return an empty vector.
    #[must_use]
    pub fn product_ratios_over_range(
        ddg: f32,
        start_k: f32,
        end_k: f32,
        step_k: f32,
    ) -> Vec<ProductRatio> {
        if !start_k.is_finite()
            || !end_k.is_finite()
            || !step_k.is_finite()
            || start_k <= 0.0
            || end_k < start_k
            || step_k <= 0.0
        {
            return Vec::new();
        }
        let count = ((end_k - start_k) / step_k).floor() as usize;
        (0..=count)
            .map(|index| Self::product_ratio(ddg, start_k + index as f32 * step_k))
            .collect()
    }
}

#[inline]
fn valid_inputs(ddg_kcal: f32, temp_k: f32) -> bool {
    ddg_kcal.is_finite() && temp_k.is_finite() && temp_k > 0.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_barrier_uses_eyring_prefactor() {
        let rate = EyringKineticLink::calculate_rate_constant(0.0, 298.15);
        assert!((rate / 6.21e12 - 1.0).abs() < 0.01);
    }

    #[test]
    fn zero_ddg_is_racemic() {
        let ratio = EyringKineticLink::product_ratio(0.0, 298.15);
        assert_eq!(ratio.percent_r, 50.0);
        assert_eq!(ratio.percent_s, 50.0);
        assert_eq!(ratio.ee_percent, 0.0);
    }

    #[test]
    fn positive_ddg_favors_r() {
        let ratio = EyringKineticLink::product_ratio(1.0, 298.15);
        assert!(ratio.percent_r > ratio.percent_s);
        assert!((ratio.percent_r + ratio.percent_s - 100.0).abs() < f32::EPSILON);
    }

    #[test]
    fn ddg_1_82_gives_about_95_percent_major_at_room_temperature() {
        let (major, minor) = EyringKineticLink::calculate_enantiomeric_ratio(1.82, 298.15);
        let ee = EyringKineticLink::calculate_enantiomeric_excess(1.82, 298.15);

        assert!((major - 95.6).abs() < 0.2);
        assert!((minor - 4.4).abs() < 0.2);
        assert!((ee - 91.2).abs() < 0.3);
    }

    #[test]
    fn approximately_2_17_kcal_gives_95_percent_ee_at_room_temperature() {
        let ee = EyringKineticLink::calculate_enantiomeric_excess(2.17, 298.15);
        assert!((ee - 95.0).abs() < 0.2);
    }

    #[test]
    fn builds_inclusive_temperature_series() {
        let ratios = EyringKineticLink::product_ratios_over_range(1.0, 280.0, 300.0, 10.0);
        assert_eq!(ratios.len(), 3);
        assert_eq!(ratios[2].temp_k, 300.0);
    }
}
