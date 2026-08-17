//! Applicability domain and predictive uncertainty for a fitted model.
//!
//! A regression can return a number for any input; whether that number deserves
//! trust is a separate question. This module carries the training-set geometry
//! needed to answer it honestly:
//!
//! * **Leverage** `h = x'(X'X)⁻¹x` measures how far a candidate sits from the
//!   centre of the training design, in the metric the fit itself defines. The
//!   conventional warning threshold is `h* = 3p/n`; beyond it a prediction is an
//!   extrapolation supported by little or no nearby training data. This is the
//!   same criterion the project's own Study 003 pre-registration applies.
//! * **Prediction interval** `ŷ ± t(0.975, n−p)·s·√(1+h)` widens automatically
//!   with leverage, so a distant ligand is reported with a correspondingly
//!   honest error bar rather than a falsely precise one.

use crate::model::MODEL_FEATURE_COUNT;
use serde::{Deserialize, Serialize};

const EPSILON: f64 = 1.0e-12;

/// Training-set geometry persisted alongside a fitted model.
///
/// The design is `[1, (x_j − mean_j)/scale_j …]` over the selected columns, in
/// the same standardized frame the fit solved in, so leverage computed here is
/// exactly the leverage of the fitted regression.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TrainingGeometry {
    /// Selected raw feature columns, in design order after the intercept.
    pub feature_indices: Vec<usize>,
    /// Per-selected-column training mean used to standardize.
    pub means: Vec<f64>,
    /// Per-selected-column training scale used to standardize.
    pub scales: Vec<f64>,
    /// `(X'X)⁻¹` in the standardized design frame, size `parameters²`.
    pub xtx_inverse: Vec<Vec<f64>>,
    /// Training observations `n`.
    pub observations: usize,
    /// Fitted parameters `p`, counting the intercept.
    pub parameters: usize,
    /// Residual standard error `s = √(RSS/(n−p))`, the irreducible scatter.
    pub residual_standard_error: f64,
    /// Conventional warning leverage `h* = 3p/n`.
    pub warning_leverage: f64,
}

impl TrainingGeometry {
    /// Standardized design vector for one expanded feature row.
    #[must_use]
    pub fn design_vector(&self, expanded: &[f32; MODEL_FEATURE_COUNT]) -> Vec<f64> {
        let mut design = Vec::with_capacity(self.parameters);
        design.push(1.0);
        for ((index, mean), scale) in self
            .feature_indices
            .iter()
            .zip(&self.means)
            .zip(&self.scales)
        {
            let value = f64::from(expanded[*index]);
            design.push((value - mean) / scale.max(EPSILON));
        }
        design
    }

    /// Leverage `h = x'(X'X)⁻¹x` of a candidate against the training design.
    #[must_use]
    pub fn leverage(&self, expanded: &[f32; MODEL_FEATURE_COUNT]) -> Option<f64> {
        let design = self.design_vector(expanded);
        if self.xtx_inverse.len() != design.len()
            || self.xtx_inverse.iter().any(|row| row.len() != design.len())
        {
            return None;
        }
        let mut leverage = 0.0;
        for (left, row) in design.iter().zip(&self.xtx_inverse) {
            for (right, value) in design.iter().zip(row) {
                leverage += left * value * right;
            }
        }
        leverage.is_finite().then_some(leverage)
    }

    /// Residual degrees of freedom `n − p`.
    #[must_use]
    pub fn degrees_of_freedom(&self) -> Option<usize> {
        self.observations
            .checked_sub(self.parameters)
            .filter(|df| *df > 0)
    }

    /// Two-sided 95 % Student-t multiplier for this fit's residual df.
    #[must_use]
    pub fn t_multiplier(&self) -> Option<f64> {
        self.degrees_of_freedom()
            .map(|df| student_t_two_sided_quantile(0.05, df as f64))
    }

    /// 95 % prediction interval for a new observation: `ŷ ± t·s·√(1+h)`.
    ///
    /// This is the interval for an individual future measurement, so it carries
    /// both the parameter uncertainty (through `h`) and the residual scatter
    /// (through `s`) — the honest band to quote for "what will this ligand do".
    #[must_use]
    pub fn prediction_interval(
        &self,
        prediction: f64,
        expanded: &[f32; MODEL_FEATURE_COUNT],
    ) -> Option<(f64, f64)> {
        let leverage = self.leverage(expanded)?;
        let multiplier = self.t_multiplier()?;
        let half_width =
            multiplier * self.residual_standard_error * (1.0 + leverage.max(0.0)).sqrt();
        half_width
            .is_finite()
            .then_some((prediction - half_width, prediction + half_width))
    }

    /// 95 % confidence interval for the fitted mean response: `ŷ ± t·s·√h`.
    #[must_use]
    pub fn confidence_interval(
        &self,
        prediction: f64,
        expanded: &[f32; MODEL_FEATURE_COUNT],
    ) -> Option<(f64, f64)> {
        let leverage = self.leverage(expanded)?;
        let multiplier = self.t_multiplier()?;
        let half_width = multiplier * self.residual_standard_error * leverage.max(0.0).sqrt();
        half_width
            .is_finite()
            .then_some((prediction - half_width, prediction + half_width))
    }
}

/// Invert a square matrix by Gauss-Jordan elimination with partial pivoting.
pub(crate) fn invert_matrix(matrix: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
    let size = matrix.len();
    if size == 0 || matrix.iter().any(|row| row.len() != size) {
        return Err("matrix to invert must be square and non-empty".into());
    }
    let mut work = matrix.to_vec();
    let mut inverse = (0..size)
        .map(|row| {
            (0..size)
                .map(|column| if row == column { 1.0 } else { 0.0 })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    for pivot in 0..size {
        let best = (pivot..size)
            .max_by(|&left, &right| work[left][pivot].abs().total_cmp(&work[right][pivot].abs()))
            .ok_or_else(|| "empty matrix".to_string())?;
        if work[best][pivot].abs() <= EPSILON {
            return Err("training design matrix is singular; leverage is undefined".into());
        }
        work.swap(pivot, best);
        inverse.swap(pivot, best);
        let divisor = work[pivot][pivot];
        for value in &mut work[pivot] {
            *value /= divisor;
        }
        for value in &mut inverse[pivot] {
            *value /= divisor;
        }
        for row in 0..size {
            if row == pivot {
                continue;
            }
            let factor = work[row][pivot];
            if factor == 0.0 {
                continue;
            }
            for column in 0..size {
                work[row][column] -= factor * work[pivot][column];
                inverse[row][column] -= factor * inverse[pivot][column];
            }
        }
    }
    Ok(inverse)
}

/// Natural log of the gamma function (Lanczos approximation, g = 7, n = 9).
fn ln_gamma(x: f64) -> f64 {
    const COEFFICIENTS: [f64; 9] = [
        0.9999999999998099,
        676.5203681218851,
        -1259.1392167224028,
        771.3234287776531,
        -176.6150291621406,
        12.507343278686905,
        -0.13857109526572012,
        9.984369578019572e-06,
        1.5056327351493116e-07,
    ];
    if x < 0.5 {
        // Reflection keeps the series in its convergent range.
        return (std::f64::consts::PI / (std::f64::consts::PI * x).sin()).ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut series = COEFFICIENTS[0];
    for (index, coefficient) in COEFFICIENTS.iter().enumerate().skip(1) {
        series += coefficient / (x + index as f64);
    }
    let t = x + 7.5;
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + series.ln()
}

/// Continued fraction for the incomplete beta function (modified Lentz).
fn beta_continued_fraction(a: f64, b: f64, x: f64) -> f64 {
    const MAX_ITERATIONS: usize = 300;
    const TINY: f64 = 1.0e-30;
    let qab = a + b;
    let qap = a + 1.0;
    let qam = a - 1.0;
    let mut c = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < TINY {
        d = TINY;
    }
    d = 1.0 / d;
    let mut result = d;
    for m in 1..=MAX_ITERATIONS {
        let m = m as f64;
        let two_m = 2.0 * m;
        // Even step.
        let numerator = m * (b - m) * x / ((qam + two_m) * (a + two_m));
        d = 1.0 + numerator * d;
        if d.abs() < TINY {
            d = TINY;
        }
        c = 1.0 + numerator / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        result *= d * c;
        // Odd step.
        let numerator = -(a + m) * (qab + m) * x / ((a + two_m) * (qap + two_m));
        d = 1.0 + numerator * d;
        if d.abs() < TINY {
            d = TINY;
        }
        c = 1.0 + numerator / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        let delta = d * c;
        result *= delta;
        if (delta - 1.0).abs() < 3.0e-16 {
            break;
        }
    }
    result
}

/// Regularized incomplete beta function `I_x(a, b)`.
fn regularized_incomplete_beta(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let front =
        (ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * x.ln() + b * (1.0 - x).ln()).exp();
    if x < (a + 1.0) / (a + b + 2.0) {
        front * beta_continued_fraction(a, b, x) / a
    } else {
        1.0 - front * beta_continued_fraction(b, a, 1.0 - x) / b
    }
}

/// Two-sided tail probability `P(|T| > t)` for Student's t with `df` degrees of
/// freedom.
fn student_t_two_sided_tail(t: f64, df: f64) -> f64 {
    if !t.is_finite() || t <= 0.0 {
        return 1.0;
    }
    regularized_incomplete_beta(0.5 * df, 0.5, df / (df + t * t))
}

/// The `t` for which `P(|T| > t) = alpha`, found by bisection on the tail.
///
/// Bisection is used rather than an inverse-beta expansion because it is short,
/// obviously correct, and this is called once per fit — never in a hot loop.
#[must_use]
pub fn student_t_two_sided_quantile(alpha: f64, df: f64) -> f64 {
    if !(0.0..1.0).contains(&alpha) || alpha <= 0.0 || df <= 0.0 {
        return f64::NAN;
    }
    let (mut low, mut high) = (0.0_f64, 1.0_f64);
    while student_t_two_sided_tail(high, df) > alpha && high < 1.0e6 {
        high *= 2.0;
    }
    for _ in 0..200 {
        let middle = 0.5 * (low + high);
        if student_t_two_sided_tail(middle, df) > alpha {
            low = middle;
        } else {
            high = middle;
        }
    }
    0.5 * (low + high)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inverts_a_known_matrix() {
        let matrix = vec![vec![4.0, 7.0], vec![2.0, 6.0]];
        let inverse = invert_matrix(&matrix).unwrap();
        // Closed form: 1/10 * [[6, -7], [-2, 4]]
        assert!((inverse[0][0] - 0.6).abs() < 1e-12);
        assert!((inverse[0][1] + 0.7).abs() < 1e-12);
        assert!((inverse[1][0] + 0.2).abs() < 1e-12);
        assert!((inverse[1][1] - 0.4).abs() < 1e-12);
    }

    #[test]
    fn rejects_a_singular_matrix() {
        let singular = vec![vec![1.0, 2.0], vec![2.0, 4.0]];
        assert!(invert_matrix(&singular).is_err());
    }

    #[test]
    fn reproduces_published_student_t_quantiles() {
        // Two-sided 95 % critical values from standard t tables.
        for (df, expected) in [
            (1.0, 12.706),
            (2.0, 4.303),
            (5.0, 2.571),
            (10.0, 2.228),
            (30.0, 2.042),
            (60.0, 2.000),
            (120.0, 1.980),
        ] {
            let actual = student_t_two_sided_quantile(0.05, df);
            assert!(
                (actual - expected).abs() < 5.0e-3,
                "t(0.975, {df}) = {actual}, expected {expected}"
            );
        }
        // Large df converges on the normal critical value.
        assert!((student_t_two_sided_quantile(0.05, 1.0e6) - 1.959_964).abs() < 1e-3);
    }

    fn geometry() -> TrainingGeometry {
        // One selected feature standardized to mean 0, scale 1, with an
        // orthonormal design: (X'X)^-1 = diag(1/n, 1/n) for centred data.
        TrainingGeometry {
            feature_indices: vec![3],
            means: vec![0.0],
            scales: vec![1.0],
            xtx_inverse: vec![vec![0.1, 0.0], vec![0.0, 0.1]],
            observations: 10,
            parameters: 2,
            residual_standard_error: 0.5,
            warning_leverage: 0.6,
        }
    }

    #[test]
    fn leverage_grows_with_distance_from_the_training_centre() {
        let geometry = geometry();
        let mut centre = [0.0_f32; MODEL_FEATURE_COUNT];
        centre[3] = 0.0;
        let mut distant = [0.0_f32; MODEL_FEATURE_COUNT];
        distant[3] = 3.0;
        let at_centre = geometry.leverage(&centre).unwrap();
        let far = geometry.leverage(&distant).unwrap();
        // At the centroid only the intercept contributes: h = 1/n.
        assert!((at_centre - 0.1).abs() < 1e-12);
        // Three standardized units out: h = 0.1 + 0.1*9.
        assert!((far - 1.0).abs() < 1e-12);
        assert!(far > geometry.warning_leverage);
    }

    #[test]
    fn prediction_interval_widens_with_leverage() {
        let geometry = geometry();
        let mut centre = [0.0_f32; MODEL_FEATURE_COUNT];
        centre[3] = 0.0;
        let mut distant = [0.0_f32; MODEL_FEATURE_COUNT];
        distant[3] = 3.0;
        let (near_low, near_high) = geometry.prediction_interval(1.0, &centre).unwrap();
        let (far_low, far_high) = geometry.prediction_interval(1.0, &distant).unwrap();
        assert!(near_high - near_low > 0.0);
        assert!(
            far_high - far_low > near_high - near_low,
            "an extrapolated ligand must carry a wider interval"
        );
        // The prediction interval always contains the mean-response interval.
        let (conf_low, conf_high) = geometry.confidence_interval(1.0, &centre).unwrap();
        assert!(near_low < conf_low && near_high > conf_high);
    }

    #[test]
    fn degrees_of_freedom_guard_against_saturated_fits() {
        let mut saturated = geometry();
        saturated.observations = 2;
        saturated.parameters = 2;
        assert_eq!(saturated.degrees_of_freedom(), None);
        assert!(saturated.t_multiplier().is_none());
    }
}
