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

use crate::model::{FeatureDomain, MODEL_FEATURE_COUNT};
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
    /// Training observations in the standardized descriptor frame, one row per
    /// observation and one column per selected descriptor (no intercept).
    ///
    /// Optional so models written before nearest-neighbour scoring existed
    /// still load; without it a distance simply cannot be reported.
    #[serde(default)]
    pub standardized_training_points: Vec<Vec<f64>>,
    /// Nearest-neighbour spacing of the training set, and the boundary derived
    /// from it. Absent when there are too few points to measure a spacing.
    #[serde(default)]
    pub neighbor_calibration: Option<NeighborCalibration>,
}

/// How densely the training set samples its own descriptor space.
///
/// Every value here is measured from the training observations: for each
/// point, the Euclidean distance in standardized descriptor space to the
/// nearest *other* training point. No constant is chosen by hand.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct NeighborCalibration {
    pub mean: f64,
    pub standard_deviation: f64,
    pub median: f64,
    pub maximum: f64,
    /// Boundary a candidate is compared against.
    ///
    /// Set to [`Self::maximum`]: a candidate is treated as sitting inside the
    /// sampled region when it is no farther from the training set than the
    /// training set's own sparsest point is from its nearest neighbour. This
    /// has no free parameter, but it is a *permissive* boundary — it is set by
    /// the loosest part of the training set, so `mean` and `standard_deviation`
    /// are recorded too for consumers who want a stricter rule.
    pub threshold: f64,
    /// Exact derivation, carried with the model so a reader need not guess.
    pub rule: String,
}

impl NeighborCalibration {
    /// Measures the training set's own nearest-neighbour spacing.
    ///
    /// Returns `None` for fewer than two points, where "distance to the nearest
    /// other point" is undefined.
    #[must_use]
    pub fn from_points(points: &[Vec<f64>]) -> Option<Self> {
        if points.len() < 2 {
            return None;
        }
        let mut distances = Vec::with_capacity(points.len());
        for (index, point) in points.iter().enumerate() {
            let nearest = points
                .iter()
                .enumerate()
                .filter(|(other, _)| *other != index)
                .filter_map(|(_, other)| euclidean(point, other))
                .fold(f64::INFINITY, f64::min);
            if !nearest.is_finite() {
                return None;
            }
            distances.push(nearest);
        }
        distances.sort_by(f64::total_cmp);
        let count = distances.len() as f64;
        let mean = distances.iter().sum::<f64>() / count;
        let variance = distances
            .iter()
            .map(|distance| (distance - mean).powi(2))
            .sum::<f64>()
            / count;
        let median = if distances.len() % 2 == 0 {
            (distances[distances.len() / 2 - 1] + distances[distances.len() / 2]) / 2.0
        } else {
            distances[distances.len() / 2]
        };
        let maximum = *distances.last()?;
        Some(Self {
            mean,
            standard_deviation: variance.sqrt(),
            median,
            maximum,
            threshold: maximum,
            rule: "threshold = max over training points of the distance to the nearest other \
                   training point, in standardized descriptor space"
                .to_owned(),
        })
    }
}

/// Where a candidate sits relative to the training set.
///
/// Each state is a stated combination of two measured quantities, so no
/// severity grade is invented: the range check answers "is this inside the box
/// the training set spans", and the neighbour distance answers "is this inside
/// a part of that box the training set actually sampled".
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum DomainVerdict {
    /// Every descriptor inside its training range, and no farther from the
    /// training set than the calibrated neighbour spacing allows.
    Interpolation,
    /// Every descriptor inside its training range, but farther from any
    /// training point than that spacing — a gap the training set did not cover.
    SparseInterpolation,
    /// At least one descriptor outside the range the model was trained on.
    Extrapolation,
    /// The model carries no training geometry, so no verdict can be reached.
    Unknown,
}

impl DomainVerdict {
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::Interpolation => "interpolation",
            Self::SparseInterpolation => "sparse_interpolation",
            Self::Extrapolation => "extrapolation",
            Self::Unknown => "unknown",
        }
    }
}

/// One descriptor that falls outside the range the model was trained on.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct DescriptorExceedance {
    pub feature: String,
    pub value: f64,
    pub training_minimum: f64,
    pub training_maximum: f64,
    /// How far outside, as a fraction of the training range width. Zero at the
    /// boundary; one means a full training range beyond it.
    pub normalized_exceedance: f64,
}

/// Structured applicability information for one candidate.
///
/// Computed from the candidate's descriptors alone. The prediction is not an
/// input, so a favourable-looking number can never make a ligand look more
/// in-domain than it is.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct ApplicabilityAssessment {
    pub verdict: DomainVerdict,
    /// Euclidean distance in standardized descriptor space to the closest
    /// training observation.
    pub nearest_training_distance: Option<f64>,
    /// The calibrated boundary that distance is compared against.
    pub nearest_training_threshold: Option<f64>,
    /// `distance / threshold`; above one is outside the sampled region.
    pub nearest_training_ratio: Option<f64>,
    /// Leverage `h = x'(X'X)⁻¹x` against the training design.
    pub leverage: Option<f64>,
    /// `h / h*` for the conventional warning leverage `h* = 3p/n`.
    pub leverage_ratio: Option<f64>,
    /// Mahalanobis distance in the training descriptor covariance, when that
    /// covariance is well enough determined to be worth reporting.
    pub mahalanobis_distance: Option<f64>,
    /// Why Mahalanobis was not reported, when it was not.
    pub mahalanobis_unavailable: Option<String>,
    /// Descriptors outside their training range, empty when all are inside.
    pub outside_range: Vec<DescriptorExceedance>,
    /// Largest `normalized_exceedance`, or zero when every descriptor is inside.
    pub maximum_extrapolation: f64,
}

fn euclidean(left: &[f64], right: &[f64]) -> Option<f64> {
    if left.len() != right.len() {
        return None;
    }
    let squared = left
        .iter()
        .zip(right)
        .map(|(left, right)| (left - right).powi(2))
        .sum::<f64>();
    squared.sqrt().is_finite().then_some(squared.sqrt())
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

    /// Standardized descriptor coordinates of a candidate, without the
    /// intercept column.
    #[must_use]
    pub fn standardized_point(&self, expanded: &[f32; MODEL_FEATURE_COUNT]) -> Vec<f64> {
        self.design_vector(expanded).split_off(1)
    }

    /// Distance in standardized descriptor space to the closest training point.
    #[must_use]
    pub fn nearest_training_distance(&self, expanded: &[f32; MODEL_FEATURE_COUNT]) -> Option<f64> {
        if self.standardized_training_points.is_empty() {
            return None;
        }
        let point = self.standardized_point(expanded);
        self.standardized_training_points
            .iter()
            .filter_map(|training| euclidean(&point, training))
            .fold(None, |best: Option<f64>, distance| {
                Some(best.map_or(distance, |best| best.min(distance)))
            })
    }

    /// Mahalanobis distance of a candidate in the training covariance.
    ///
    /// The standardized descriptor columns are centred on the training means by
    /// construction, so `Z'Z` is block diagonal and leverage decomposes as
    /// `h = 1/n + z'S⁻¹z` with `S = (n−1)·Cov`. The Mahalanobis distance is
    /// therefore `√((n−1)(h − 1/n))`, recovered exactly from the geometry
    /// already stored — no second covariance estimate, and nothing that can
    /// disagree with the leverage the same model reports.
    ///
    /// Returns the reason instead of a number when that decomposition cannot be
    /// trusted: a covariance estimated from too few observations, or a stored
    /// matrix whose block structure has been disturbed, gives a figure that
    /// looks authoritative and is not.
    pub fn mahalanobis_distance(
        &self,
        expanded: &[f32; MODEL_FEATURE_COUNT],
    ) -> Result<f64, String> {
        let descriptors = self.parameters.saturating_sub(1);
        if descriptors == 0 {
            return Err("model has no descriptor columns".into());
        }
        let observations = self.observations;
        // The sample covariance of k descriptors is singular unless n − 1 ≥ k.
        if observations < descriptors + 2 {
            return Err(format!(
                "covariance of {descriptors} descriptor(s) is not estimable from {observations} \
                 observations"
            ));
        }
        if self.xtx_inverse.len() != self.parameters {
            return Err("stored (X'X)^-1 does not match the parameter count".into());
        }
        // Verify the block structure the decomposition relies on rather than
        // assuming it: a hand-edited or differently-centred matrix must not be
        // read as if it were centred.
        let expected_intercept = 1.0 / observations as f64;
        if (self.xtx_inverse[0][0] - expected_intercept).abs()
            > 1.0e-6 * expected_intercept.max(1.0)
        {
            return Err("training design is not centred, so leverage does not decompose".into());
        }
        if self.xtx_inverse[0]
            .iter()
            .skip(1)
            .any(|value| value.abs() > 1.0e-6)
        {
            return Err("training design is not centred, so leverage does not decompose".into());
        }
        let leverage = self
            .leverage(expanded)
            .ok_or_else(|| "leverage could not be computed".to_string())?;
        let squared = (observations - 1) as f64 * (leverage - expected_intercept);
        if !squared.is_finite() || squared < -1.0e-9 {
            return Err("leverage decomposition produced a negative squared distance".into());
        }
        Ok(squared.max(0.0).sqrt())
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

/// Scores a candidate against a model's training distribution.
///
/// Takes only descriptors. There is deliberately no way to pass a prediction
/// in, so applicability can never be influenced by whether the predicted value
/// happens to look good.
#[must_use]
pub fn assess_applicability(
    geometry: Option<&TrainingGeometry>,
    ranges: &[FeatureDomain],
    selected: &[usize],
    expanded: &[f32; MODEL_FEATURE_COUNT],
) -> ApplicabilityAssessment {
    let outside_range = selected
        .iter()
        .zip(ranges)
        .filter_map(|(column, range)| {
            let value = f64::from(expanded[*column]);
            let width = range.maximum - range.minimum;
            let overshoot = if value > range.maximum {
                value - range.maximum
            } else if value < range.minimum {
                range.minimum - value
            } else {
                return None;
            };
            // A zero-width training range means every training value was
            // identical; any departure is unbounded rather than a fraction.
            let normalized = if width > 0.0 {
                overshoot / width
            } else {
                f64::INFINITY
            };
            Some(DescriptorExceedance {
                feature: range.feature.clone(),
                value,
                training_minimum: range.minimum,
                training_maximum: range.maximum,
                normalized_exceedance: normalized,
            })
        })
        .collect::<Vec<_>>();
    let maximum_extrapolation = outside_range
        .iter()
        .map(|exceedance| exceedance.normalized_exceedance)
        .fold(0.0_f64, f64::max);

    let Some(geometry) = geometry else {
        return ApplicabilityAssessment {
            verdict: if outside_range.is_empty() {
                DomainVerdict::Unknown
            } else {
                // A range violation is decidable without any geometry.
                DomainVerdict::Extrapolation
            },
            nearest_training_distance: None,
            nearest_training_threshold: None,
            nearest_training_ratio: None,
            leverage: None,
            leverage_ratio: None,
            mahalanobis_distance: None,
            mahalanobis_unavailable: Some("model carries no training geometry".to_owned()),
            outside_range,
            maximum_extrapolation,
        };
    };

    let distance = geometry.nearest_training_distance(expanded);
    let threshold = geometry
        .neighbor_calibration
        .as_ref()
        .map(|calibration| calibration.threshold);
    let ratio = match (distance, threshold) {
        (Some(distance), Some(threshold)) if threshold > 0.0 => Some(distance / threshold),
        _ => None,
    };
    let leverage = geometry.leverage(expanded);
    let leverage_ratio = match leverage {
        Some(leverage) if geometry.warning_leverage > 0.0 => {
            Some(leverage / geometry.warning_leverage)
        }
        _ => None,
    };
    let (mahalanobis_distance, mahalanobis_unavailable) =
        match geometry.mahalanobis_distance(expanded) {
            Ok(distance) => (Some(distance), None),
            Err(reason) => (None, Some(reason)),
        };

    let verdict = if !outside_range.is_empty() {
        DomainVerdict::Extrapolation
    } else {
        match (distance, threshold) {
            (Some(distance), Some(threshold)) if distance > threshold => {
                DomainVerdict::SparseInterpolation
            }
            (Some(_), Some(_)) => DomainVerdict::Interpolation,
            // Inside every range, but the model records no neighbour spacing to
            // check against; say so rather than implying a full verdict.
            _ => DomainVerdict::Unknown,
        }
    };

    ApplicabilityAssessment {
        verdict,
        nearest_training_distance: distance,
        nearest_training_threshold: threshold,
        nearest_training_ratio: ratio,
        leverage,
        leverage_ratio,
        mahalanobis_distance,
        mahalanobis_unavailable,
        outside_range,
        maximum_extrapolation,
    }
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

    /// A deliberately obvious 2-D descriptor space: nine training points on the
    /// integer grid −1..1 in each of two standardized descriptors. Nearest
    /// neighbours are all exactly 1.0 apart, so the calibrated boundary is 1.0
    /// and every expectation below can be checked by eye.
    fn grid_geometry() -> TrainingGeometry {
        let mut points = Vec::new();
        for x in [-1.0_f64, 0.0, 1.0] {
            for y in [-1.0_f64, 0.0, 1.0] {
                points.push(vec![x, y]);
            }
        }
        let observations = points.len();
        let calibration = NeighborCalibration::from_points(&points);
        // Columns 2 and 3 are B1 and B5; centred already, so (X'X) is
        // diag(n, Σx², Σy²) = diag(9, 6, 6).
        TrainingGeometry {
            feature_indices: vec![2, 3],
            means: vec![0.0, 0.0],
            scales: vec![1.0, 1.0],
            xtx_inverse: vec![
                vec![1.0 / observations as f64, 0.0, 0.0],
                vec![0.0, 1.0 / 6.0, 0.0],
                vec![0.0, 0.0, 1.0 / 6.0],
            ],
            observations,
            parameters: 3,
            residual_standard_error: 0.1,
            warning_leverage: 3.0 * 3.0 / observations as f64,
            standardized_training_points: points,
            neighbor_calibration: calibration,
        }
    }

    fn grid_ranges() -> Vec<FeatureDomain> {
        vec![
            FeatureDomain {
                feature: "B1_boltz".into(),
                minimum: -1.0,
                maximum: 1.0,
            },
            FeatureDomain {
                feature: "B5_boltz".into(),
                minimum: -1.0,
                maximum: 1.0,
            },
        ]
    }

    fn at(b1: f32, b5: f32) -> [f32; MODEL_FEATURE_COUNT] {
        let mut expanded = [0.0_f32; MODEL_FEATURE_COUNT];
        expanded[0] = 1.0;
        expanded[2] = b1;
        expanded[3] = b5;
        expanded
    }

    fn assess(b1: f32, b5: f32) -> ApplicabilityAssessment {
        assess_applicability(Some(&grid_geometry()), &grid_ranges(), &[2, 3], &at(b1, b5))
    }

    #[test]
    fn calibration_measures_the_training_spacing() {
        let calibration = grid_geometry().neighbor_calibration.unwrap();

        // Every grid point has an orthogonal neighbour exactly 1.0 away.
        assert!((calibration.mean - 1.0).abs() < 1.0e-12);
        assert!((calibration.maximum - 1.0).abs() < 1.0e-12);
        assert!((calibration.median - 1.0).abs() < 1.0e-12);
        assert!(calibration.standard_deviation < 1.0e-12);
        assert!((calibration.threshold - calibration.maximum).abs() < 1.0e-12);
        assert!(calibration.rule.contains("nearest other"));
    }

    #[test]
    fn a_training_point_is_zero_distance_and_interpolating() {
        let assessment = assess(0.0, 0.0);

        assert_eq!(assessment.verdict, DomainVerdict::Interpolation);
        assert_eq!(assessment.nearest_training_distance, Some(0.0));
        assert!(assessment.outside_range.is_empty());
        assert_eq!(assessment.maximum_extrapolation, 0.0);
    }

    #[test]
    fn a_point_between_training_points_interpolates() {
        // Midway between (0,0) and (1,0): 0.5 away, inside the 1.0 boundary.
        let assessment = assess(0.5, 0.0);

        assert_eq!(assessment.verdict, DomainVerdict::Interpolation);
        assert!((assessment.nearest_training_distance.unwrap() - 0.5).abs() < 1.0e-12);
        assert!(assessment.nearest_training_ratio.unwrap() < 1.0);
    }

    #[test]
    fn a_point_outside_a_descriptor_range_extrapolates() {
        // Half a training range beyond the B1 maximum of 1.0.
        let assessment = assess(2.0, 0.0);

        assert_eq!(assessment.verdict, DomainVerdict::Extrapolation);
        assert_eq!(assessment.outside_range.len(), 1);
        let exceedance = &assessment.outside_range[0];
        assert_eq!(exceedance.feature, "B1_boltz");
        assert_eq!(exceedance.value, 2.0);
        // One full training-range width past the boundary: (2 − 1) / (1 − −1).
        assert!((exceedance.normalized_exceedance - 0.5).abs() < 1.0e-12);
        assert!((assessment.maximum_extrapolation - 0.5).abs() < 1.0e-12);
    }

    #[test]
    fn both_descriptors_outside_are_both_reported() {
        let assessment = assess(3.0, -2.0);

        assert_eq!(assessment.verdict, DomainVerdict::Extrapolation);
        assert_eq!(assessment.outside_range.len(), 2);
        // B1 is 2.0 past the maximum, B5 is 1.0 past the minimum, over a width
        // of 2.0 each.
        assert!((assessment.maximum_extrapolation - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn a_gap_inside_the_range_is_sparse_interpolation() {
        // A 5x5 grid with the middle 3x3 removed: the ring is still spaced 1.0
        // apart, but the centre is now 2.0 from the nearest remaining point —
        // unambiguously a gap rather than a rounding argument.
        let mut points = Vec::new();
        for x in [-2.0_f64, -1.0, 0.0, 1.0, 2.0] {
            for y in [-2.0_f64, -1.0, 0.0, 1.0, 2.0] {
                if x.abs() <= 1.0 && y.abs() <= 1.0 {
                    continue;
                }
                points.push(vec![x, y]);
            }
        }
        let mut geometry = grid_geometry();
        geometry.observations = points.len();
        geometry.neighbor_calibration = NeighborCalibration::from_points(&points);
        geometry.standardized_training_points = points;
        let threshold = geometry.neighbor_calibration.as_ref().unwrap().threshold;
        let ranges = vec![
            FeatureDomain {
                feature: "B1_boltz".into(),
                minimum: -2.0,
                maximum: 2.0,
            },
            FeatureDomain {
                feature: "B5_boltz".into(),
                minimum: -2.0,
                maximum: 2.0,
            },
        ];

        let assessment = assess_applicability(Some(&geometry), &ranges, &[2, 3], &at(0.0, 0.0));

        assert!(assessment.outside_range.is_empty(), "still inside the box");
        let distance = assessment.nearest_training_distance.unwrap();
        assert!(
            distance > threshold,
            "the hole must be wider than the training spacing: {distance} vs {threshold}"
        );
        assert_eq!(assessment.verdict, DomainVerdict::SparseInterpolation);
        assert!(assessment.nearest_training_ratio.unwrap() > 1.0);
    }

    #[test]
    fn mahalanobis_matches_a_hand_computed_value() {
        // Columns have Σx² = 6 over n = 9, so the sample variance is
        // 6/(9−1) = 0.75 and the Mahalanobis distance of (1, 0) is
        // √(1²/0.75) = 1.1547.
        let assessment = assess(1.0, 0.0);
        let expected = (1.0_f64 / 0.75).sqrt();

        let actual = assessment.mahalanobis_distance.expect("estimable here");
        assert!(
            (actual - expected).abs() < 1.0e-9,
            "Mahalanobis {actual}, expected {expected}"
        );
        assert!(assessment.mahalanobis_unavailable.is_none());
    }

    #[test]
    fn mahalanobis_is_declined_when_the_covariance_is_not_estimable() {
        let mut geometry = grid_geometry();
        // Two descriptors need at least four observations for a covariance.
        geometry.observations = 3;
        geometry.xtx_inverse[0][0] = 1.0 / 3.0;

        let error = geometry.mahalanobis_distance(&at(1.0, 0.0)).unwrap_err();

        assert!(error.contains("not estimable"), "{error}");
    }

    #[test]
    fn mahalanobis_is_declined_when_the_design_is_not_centred() {
        let mut geometry = grid_geometry();
        // Break the block structure the decomposition depends on.
        geometry.xtx_inverse[0][1] = 0.4;

        let error = geometry.mahalanobis_distance(&at(1.0, 0.0)).unwrap_err();

        assert!(error.contains("not centred"), "{error}");
    }

    #[test]
    fn a_model_without_geometry_reports_unknown_but_still_checks_ranges() {
        let inside = assess_applicability(None, &grid_ranges(), &[2, 3], &at(0.0, 0.0));
        assert_eq!(inside.verdict, DomainVerdict::Unknown);
        assert!(inside.nearest_training_distance.is_none());
        assert!(inside.mahalanobis_unavailable.is_some());

        // A range violation needs no geometry to be decidable.
        let outside = assess_applicability(None, &grid_ranges(), &[2, 3], &at(5.0, 0.0));
        assert_eq!(outside.verdict, DomainVerdict::Extrapolation);
        assert_eq!(outside.outside_range.len(), 1);
    }

    #[test]
    fn a_model_without_calibration_cannot_claim_interpolation() {
        let mut geometry = grid_geometry();
        geometry.neighbor_calibration = None;

        let assessment =
            assess_applicability(Some(&geometry), &grid_ranges(), &[2, 3], &at(0.0, 0.0));

        assert_eq!(
            assessment.verdict,
            DomainVerdict::Unknown,
            "no measured spacing means no interpolation claim"
        );
        assert!(assessment.nearest_training_distance.is_some());
        assert!(assessment.nearest_training_threshold.is_none());
    }

    #[test]
    fn calibration_needs_at_least_two_points() {
        assert!(NeighborCalibration::from_points(&[]).is_none());
        assert!(NeighborCalibration::from_points(&[vec![0.0, 0.0]]).is_none());
        assert!(NeighborCalibration::from_points(&[vec![0.0], vec![3.0]]).is_some());
    }

    #[test]
    fn a_zero_width_training_range_never_reports_a_finite_fraction() {
        // Every training value identical: any departure is unbounded, not a
        // fraction of a range that does not exist.
        let ranges = vec![FeatureDomain {
            feature: "B1_boltz".into(),
            minimum: 2.0,
            maximum: 2.0,
        }];
        let assessment = assess_applicability(None, &ranges, &[2], &at(2.5, 0.0));

        assert_eq!(assessment.verdict, DomainVerdict::Extrapolation);
        assert!(
            assessment.outside_range[0]
                .normalized_exceedance
                .is_infinite()
        );
    }

    #[test]
    fn applicability_ignores_the_prediction_entirely() {
        // Two candidates at the same descriptor point must score identically no
        // matter what a model would predict for them; the assessment has no
        // access to a prediction at all, and this pins that property.
        let first = assess(0.25, 0.25);
        let second = assess(0.25, 0.25);

        assert_eq!(first, second);

        // Distance tracks position in descriptor space and nothing else: a
        // point sitting on a training observation scores zero, and one in the
        // middle of a cell scores the half-diagonal, whatever a model would
        // predict at either.
        let on_a_point = assess(0.0, 0.0);
        let mid_cell = assess(0.5, 0.5);
        assert_eq!(on_a_point.nearest_training_distance, Some(0.0));
        assert!(
            (mid_cell.nearest_training_distance.unwrap() - 0.5_f64.hypot(0.5)).abs() < 1.0e-12,
            "distance must track geometry, not desirability: got {:?}",
            mid_cell.nearest_training_distance
        );
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
            standardized_training_points: Vec::new(),
            neighbor_calibration: None,
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
