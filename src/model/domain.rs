//! Applicability domain and predictive uncertainty for a fitted model.
//!
//! A regression can return a number for any input. This module reports
//! descriptor-space diagnostics and conditional model uncertainty; neither
//! establishes chemical reliability:
//!
//! * **Leverage** `h = x'(X'X)⁻¹x` measures how far a candidate sits from the
//!   centre of the training design, in the metric the fit itself defines. The
//!   conventional warning threshold is `h* = 3p/n`; beyond it a prediction is an
//!   extrapolation supported by little or no nearby training data. This is the
//!   same criterion the project's own Study 003 pre-registration applies.
//! * **Prediction interval** `ŷ ± t(0.975, n−p)·s·√(1+h)` widens automatically
//!   with leverage under the fixed linear model and IID homoscedastic normal
//!   error assumptions. It does not guarantee empirical or chemical coverage.

use crate::model::{FeatureDomain, MODEL_FEATURE_COUNT};
use serde::{Deserialize, Serialize};
use std::cell::OnceCell;

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
    /// Identifier of each training observation, positionally aligned with
    /// `standardized_training_points`, so a screened ligand's nearest training
    /// neighbour can be named rather than only measured.
    ///
    /// Identifiers only. No experimental response is recorded here: a model
    /// document must never become a channel for a blinded target value.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub training_labels: Vec<String>,
    /// Nearest-neighbour spacing of the training set, and the boundary derived
    /// from it. Absent when there are too few points to measure a spacing.
    #[serde(default)]
    pub neighbor_calibration: Option<NeighborCalibration>,
}

/// Reuses one fit's Student-t multiplier across prediction interval queries.
///
/// The geometry remains immutably borrowed for this evaluator's lifetime.
/// Construction does no numerical work: every query first performs the same
/// leverage calculation as [`TrainingGeometry::prediction_interval`], then
/// initializes or reuses the multiplier, including an unavailable value.
/// Candidate-specific interval availability is never cached.
pub struct PredictionIntervalEvaluator<'a> {
    geometry: &'a TrainingGeometry,
    multiplier: OnceCell<Option<f64>>,
}

impl<'a> PredictionIntervalEvaluator<'a> {
    /// Borrow a fit without validating it or calculating its multiplier.
    #[must_use]
    pub fn new(geometry: &'a TrainingGeometry) -> Self {
        Self {
            geometry,
            multiplier: OnceCell::new(),
        }
    }

    /// The same 95% prediction interval as the uncached geometry method.
    #[must_use]
    pub fn prediction_interval(
        &self,
        prediction: f64,
        expanded: &[f32; MODEL_FEATURE_COUNT],
    ) -> Option<(f64, f64)> {
        crate::profile_scope!(
            "uncertainty",
            "model::TrainingGeometry::prediction_interval"
        );
        let leverage = self.geometry.leverage(expanded)?;
        let multiplier = (*self.multiplier.get_or_init(|| {
            crate::profile_scope!(
                "uncertainty_cache",
                "model::PredictionIntervalEvaluator::initialize_multiplier"
            );
            self.geometry.t_multiplier()
        }))?;
        let half_width =
            multiplier * self.geometry.residual_standard_error * (1.0 + leverage.max(0.0)).sqrt();
        (half_width.is_finite()
            && (prediction - half_width).is_finite()
            && (prediction + half_width).is_finite())
        .then_some((prediction - half_width, prediction + half_width))
    }
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

/// Which statistic of the training set's own nearest-neighbour spacing a
/// candidate's distance is compared against.
///
/// Every variant is a statistic *of that distribution*, computed from the
/// training data at fit time and serialized into the model, so none of them
/// introduces a severity grade or a cutoff chosen to make results look a
/// particular way. The distribution is the n nearest-neighbour distances
/// `d_i = min_{j != i} ||z_i - z_j||` over the standardized training points.
///
/// [`Self::MaxNeighbor`] is the default because it is the only variant with no
/// distributional assumption at all. The two `mean + kσ` variants are stricter
/// and assume the neighbour distances are roughly symmetric; under approximate
/// normality they sit near the 84th and 97.7th percentiles of that
/// distribution. They exist because the maximum is set by the single sparsest
/// training point, so one isolated observation widens the boundary for every
/// candidate.
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum DomainRule {
    /// `threshold = maximum`. Permissive: the sparsest training point sets it.
    #[default]
    MaxNeighbor,
    /// `threshold = mean + σ`.
    MeanPlusSd,
    /// `threshold = mean + 2σ`.
    MeanPlusTwoSd,
}

impl DomainRule {
    /// The boundary this rule derives from a training set's spacing.
    ///
    /// Never negative: a degenerate training set whose neighbour distances are
    /// all identical has `σ = 0`, and every rule then collapses onto the mean.
    #[must_use]
    pub fn threshold(self, calibration: &NeighborCalibration) -> f64 {
        let derived = match self {
            Self::MaxNeighbor => calibration.threshold,
            Self::MeanPlusSd => calibration.mean + calibration.standard_deviation,
            Self::MeanPlusTwoSd => calibration.mean + 2.0 * calibration.standard_deviation,
        };
        derived.max(0.0)
    }

    /// A short label for reports.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::MaxNeighbor => "max_neighbor",
            Self::MeanPlusSd => "mean_plus_sd",
            Self::MeanPlusTwoSd => "mean_plus_2sd",
        }
    }

    /// The exact derivation, so a report never leaves a reader guessing.
    #[must_use]
    pub fn rule(self) -> &'static str {
        match self {
            Self::MaxNeighbor => {
                "threshold = max over training points of the distance to the nearest other training point, in standardized descriptor space"
            }
            Self::MeanPlusSd => {
                "threshold = mean + 1 standard deviation of the training nearest-neighbour distances, in standardized descriptor space"
            }
            Self::MeanPlusTwoSd => {
                "threshold = mean + 2 standard deviations of the training nearest-neighbour distances, in standardized descriptor space"
            }
        }
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
    /// The entire assessment is unavailable when its input metadata is invalid.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub unavailable: Option<String>,
    /// Euclidean distance in standardized descriptor space to the closest
    /// training observation.
    pub nearest_training_distance: Option<f64>,
    /// Identifier of that nearest training observation, when recorded.
    pub nearest_training_label: Option<String>,
    /// The calibrated boundary that distance is compared against.
    pub nearest_training_threshold: Option<f64>,
    /// Which statistic of the training spacing produced that boundary.
    pub nearest_training_rule: DomainRule,
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
    /// Descriptors outside their training range; empty also when assessment is unavailable.
    pub outside_range: Vec<DescriptorExceedance>,
    /// Largest `normalized_exceedance`, or zero when every descriptor is inside.
    pub maximum_extrapolation: f64,
}

fn euclidean(left: &[f64], right: &[f64]) -> Option<f64> {
    if left.is_empty() || left.len() != right.len() {
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
    /// Validate persisted numerical geometry before it can be used for inference.
    /// This checks structure and finiteness, not experimental applicability.
    pub fn validate(&self) -> Result<(), String> {
        let columns = self.feature_indices.len();
        if self.observations == 0
            || self.parameters != columns + 1
            || self.means.len() != columns
            || self.scales.len() != columns
            || self
                .feature_indices
                .iter()
                .any(|&i| i == 0 || i >= MODEL_FEATURE_COUNT)
            || self
                .feature_indices
                .iter()
                .enumerate()
                .any(|(i, value)| self.feature_indices[..i].contains(value))
        {
            return Err("training geometry has inconsistent dimensions or feature indices".into());
        }
        if self.means.iter().any(|value| !value.is_finite())
            || self
                .scales
                .iter()
                .any(|value| !value.is_finite() || *value <= 0.0)
            || !self.residual_standard_error.is_finite()
            || self.residual_standard_error < 0.0
            || !self.warning_leverage.is_finite()
            || self.warning_leverage < 0.0
        {
            return Err("training geometry contains invalid scales or nonfinite statistics".into());
        }
        if self.xtx_inverse.len() != self.parameters
            || self.xtx_inverse.iter().any(|row| {
                row.len() != self.parameters || row.iter().any(|value| !value.is_finite())
            })
        {
            return Err("stored design inverse is not a finite square parameter matrix".into());
        }
        for i in 0..self.parameters {
            if self.xtx_inverse[i][i] <= 0.0 {
                return Err("stored design inverse has a nonpositive diagonal".into());
            }
            for j in 0..i {
                let (a, b) = (self.xtx_inverse[i][j], self.xtx_inverse[j][i]);
                if (a - b).abs() > 1.0e-10 * a.abs().max(b.abs()).max(1.0) {
                    return Err("stored design inverse is not symmetric".into());
                }
            }
        }
        if !self.standardized_training_points.is_empty()
            && (self.standardized_training_points.len() != self.observations
                || self
                    .standardized_training_points
                    .iter()
                    .any(|row| row.len() != columns || row.iter().any(|value| !value.is_finite())))
        {
            return Err(
                "training covariance is not estimable from inconsistent stored points".into(),
            );
        }
        if !self.training_labels.is_empty() && self.training_labels.len() != self.observations {
            return Err("training geometry labels do not match observations".into());
        }
        if let Some(calibration) = &self.neighbor_calibration
            && [
                calibration.mean,
                calibration.standard_deviation,
                calibration.median,
                calibration.maximum,
                calibration.threshold,
            ]
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err("nearest-neighbor calibration must be finite and nonnegative".into());
        }
        Ok(())
    }

    /// Standardized design vector for one expanded feature row.
    #[must_use]
    pub fn design_vector(&self, expanded: &[f32; MODEL_FEATURE_COUNT]) -> Vec<f64> {
        if self.validate().is_err() {
            return Vec::new();
        }
        let mut design = Vec::with_capacity(self.parameters);
        design.push(1.0);
        for ((index, mean), scale) in self
            .feature_indices
            .iter()
            .zip(&self.means)
            .zip(&self.scales)
        {
            let value = f64::from(expanded[*index]);
            design.push((value - mean) / scale);
        }
        design
    }

    /// Leverage `h = x'(X'X)⁻¹x` of a candidate against the training design.
    #[must_use]
    pub fn leverage(&self, expanded: &[f32; MODEL_FEATURE_COUNT]) -> Option<f64> {
        let design = self.design_vector(expanded);
        if design.len() != self.parameters
            || self.xtx_inverse.len() != design.len()
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
        (leverage.is_finite() && leverage >= 0.0).then_some(leverage)
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
            .filter(|value| value.is_finite())
    }

    /// 95 % prediction interval for a new observation: `ŷ ± t·s·√(1+h)`.
    ///
    /// This is the interval for an individual future measurement, so it carries
    /// both the parameter uncertainty (through `h`) and the residual scatter
    /// (through `s`), conditional on the fixed linear model and IID homoscedastic
    /// normal errors. This does not establish chemical predictive coverage.
    #[must_use]
    pub fn prediction_interval(
        &self,
        prediction: f64,
        expanded: &[f32; MODEL_FEATURE_COUNT],
    ) -> Option<(f64, f64)> {
        crate::profile_scope!(
            "uncertainty",
            "model::TrainingGeometry::prediction_interval"
        );
        let leverage = self.leverage(expanded)?;
        let multiplier = self.t_multiplier()?;
        let half_width =
            multiplier * self.residual_standard_error * (1.0 + leverage.max(0.0)).sqrt();
        (half_width.is_finite()
            && (prediction - half_width).is_finite()
            && (prediction + half_width).is_finite())
        .then_some((prediction - half_width, prediction + half_width))
    }

    /// Standardized descriptor coordinates of a candidate, without the
    /// intercept column.
    #[must_use]
    pub fn standardized_point(&self, expanded: &[f32; MODEL_FEATURE_COUNT]) -> Vec<f64> {
        self.design_vector(expanded).into_iter().skip(1).collect()
    }

    /// Nearest training observation as `(index, distance)`.
    #[must_use]
    pub fn nearest_training_neighbour(
        &self,
        expanded: &[f32; MODEL_FEATURE_COUNT],
    ) -> Option<(usize, f64)> {
        crate::profile_scope!(
            "applicability",
            "model::TrainingGeometry::nearest_training_neighbour"
        );
        if self.standardized_training_points.is_empty() {
            return None;
        }
        let point = self.standardized_point(expanded);
        self.standardized_training_points
            .iter()
            .enumerate()
            .filter_map(|(index, training)| Some((index, euclidean(&point, training)?)))
            .min_by(|left, right| left.1.total_cmp(&right.1))
    }

    /// Identifier of the nearest training observation, when the model records
    /// training labels.
    #[must_use]
    pub fn nearest_training_label(&self, expanded: &[f32; MODEL_FEATURE_COUNT]) -> Option<String> {
        let (index, _) = self.nearest_training_neighbour(expanded)?;
        self.training_labels.get(index).cloned()
    }

    /// Distance in standardized descriptor space to the closest training point.
    #[must_use]
    pub fn nearest_training_distance(&self, expanded: &[f32; MODEL_FEATURE_COUNT]) -> Option<f64> {
        crate::profile_scope!(
            "applicability",
            "model::TrainingGeometry::nearest_training_distance"
        );
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

    /// Ordinary Mahalanobis distance from the stored points' sample covariance.
    ///
    /// Reconstruct the unregularized covariance instead of deriving it from the
    /// fitted design inverse, which can contain a ridge floor. A singular or
    /// numerically non-invertible covariance has no ordinary distance here.
    /// Legacy models without the training points report unavailability.
    pub fn mahalanobis_distance(
        &self,
        expanded: &[f32; MODEL_FEATURE_COUNT],
    ) -> Result<f64, String> {
        self.validate()?;
        let dimensions = self.feature_indices.len();
        let points = &self.standardized_training_points;
        if dimensions == 0 || self.observations <= dimensions || self.observations < 2 {
            return Err("training covariance is not estimable from these dimensions".into());
        }
        if points.is_empty() {
            return Err("ordinary covariance unavailable: model carries no training points".into());
        }
        let mut means = vec![0.0; dimensions];
        for row in points {
            for (mean, value) in means.iter_mut().zip(row) {
                *mean += value / points.len() as f64;
            }
        }
        let mut covariance = vec![vec![0.0; dimensions]; dimensions];
        for row in points {
            for i in 0..dimensions {
                for j in 0..dimensions {
                    covariance[i][j] +=
                        (row[i] - means[i]) * (row[j] - means[j]) / (points.len() - 1) as f64;
                }
            }
        }
        let inverse = invert_matrix(&covariance).map_err(|_| {
            "ordinary covariance is singular or numerically non-invertible".to_owned()
        })?;
        let point = self.standardized_point(expanded);
        let delta = point
            .iter()
            .zip(&means)
            .map(|(x, mean)| x - mean)
            .collect::<Vec<_>>();
        let mut squared = 0.0;
        for i in 0..dimensions {
            for j in 0..dimensions {
                squared += delta[i] * inverse[i][j] * delta[j];
            }
        }
        if !squared.is_finite() || squared < 0.0 {
            return Err("ordinary covariance distance is not finite and nonnegative".into());
        }
        Ok(squared.sqrt())
    }

    /// 95 % confidence interval for the fitted mean response: `ŷ ± t·s·√h`.
    #[must_use]
    pub fn confidence_interval(
        &self,
        prediction: f64,
        expanded: &[f32; MODEL_FEATURE_COUNT],
    ) -> Option<(f64, f64)> {
        crate::profile_scope!(
            "uncertainty",
            "model::TrainingGeometry::confidence_interval"
        );
        let leverage = self.leverage(expanded)?;
        let multiplier = self.t_multiplier()?;
        let half_width = multiplier * self.residual_standard_error * leverage.max(0.0).sqrt();
        (half_width.is_finite()
            && (prediction - half_width).is_finite()
            && (prediction + half_width).is_finite())
        .then_some((prediction - half_width, prediction + half_width))
    }
}

/// Invert a square matrix by Gauss-Jordan elimination with partial pivoting.
pub(crate) fn invert_matrix(matrix: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
    let size = matrix.len();
    if size == 0
        || matrix
            .iter()
            .any(|row| row.len() != size || row.iter().any(|value| !value.is_finite()))
    {
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
    if inverse.iter().flatten().any(|value| !value.is_finite()) {
        return Err("matrix inverse is not finite".into());
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
/// A finite upper bracket is established before bisection. Invalid arguments,
/// a nonfinite tail evaluation, or a quantile beyond the squared-argument range
/// of the incomplete-beta calculation return `NaN` (numerically unavailable),
/// never a silently truncated bracket endpoint.
#[must_use]
pub fn student_t_two_sided_quantile(alpha: f64, df: f64) -> f64 {
    crate::profile_scope!("uncertainty", "model::student_t_two_sided_quantile");
    if !(0.0..1.0).contains(&alpha) || alpha <= 0.0 || !df.is_finite() || df <= 0.0 {
        return f64::NAN;
    }
    let (mut low, mut high) = (0.0_f64, 1.0_f64);
    // Near alpha=1 the desired central mass is small. Evaluate its beta
    // representation directly; df/(df+t²) rounds to 1 near t=0 and would
    // destroy this information before the tail probability is calculated.
    // The complementary branch preserves the usual small-alpha tail path.
    let below_quantile = |t: f64| -> Option<bool> {
        if alpha > 0.5 {
            let squared = t * t;
            let central = regularized_incomplete_beta(0.5, 0.5 * df, squared / (df + squared));
            central.is_finite().then_some(central < 1.0 - alpha)
        } else {
            let tail = student_t_two_sided_tail(t, df);
            tail.is_finite().then_some(tail > alpha)
        }
    };
    loop {
        // df + high² must remain finite: overflow must not masquerade as a
        // zero tail and therefore a valid bracket for an extreme quantile.
        if !(df + high * high).is_finite() {
            return f64::NAN;
        }
        let Some(below) = below_quantile(high) else {
            return f64::NAN;
        };
        if !below {
            break;
        }
        high *= 2.0;
    }
    for _ in 0..200 {
        let middle = 0.5 * (low + high);
        let Some(below) = below_quantile(middle) else {
            return f64::NAN;
        };
        if below {
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
    rule: DomainRule,
) -> ApplicabilityAssessment {
    crate::profile_scope!("applicability", "model::assess_applicability");
    let unavailable = |reason: String| ApplicabilityAssessment {
        verdict: DomainVerdict::Unknown,
        unavailable: Some(reason.clone()),
        nearest_training_distance: None,
        nearest_training_label: None,
        nearest_training_threshold: None,
        nearest_training_rule: rule,
        nearest_training_ratio: None,
        leverage: None,
        leverage_ratio: None,
        mahalanobis_distance: None,
        mahalanobis_unavailable: Some(reason),
        outside_range: Vec::new(),
        maximum_extrapolation: 0.0,
    };
    if selected.len() != ranges.len()
        || selected.iter().enumerate().any(|(i, &column)| {
            column == 0 || column >= MODEL_FEATURE_COUNT || selected[..i].contains(&column)
        })
    {
        return unavailable("invalid applicability range dimensions or feature indices".into());
    }
    if ranges.iter().any(|range| {
        !range.minimum.is_finite() || !range.maximum.is_finite() || range.minimum > range.maximum
    }) {
        return unavailable("invalid applicability descriptor ranges".into());
    }
    if let Some(geometry) = geometry {
        if let Err(reason) = geometry.validate() {
            return unavailable(format!("invalid training geometry: {reason}"));
        }
        if geometry
            .feature_indices
            .iter()
            .any(|&column| !expanded[column].is_finite())
        {
            return unavailable("applicability requires finite selected descriptors".into());
        }
    }
    if selected.iter().any(|&column| !expanded[column].is_finite()) {
        return unavailable("applicability requires finite selected descriptors".into());
    }
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
            unavailable: None,
            verdict: if outside_range.is_empty() {
                DomainVerdict::Unknown
            } else {
                // A range violation is decidable without any geometry.
                DomainVerdict::Extrapolation
            },
            nearest_training_distance: None,
            nearest_training_label: None,
            nearest_training_threshold: None,
            nearest_training_rule: rule,
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
        .map(|calibration| rule.threshold(calibration));
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
        unavailable: None,
        nearest_training_distance: distance,
        nearest_training_label: geometry.nearest_training_label(expanded),
        nearest_training_threshold: threshold,
        nearest_training_rule: rule,
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
#[path = "domain_interval_cache_tests.rs"]
mod interval_cache_tests;

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

    #[test]
    fn brackets_extreme_cauchy_tail_without_silent_clipping() {
        // df=1 is exactly the Cauchy distribution: t=cot(pi*alpha/2).
        for alpha in [0.05, 1.0e-8, 1.0e-20, 1.0e-100] {
            let expected = 1.0 / (std::f64::consts::PI * alpha / 2.0).tan();
            let actual = student_t_two_sided_quantile(alpha, 1.0);
            assert!(
                (actual / expected - 1.0).abs() < 2.0e-12,
                "{alpha}: {actual} vs {expected}"
            );
        }
        assert!(student_t_two_sided_quantile(1.0e-300, 1.0).is_nan());
        for df in [f64::NAN, f64::INFINITY, 0.0, -1.0] {
            assert!(student_t_two_sided_quantile(0.05, df).is_nan());
        }
    }

    #[test]
    fn near_unit_alpha_preserves_small_central_mass() {
        for alpha in [0.75, 0.999_999_99, f64::from_bits(1.0_f64.to_bits() - 1)] {
            // Cauchy complement form avoids subtracting nearly equal angles.
            let expected = (std::f64::consts::PI * (1.0 - alpha) / 2.0).tan();
            let actual = student_t_two_sided_quantile(alpha, 1.0);
            assert!(
                (actual / expected - 1.0).abs() < 2.0e-12,
                "alpha {alpha}: {actual} versus analytic Cauchy {expected}"
            );
        }
        for alpha in [0.999_999_99, f64::from_bits(1.0_f64.to_bits() - 1)] {
            let mass = 1.0 - alpha;
            // Analytic t densities at zero. The central inverse is
            // mass/(2*density) + O(mass³); for these dfs and masses the
            // relative cubic contribution is below 2e-16.
            for (df, density) in [
                (2.0, 1.0 / (2.0 * 2.0_f64.sqrt())),
                (3.0, 2.0 / (std::f64::consts::PI * 3.0_f64.sqrt())),
                (4.0, 3.0 / 8.0),
                (10.0, 315.0 / (256.0 * 10.0_f64.sqrt())),
            ] {
                let expected = mass / (2.0 * density);
                let actual = student_t_two_sided_quantile(alpha, df);
                assert!(
                    (actual / expected - 1.0).abs() < 2.0e-12,
                    "df {df}, alpha {alpha}: {actual} versus central limit {expected}"
                );
            }
        }
    }

    #[test]
    fn public_assessment_reports_malformed_metadata_without_panicking() {
        for case in 0..4 {
            let mut geometry = grid_geometry();
            match case {
                0 => geometry.feature_indices[0] = MODEL_FEATURE_COUNT,
                1 => geometry.xtx_inverse[0].clear(),
                2 => geometry.scales[0] = 0.0,
                _ => geometry.standardized_training_points[0].clear(),
            }
            let assessment = assess_applicability(
                Some(&geometry),
                &grid_ranges(),
                &geometry.feature_indices,
                &at(0.0, 0.0),
                DomainRule::default(),
            );
            assert_eq!(assessment.verdict, DomainVerdict::Unknown);
            assert!(assessment.unavailable.is_some());
            assert!(assessment.leverage.is_none());
            assert!(assessment.nearest_training_distance.is_none());
            assert!(assessment.mahalanobis_distance.is_none());
        }
        let assessment = assess_applicability(
            None,
            &grid_ranges(),
            &[MODEL_FEATURE_COUNT, 3],
            &at(0.0, 0.0),
            DomainRule::default(),
        );
        assert_eq!(assessment.verdict, DomainVerdict::Unknown);
        assert!(assessment.unavailable.is_some());
        let assessment = assess_applicability(
            Some(&grid_geometry()),
            &grid_ranges(),
            &[2, 3],
            &at(f32::NAN, 0.0),
            DomainRule::default(),
        );
        assert_eq!(assessment.verdict, DomainVerdict::Unknown);
        assert!(assessment.unavailable.is_some());
    }

    #[test]
    fn singular_covariance_is_not_a_regularized_mahalanobis_distance() {
        let mut geometry = grid_geometry();
        geometry.standardized_training_points =
            (-4..=4).map(|x| vec![f64::from(x), f64::from(x)]).collect();
        // A finite, positive regularized inverse must not hide rank-one points.
        geometry.xtx_inverse = vec![
            vec![1.0 / 9.0, 0.0, 0.0],
            vec![0.0, 5.0e9, -5.0e9],
            vec![0.0, -5.0e9, 5.0e9],
        ];
        let assessment = assess_applicability(
            Some(&geometry),
            &[],
            &[],
            &at(1.0, -1.0),
            DomainRule::default(),
        );
        assert!(assessment.mahalanobis_distance.is_none());
        assert!(
            assessment
                .mahalanobis_unavailable
                .unwrap()
                .contains("singular")
        );
    }

    #[test]
    fn malformed_geometry_declines_inference_without_panicking() {
        for kind in 0..4 {
            let mut geometry = grid_geometry();
            match kind {
                0 => geometry.feature_indices[0] = MODEL_FEATURE_COUNT,
                1 => geometry.xtx_inverse[0].clear(),
                2 => geometry.scales[0] = f64::NAN,
                _ => geometry.standardized_training_points[0].clear(),
            }
            assert!(geometry.validate().is_err());
            assert!(geometry.leverage(&at(0.0, 0.0)).is_none());
            assert!(geometry.mahalanobis_distance(&at(0.0, 0.0)).is_err());
            assert!(geometry.nearest_training_distance(&at(0.0, 0.0)).is_none());
        }
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
            training_labels: (0..points.len()).map(|i| format!("T{i}")).collect(),
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
        assess_applicability(
            Some(&grid_geometry()),
            &grid_ranges(),
            &[2, 3],
            &at(b1, b5),
            DomainRule::default(),
        )
    }

    /// A one-dimensional training set with one isolated observation, so the
    /// maximum neighbour distance is set by a single outlying point and the
    /// permissive default is visibly looser than the mean-based rules.
    ///
    /// Points at 0..=8 and one isolated at 20. Nine neighbour distances of 1
    /// and one of 12, so mean = 2.1, variance = (9·1.1² + 9.9²)/10 = 10.89 and
    /// σ = 3.3 exactly. The three rules therefore land far apart: 12, 5.4, 8.7.
    fn isolated_point_geometry() -> TrainingGeometry {
        let points: Vec<Vec<f64>> = [0.0_f64, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 20.0]
            .into_iter()
            .map(|x| vec![x])
            .collect();
        let observations = points.len();
        TrainingGeometry {
            feature_indices: vec![2],
            means: vec![0.0],
            scales: vec![1.0],
            xtx_inverse: vec![vec![1.0 / observations as f64, 0.0], vec![0.0, 1.0 / 50.0]],
            observations,
            parameters: 2,
            residual_standard_error: 0.1,
            warning_leverage: 3.0 * 2.0 / observations as f64,
            neighbor_calibration: NeighborCalibration::from_points(&points),
            training_labels: (0..points.len()).map(|i| format!("T{i}")).collect(),
            standardized_training_points: points,
        }
    }

    #[test]
    fn each_rule_derives_its_threshold_from_the_calibration() {
        let calibration = isolated_point_geometry().neighbor_calibration.unwrap();
        assert!((calibration.mean - 2.1).abs() < 1.0e-12);
        assert!((calibration.standard_deviation - 3.3).abs() < 1.0e-12);
        assert!((calibration.maximum - 12.0).abs() < 1.0e-12);

        // Each is exactly the documented statistic, nothing tuned.
        assert!((DomainRule::MaxNeighbor.threshold(&calibration) - 12.0).abs() < 1.0e-12);
        assert!((DomainRule::MeanPlusSd.threshold(&calibration) - 5.4).abs() < 1.0e-12);
        assert!((DomainRule::MeanPlusTwoSd.threshold(&calibration) - 8.7).abs() < 1.0e-12);
    }

    #[test]
    fn a_uniform_training_set_collapses_every_rule_onto_one_boundary() {
        // The grid's neighbour distances are all exactly 1.0, so σ = 0 and the
        // mean-based rules have nothing to widen by.
        let calibration = grid_geometry().neighbor_calibration.unwrap();
        assert!(calibration.standard_deviation.abs() < 1.0e-12);
        for rule in [
            DomainRule::MaxNeighbor,
            DomainRule::MeanPlusSd,
            DomainRule::MeanPlusTwoSd,
        ] {
            assert!(
                (rule.threshold(&calibration) - 1.0).abs() < 1.0e-12,
                "{} moved a boundary a uniform training set does not support",
                rule.label()
            );
        }
    }

    #[test]
    fn a_stricter_rule_reclassifies_a_point_the_default_accepts() {
        let geometry = isolated_point_geometry();
        let ranges = vec![FeatureDomain {
            feature: "B1_boltz".to_owned(),
            minimum: 0.0,
            maximum: 20.0,
        }];
        // 14.0 sits in the 8 -> 20 gap, 6.0 from the nearest training point:
        // inside the permissive maximum of 12.0, outside mean + σ = 5.4.
        let point = at(14.0, 0.0);

        let permissive = assess_applicability(
            Some(&geometry),
            &ranges,
            &[2],
            &point,
            DomainRule::MaxNeighbor,
        );
        assert_eq!(permissive.verdict, DomainVerdict::Interpolation);
        assert_eq!(permissive.nearest_training_rule, DomainRule::MaxNeighbor);

        let strict = assess_applicability(
            Some(&geometry),
            &ranges,
            &[2],
            &point,
            DomainRule::MeanPlusSd,
        );
        assert_eq!(strict.verdict, DomainVerdict::SparseInterpolation);
        assert_eq!(strict.nearest_training_rule, DomainRule::MeanPlusSd);

        // The measurement is identical either way; only the boundary moved, so
        // a rule change can never be mistaken for a different distance.
        assert_eq!(
            permissive.nearest_training_distance,
            strict.nearest_training_distance
        );
        assert!(strict.nearest_training_threshold < permissive.nearest_training_threshold);

        // mean + 2σ = 8.7 still accepts it: the rules are ordered by the
        // statistic they name, not by a severity grade layered on top.
        let two_sd = assess_applicability(
            Some(&geometry),
            &ranges,
            &[2],
            &point,
            DomainRule::MeanPlusTwoSd,
        );
        assert_eq!(two_sd.verdict, DomainVerdict::Interpolation);
    }

    #[test]
    fn a_rule_never_yields_a_negative_boundary() {
        // A degenerate calibration cannot push mean - nothing below zero, but
        // guard the invariant explicitly: a negative boundary would make every
        // candidate sparse, including the training points themselves.
        let calibration = NeighborCalibration {
            mean: 0.0,
            standard_deviation: 0.0,
            median: 0.0,
            maximum: 0.0,
            threshold: -1.0,
            rule: "hand-built degenerate calibration".to_owned(),
        };
        for rule in [
            DomainRule::MaxNeighbor,
            DomainRule::MeanPlusSd,
            DomainRule::MeanPlusTwoSd,
        ] {
            assert!(rule.threshold(&calibration) >= 0.0, "{}", rule.label());
        }
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
        geometry.training_labels.clear();
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

        let assessment = assess_applicability(
            Some(&geometry),
            &ranges,
            &[2, 3],
            &at(0.0, 0.0),
            DomainRule::default(),
        );

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

        assert!(error.contains("not symmetric"), "{error}");
    }

    #[test]
    fn a_model_without_geometry_reports_unknown_but_still_checks_ranges() {
        let inside = assess_applicability(
            None,
            &grid_ranges(),
            &[2, 3],
            &at(0.0, 0.0),
            DomainRule::default(),
        );
        assert_eq!(inside.verdict, DomainVerdict::Unknown);
        assert!(inside.nearest_training_distance.is_none());
        assert!(inside.mahalanobis_unavailable.is_some());

        // A range violation needs no geometry to be decidable.
        let outside = assess_applicability(
            None,
            &grid_ranges(),
            &[2, 3],
            &at(5.0, 0.0),
            DomainRule::default(),
        );
        assert_eq!(outside.verdict, DomainVerdict::Extrapolation);
        assert_eq!(outside.outside_range.len(), 1);
    }

    #[test]
    fn a_model_without_calibration_cannot_claim_interpolation() {
        let mut geometry = grid_geometry();
        geometry.neighbor_calibration = None;

        let assessment = assess_applicability(
            Some(&geometry),
            &grid_ranges(),
            &[2, 3],
            &at(0.0, 0.0),
            DomainRule::default(),
        );

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
        let assessment =
            assess_applicability(None, &ranges, &[2], &at(2.5, 0.0), DomainRule::default());

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
            training_labels: Vec::new(),
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
