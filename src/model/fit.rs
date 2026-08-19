use super::domain::{NeighborCalibration, TrainingGeometry, invert_matrix};
use super::{MODEL_FEATURE_COUNT, MODEL_FEATURE_NAMES, expand_features};
use crate::storage::PackedReactionRecord;
use serde::{Deserialize, Serialize};

const EPSILON: f64 = 1.0e-12;

/// Controls deterministic physical-organic model development.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct FitOptions {
    /// Maximum non-intercept terms. The fitter additionally enforces `p < n/3`.
    pub max_terms: usize,
    /// Fixed-feature bootstrap replicates used for coefficient intervals.
    pub bootstrap_samples: usize,
    /// Fixed-feature response permutations used for a null-model test.
    pub permutation_samples: usize,
    /// Deterministic resampling seed.
    pub seed: u64,
}

impl Default for FitOptions {
    fn default() -> Self {
        Self {
            max_terms: 3,
            bootstrap_samples: 1_000,
            permutation_samples: 500,
            seed: 20_260_725,
        }
    }
}

/// Regression accuracy for one partition or validation procedure.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ModelMetrics {
    pub count: usize,
    pub r2: Option<f64>,
    pub mae: f64,
    pub rmse: f64,
}

/// One fitted regularized baseline.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BaselineReport {
    pub model: String,
    pub regularization: f64,
    pub weights: [f32; MODEL_FEATURE_COUNT],
    pub training: ModelMetrics,
    pub nested_loo: ModelMetrics,
}

/// A percentile confidence interval for one raw-scale model coefficient.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CoefficientInterval {
    pub feature: String,
    pub estimate: f64,
    pub lower_95: f64,
    pub upper_95: f64,
}

/// Raw feature interval used for a simple applicability-domain check.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FeatureDomain {
    pub feature: String,
    pub minimum: f64,
    pub maximum: f64,
}

/// Reproducible, interpretable physical-organic model artifact.
///
/// `weights` are always expressed on the raw eight-feature StericX scale and
/// can therefore be passed directly to `RegressXPredictor`.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ScientificFitReport {
    pub schema_version: u32,
    pub model: String,
    pub training_count: usize,
    pub training_group_count: usize,
    pub feature_names: Vec<String>,
    pub selected_feature_indices: Vec<usize>,
    pub selected_features: Vec<String>,
    pub weights: [f32; MODEL_FEATURE_COUNT],
    pub standardized_means: [f64; MODEL_FEATURE_COUNT],
    pub standardized_scales: [f64; MODEL_FEATURE_COUNT],
    pub training: ModelMetrics,
    pub fixed_feature_loo: ModelMetrics,
    pub fixed_feature_group_loo: ModelMetrics,
    pub ridge_baseline: BaselineReport,
    pub lasso_baseline: BaselineReport,
    pub coefficient_intervals: Vec<CoefficientInterval>,
    pub response_permutation_p_value: f64,
    pub correlation_matrix: Vec<Vec<f64>>,
    pub variance_inflation_factors: Vec<Option<f64>>,
    pub applicability_domain: Vec<FeatureDomain>,
    /// Training-set geometry enabling leverage-based domain checks and true
    /// prediction intervals. Optional so model files written before this was
    /// recorded still deserialize; consumers must degrade honestly when absent.
    #[serde(default)]
    pub training_geometry: Option<TrainingGeometry>,
    pub notes: Vec<String>,
}

#[derive(Clone, Debug)]
struct Fit {
    weights: [f32; MODEL_FEATURE_COUNT],
    means: [f64; MODEL_FEATURE_COUNT],
    scales: [f64; MODEL_FEATURE_COUNT],
    predictions: Vec<f64>,
}

/// Fits a compact, mechanistically interpretable model to selected records.
///
/// Forward selection uses BIC, rejects candidate pairs with absolute
/// correlation above 0.95, and caps the number of terms below one third of the
/// training count. All means and scales are learned exclusively from the
/// supplied training rows.
pub fn fit_scientific_model(
    records: &[PackedReactionRecord],
    training_indices: &[usize],
    options: FitOptions,
) -> Result<ScientificFitReport, String> {
    let unique_groups = (0..training_indices.len())
        .map(|index| format!("row_{index}"))
        .collect::<Vec<_>>();
    fit_scientific_model_grouped(records, training_indices, &unique_groups, options)
}

/// Fits a model and additionally validates by holding out whole ligand groups.
pub fn fit_scientific_model_grouped(
    records: &[PackedReactionRecord],
    training_indices: &[usize],
    training_groups: &[String],
    options: FitOptions,
) -> Result<ScientificFitReport, String> {
    if training_indices.len() < 4 {
        return Err("scientific fitting requires at least four training records".into());
    }
    if training_groups.len() != training_indices.len()
        || training_groups.iter().any(|group| group.trim().is_empty())
    {
        return Err("training groups must contain one non-empty value per record".into());
    }
    if options.max_terms == 0 {
        return Err("max_terms must be positive".into());
    }
    let mut features = Vec::with_capacity(training_indices.len());
    let mut targets = Vec::with_capacity(training_indices.len());
    for &index in training_indices {
        let record = records
            .get(index)
            .ok_or_else(|| format!("training index {index} is out of bounds"))?;
        if !record.exp_ddg.is_finite() {
            return Err(format!("training record {index} has a non-finite target"));
        }
        let row = expand_features(record).map(f64::from);
        if row.iter().any(|value| !value.is_finite()) {
            return Err(format!("training record {index} has non-finite features"));
        }
        features.push(row);
        targets.push(f64::from(record.exp_ddg));
    }

    let term_limit = options
        .max_terms
        .min((training_indices.len().saturating_sub(1) / 3).max(1));
    let selected = select_features(&features, &targets, term_limit)?;
    if selected.is_empty() {
        return Err("no non-constant descriptor improved the intercept model".into());
    }

    let ols = fit_linear(&features, &targets, &selected, 0.0)?;
    let training = metrics(&targets, &ols.predictions);
    let loo_predictions = leave_one_out(&features, &targets, &selected, FitKind::Ols(0.0))?;
    let fixed_feature_loo = metrics(&targets, &loo_predictions);
    let group_loo_predictions = leave_groups_out(&features, &targets, training_groups, &selected)?;
    let fixed_feature_group_loo = metrics(&targets, &group_loo_predictions);
    let training_group_count = training_groups
        .iter()
        .collect::<std::collections::BTreeSet<_>>()
        .len();

    let ridge_alpha = tune_regularization(
        &features,
        &targets,
        &selected,
        &[1.0e-4, 1.0e-3, 1.0e-2, 0.1, 1.0, 10.0],
        false,
    )?;
    let ridge = fit_linear(&features, &targets, &selected, ridge_alpha)?;
    let ridge_loo = nested_regularized_predictions(
        &features,
        &targets,
        &selected,
        &[1.0e-4, 1.0e-3, 1.0e-2, 0.1, 1.0, 10.0],
        false,
    )?;

    let lasso_alpha = tune_regularization(
        &features,
        &targets,
        &selected,
        &[1.0e-4, 1.0e-3, 1.0e-2, 0.05, 0.1, 0.5],
        true,
    )?;
    let lasso = fit_lasso(&features, &targets, &selected, lasso_alpha)?;
    let lasso_loo = nested_regularized_predictions(
        &features,
        &targets,
        &selected,
        &[1.0e-4, 1.0e-3, 1.0e-2, 0.05, 0.1, 0.5],
        true,
    )?;

    let coefficient_intervals = bootstrap_intervals(
        &features,
        &targets,
        &selected,
        &ols.weights,
        options.bootstrap_samples,
        options.seed,
    );
    let permutation_p_value = permutation_test(
        &features,
        &targets,
        &selected,
        training.r2.unwrap_or(0.0),
        options.permutation_samples,
        options.seed ^ 0x9e37_79b9_7f4a_7c15,
    );
    let correlation_matrix = correlation_matrix(&features);
    let variance_inflation_factors = vif(&features, &selected);
    let applicability_domain = selected
        .iter()
        .map(|&column| {
            let (minimum, maximum) = feature_bounds(&features, column);
            FeatureDomain {
                feature: MODEL_FEATURE_NAMES[column].into(),
                minimum,
                maximum,
            }
        })
        .collect();

    Ok(ScientificFitReport {
        schema_version: 1,
        model: "mechanistically_constrained_ols".into(),
        training_count: training_indices.len(),
        training_group_count,
        feature_names: MODEL_FEATURE_NAMES
            .iter()
            .map(|name| (*name).into())
            .collect(),
        selected_feature_indices: selected.clone(),
        selected_features: selected
            .iter()
            .map(|&index| MODEL_FEATURE_NAMES[index].into())
            .collect(),
        weights: ols.weights,
        standardized_means: ols.means,
        standardized_scales: ols.scales,
        training,
        fixed_feature_loo,
        fixed_feature_group_loo,
        ridge_baseline: BaselineReport {
            model: "ridge".into(),
            regularization: ridge_alpha,
            weights: ridge.weights,
            training: metrics(&targets, &ridge.predictions),
            nested_loo: metrics(&targets, &ridge_loo),
        },
        lasso_baseline: BaselineReport {
            model: "lasso".into(),
            regularization: lasso_alpha,
            weights: lasso.weights,
            training: metrics(&targets, &lasso.predictions),
            nested_loo: metrics(&targets, &lasso_loo),
        },
        coefficient_intervals,
        response_permutation_p_value: permutation_p_value,
        correlation_matrix,
        variance_inflation_factors,
        applicability_domain,
        training_geometry: training_geometry(&features, &targets, &selected, &ols.predictions),
        notes: vec![
            "Descriptors were standardized using training rows only.".into(),
            "Leverage h = x'(X'X)^-1 x is reported against the warning leverage h* = 3p/n; \
             beyond it a prediction is an extrapolation."
                .into(),
            "Forward selection used BIC and rejected |r| > 0.95 descriptor pairs.".into(),
            "LOO, bootstrap, and permutation diagnostics keep the selected descriptor set fixed."
                .into(),
            "Group LOO holds out every row sharing the same ligand-group label.".into(),
            "A prospective claim still requires an experimentally untouched test set.".into(),
        ],
    })
}

fn select_features(
    features: &[[f64; MODEL_FEATURE_COUNT]],
    targets: &[f64],
    term_limit: usize,
) -> Result<Vec<usize>, String> {
    let mut selected = Vec::new();
    let intercept_mean = targets.iter().sum::<f64>() / targets.len() as f64;
    let intercept_rss = targets
        .iter()
        .map(|target| (target - intercept_mean).powi(2))
        .sum::<f64>();
    let mut current_bic = bic(intercept_rss, targets.len(), 1);

    while selected.len() < term_limit {
        let mut best: Option<(usize, f64)> = None;
        for candidate in 1..MODEL_FEATURE_COUNT {
            if selected.contains(&candidate) || feature_scale(features, candidate) <= EPSILON {
                continue;
            }
            if selected
                .iter()
                .any(|&existing| correlation(features, candidate, existing).abs() > 0.95)
            {
                continue;
            }
            let mut trial = selected.clone();
            trial.push(candidate);
            let fit = fit_linear(features, targets, &trial, 0.0)?;
            let rss = residual_sum_of_squares(targets, &fit.predictions);
            let trial_bic = bic(rss, targets.len(), trial.len() + 1);
            if best.is_none_or(|(_, best_bic)| trial_bic < best_bic) {
                best = Some((candidate, trial_bic));
            }
        }
        let Some((candidate, candidate_bic)) = best else {
            break;
        };
        if candidate_bic >= current_bic - 2.0 {
            break;
        }
        selected.push(candidate);
        current_bic = candidate_bic;
    }
    Ok(selected)
}

fn fit_linear(
    features: &[[f64; MODEL_FEATURE_COUNT]],
    targets: &[f64],
    selected: &[usize],
    ridge: f64,
) -> Result<Fit, String> {
    let (means, scales) = standardization(features, selected);
    let columns = selected.len() + 1;
    let mut normal = vec![vec![0.0; columns]; columns];
    let mut rhs = vec![0.0; columns];

    for (row, &target) in features.iter().zip(targets) {
        let mut design = Vec::with_capacity(columns);
        design.push(1.0);
        design.extend(
            selected
                .iter()
                .map(|&column| (row[column] - means[column]) / scales[column]),
        );
        for left in 0..columns {
            rhs[left] += design[left] * target;
            for right in 0..columns {
                normal[left][right] += design[left] * design[right];
            }
        }
    }
    for (index, row) in normal.iter_mut().enumerate().skip(1) {
        row[index] += ridge.max(1.0e-10);
    }
    let beta = solve_linear_system(normal, rhs)?;
    let weights = raw_weights(&beta, selected, &means, &scales);
    let predictions = features.iter().map(|row| predict(row, &weights)).collect();
    Ok(Fit {
        weights,
        means,
        scales,
        predictions,
    })
}

fn fit_lasso(
    features: &[[f64; MODEL_FEATURE_COUNT]],
    targets: &[f64],
    selected: &[usize],
    alpha: f64,
) -> Result<Fit, String> {
    let (means, scales) = standardization(features, selected);
    let target_mean = targets.iter().sum::<f64>() / targets.len() as f64;
    let mut beta = vec![0.0; selected.len()];
    for _ in 0..2_000 {
        let mut maximum_change = 0.0_f64;
        for (coefficient, &column) in selected.iter().enumerate() {
            let mut numerator = 0.0;
            let mut denominator = 0.0;
            for (row, &target) in features.iter().zip(targets) {
                let x = (row[column] - means[column]) / scales[column];
                let other = selected
                    .iter()
                    .enumerate()
                    .filter(|(index, _)| *index != coefficient)
                    .map(|(index, &other_column)| {
                        beta[index] * (row[other_column] - means[other_column])
                            / scales[other_column]
                    })
                    .sum::<f64>();
                numerator += x * (target - target_mean - other);
                denominator += x * x;
            }
            let updated = soft_threshold(numerator / features.len() as f64, alpha)
                / (denominator / features.len() as f64).max(EPSILON);
            maximum_change = maximum_change.max((updated - beta[coefficient]).abs());
            beta[coefficient] = updated;
        }
        if maximum_change < 1.0e-10 {
            break;
        }
    }
    let mut standardized = Vec::with_capacity(selected.len() + 1);
    standardized.push(target_mean);
    standardized.extend(beta);
    let weights = raw_weights(&standardized, selected, &means, &scales);
    let predictions = features.iter().map(|row| predict(row, &weights)).collect();
    Ok(Fit {
        weights,
        means,
        scales,
        predictions,
    })
}

#[derive(Clone, Copy)]
enum FitKind {
    Ols(f64),
    Lasso(f64),
}

fn leave_one_out(
    features: &[[f64; MODEL_FEATURE_COUNT]],
    targets: &[f64],
    selected: &[usize],
    kind: FitKind,
) -> Result<Vec<f64>, String> {
    let mut predictions = Vec::with_capacity(targets.len());
    for held_out in 0..targets.len() {
        let train_features = features
            .iter()
            .enumerate()
            .filter(|(index, _)| *index != held_out)
            .map(|(_, row)| *row)
            .collect::<Vec<_>>();
        let train_targets = targets
            .iter()
            .enumerate()
            .filter(|(index, _)| *index != held_out)
            .map(|(_, target)| *target)
            .collect::<Vec<_>>();
        let fit = match kind {
            FitKind::Ols(alpha) => fit_linear(&train_features, &train_targets, selected, alpha)?,
            FitKind::Lasso(alpha) => fit_lasso(&train_features, &train_targets, selected, alpha)?,
        };
        predictions.push(predict(&features[held_out], &fit.weights));
    }
    Ok(predictions)
}

fn leave_groups_out(
    features: &[[f64; MODEL_FEATURE_COUNT]],
    targets: &[f64],
    groups: &[String],
    selected: &[usize],
) -> Result<Vec<f64>, String> {
    let unique_groups = groups
        .iter()
        .map(String::as_str)
        .collect::<std::collections::BTreeSet<_>>();
    if unique_groups.len() < 2 {
        return Err("group validation requires at least two ligand groups".into());
    }
    let mut predictions = vec![f64::NAN; targets.len()];
    for held_out_group in unique_groups {
        let training_indices = groups
            .iter()
            .enumerate()
            .filter(|(_, group)| group.as_str() != held_out_group)
            .map(|(index, _)| index)
            .collect::<Vec<_>>();
        if training_indices.len() <= selected.len() {
            return Err(format!(
                "holding out group `{held_out_group}` leaves too few training rows"
            ));
        }
        let train_features = training_indices
            .iter()
            .map(|&index| features[index])
            .collect::<Vec<_>>();
        let train_targets = training_indices
            .iter()
            .map(|&index| targets[index])
            .collect::<Vec<_>>();
        let fit = fit_linear(&train_features, &train_targets, selected, 0.0)?;
        for (index, group) in groups.iter().enumerate() {
            if group == held_out_group {
                predictions[index] = predict(&features[index], &fit.weights);
            }
        }
    }
    if predictions.iter().any(|prediction| !prediction.is_finite()) {
        return Err("group validation failed to predict every training row".into());
    }
    Ok(predictions)
}

fn nested_regularized_predictions(
    features: &[[f64; MODEL_FEATURE_COUNT]],
    targets: &[f64],
    selected: &[usize],
    grid: &[f64],
    lasso: bool,
) -> Result<Vec<f64>, String> {
    let mut predictions = Vec::with_capacity(targets.len());
    for held_out in 0..targets.len() {
        let train_features = features
            .iter()
            .enumerate()
            .filter(|(index, _)| *index != held_out)
            .map(|(_, row)| *row)
            .collect::<Vec<_>>();
        let train_targets = targets
            .iter()
            .enumerate()
            .filter(|(index, _)| *index != held_out)
            .map(|(_, target)| *target)
            .collect::<Vec<_>>();
        let alpha = tune_regularization(&train_features, &train_targets, selected, grid, lasso)?;
        let fit = if lasso {
            fit_lasso(&train_features, &train_targets, selected, alpha)?
        } else {
            fit_linear(&train_features, &train_targets, selected, alpha)?
        };
        predictions.push(predict(&features[held_out], &fit.weights));
    }
    Ok(predictions)
}

fn tune_regularization(
    features: &[[f64; MODEL_FEATURE_COUNT]],
    targets: &[f64],
    selected: &[usize],
    grid: &[f64],
    lasso: bool,
) -> Result<f64, String> {
    let mut best = None;
    for &alpha in grid {
        let kind = if lasso {
            FitKind::Lasso(alpha)
        } else {
            FitKind::Ols(alpha)
        };
        let predictions = leave_one_out(features, targets, selected, kind)?;
        let mse = residual_sum_of_squares(targets, &predictions) / targets.len() as f64;
        if best.is_none_or(|(_, best_mse)| mse < best_mse) {
            best = Some((alpha, mse));
        }
    }
    best.map(|(alpha, _)| alpha)
        .ok_or_else(|| "regularization grid is empty".into())
}

fn standardization(
    features: &[[f64; MODEL_FEATURE_COUNT]],
    selected: &[usize],
) -> ([f64; MODEL_FEATURE_COUNT], [f64; MODEL_FEATURE_COUNT]) {
    let mut means = [0.0; MODEL_FEATURE_COUNT];
    let mut scales = [1.0; MODEL_FEATURE_COUNT];
    for &column in selected {
        means[column] = features.iter().map(|row| row[column]).sum::<f64>() / features.len() as f64;
        scales[column] = feature_scale(features, column).max(EPSILON);
    }
    (means, scales)
}

fn feature_scale(features: &[[f64; MODEL_FEATURE_COUNT]], column: usize) -> f64 {
    let mean = features.iter().map(|row| row[column]).sum::<f64>() / features.len() as f64;
    (features
        .iter()
        .map(|row| (row[column] - mean).powi(2))
        .sum::<f64>()
        / features.len() as f64)
        .sqrt()
}

fn raw_weights(
    beta: &[f64],
    selected: &[usize],
    means: &[f64; MODEL_FEATURE_COUNT],
    scales: &[f64; MODEL_FEATURE_COUNT],
) -> [f32; MODEL_FEATURE_COUNT] {
    let mut weights = [0.0_f32; MODEL_FEATURE_COUNT];
    let mut intercept = beta[0];
    for (position, &column) in selected.iter().enumerate() {
        let raw = beta[position + 1] / scales[column];
        weights[column] = raw as f32;
        intercept -= raw * means[column];
    }
    weights[0] = intercept as f32;
    weights
}

fn predict(row: &[f64; MODEL_FEATURE_COUNT], weights: &[f32; MODEL_FEATURE_COUNT]) -> f64 {
    row.iter()
        .zip(weights)
        .map(|(feature, weight)| feature * f64::from(*weight))
        .sum()
}

/// Recompute the standardized normal matrix for the final model and invert it,
/// capturing everything needed to judge a new ligand's distance from the
/// training set and to form a real prediction interval.
fn training_geometry(
    features: &[[f64; MODEL_FEATURE_COUNT]],
    targets: &[f64],
    selected: &[usize],
    predictions: &[f64],
) -> Option<TrainingGeometry> {
    let (means, scales) = standardization(features, selected);
    let columns = selected.len() + 1;
    let observations = features.len();
    if observations <= columns {
        // A saturated fit has no residual degrees of freedom, so neither the
        // residual standard error nor an interval derived from it is defined.
        return None;
    }
    let mut normal = vec![vec![0.0; columns]; columns];
    for row in features {
        let mut design = Vec::with_capacity(columns);
        design.push(1.0);
        design.extend(
            selected
                .iter()
                .map(|&column| (row[column] - means[column]) / scales[column]),
        );
        for left in 0..columns {
            for right in 0..columns {
                normal[left][right] += design[left] * design[right];
            }
        }
    }
    // Mirrors the ridge floor `fit_linear` applies, so leverage is the leverage
    // of the model actually fitted rather than of an idealized one.
    for (index, row) in normal.iter_mut().enumerate().skip(1) {
        row[index] += EPSILON;
    }
    let xtx_inverse = invert_matrix(&normal).ok()?;
    // The same standardized coordinates the fit solved in, kept so a screened
    // ligand's distance to the training set can be measured in that frame.
    let standardized_training_points = features
        .iter()
        .map(|row| {
            selected
                .iter()
                .map(|&column| (row[column] - means[column]) / scales[column])
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let neighbor_calibration = NeighborCalibration::from_points(&standardized_training_points);
    let rss = residual_sum_of_squares(targets, predictions);
    let residual_standard_error = (rss / (observations - columns) as f64).sqrt();
    if !residual_standard_error.is_finite() {
        return None;
    }
    Some(TrainingGeometry {
        feature_indices: selected.to_vec(),
        means: selected.iter().map(|&column| means[column]).collect(),
        scales: selected.iter().map(|&column| scales[column]).collect(),
        xtx_inverse,
        observations,
        parameters: columns,
        residual_standard_error,
        warning_leverage: 3.0 * columns as f64 / observations as f64,
        standardized_training_points,
        neighbor_calibration,
    })
}

fn solve_linear_system(mut matrix: Vec<Vec<f64>>, mut rhs: Vec<f64>) -> Result<Vec<f64>, String> {
    let size = rhs.len();
    for pivot in 0..size {
        let best = (pivot..size)
            .max_by(|&left, &right| {
                matrix[left][pivot]
                    .abs()
                    .total_cmp(&matrix[right][pivot].abs())
            })
            .ok_or_else(|| "empty normal equation".to_string())?;
        if matrix[best][pivot].abs() <= EPSILON {
            return Err("descriptor matrix is singular".into());
        }
        matrix.swap(pivot, best);
        rhs.swap(pivot, best);
        let divisor = matrix[pivot][pivot];
        for value in &mut matrix[pivot][pivot..] {
            *value /= divisor;
        }
        rhs[pivot] /= divisor;
        let pivot_row = matrix[pivot].clone();
        for row in 0..size {
            if row == pivot {
                continue;
            }
            let factor = matrix[row][pivot];
            for (value, pivot_value) in matrix[row][pivot..].iter_mut().zip(&pivot_row[pivot..]) {
                *value -= factor * pivot_value;
            }
            rhs[row] -= factor * rhs[pivot];
        }
    }
    Ok(rhs)
}

fn metrics(actual: &[f64], predicted: &[f64]) -> ModelMetrics {
    let count = actual.len();
    let mean = actual.iter().sum::<f64>() / count as f64;
    let rss = residual_sum_of_squares(actual, predicted);
    let total = actual
        .iter()
        .map(|value| (value - mean).powi(2))
        .sum::<f64>();
    ModelMetrics {
        count,
        r2: (total > EPSILON).then_some(1.0 - rss / total),
        mae: actual
            .iter()
            .zip(predicted)
            .map(|(actual, predicted)| (actual - predicted).abs())
            .sum::<f64>()
            / count as f64,
        rmse: (rss / count as f64).sqrt(),
    }
}

fn residual_sum_of_squares(actual: &[f64], predicted: &[f64]) -> f64 {
    actual
        .iter()
        .zip(predicted)
        .map(|(actual, predicted)| (actual - predicted).powi(2))
        .sum()
}

fn bic(rss: f64, count: usize, parameter_count: usize) -> f64 {
    count as f64 * (rss.max(EPSILON) / count as f64).ln()
        + parameter_count as f64 * (count as f64).ln()
}

fn correlation(features: &[[f64; MODEL_FEATURE_COUNT]], left: usize, right: usize) -> f64 {
    let left_mean = features.iter().map(|row| row[left]).sum::<f64>() / features.len() as f64;
    let right_mean = features.iter().map(|row| row[right]).sum::<f64>() / features.len() as f64;
    let numerator = features
        .iter()
        .map(|row| (row[left] - left_mean) * (row[right] - right_mean))
        .sum::<f64>();
    let left_norm = features
        .iter()
        .map(|row| (row[left] - left_mean).powi(2))
        .sum::<f64>()
        .sqrt();
    let right_norm = features
        .iter()
        .map(|row| (row[right] - right_mean).powi(2))
        .sum::<f64>()
        .sqrt();
    if left_norm <= EPSILON || right_norm <= EPSILON {
        0.0
    } else {
        numerator / (left_norm * right_norm)
    }
}

fn correlation_matrix(features: &[[f64; MODEL_FEATURE_COUNT]]) -> Vec<Vec<f64>> {
    (0..MODEL_FEATURE_COUNT)
        .map(|left| {
            (0..MODEL_FEATURE_COUNT)
                .map(|right| {
                    if left == right && feature_scale(features, left) > EPSILON {
                        1.0
                    } else {
                        correlation(features, left, right)
                    }
                })
                .collect()
        })
        .collect()
}

fn vif(features: &[[f64; MODEL_FEATURE_COUNT]], selected: &[usize]) -> Vec<Option<f64>> {
    let mut result = vec![None; MODEL_FEATURE_COUNT];
    for &target_column in selected {
        let predictors = selected
            .iter()
            .copied()
            .filter(|column| *column != target_column)
            .collect::<Vec<_>>();
        if predictors.is_empty() {
            result[target_column] = Some(1.0);
            continue;
        }
        let target = features
            .iter()
            .map(|row| row[target_column])
            .collect::<Vec<_>>();
        if let Ok(fit) = fit_linear(features, &target, &predictors, 0.0) {
            let r2 = metrics(&target, &fit.predictions).r2.unwrap_or(0.0);
            result[target_column] = Some(1.0 / (1.0 - r2).max(1.0e-9));
        }
    }
    result
}

fn feature_bounds(features: &[[f64; MODEL_FEATURE_COUNT]], column: usize) -> (f64, f64) {
    features.iter().map(|row| row[column]).fold(
        (f64::INFINITY, f64::NEG_INFINITY),
        |(minimum, maximum), value| (minimum.min(value), maximum.max(value)),
    )
}

fn bootstrap_intervals(
    features: &[[f64; MODEL_FEATURE_COUNT]],
    targets: &[f64],
    selected: &[usize],
    estimates: &[f32; MODEL_FEATURE_COUNT],
    samples: usize,
    seed: u64,
) -> Vec<CoefficientInterval> {
    let mut rng = DeterministicRng::new(seed);
    let mut distributions = (0..MODEL_FEATURE_COUNT)
        .map(|_| Vec::with_capacity(samples))
        .collect::<Vec<_>>();
    for _ in 0..samples {
        let indices = (0..features.len())
            .map(|_| rng.index(features.len()))
            .collect::<Vec<_>>();
        let sampled_features = indices
            .iter()
            .map(|&index| features[index])
            .collect::<Vec<_>>();
        let sampled_targets = indices
            .iter()
            .map(|&index| targets[index])
            .collect::<Vec<_>>();
        if let Ok(fit) = fit_linear(&sampled_features, &sampled_targets, selected, 1.0e-8) {
            for (column, values) in distributions.iter_mut().enumerate() {
                values.push(f64::from(fit.weights[column]));
            }
        }
    }
    let reported = std::iter::once(0)
        .chain(selected.iter().copied())
        .collect::<Vec<_>>();
    reported
        .into_iter()
        .map(|column| {
            let values = &mut distributions[column];
            values.sort_by(f64::total_cmp);
            CoefficientInterval {
                feature: MODEL_FEATURE_NAMES[column].into(),
                estimate: f64::from(estimates[column]),
                lower_95: percentile(values, 0.025).unwrap_or(f64::from(estimates[column])),
                upper_95: percentile(values, 0.975).unwrap_or(f64::from(estimates[column])),
            }
        })
        .collect()
}

fn permutation_test(
    features: &[[f64; MODEL_FEATURE_COUNT]],
    targets: &[f64],
    selected: &[usize],
    observed_r2: f64,
    samples: usize,
    seed: u64,
) -> f64 {
    let mut rng = DeterministicRng::new(seed);
    let mut shuffled = targets.to_vec();
    let mut extreme = 0_usize;
    for _ in 0..samples {
        for index in (1..shuffled.len()).rev() {
            let other = rng.index(index + 1);
            shuffled.swap(index, other);
        }
        if let Ok(fit) = fit_linear(features, &shuffled, selected, 0.0)
            && metrics(&shuffled, &fit.predictions).r2.unwrap_or(0.0) >= observed_r2
        {
            extreme += 1;
        }
    }
    (extreme + 1) as f64 / (samples + 1) as f64
}

fn percentile(values: &[f64], probability: f64) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    let position = probability * (values.len() - 1) as f64;
    let lower = position.floor() as usize;
    let upper = position.ceil() as usize;
    let fraction = position - lower as f64;
    Some(values[lower] * (1.0 - fraction) + values[upper] * fraction)
}

fn soft_threshold(value: f64, threshold: f64) -> f64 {
    if value > threshold {
        value - threshold
    } else if value < -threshold {
        value + threshold
    } else {
        0.0
    }
}

struct DeterministicRng(u64);

impl DeterministicRng {
    fn new(seed: u64) -> Self {
        Self(seed.max(1))
    }

    fn next(&mut self) -> u64 {
        let mut value = self.0;
        value ^= value << 13;
        value ^= value >> 7;
        value ^= value << 17;
        self.0 = value;
        value
    }

    fn index(&mut self, upper: usize) -> usize {
        (self.next() % upper as u64) as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn recovers_a_sparse_interpretable_model() {
        let records = (0..30)
            .map(|index| {
                let l = 1.0 + index as f32 * 0.1;
                PackedReactionRecord {
                    l,
                    b1: 2.0 + (index % 3) as f32 * 0.1,
                    b5: 4.0,
                    nbo_charge: -0.5,
                    ir_freq: 1_650.0,
                    exp_ddg: 0.75 + 1.5 * l,
                    ..PackedReactionRecord::default()
                }
            })
            .collect::<Vec<_>>();
        let options = FitOptions {
            bootstrap_samples: 50,
            permutation_samples: 50,
            ..FitOptions::default()
        };
        let report =
            fit_scientific_model(&records, &(0..records.len()).collect::<Vec<_>>(), options)
                .unwrap();
        assert!(report.selected_features.contains(&"L_boltz".to_string()));
        assert!(report.training.r2.unwrap() > 0.999);
        assert!(report.fixed_feature_loo.r2.unwrap() > 0.999);
        assert!(report.fixed_feature_group_loo.r2.unwrap() > 0.999);
        assert_eq!(report.training_group_count, records.len());
        assert!(report.weights.iter().all(|weight| weight.is_finite()));
    }
}
