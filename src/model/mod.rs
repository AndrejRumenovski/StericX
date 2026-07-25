//! Feature construction and vectorized multivariate linear regression.

mod features;
mod fit;
mod regress;

pub use features::{MODEL_FEATURE_COUNT, MODEL_FEATURE_NAMES, expand_features};
pub use fit::{
    BaselineReport, CoefficientInterval, FeatureDomain, FitOptions, ModelMetrics,
    ScientificFitReport, fit_scientific_model, fit_scientific_model_grouped,
};
pub use regress::RegressXPredictor;
