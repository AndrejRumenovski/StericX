//! Feature construction, model training, inference, and evaluation.
//!
//! The layering is deliberate:
//!
//! * [`fit`] owns the scientific methodology — descriptor selection,
//!   regularized baselines, and validation statistics;
//! * [`dataset`] turns row provenance into a train/frozen partition;
//! * [`training`] composes those two into one reusable entry point,
//!   [`training::train_scientific_model`];
//! * [`evaluation`] scores frozen predictions once targets are revealed;
//! * [`regress`] evaluates a fitted model at high throughput.
//!
//! Front ends such as the `stericx` command line are expected to handle only
//! I/O and presentation, and to call [`training::train_scientific_model`] and
//! [`evaluation::score_frozen_predictions`] for everything else.

pub mod dataset;
pub mod evaluation;
mod features;
mod fit;
mod regress;
pub mod training;

pub use dataset::{ReactionLabel, TrainingSplit, is_supported_split};
pub use evaluation::{EvaluationSummary, ScoredPrediction, score_frozen_predictions};
pub use features::{MODEL_FEATURE_COUNT, MODEL_FEATURE_NAMES, expand_features};
pub use fit::{
    BaselineReport, CoefficientInterval, FeatureDomain, FitOptions, ModelMetrics,
    ScientificFitReport, fit_scientific_model, fit_scientific_model_grouped,
};
pub use regress::RegressXPredictor;
pub use training::{FrozenPrediction, TrainedModel, applicability_status, train_scientific_model};
