//! High-performance molecular steric featurization, packed storage, regression,
//! and kinetic simulation.

pub mod geometry;
pub mod kinetics;
pub mod model;
#[cfg(feature = "profiling")]
pub mod profiling;
pub mod storage;

/// Time a function or lexical scope in diagnostic builds only.
///
/// Stage totals use exclusive time, so nested scopes are not double counted.
#[macro_export]
macro_rules! profile_scope {
    ($guard:ident, $stage:literal, $function:literal) => {
        #[cfg(feature = "profiling")]
        let $guard = {
            static SITE: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
            $crate::profiling::Span::enter(&SITE, $stage, $function)
        };
    };
    ($stage:literal, $function:literal) => {
        #[cfg(feature = "profiling")]
        let _stericx_profile_span = {
            static SITE: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
            $crate::profiling::Span::enter(&SITE, $stage, $function)
        };
    };
}

/// End a named diagnostic phase before its enclosing function returns.
/// Use with `profile_scope!(guard_name, "stage", "function")`.
#[macro_export]
macro_rules! profile_end {
    ($guard:ident) => {
        #[cfg(feature = "profiling")]
        drop($guard);
    };
}

pub use geometry::{
    Atom, BuriedVolumeCalculator, BuriedVolumeConfig, BuriedVolumeEnsembleParams,
    BuriedVolumeError, BuriedVolumeParams, Molecule, PyramidalizationCalculator,
    PyramidalizationParams, SterimolCalculator, SterimolParams, bonded_neighbors,
    coordination_center, covalent_radius, parse_coordinate_file,
};
pub use kinetics::{EyringKineticLink, ProductRatio};
pub use model::{
    EvaluationSummary, FitOptions, FrozenPrediction, ReactionLabel, RegressXPredictor,
    ScientificFitReport, ScoredPrediction, TrainedModel, TrainingSplit, fit_scientific_model,
    fit_scientific_model_grouped, score_frozen_predictions, train_scientific_model,
};
pub use storage::{
    PackedBuriedVolumeRecord, PackedReactionRecord, PackedReactionRecordV2, SigPackHeaderV2,
    SigPackReader, SigPackV2Reader, SigPackV2Writer, SigPackWriter,
};
