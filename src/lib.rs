//! High-performance molecular steric featurization, packed storage, regression,
//! and kinetic simulation.

pub mod geometry;
pub mod kinetics;
pub mod model;
pub mod storage;

pub use geometry::{
    Atom, BuriedVolumeCalculator, BuriedVolumeConfig, BuriedVolumeEnsembleParams,
    BuriedVolumeError, BuriedVolumeParams, Molecule, SterimolCalculator, SterimolParams,
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
