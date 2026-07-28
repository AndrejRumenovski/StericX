//! High-performance molecular steric featurization, packed storage, regression,
//! and kinetic simulation.

pub mod geometry;
pub mod kinetics;
pub mod model;
pub mod storage;

pub use geometry::{
    Atom, BuriedVolumeCalculator, BuriedVolumeConfig, BuriedVolumeEnsembleParams,
    BuriedVolumeError, BuriedVolumeParams, Molecule, SterimolCalculator, SterimolParams,
    coordination_center, covalent_radius, parse_coordinate_file,
};
pub use kinetics::{EyringKineticLink, ProductRatio};
pub use model::{
    FitOptions, RegressXPredictor, ScientificFitReport, fit_scientific_model,
    fit_scientific_model_grouped,
};
pub use storage::{
    PackedBuriedVolumeRecord, PackedReactionRecord, PackedReactionRecordV2, SigPackHeaderV2,
    SigPackReader, SigPackV2Reader, SigPackV2Writer, SigPackWriter,
};
