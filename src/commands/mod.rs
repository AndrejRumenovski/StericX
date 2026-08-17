//! One module per `stericx` subcommand handler.

pub(crate) mod buried_volume;
pub(crate) mod evaluate;
pub(crate) mod fit;
pub(crate) mod parse;
pub(crate) mod predict;
pub(crate) mod search;
pub(crate) mod simulate;

use serde::{Deserialize, Serialize};

/// A frozen non-training prediction: written by `fit` and re-read by `evaluate`.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub(crate) struct FrozenPredictionRow {
    #[serde(rename = "Reaction_ID")]
    pub(crate) reaction_id: String,
    #[serde(rename = "Ligand_Group")]
    pub(crate) ligand_group: String,
    #[serde(rename = "Dataset_Split")]
    pub(crate) dataset_split: String,
    #[serde(rename = "Predicted_ddG_kcal_mol")]
    pub(crate) predicted_ddg: f32,
    #[serde(rename = "Applicability_Domain")]
    pub(crate) applicability_domain: String,
}
