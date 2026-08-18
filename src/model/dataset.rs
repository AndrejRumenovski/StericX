//! Row labels and the train/frozen partition consumed by model training.
//!
//! The packed `.sigpack` matrix carries descriptors and targets but no
//! provenance. Callers supply one [`ReactionLabel`] per record, in the same
//! order, and [`TrainingSplit`] derives the partition used for fitting and for
//! freezing predictions.

use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

/// Dataset split values accepted by the training pipeline.
pub const SUPPORTED_DATASET_SPLITS: [&str; 4] = ["train", "external", "blind", "test"];

/// Returns whether `value` names a supported dataset split, ignoring case.
#[must_use]
pub fn is_supported_split(value: &str) -> bool {
    SUPPORTED_DATASET_SPLITS
        .iter()
        .any(|split| value.eq_ignore_ascii_case(split))
}

/// Provenance for one reaction row, aligned by position with a record matrix.
///
/// Field names and aliases match the reaction/metadata CSV schema so the same
/// type deserializes both the raw preparation table and the split file.
#[derive(Clone, Debug, Default, Deserialize, Serialize, PartialEq, Eq)]
pub struct ReactionLabel {
    /// Stable reaction identifier.
    #[serde(rename = "Reaction_ID", alias = "reaction_id")]
    pub reaction_id: String,
    /// Dataset split; only `train` contributes to model fitting.
    #[serde(rename = "Dataset_Split", alias = "dataset_split")]
    pub dataset_split: String,
    /// Scaffold-group label used for leave-one-scaffold-group-out validation.
    #[serde(rename = "Ligand_Group", alias = "ligand_group", default)]
    pub ligand_group: String,
}

impl ReactionLabel {
    /// Creates a label from its three components.
    #[must_use]
    pub fn new(
        reaction_id: impl Into<String>,
        dataset_split: impl Into<String>,
        ligand_group: impl Into<String>,
    ) -> Self {
        Self {
            reaction_id: reaction_id.into(),
            dataset_split: dataset_split.into(),
            ligand_group: ligand_group.into(),
        }
    }

    /// Returns whether this row belongs to the training partition.
    #[must_use]
    pub fn is_training(&self) -> bool {
        self.dataset_split.eq_ignore_ascii_case("train")
    }

    /// Returns the grouping label, falling back to the reaction identifier.
    ///
    /// An absent or whitespace-only `Ligand_Group` makes the row its own group,
    /// which degrades group validation to plain leave-one-out for that row
    /// rather than silently merging unrelated scaffolds.
    #[must_use]
    pub fn training_group(&self) -> &str {
        let group = self.ligand_group.trim();
        if group.is_empty() {
            self.reaction_id.as_str()
        } else {
            group
        }
    }
}

/// Partition of a record matrix into training rows and frozen prediction rows.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TrainingSplit {
    training_indices: Vec<usize>,
    frozen_indices: Vec<usize>,
    training_groups: Vec<String>,
}

impl TrainingSplit {
    /// Derives the partition from positionally aligned row labels.
    ///
    /// Rows split as `train` become training rows; every other supported split
    /// is frozen for later revelation. Fails when no row is left to freeze,
    /// because a run with nothing held out cannot produce a blind artifact.
    pub fn from_labels(labels: &[ReactionLabel]) -> Result<Self, String> {
        let training_indices = labels
            .iter()
            .enumerate()
            .filter(|(_, label)| label.is_training())
            .map(|(index, _)| index)
            .collect::<Vec<_>>();
        let frozen_indices = labels
            .iter()
            .enumerate()
            .filter(|(_, label)| !label.is_training())
            .map(|(index, _)| index)
            .collect::<Vec<_>>();
        if frozen_indices.is_empty() {
            return Err("metadata contains no non-training rows to freeze".into());
        }
        let training_groups = training_indices
            .iter()
            .map(|&index| labels[index].training_group().to_owned())
            .collect::<Vec<_>>();
        Ok(Self {
            training_indices,
            frozen_indices,
            training_groups,
        })
    }

    /// Ascending record indices used to fit the model.
    #[must_use]
    pub fn training_indices(&self) -> &[usize] {
        &self.training_indices
    }

    /// Ascending record indices whose predictions are frozen before revelation.
    #[must_use]
    pub fn frozen_indices(&self) -> &[usize] {
        &self.frozen_indices
    }

    /// Group labels aligned with [`Self::training_indices`].
    #[must_use]
    pub fn training_groups(&self) -> &[String] {
        &self.training_groups
    }

    /// Number of distinct training group labels.
    #[must_use]
    pub fn training_group_count(&self) -> usize {
        self.training_groups.iter().collect::<BTreeSet<_>>().len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn labels() -> Vec<ReactionLabel> {
        vec![
            ReactionLabel::new("A", "train", "scaffold_1"),
            ReactionLabel::new("B", "TRAIN", "scaffold_1"),
            ReactionLabel::new("C", "train", "   "),
            ReactionLabel::new("D", "blind", "scaffold_2"),
        ]
    }

    #[test]
    fn partitions_rows_and_falls_back_to_reaction_id_groups() {
        let split = TrainingSplit::from_labels(&labels()).unwrap();

        assert_eq!(split.training_indices(), [0, 1, 2]);
        assert_eq!(split.frozen_indices(), [3]);
        assert_eq!(split.training_groups(), ["scaffold_1", "scaffold_1", "C"]);
        assert_eq!(split.training_group_count(), 2);
    }

    #[test]
    fn rejects_a_dataset_without_frozen_rows() {
        let labels = vec![
            ReactionLabel::new("A", "train", "g"),
            ReactionLabel::new("B", "train", "g"),
        ];

        assert_eq!(
            TrainingSplit::from_labels(&labels).unwrap_err(),
            "metadata contains no non-training rows to freeze"
        );
    }

    #[test]
    fn recognizes_supported_splits_case_insensitively() {
        assert!(is_supported_split("Train"));
        assert!(is_supported_split("BLIND"));
        assert!(is_supported_split("external"));
        assert!(is_supported_split("test"));
        assert!(!is_supported_split("validation"));
    }
}
