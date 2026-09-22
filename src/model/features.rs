use crate::storage::PackedReactionRecord;
use serde::{Deserialize, Serialize};

/// Meaning declared for the stored Sterimol inputs, not inferred from names.
/// A declaration records the data provider's method; it does not validate a
/// conformer population or establish thermodynamic equilibrium.
#[derive(Clone, Copy, Debug, Default, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum DescriptorAggregation {
    /// Older model files did not record the method.
    #[default]
    Unknown,
    /// Use supplied record values; their source aggregation is not established.
    SuppliedRecordValues,
    /// Each supplied training descriptor was computed from one geometry.
    SingleGeometry,
    /// Descriptors are means under explicitly supplied conformer weights.
    SuppliedWeightMean,
}

impl DescriptorAggregation {
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Unknown => "unknown",
            Self::SuppliedRecordValues => "supplied_record_values",
            Self::SingleGeometry => "single_geometry",
            Self::SuppliedWeightMean => "supplied_weight_mean",
        }
    }

    pub fn from_label(value: &str) -> Result<Self, String> {
        match value {
            "unknown" => Ok(Self::Unknown),
            "supplied_record_values" => Ok(Self::SuppliedRecordValues),
            "single_geometry" => Ok(Self::SingleGeometry),
            "supplied_weight_mean" => Ok(Self::SuppliedWeightMean),
            _ => Err(format!("unsupported descriptor aggregation `{value}`")),
        }
    }
}

/// Number of coefficients in the physical-organic regression model.
pub const MODEL_FEATURE_COUNT: usize = 8;

/// Stable names for the physical-organic regression columns.
pub const MODEL_FEATURE_NAMES: [&str; MODEL_FEATURE_COUNT] = [
    "intercept",
    "L_boltz",
    "B1_boltz",
    "B5_boltz",
    "nbo_charge",
    "B1_x_nbo_charge",
    "B5_x_nbo_charge",
    "ir_frequency",
];

/// Expands one packed observation into the fixed regression feature vector.
///
/// The leading constant carries the model intercept. The interaction terms
/// couple Sterimol widths to donor electronics in the style of interpretable
/// physical-organic linear free-energy relationships.
#[must_use]
#[inline]
pub fn expand_features(record: &PackedReactionRecord) -> [f32; 8] {
    [
        1.0,
        record.l,
        record.b1,
        record.b5,
        record.nbo_charge,
        record.b1 * record.nbo_charge,
        record.b5 * record.nbo_charge,
        record.ir_freq,
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn creates_exact_physical_organic_interactions() {
        let record = PackedReactionRecord {
            l: 2.0,
            b1: 3.0,
            b5: 5.0,
            nbo_charge: -0.2,
            ir_freq: 1_600.0,
            ..PackedReactionRecord::default()
        };

        assert_eq!(
            expand_features(&record),
            [1.0, 2.0, 3.0, 5.0, -0.2, -0.6, -1.0, 1_600.0]
        );
    }
}
