//! Versioned portable model format.
//!
//! A portable model is the model artifact written by `stericx fit` **plus** the
//! provenance and inference metadata needed to score new structures on another
//! machine, without the training data, the training code, or the study driver.
//!
//! # Versioning
//!
//! The format extends the existing artifact instead of replacing it. A
//! `schema_version: 2` document contains every key a `schema_version: 1`
//! document contains, in the same place, and adds three sections:
//!
//! ```text
//! {
//!   "schema_version": 2,
//!   ... every schema_version 1 field, unchanged ...
//!   "inference":  { response, feature space, terms, standardization, ranges },
//!   "provenance": { model id, StericX version, training data digests, reaction },
//!   "created":    { timestamp, producer }
//! }
//! ```
//!
//! Consequences of that layout:
//!
//! * any reader of the version 1 artifact keeps working against version 2;
//! * a version 1 artifact loads here as a [`PortableModel`] that reports itself
//!   as legacy, with the three sections absent rather than invented.
//!
//! [`PortableModel::from_json`] rejects any `schema_version` above
//! [`PORTABLE_SCHEMA_VERSION`], so a future format is never partially read.
//!
//! # Numeric fidelity
//!
//! Reading a model and writing it again must not move a single bit of a
//! validation statistic. That requires `serde_json`'s `float_roundtrip`
//! feature: the default parser is not correctly rounded and shifts some
//! seventeen-digit `f64` values by one unit in the last place. The feature is
//! enabled in `Cargo.toml` and guarded by a test.
//!
//! # Trust boundary
//!
//! The three added sections restate values that also live in the flattened fit
//! report — coefficients, standardization, and descriptor ranges.
//! [`PortableModel::validate`] requires the two copies to agree, which turns
//! that redundancy into an integrity check: a document whose `inference` block
//! was edited without editing `weights` is rejected as malformed rather than
//! silently scoring with one of the two.

use super::{
    FitOptions, MODEL_FEATURE_COUNT, MODEL_FEATURE_NAMES, RegressXPredictor, ScientificFitReport,
};
use serde::{Deserialize, Serialize};
use std::error::Error;
use std::fmt::{Display, Formatter};
use std::time::{SystemTime, UNIX_EPOCH};

/// Highest document schema version this build writes and accepts.
pub const PORTABLE_SCHEMA_VERSION: u32 = 2;

/// Lowest document schema version this build accepts.
///
/// Version 1 is the pre-portable artifact: readable, but not self-sufficient
/// for inference on another machine.
pub const LEGACY_SCHEMA_VERSION: u32 = 1;

/// Relative slack allowed when cross-checking the `inference` section against
/// the flattened fit report. Both copies are machine generated from one value,
/// so this tolerates re-serialization only, not edited numbers.
const CROSS_CHECK_TOLERANCE: f64 = 1.0e-12;

/// Failure to read, write, or validate a portable model document.
#[derive(Debug)]
pub enum ModelFormatError {
    /// The document declares a schema this build cannot interpret.
    UnsupportedSchemaVersion { found: u32, maximum: u32 },
    /// A schema version 2 document is missing a section it must carry.
    MissingSection { section: &'static str },
    /// The document is structurally readable but scientifically invalid.
    Malformed(String),
    /// The document is not well-formed JSON, or a required field is absent.
    Json(serde_json::Error),
}

impl Display for ModelFormatError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnsupportedSchemaVersion { found, maximum } => write!(
                formatter,
                "model schema version {found} is newer than the supported maximum {maximum}"
            ),
            Self::MissingSection { section } => write!(
                formatter,
                "portable model is missing the required `{section}` section"
            ),
            Self::Malformed(message) => write!(formatter, "malformed model: {message}"),
            Self::Json(error) => write!(formatter, "model document could not be parsed: {error}"),
        }
    }
}

impl Error for ModelFormatError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Json(error) => Some(error),
            _ => None,
        }
    }
}

impl From<serde_json::Error> for ModelFormatError {
    fn from(value: serde_json::Error) -> Self {
        Self::Json(value)
    }
}

/// What the model predicts, in enough detail to interpret a number.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
pub struct ResponseSpec {
    /// Machine-readable response name, for example `ddg_double_dagger`.
    pub name: String,
    /// Physical units of the predicted value.
    pub units: String,
    /// Human-readable definition of the quantity.
    pub description: String,
    /// Sign convention, so a consumer cannot invert the selectivity.
    pub sign_convention: String,
    /// Temperature the response refers to, when the study fixes one.
    pub temperature_k: Option<f32>,
}

impl ResponseSpec {
    /// The response every current StericX model is fitted against.
    #[must_use]
    pub fn transition_state_energy_difference(temperature_k: Option<f32>) -> Self {
        Self {
            name: "ddg_double_dagger".into(),
            units: "kcal/mol".into(),
            description: "Transition-state free-energy difference between competing \
                          enantiomeric pathways."
                .into(),
            sign_convention: "ddG = G(S) - G(R); positive values favor the R product.".into(),
            temperature_k,
        }
    }
}

/// How one column of the model feature vector is built from a record.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum FeatureTransform {
    /// The literal `1.0` column carrying the intercept.
    Constant,
    /// A descriptor copied straight from the packed record.
    Descriptor { descriptor: String },
    /// A product of two or more record descriptors.
    Interaction { factors: Vec<String> },
}

/// The complete feature space a model is defined over.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
pub struct FeatureSpaceSpec {
    /// Stable identifier for this feature construction.
    pub definition: String,
    /// Column names, in model order.
    pub feature_names: Vec<String>,
    /// Column construction rules, aligned with `feature_names`.
    pub transformations: Vec<FeatureTransform>,
}

impl FeatureSpaceSpec {
    /// The eight-column physical-organic feature space used by StericX.
    #[must_use]
    pub fn physical_organic_v1() -> Self {
        Self {
            definition: "stericx.physical_organic.v1".into(),
            feature_names: MODEL_FEATURE_NAMES
                .iter()
                .map(|name| (*name).to_owned())
                .collect(),
            transformations: vec![
                FeatureTransform::Constant,
                FeatureTransform::Descriptor {
                    descriptor: "sterimol_l".into(),
                },
                FeatureTransform::Descriptor {
                    descriptor: "sterimol_b1".into(),
                },
                FeatureTransform::Descriptor {
                    descriptor: "sterimol_b5".into(),
                },
                FeatureTransform::Descriptor {
                    descriptor: "nbo_charge".into(),
                },
                FeatureTransform::Interaction {
                    factors: vec!["sterimol_b1".into(), "nbo_charge".into()],
                },
                FeatureTransform::Interaction {
                    factors: vec!["sterimol_b5".into(), "nbo_charge".into()],
                },
                FeatureTransform::Descriptor {
                    descriptor: "ir_frequency".into(),
                },
            ],
        }
    }
}

/// One selected descriptor and its raw-scale coefficient.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
pub struct ModelTerm {
    pub feature_index: usize,
    pub feature_name: String,
    /// Coefficient on the raw descriptor scale, directly usable for inference.
    pub coefficient: f64,
    /// Training mean subtracted before the model was fitted.
    pub training_mean: f64,
    /// Training standard deviation used to scale the descriptor.
    pub training_standard_deviation: f64,
    /// Smallest training value, for the applicability-domain check.
    pub training_minimum: f64,
    /// Largest training value, for the applicability-domain check.
    pub training_maximum: f64,
}

/// Everything required to evaluate the model on a new structure.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
pub struct InferenceSpec {
    pub response: ResponseSpec,
    pub feature_space: FeatureSpaceSpec,
    /// Raw-scale intercept.
    pub intercept: f64,
    /// Selected descriptors, in selection order.
    pub terms: Vec<ModelTerm>,
}

/// A content digest for one training input artifact.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
pub struct DatasetDigest {
    /// Artifact name, for example `reactions.sigpack`.
    pub artifact: String,
    /// Digest algorithm, for example `sha256` or `fnv1a64`.
    pub algorithm: String,
    /// Lower-case hexadecimal digest.
    pub digest: String,
    /// Size of the digested artifact in bytes.
    pub byte_count: u64,
}

/// Training-set provenance and the configuration that produced the fit.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
pub struct TrainingProvenance {
    pub record_count: usize,
    pub group_count: usize,
    /// One digest per training input. Never empty in a valid document.
    pub dataset_digests: Vec<DatasetDigest>,
    /// Fitting configuration, so the model can be re-derived.
    pub fit_options: FitOptions,
}

/// Chemistry context that cannot be derived from the fit itself.
///
/// Every field is optional and is serialized explicitly as `null` when unknown.
/// Nothing here is ever defaulted to a plausible value; see
/// [`PortableModel::missing_provenance`].
#[derive(Clone, Debug, Default, Deserialize, Serialize, PartialEq)]
pub struct ReactionProvenance {
    pub reaction_family: Option<String>,
    pub catalyst_metal: Option<String>,
    pub ligand_class: Option<String>,
    pub source_url: Option<String>,
    #[serde(default)]
    pub notes: Vec<String>,
}

/// Model identity and where it came from.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
pub struct ModelProvenance {
    /// Stable identifier for this fitted model.
    pub model_id: String,
    /// Version of the StericX crate that produced the fit.
    pub stericx_version: String,
    /// Record layout the descriptors were read from.
    pub record_format: String,
    pub training: TrainingProvenance,
    pub reaction: ReactionProvenance,
}

/// When and by what the document was written.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
pub struct CreationMetadata {
    /// RFC 3339 UTC timestamp.
    pub created_utc: String,
    /// Producing command or program.
    pub produced_by: String,
}

impl CreationMetadata {
    /// Stamps the current UTC time.
    #[must_use]
    pub fn now(produced_by: impl Into<String>) -> Self {
        Self {
            created_utc: rfc3339_utc(SystemTime::now()),
            produced_by: produced_by.into(),
        }
    }
}

/// Whether a finding blocks use of the model or merely limits it.
#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Severity {
    /// The model cannot be trusted for inference.
    Error,
    /// The model is usable but records less than it should.
    Warning,
}

/// One problem found in a model document.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
pub struct ModelIssue {
    pub severity: Severity,
    /// Stable machine-readable slug, e.g. `invalid_scale`.
    pub code: String,
    /// Dotted path to the offending field, e.g. `inference.terms[0]`.
    pub location: String,
    /// What is wrong, in terms a reader can act on.
    pub message: String,
}

impl ModelIssue {
    fn error(code: &str, location: impl Into<String>, message: String) -> Self {
        Self {
            severity: Severity::Error,
            code: code.to_owned(),
            location: location.into(),
            message,
        }
    }

    fn warning(code: &str, location: impl Into<String>, message: String) -> Self {
        Self {
            severity: Severity::Warning,
            code: code.to_owned(),
            location: location.into(),
            message,
        }
    }
}

/// One selected descriptor as reported by [`PortableModel::summary`].
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
pub struct DescriptorSummary {
    pub name: String,
    pub coefficient: f64,
    pub training_mean: f64,
    pub training_standard_deviation: f64,
    pub training_minimum: f64,
    pub training_maximum: f64,
}

/// A concise scientific description of a saved model.
///
/// Fields the document does not record are `None` or empty rather than filled
/// with a plausible value; a legacy model reports most of the provenance as
/// absent.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
pub struct ModelSummary {
    pub schema_version: u32,
    pub portable: bool,
    pub model: String,
    pub reaction_family: Option<String>,
    pub target: Option<ResponseSpec>,
    pub training_observations: usize,
    pub training_groups: usize,
    pub intercept: f64,
    pub descriptors: Vec<DescriptorSummary>,
    pub training_r2: Option<f64>,
    pub training_rmse: f64,
    /// Leave-one-out Q², the cross-validated coefficient of determination.
    pub loo_q2: Option<f64>,
    pub loo_rmse: f64,
    pub group_loo_q2: Option<f64>,
    pub group_loo_rmse: f64,
    pub dataset_digests: Vec<DatasetDigest>,
    pub model_id: Option<String>,
    pub stericx_version: Option<String>,
    pub created_utc: Option<String>,
    pub missing_provenance: Vec<String>,
}

/// A fitted model plus the metadata that makes it portable.
///
/// Serializes as a strict superset of the `schema_version: 1` artifact: the fit
/// report is flattened to the top level and the added sections follow it.
///
/// Deliberately not `Deserialize`. Reading goes through
/// [`PortableModel::from_json`], so a document can never reach a caller without
/// passing [`PortableModel::validate`], and each section is parsed with its own
/// concrete type instead of through the buffered `#[serde(flatten)]` path.
#[derive(Clone, Debug, Serialize)]
pub struct PortableModel {
    /// The unchanged `stericx fit` artifact, flattened to the document root.
    #[serde(flatten)]
    pub fit: ScientificFitReport,
    /// Present from schema version 2 onward.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub inference: Option<InferenceSpec>,
    /// Present from schema version 2 onward.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub provenance: Option<ModelProvenance>,
    /// Present from schema version 2 onward.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub created: Option<CreationMetadata>,
}

/// The sections a schema version 2 document adds, parsed on their own so the
/// flattened fit report can be deserialized with its exact field types.
#[derive(Deserialize)]
struct PortableSections {
    #[serde(default)]
    inference: Option<InferenceSpec>,
    #[serde(default)]
    provenance: Option<ModelProvenance>,
    #[serde(default)]
    created: Option<CreationMetadata>,
}

impl PortableModel {
    /// Builds a portable model from a fit report and explicit provenance.
    ///
    /// The caller supplies provenance rather than the builder guessing it: a
    /// model that cannot say which data trained it must not claim to know.
    pub fn from_fit_report(
        mut fit: ScientificFitReport,
        response: ResponseSpec,
        provenance: ModelProvenance,
        created: CreationMetadata,
    ) -> Result<Self, ModelFormatError> {
        let terms = fit
            .selected_feature_indices
            .iter()
            .map(|&column| {
                let domain = fit
                    .applicability_domain
                    .iter()
                    .find(|domain| domain.feature == MODEL_FEATURE_NAMES[column])
                    .ok_or_else(|| {
                        ModelFormatError::Malformed(format!(
                            "selected feature {} has no applicability domain entry",
                            MODEL_FEATURE_NAMES[column]
                        ))
                    })?;
                Ok(ModelTerm {
                    feature_index: column,
                    feature_name: MODEL_FEATURE_NAMES[column].to_owned(),
                    coefficient: f64::from(fit.weights[column]),
                    training_mean: fit.standardized_means[column],
                    training_standard_deviation: fit.standardized_scales[column],
                    training_minimum: domain.minimum,
                    training_maximum: domain.maximum,
                })
            })
            .collect::<Result<Vec<_>, ModelFormatError>>()?;

        fit.schema_version = PORTABLE_SCHEMA_VERSION;
        let model = Self {
            inference: Some(InferenceSpec {
                response,
                feature_space: FeatureSpaceSpec::physical_organic_v1(),
                intercept: f64::from(fit.weights[0]),
                terms,
            }),
            provenance: Some(provenance),
            created: Some(created),
            fit,
        };
        model.validate()?;
        Ok(model)
    }

    /// Parses and validates a model document.
    ///
    /// Rejects schema versions above [`PORTABLE_SCHEMA_VERSION`] before looking
    /// at anything else, so a future document is never partially interpreted.
    pub fn from_json(text: &str) -> Result<Self, ModelFormatError> {
        let probe: SchemaProbe = serde_json::from_str(text)?;
        if probe.schema_version > PORTABLE_SCHEMA_VERSION {
            return Err(ModelFormatError::UnsupportedSchemaVersion {
                found: probe.schema_version,
                maximum: PORTABLE_SCHEMA_VERSION,
            });
        }
        // Two typed passes over the same text. Each ignores the keys it does
        // not own, which keeps every float on its declared type.
        let fit: ScientificFitReport = serde_json::from_str(text)?;
        let sections: PortableSections = serde_json::from_str(text)?;
        let model = Self {
            fit,
            inference: sections.inference,
            provenance: sections.provenance,
            created: sections.created,
        };
        model.validate()?;
        Ok(model)
    }

    /// Serializes the document deterministically as pretty-printed JSON.
    ///
    /// Field order follows the Rust declarations, so the same model always
    /// produces byte-identical output.
    pub fn to_json(&self) -> Result<String, ModelFormatError> {
        self.validate()?;
        Ok(serde_json::to_string_pretty(self)?)
    }

    /// Returns the document schema version.
    #[must_use]
    pub fn schema_version(&self) -> u32 {
        self.fit.schema_version
    }

    /// Returns whether this document carries the portable sections.
    ///
    /// A `schema_version: 1` artifact is readable but not portable: it cannot
    /// state which data trained it or what its response means.
    #[must_use]
    pub fn is_portable(&self) -> bool {
        self.schema_version() >= PORTABLE_SCHEMA_VERSION
    }

    /// Returns the inference section, or an error for a legacy document.
    pub fn inference(&self) -> Result<&InferenceSpec, ModelFormatError> {
        self.inference
            .as_ref()
            .ok_or(ModelFormatError::MissingSection {
                section: "inference",
            })
    }

    /// Returns the provenance section, or an error for a legacy document.
    pub fn provenance(&self) -> Result<&ModelProvenance, ModelFormatError> {
        self.provenance
            .as_ref()
            .ok_or(ModelFormatError::MissingSection {
                section: "provenance",
            })
    }

    /// Builds a predictor for the model's raw-scale weights.
    #[must_use]
    pub fn predictor(&self) -> RegressXPredictor {
        RegressXPredictor::new(self.fit.weights)
    }

    /// Names the optional chemistry-context fields this model does not record.
    ///
    /// Absent context is reported, never filled in. Callers that require a
    /// complete record can refuse a model whose list is non-empty.
    #[must_use]
    pub fn missing_provenance(&self) -> Vec<&'static str> {
        let Some(provenance) = self.provenance.as_ref() else {
            return vec![
                "reaction_family",
                "catalyst_metal",
                "ligand_class",
                "source_url",
            ];
        };
        let reaction = &provenance.reaction;
        [
            ("reaction_family", reaction.reaction_family.is_none()),
            ("catalyst_metal", reaction.catalyst_metal.is_none()),
            ("ligand_class", reaction.ligand_class.is_none()),
            ("source_url", reaction.source_url.is_none()),
        ]
        .into_iter()
        .filter_map(|(name, missing)| missing.then_some(name))
        .collect()
    }

    /// Builds a concise scientific description of the model.
    ///
    /// Everything comes from the document. Values a legacy model cannot state —
    /// its response, training data, and origin — are reported as absent.
    #[must_use]
    pub fn summary(&self) -> ModelSummary {
        let fit = &self.fit;
        let inference = self.inference.as_ref();
        let provenance = self.provenance.as_ref();
        let descriptors = fit
            .selected_feature_indices
            .iter()
            .enumerate()
            .filter(|(_, column)| **column < MODEL_FEATURE_COUNT)
            .map(|(position, &column)| {
                let domain = fit.applicability_domain.get(position);
                DescriptorSummary {
                    name: MODEL_FEATURE_NAMES[column].to_owned(),
                    coefficient: f64::from(fit.weights[column]),
                    training_mean: fit.standardized_means[column],
                    training_standard_deviation: fit.standardized_scales[column],
                    training_minimum: domain.map_or(f64::NAN, |domain| domain.minimum),
                    training_maximum: domain.map_or(f64::NAN, |domain| domain.maximum),
                }
            })
            .collect();
        ModelSummary {
            schema_version: fit.schema_version,
            portable: self.is_portable(),
            model: fit.model.clone(),
            reaction_family: provenance
                .and_then(|provenance| provenance.reaction.reaction_family.clone()),
            target: inference.map(|inference| inference.response.clone()),
            training_observations: fit.training_count,
            training_groups: fit.training_group_count,
            intercept: f64::from(fit.weights[0]),
            descriptors,
            training_r2: fit.training.r2,
            training_rmse: fit.training.rmse,
            loo_q2: fit.fixed_feature_loo.r2,
            loo_rmse: fit.fixed_feature_loo.rmse,
            group_loo_q2: fit.fixed_feature_group_loo.r2,
            group_loo_rmse: fit.fixed_feature_group_loo.rmse,
            dataset_digests: provenance
                .map(|provenance| provenance.training.dataset_digests.clone())
                .unwrap_or_default(),
            model_id: provenance.map(|provenance| provenance.model_id.clone()),
            stericx_version: provenance.map(|provenance| provenance.stericx_version.clone()),
            created_utc: self
                .created
                .as_ref()
                .map(|created| created.created_utc.clone()),
            missing_provenance: self
                .missing_provenance()
                .into_iter()
                .map(str::to_owned)
                .collect(),
        }
    }

    /// Checks schema support, structural consistency, and numeric sanity.
    ///
    /// Returns the first error, preserving the typed variants callers match on.
    /// [`PortableModel::issues`] reports every problem at once instead.
    pub fn validate(&self) -> Result<(), ModelFormatError> {
        // Structural failures keep their typed variants so callers can branch
        // on them; everything else is reported through the shared collector.
        let version = self.schema_version();
        if version > PORTABLE_SCHEMA_VERSION {
            return Err(ModelFormatError::UnsupportedSchemaVersion {
                found: version,
                maximum: PORTABLE_SCHEMA_VERSION,
            });
        }
        if version >= PORTABLE_SCHEMA_VERSION {
            for (section, present) in [
                ("inference", self.inference.is_some()),
                ("provenance", self.provenance.is_some()),
                ("created", self.created.is_some()),
            ] {
                if !present {
                    return Err(ModelFormatError::MissingSection { section });
                }
            }
        }
        match self
            .issues()
            .into_iter()
            .find(|issue| issue.severity == Severity::Error)
        {
            Some(issue) => Err(ModelFormatError::Malformed(issue.message)),
            None => Ok(()),
        }
    }

    /// Reports every problem with the document, most structural first.
    ///
    /// Unlike [`PortableModel::validate`] this does not stop at the first
    /// failure, so `stericx model validate` can show a complete list.
    #[must_use]
    pub fn issues(&self) -> Vec<ModelIssue> {
        let mut issues = Vec::new();
        let version = self.schema_version();
        if version > PORTABLE_SCHEMA_VERSION {
            issues.push(ModelIssue::error(
                "unsupported_schema_version",
                "schema_version",
                format!(
                    "model schema version {version} is newer than the supported maximum \
                     {PORTABLE_SCHEMA_VERSION}"
                ),
            ));
            // Nothing below can be trusted against an unknown schema.
            return issues;
        }
        if version < LEGACY_SCHEMA_VERSION {
            issues.push(ModelIssue::error(
                "schema_version_too_old",
                "schema_version",
                format!(
                    "schema version {version} is below the minimum supported version \
                     {LEGACY_SCHEMA_VERSION}"
                ),
            ));
        }
        self.collect_fit_issues(&mut issues);
        if version >= PORTABLE_SCHEMA_VERSION {
            self.collect_section_issues(&mut issues);
        } else {
            issues.push(ModelIssue::warning(
                "legacy_schema",
                "schema_version",
                format!(
                    "schema version {version} predates the portable format: this model cannot \
                     state its response, training data, or origin"
                ),
            ));
        }
        issues
    }

    /// Checks the fields every schema version carries.
    fn collect_fit_issues(&self, issues: &mut Vec<ModelIssue>) {
        let fit = &self.fit;
        if fit.feature_names.len() != MODEL_FEATURE_COUNT {
            issues.push(ModelIssue::error(
                "feature_name_count",
                "feature_names",
                format!(
                    "expected {MODEL_FEATURE_COUNT} feature names, found {}",
                    fit.feature_names.len()
                ),
            ));
        }
        if fit.selected_feature_indices.is_empty() {
            issues.push(ModelIssue::error(
                "no_descriptors",
                "selected_feature_indices",
                "model selects no descriptors".to_owned(),
            ));
        }
        let selection_aligned = fit.selected_feature_indices.len() == fit.selected_features.len();
        if !selection_aligned {
            issues.push(ModelIssue::error(
                "selection_length_mismatch",
                "selected_features",
                "selected_feature_indices and selected_features have different lengths".to_owned(),
            ));
        }
        let domain_aligned = fit.applicability_domain.len() == fit.selected_feature_indices.len();
        if !domain_aligned {
            issues.push(ModelIssue::error(
                "domain_length_mismatch",
                "applicability_domain",
                "applicability_domain does not cover every selected descriptor".to_owned(),
            ));
        }
        if let Some(weight) = fit.weights.iter().find(|weight| !weight.is_finite()) {
            issues.push(ModelIssue::error(
                "non_finite_weight",
                "weights",
                format!("model weight {weight} is not finite"),
            ));
        }

        let mut seen: Vec<usize> = Vec::with_capacity(fit.selected_feature_indices.len());
        for (position, &column) in fit.selected_feature_indices.iter().enumerate() {
            if column == 0 || column >= MODEL_FEATURE_COUNT {
                issues.push(ModelIssue::error(
                    "invalid_feature_index",
                    format!("selected_feature_indices[{position}]"),
                    format!("selected feature index {column} is not a valid descriptor column"),
                ));
                // Every check below would index past the feature vector.
                continue;
            }
            if seen.contains(&column) {
                issues.push(ModelIssue::error(
                    "duplicate_feature_index",
                    format!("selected_feature_indices[{position}]"),
                    format!("selected feature index {column} appears more than once"),
                ));
            }
            seen.push(column);
            if selection_aligned && fit.selected_features[position] != MODEL_FEATURE_NAMES[column] {
                issues.push(ModelIssue::error(
                    "descriptor_name_mismatch",
                    format!("selected_features[{position}]"),
                    format!(
                        "selected feature {} does not match column {column} ({})",
                        fit.selected_features[position], MODEL_FEATURE_NAMES[column]
                    ),
                ));
            }
            let scale = fit.standardized_scales[column];
            if !scale.is_finite() || scale <= 0.0 {
                issues.push(ModelIssue::error(
                    "invalid_scale",
                    format!("standardized_scales[{column}]"),
                    format!(
                        "standardized scale {scale} for {} is not a positive finite number",
                        MODEL_FEATURE_NAMES[column]
                    ),
                ));
            }
            if !fit.standardized_means[column].is_finite() {
                issues.push(ModelIssue::error(
                    "non_finite_mean",
                    format!("standardized_means[{column}]"),
                    format!(
                        "standardized mean for {} is not finite",
                        MODEL_FEATURE_NAMES[column]
                    ),
                ));
            }
            if !domain_aligned {
                continue;
            }
            let domain = &fit.applicability_domain[position];
            if !domain.minimum.is_finite() || !domain.maximum.is_finite() {
                issues.push(ModelIssue::error(
                    "non_finite_range",
                    format!("applicability_domain[{position}]"),
                    format!("applicability range for {} is not finite", domain.feature),
                ));
            } else if domain.minimum > domain.maximum {
                issues.push(ModelIssue::error(
                    "inverted_range",
                    format!("applicability_domain[{position}]"),
                    format!(
                        "applicability range for {} is inverted: [{}, {}]",
                        domain.feature, domain.minimum, domain.maximum
                    ),
                ));
            }
        }
        if fit.training_count == 0 {
            issues.push(ModelIssue::error(
                "empty_training_set",
                "training_count",
                "training_count is zero".to_owned(),
            ));
        }
    }

    /// Checks the sections a schema version 2 document must carry.
    fn collect_section_issues(&self, issues: &mut Vec<ModelIssue>) {
        for (section, present) in [
            ("inference", self.inference.is_some()),
            ("provenance", self.provenance.is_some()),
            ("created", self.created.is_some()),
        ] {
            if !present {
                issues.push(ModelIssue::error(
                    "missing_section",
                    section,
                    format!("portable model is missing the required `{section}` section"),
                ));
            }
        }
        let (Some(inference), Some(provenance), Some(created)) = (
            self.inference.as_ref(),
            self.provenance.as_ref(),
            self.created.as_ref(),
        ) else {
            return;
        };

        if created.created_utc.trim().is_empty() {
            issues.push(ModelIssue::error(
                "empty_field",
                "created.created_utc",
                "created.created_utc is empty".to_owned(),
            ));
        }
        if provenance.model_id.trim().is_empty() {
            issues.push(ModelIssue::error(
                "empty_field",
                "provenance.model_id",
                "provenance.model_id is empty".to_owned(),
            ));
        }
        if provenance.stericx_version.trim().is_empty() {
            issues.push(ModelIssue::error(
                "empty_field",
                "provenance.stericx_version",
                "provenance.stericx_version is empty".to_owned(),
            ));
        }

        let feature_space = &inference.feature_space;
        if feature_space.feature_names != self.fit.feature_names {
            issues.push(ModelIssue::error(
                "feature_space_mismatch",
                "inference.feature_space.feature_names",
                "inference.feature_space.feature_names disagrees with feature_names".to_owned(),
            ));
        }
        if feature_space.transformations.len() != MODEL_FEATURE_COUNT {
            issues.push(ModelIssue::error(
                "transformation_count",
                "inference.feature_space.transformations",
                format!(
                    "expected {MODEL_FEATURE_COUNT} feature transformations, found {}",
                    feature_space.transformations.len()
                ),
            ));
        } else if !matches!(feature_space.transformations[0], FeatureTransform::Constant) {
            issues.push(ModelIssue::error(
                "missing_intercept_column",
                "inference.feature_space.transformations[0]",
                "feature column 0 must be the constant intercept column".to_owned(),
            ));
        }

        push_cross_check(
            issues,
            "inference.intercept",
            inference.intercept,
            f64::from(self.fit.weights[0]),
            "inference.intercept",
        );
        if inference.terms.len() != self.fit.selected_feature_indices.len() {
            issues.push(ModelIssue::error(
                "term_count_mismatch",
                "inference.terms",
                "inference.terms does not cover every selected descriptor".to_owned(),
            ));
        }
        for (position, term) in inference.terms.iter().enumerate() {
            let Some(&column) = self.fit.selected_feature_indices.get(position) else {
                break;
            };
            let location = format!("inference.terms[{position}]");
            if term.feature_index != column {
                issues.push(ModelIssue::error(
                    "term_index_mismatch",
                    &location,
                    format!(
                        "inference term {position} refers to column {} but the model selected \
                         {column}",
                        term.feature_index
                    ),
                ));
            }
            if column >= MODEL_FEATURE_COUNT {
                continue;
            }
            if term.feature_name != MODEL_FEATURE_NAMES[column] {
                issues.push(ModelIssue::error(
                    "term_name_mismatch",
                    &location,
                    format!(
                        "inference term {position} is named {} but column {column} is {}",
                        term.feature_name, MODEL_FEATURE_NAMES[column]
                    ),
                ));
            }
            if !term.coefficient.is_finite() {
                issues.push(ModelIssue::error(
                    "non_finite_coefficient",
                    &location,
                    format!("coefficient for {} is not finite", term.feature_name),
                ));
            }
            if !term.training_standard_deviation.is_finite()
                || term.training_standard_deviation <= 0.0
            {
                issues.push(ModelIssue::error(
                    "invalid_scale",
                    &location,
                    format!(
                        "training standard deviation for {} is not a positive finite number",
                        term.feature_name
                    ),
                ));
            }
            if term.training_minimum > term.training_maximum {
                issues.push(ModelIssue::error(
                    "inverted_range",
                    &location,
                    format!("training range for {} is inverted", term.feature_name),
                ));
            }
            push_cross_check(
                issues,
                &location,
                term.coefficient,
                f64::from(self.fit.weights[column]),
                &format!("coefficient for {}", term.feature_name),
            );
            push_cross_check(
                issues,
                &location,
                term.training_mean,
                self.fit.standardized_means[column],
                &format!("training mean for {}", term.feature_name),
            );
            push_cross_check(
                issues,
                &location,
                term.training_standard_deviation,
                self.fit.standardized_scales[column],
                &format!("training standard deviation for {}", term.feature_name),
            );
            let Some(domain) = self.fit.applicability_domain.get(position) else {
                continue;
            };
            push_cross_check(
                issues,
                &location,
                term.training_minimum,
                domain.minimum,
                &format!("training minimum for {}", term.feature_name),
            );
            push_cross_check(
                issues,
                &location,
                term.training_maximum,
                domain.maximum,
                &format!("training maximum for {}", term.feature_name),
            );
        }

        if inference.response.name.trim().is_empty() || inference.response.units.trim().is_empty() {
            issues.push(ModelIssue::error(
                "incomplete_response",
                "inference.response",
                "inference.response must name the predicted quantity and its units".to_owned(),
            ));
        }

        let training = &provenance.training;
        if training.record_count != self.fit.training_count {
            issues.push(ModelIssue::error(
                "training_count_mismatch",
                "provenance.training.record_count",
                format!(
                    "provenance records {} training rows but the fit reports {}",
                    training.record_count, self.fit.training_count
                ),
            ));
        }
        if training.group_count != self.fit.training_group_count {
            issues.push(ModelIssue::error(
                "training_group_mismatch",
                "provenance.training.group_count",
                format!(
                    "provenance records {} training groups but the fit reports {}",
                    training.group_count, self.fit.training_group_count
                ),
            ));
        }
        if training.dataset_digests.is_empty() {
            issues.push(ModelIssue::error(
                "missing_dataset_digest",
                "provenance.training.dataset_digests",
                "provenance.training.dataset_digests is empty; a portable model must identify \
                 its training data"
                    .to_owned(),
            ));
        }
        for (position, digest) in training.dataset_digests.iter().enumerate() {
            let location = format!("provenance.training.dataset_digests[{position}]");
            if digest.artifact.trim().is_empty()
                || digest.algorithm.trim().is_empty()
                || digest.digest.trim().is_empty()
            {
                issues.push(ModelIssue::error(
                    "incomplete_digest",
                    &location,
                    "every dataset digest needs an artifact, algorithm, and digest".to_owned(),
                ));
                continue;
            }
            if !digest.digest.chars().all(|c| c.is_ascii_hexdigit()) {
                issues.push(ModelIssue::error(
                    "malformed_digest",
                    &location,
                    format!("dataset digest for {} is not hexadecimal", digest.artifact),
                ));
            }
        }

        for missing in self.missing_provenance() {
            issues.push(ModelIssue::warning(
                "unrecorded_context",
                format!("provenance.reaction.{missing}"),
                format!("{missing} is not recorded; the model does not state it"),
            ));
        }
    }
}

/// Outcome of inspecting a model document that may not even parse.
#[derive(Debug)]
pub struct Diagnosis {
    /// The parsed model, absent when the document could not be read at all.
    pub model: Option<PortableModel>,
    /// Every problem found, most structural first.
    pub issues: Vec<ModelIssue>,
}

impl Diagnosis {
    /// Returns whether any finding blocks use of the model.
    #[must_use]
    pub fn has_errors(&self) -> bool {
        self.issues
            .iter()
            .any(|issue| issue.severity == Severity::Error)
    }

    /// Counts findings at one severity.
    #[must_use]
    pub fn count(&self, severity: Severity) -> usize {
        self.issues
            .iter()
            .filter(|issue| issue.severity == severity)
            .count()
    }
}

/// Diagnoses a model document, reporting specifics rather than a parse error.
///
/// The document is read in stages — JSON syntax, then `schema_version`, then
/// the typed fields, then the semantic checks — so a reader is told which field
/// is wrong instead of being handed a decoder message about the whole file.
#[must_use]
pub fn diagnose(text: &str) -> Diagnosis {
    let value: serde_json::Value = match serde_json::from_str(text) {
        Ok(value) => value,
        Err(error) => {
            return Diagnosis {
                model: None,
                issues: vec![ModelIssue::error(
                    "invalid_json",
                    format!("line {} column {}", error.line(), error.column()),
                    format!("document is not valid JSON: {error}"),
                )],
            };
        }
    };
    let Some(object) = value.as_object() else {
        return Diagnosis {
            model: None,
            issues: vec![ModelIssue::error(
                "not_an_object",
                "$",
                "a model document must be a JSON object".to_owned(),
            )],
        };
    };
    match object.get("schema_version") {
        None => {
            return Diagnosis {
                model: None,
                issues: vec![ModelIssue::error(
                    "missing_schema_version",
                    "schema_version",
                    "document has no `schema_version`, so its format cannot be determined"
                        .to_owned(),
                )],
            };
        }
        Some(version) if version.as_u64().is_none() => {
            return Diagnosis {
                model: None,
                issues: vec![ModelIssue::error(
                    "invalid_schema_version",
                    "schema_version",
                    format!("`schema_version` must be a non-negative integer, found {version}"),
                )],
            };
        }
        Some(version) => {
            let version = version.as_u64().unwrap_or_default();
            if version > u64::from(PORTABLE_SCHEMA_VERSION) {
                return Diagnosis {
                    model: None,
                    issues: vec![ModelIssue::error(
                        "unsupported_schema_version",
                        "schema_version",
                        format!(
                            "model schema version {version} is newer than the supported maximum \
                             {PORTABLE_SCHEMA_VERSION}; upgrade StericX to read it"
                        ),
                    )],
                };
            }
        }
    }

    // Fixed-width numeric arrays are checked here, by name, because the
    // decoder would only report a length or type mismatch against the whole
    // document without saying which field carried it.
    let mut issues = Vec::new();
    check_fixed_width_arrays(object, &mut issues);
    if !issues.is_empty() {
        return Diagnosis {
            model: None,
            issues,
        };
    }

    let fit: ScientificFitReport = match serde_json::from_str(text) {
        Ok(fit) => fit,
        Err(error) => {
            return Diagnosis {
                model: None,
                issues: vec![field_issue(&error)],
            };
        }
    };
    let sections: PortableSections = match serde_json::from_str(text) {
        Ok(sections) => sections,
        Err(error) => {
            return Diagnosis {
                model: None,
                issues: vec![field_issue(&error)],
            };
        }
    };
    let model = PortableModel {
        fit,
        inference: sections.inference,
        provenance: sections.provenance,
        created: sections.created,
    };
    let issues = model.issues();
    Diagnosis {
        model: Some(model),
        issues,
    }
}

/// Checks the arrays whose width is fixed by the feature space.
///
/// Reports the field, the expected width, and any element that is not a finite
/// number, so a truncated or null-bearing array is named rather than surfacing
/// as a decoder complaint about the document.
fn check_fixed_width_arrays(
    object: &serde_json::Map<String, serde_json::Value>,
    issues: &mut Vec<ModelIssue>,
) {
    for field in ["weights", "standardized_means", "standardized_scales"] {
        let Some(value) = object.get(field) else {
            continue;
        };
        let Some(array) = value.as_array() else {
            issues.push(ModelIssue::error(
                "invalid_field",
                field,
                format!("`{field}` must be an array of {MODEL_FEATURE_COUNT} numbers"),
            ));
            continue;
        };
        if array.len() != MODEL_FEATURE_COUNT {
            issues.push(ModelIssue::error(
                "dimension_mismatch",
                field,
                format!(
                    "`{field}` has {} entries but the feature space has {MODEL_FEATURE_COUNT} \
                     columns",
                    array.len()
                ),
            ));
        }
        for (index, entry) in array.iter().enumerate() {
            match entry.as_f64() {
                Some(number) if number.is_finite() => {}
                Some(number) => issues.push(ModelIssue::error(
                    "non_finite_value",
                    format!("{field}[{index}]"),
                    format!("`{field}[{index}]` is {number}, which cannot be used for inference"),
                )),
                None => issues.push(ModelIssue::error(
                    "non_numeric_value",
                    format!("{field}[{index}]"),
                    format!("`{field}[{index}]` is {entry}, expected a number"),
                )),
            }
        }
    }
    if let Some(names) = object
        .get("feature_names")
        .and_then(|value| value.as_array())
        && names.len() != MODEL_FEATURE_COUNT
    {
        issues.push(ModelIssue::error(
            "dimension_mismatch",
            "feature_names",
            format!(
                "`feature_names` has {} entries but the feature space has \
                 {MODEL_FEATURE_COUNT} columns",
                names.len()
            ),
        ));
    }
}

/// Turns a decoder failure into a finding that names the offending field.
fn field_issue(error: &serde_json::Error) -> ModelIssue {
    let text = error.to_string();
    let field = text
        .split_once('`')
        .and_then(|(_, rest)| rest.split_once('`'))
        .map(|(name, _)| name.to_owned());
    match (field, text.starts_with("missing field")) {
        (Some(field), true) => ModelIssue::error(
            "missing_field",
            field.clone(),
            format!("required field `{field}` is absent"),
        ),
        (Some(field), false) => ModelIssue::error(
            "invalid_field",
            field.clone(),
            format!("field `{field}` could not be read: {text}"),
        ),
        _ => ModelIssue::error(
            "invalid_document",
            "$",
            format!("model is unreadable: {text}"),
        ),
    }
}

/// Minimal view used to read `schema_version` before full deserialization.
#[derive(Deserialize)]
struct SchemaProbe {
    schema_version: u32,
}

/// Records a disagreement between the `inference` section and the fit report.
///
/// The two copies are generated from one value, so this tolerates
/// re-serialization only. Any real edit to either shows up here.
fn push_cross_check(
    issues: &mut Vec<ModelIssue>,
    location: &str,
    actual: f64,
    expected: f64,
    label: &str,
) {
    let slack = CROSS_CHECK_TOLERANCE * expected.abs().max(1.0);
    if (actual - expected).abs() > slack {
        issues.push(ModelIssue::error(
            "inference_disagrees_with_fit",
            location,
            format!("{label} is {actual} but the fit report says {expected}"),
        ));
    }
}

/// Formats a system time as an RFC 3339 UTC timestamp.
fn rfc3339_utc(time: SystemTime) -> String {
    let seconds = time
        .duration_since(UNIX_EPOCH)
        .map(|elapsed| elapsed.as_secs())
        .unwrap_or_default();
    let (year, month, day) = civil_from_days((seconds / 86_400) as i64);
    let time_of_day = seconds % 86_400;
    format!(
        "{year:04}-{month:02}-{day:02}T{:02}:{:02}:{:02}Z",
        time_of_day / 3_600,
        (time_of_day % 3_600) / 60,
        time_of_day % 60
    )
}

/// Converts days since the Unix epoch into a proleptic Gregorian date.
fn civil_from_days(days: i64) -> (i64, u32, u32) {
    let shifted = days + 719_468;
    let era = if shifted >= 0 {
        shifted
    } else {
        shifted - 146_096
    } / 146_097;
    let day_of_era = shifted - era * 146_097;
    let year_of_era =
        (day_of_era - day_of_era / 1_460 + day_of_era / 36_524 - day_of_era / 146_096) / 365;
    let year = year_of_era + era * 400;
    let day_of_year = day_of_era - (365 * year_of_era + year_of_era / 4 - year_of_era / 100);
    let month_position = (5 * day_of_year + 2) / 153;
    let day = (day_of_year - (153 * month_position + 2) / 5 + 1) as u32;
    let month = if month_position < 10 {
        month_position + 3
    } else {
        month_position - 9
    } as u32;
    (if month <= 2 { year + 1 } else { year }, month, day)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn formats_known_utc_timestamps() {
        assert_eq!(
            rfc3339_utc(UNIX_EPOCH),
            "1970-01-01T00:00:00Z",
            "the epoch itself"
        );
        assert_eq!(
            rfc3339_utc(UNIX_EPOCH + std::time::Duration::from_secs(1_774_000_000)),
            "2026-03-20T09:46:40Z"
        );
        // 2024 was a leap year; day 60 must be 29 February.
        assert_eq!(
            rfc3339_utc(UNIX_EPOCH + std::time::Duration::from_secs(1_709_164_800)),
            "2024-02-29T00:00:00Z"
        );
    }
}
