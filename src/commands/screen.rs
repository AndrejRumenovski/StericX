//! `stericx screen`: rank a ligand library with a fitted reaction model,
//! reporting predicted performance, coefficient uncertainty, and
//! applicability-domain warnings.
//!
//! The model decides what the library must provide. StericX's regression space
//! mixes geometry (`L`, `B1`, `B5`) with donor electronics (`nbo_charge`,
//! `ir_frequency`) and their interactions, so a model that selected an
//! electronic term cannot be screened from geometry alone. Rather than invent
//! the missing quantity, `screen` inspects the fitted weights, works out which
//! inputs actually carry a nonzero coefficient, and refuses to run when the
//! library cannot supply one.

use crate::cli::{DescriptorFormat, SterimolAxis};
use crate::descriptors::descriptors_for_file;
use rayon::prelude::*;
use serde::Serialize;
use std::collections::HashMap;
use std::error::Error;
use std::path::{Path, PathBuf};
use steric_x::model::{
    FeatureTransform, InferenceSpec, MODEL_FEATURE_NAMES, Optimization, PortableModel,
    assess_applicability, expand_features,
};
use steric_x::{BuriedVolumeConfig, EyringKineticLink, PackedReactionRecord, ScientificFitReport};

/// The only feature construction this build knows how to build descriptors for.
const KNOWN_FEATURE_SPACE: &str = "stericx.physical_organic.v1";

/// Model feature positions, mirroring `MODEL_FEATURE_NAMES`.
const F_L: usize = 1;
const F_B1: usize = 2;
const F_B5: usize = 3;
const F_NBO: usize = 4;
const F_B1_NBO: usize = 5;
const F_B5_NBO: usize = 6;
const F_IR: usize = 7;

/// Which raw inputs a fitted model actually needs, derived from its weights.
#[derive(Clone, Copy, Debug, Default)]
struct RequiredInputs {
    l: bool,
    b1: bool,
    b5: bool,
    nbo_charge: bool,
    ir_frequency: bool,
}

impl RequiredInputs {
    /// A feature with a zero coefficient cannot influence the prediction, so
    /// its inputs are genuinely not required — this keeps a geometry-only model
    /// screenable from a geometry-only library.
    fn from_weights(weights: &[f32; 8]) -> Self {
        let used = |index: usize| weights[index] != 0.0;
        Self {
            l: used(F_L),
            b1: used(F_B1) || used(F_B1_NBO),
            b5: used(F_B5) || used(F_B5_NBO),
            nbo_charge: used(F_NBO) || used(F_B1_NBO) || used(F_B5_NBO),
            ir_frequency: used(F_IR),
        }
    }

    /// Derives the required descriptors from the feature space the model
    /// itself records, mapping each selected term through its stored
    /// transformation rather than assuming this build's layout.
    ///
    /// An unrecognised feature space or descriptor name is refused: screening a
    /// model against a feature construction StericX does not implement would
    /// silently substitute the wrong quantity.
    fn from_inference(spec: &InferenceSpec) -> Result<Self, String> {
        if spec.feature_space.definition != KNOWN_FEATURE_SPACE {
            return Err(format!(
                "model declares feature space `{}`, but this build implements `{}`; \
                 StericX will not guess how to build its descriptors",
                spec.feature_space.definition, KNOWN_FEATURE_SPACE
            ));
        }
        let mut required = Self::default();
        for term in &spec.terms {
            if term.coefficient == 0.0 {
                continue;
            }
            let transform = spec
                .feature_space
                .transformations
                .get(term.feature_index)
                .ok_or_else(|| {
                    format!(
                        "model term `{}` refers to column {} but the feature space defines {} \
                         columns",
                        term.feature_name,
                        term.feature_index,
                        spec.feature_space.transformations.len()
                    )
                })?;
            let descriptors: Vec<&str> = match transform {
                FeatureTransform::Constant => Vec::new(),
                FeatureTransform::Descriptor { descriptor } => vec![descriptor.as_str()],
                FeatureTransform::Interaction { factors } => {
                    factors.iter().map(String::as_str).collect()
                }
            };
            for descriptor in descriptors {
                match descriptor {
                    "sterimol_l" => required.l = true,
                    "sterimol_b1" => required.b1 = true,
                    "sterimol_b5" => required.b5 = true,
                    "nbo_charge" => required.nbo_charge = true,
                    "ir_frequency" => required.ir_frequency = true,
                    other => {
                        return Err(format!(
                            "model term `{}` needs descriptor `{other}`, which this build cannot \
                             calculate",
                            term.feature_name
                        ));
                    }
                }
            }
        }
        Ok(required)
    }

    fn missing_from(&self, available: &Available) -> Vec<&'static str> {
        let mut missing = Vec::new();
        if self.l && !available.l {
            missing.push("sterimol_l");
        }
        if self.b1 && !available.b1 {
            missing.push("sterimol_b1");
        }
        if self.b5 && !available.b5 {
            missing.push("sterimol_b5");
        }
        if self.nbo_charge && !available.nbo_charge {
            missing.push("nbo_charge");
        }
        if self.ir_frequency && !available.ir_frequency {
            missing.push("ir_frequency");
        }
        missing
    }
}

/// A loaded library: its members, what the source can supply, and the members
/// that could not be loaded at all.
type LoadedLibrary = (Vec<Candidate>, Available, Vec<Exclusion>);

/// Which inputs a library source actually carries.
#[derive(Clone, Copy, Debug, Default)]
struct Available {
    l: bool,
    b1: bool,
    b5: bool,
    nbo_charge: bool,
    ir_frequency: bool,
}

/// One library member awaiting prediction.
#[derive(Clone, Debug, Default)]
struct Candidate {
    /// Position in the library as read, the final tiebreak when a prediction
    /// and an identifier are both shared.
    library_index: usize,
    label: String,
    /// Human-readable name when the library supplies one, e.g. a SMILES.
    name: Option<String>,
    /// Geometry referenced by the row, used to fill Sterimol terms a CSV lacks.
    geometry: Option<PathBuf>,
    l: Option<f32>,
    b1: Option<f32>,
    b5: Option<f32>,
    nbo_charge: Option<f32>,
    ir_frequency: Option<f32>,
}

impl Candidate {
    fn to_record(&self) -> PackedReactionRecord {
        PackedReactionRecord {
            l: self.l.unwrap_or(0.0),
            b1: self.b1.unwrap_or(0.0),
            b5: self.b5.unwrap_or(0.0),
            nbo_charge: self.nbo_charge.unwrap_or(0.0),
            ir_freq: self.ir_frequency.unwrap_or(0.0),
            ..PackedReactionRecord::default()
        }
    }

    /// Names the required descriptors this ligand does not carry.
    fn missing_values(&self, required: RequiredInputs) -> Vec<String> {
        [
            (required.l, self.l.is_some(), "sterimol_l"),
            (required.b1, self.b1.is_some(), "sterimol_b1"),
            (required.b5, self.b5.is_some(), "sterimol_b5"),
            (required.nbo_charge, self.nbo_charge.is_some(), "nbo_charge"),
            (
                required.ir_frequency,
                self.ir_frequency.is_some(),
                "ir_frequency",
            ),
        ]
        .into_iter()
        .filter(|(needed, present, _)| *needed && !*present)
        .map(|(_, _, name)| name.to_owned())
        .collect()
    }

    /// The descriptor values the model actually consumed, in model order.
    fn descriptor_values(&self, required: RequiredInputs) -> Vec<DescriptorValue> {
        [
            (required.l, self.l, "sterimol_l"),
            (required.b1, self.b1, "sterimol_b1"),
            (required.b5, self.b5, "sterimol_b5"),
            (required.nbo_charge, self.nbo_charge, "nbo_charge"),
            (required.ir_frequency, self.ir_frequency, "ir_frequency"),
        ]
        .into_iter()
        .filter(|(needed, _, _)| *needed)
        .filter_map(|(_, value, name)| {
            value.map(|value| DescriptorValue {
                name: name.to_owned(),
                value: f64::from(value),
            })
        })
        .collect()
    }

    fn has(&self, required: RequiredInputs) -> bool {
        (!required.l || self.l.is_some())
            && (!required.b1 || self.b1.is_some())
            && (!required.b5 || self.b5.is_some())
            && (!required.nbo_charge || self.nbo_charge.is_some())
            && (!required.ir_frequency || self.ir_frequency.is_some())
    }
}

/// One feature that fell outside the training range, and by how much.
#[derive(Clone, Debug, Serialize)]
struct DomainExceedance {
    feature: String,
    value: f64,
    training_minimum: f64,
    training_maximum: f64,
    /// How far outside, as a fraction of the training range width.
    exceedance_fraction: f64,
}

/// How the screened candidates were ordered.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RankingOrder {
    /// Best is the largest prediction.
    Descending,
    /// Best is the smallest prediction.
    Ascending,
    /// Best is the largest absolute prediction, whatever the sign.
    Magnitude,
}

impl RankingOrder {
    fn label(self) -> &'static str {
        match self {
            Self::Descending => "descending",
            Self::Ascending => "ascending",
            Self::Magnitude => "magnitude_descending",
        }
    }

    /// Orders two predictions best-first.
    ///
    /// `total_cmp` gives a total order over every f64 including NaN, so the
    /// comparison never depends on operand order.
    fn compare(self, left: f64, right: f64) -> std::cmp::Ordering {
        match self {
            Self::Descending => right.total_cmp(&left),
            Self::Ascending => left.total_cmp(&right),
            Self::Magnitude => right.abs().total_cmp(&left.abs()),
        }
    }
}

/// Resolves the ranking direction from the model, honouring an explicit
/// override.
///
/// A model that does not state a direction is not ranked on a guess: the caller
/// is asked to say which way is better, because for a signed selectivity
/// response neither direction is universally right.
fn resolve_order(
    optimization: Option<Optimization>,
    ascending: bool,
    descending: bool,
) -> Result<(RankingOrder, bool), Box<dyn Error>> {
    let from_model = match optimization {
        Some(Optimization::Maximize) => Some(RankingOrder::Descending),
        Some(Optimization::Minimize) => Some(RankingOrder::Ascending),
        Some(Optimization::MaximizeMagnitude) => Some(RankingOrder::Magnitude),
        Some(Optimization::Unspecified) | None => None,
    };
    let requested = match (ascending, descending) {
        (true, false) => Some(RankingOrder::Ascending),
        (false, true) => Some(RankingOrder::Descending),
        _ => None,
    };
    match (from_model, requested) {
        (_, Some(requested)) => {
            let overridden = from_model.is_some_and(|model| model != requested);
            Ok((requested, overridden))
        }
        (Some(model), None) => Ok((model, false)),
        (None, None) => Err(
            "this model does not record which direction of its response is \
             better, so `screen` will not rank on a guess. Re-fit with \
             `--optimize maximize|minimize|maximize-magnitude`, or state the \
             order for this run with `--ascending` or `--descending`."
                .into(),
        ),
    }
}

/// One descriptor the model actually consumed for a ligand.
#[derive(Clone, Debug, Serialize)]
struct DescriptorValue {
    name: String,
    value: f64,
}

/// A library member that could not be screened, and why.
///
/// Every exclusion is recorded rather than dropped: a screen that silently
/// shrank its candidate set would misrepresent the search that was performed.
#[derive(Clone, Debug, Serialize)]
struct Exclusion {
    ligand: String,
    /// Stable slug, e.g. `missing_descriptors`.
    reason: String,
    /// Descriptors the model needs that this ligand does not supply.
    missing_descriptors: Vec<String>,
    detail: String,
}

#[derive(Clone, Debug, Serialize)]
struct ScreenHit {
    rank: usize,
    #[serde(skip)]
    library_index: usize,
    ligand: String,
    /// Present only when the library supplies a name distinct from the id.
    ligand_name: Option<String>,
    /// The raw model output, always retained even when a transformed value is
    /// also reported.
    predicted_ddg_kcal_mol: f64,
    /// The descriptor values this prediction consumed, in model order.
    descriptors: Vec<DescriptorValue>,
    /// Conservative bounds from the bootstrap coefficient intervals. NOT an
    /// OLS prediction interval — see the note emitted with the report.
    coefficient_band_low: Option<f64>,
    coefficient_band_high: Option<f64>,
    /// True 95 % prediction interval `ŷ ± t·s·√(1+h)`, when the model carries
    /// the training geometry needed to compute one.
    prediction_interval_low: Option<f64>,
    prediction_interval_high: Option<f64>,
    /// Leverage of this ligand against the training design.
    leverage: Option<f64>,
    /// `leverage / warning_leverage`; above 1.0 the ligand is an extrapolation.
    leverage_ratio: Option<f64>,
    /// Graded verdict combining the range check and the leverage check.
    trust: String,
    /// Where the candidate sits relative to the training set: `interpolation`,
    /// `sparse_interpolation`, `extrapolation`, or `unknown`.
    domain_verdict: String,
    /// Standardized distance to the closest training observation.
    nearest_training_distance: Option<f64>,
    /// The training set's own nearest-neighbour spacing, the boundary above.
    nearest_training_threshold: Option<f64>,
    /// `distance / threshold`; above one is outside the sampled region.
    nearest_training_ratio: Option<f64>,
    /// Mahalanobis distance, when the training covariance supports one.
    mahalanobis_distance: Option<f64>,
    /// Largest range overshoot as a fraction of the training range width.
    maximum_extrapolation: f64,
    /// Signed by the ΔΔG‡ convention: positive favours R, negative the same
    /// excess of the opposite enantiomer.
    predicted_ee_percent: f64,
    applicability: String,
    outside_domain: Vec<DomainExceedance>,
}

#[derive(Clone, Debug, Serialize)]
struct ScreenReport {
    model_path: String,
    model: String,
    selected_features: Vec<String>,
    required_inputs: Vec<String>,
    training_count: usize,
    training_r2: Option<f64>,
    /// Residual scatter of the fit: the irreducible spread these predictions
    /// inherit, reported separately from the coefficient band.
    training_rmse_kcal_mol: f64,
    temperature_k: f32,
    library_size: usize,
    /// Candidates that produced a prediction, before `--top` was applied.
    screened: usize,
    /// Candidates actually listed after `--top`.
    returned: usize,
    /// The order the table is in.
    ranking_order: String,
    /// The direction the model records, when it records one.
    model_optimization: String,
    /// True when the caller asked for an order the model disagrees with.
    ranking_overridden: bool,
    skipped: usize,
    /// Why each excluded candidate was excluded.
    excluded: Vec<Exclusion>,
    /// Exclusion counts by reason, largest first.
    exclusion_summary: Vec<(String, usize)>,
    inside_domain: usize,
    /// Warning leverage h* = 3p/n for this fit, when available.
    warning_leverage: Option<f64>,
    /// Count of ligands whose leverage exceeds h*.
    high_leverage: usize,
    /// How many ligands each trust grade covers.
    trust_summary: Vec<(String, usize)>,
    /// How many ligands each applicability verdict covers.
    domain_summary: Vec<(String, usize)>,
    /// How the neighbour boundary was derived, carried from the model.
    neighbor_rule: Option<String>,
    uncertainty_note: String,
    hits: Vec<ScreenHit>,
}

pub(crate) struct ScreenArgs<'a> {
    pub(crate) model: &'a Path,
    pub(crate) library: &'a Path,
    pub(crate) top: Option<usize>,
    pub(crate) temperature: f32,
    pub(crate) inside_domain_only: bool,
    pub(crate) ascending: bool,
    pub(crate) descending: bool,
    pub(crate) donor_element: &'a str,
    pub(crate) sterimol_axis: SterimolAxis,
    pub(crate) format: DescriptorFormat,
    pub(crate) config: BuriedVolumeConfig,
}

/// Loads a model through the portable-format reader, so both the legacy
/// artifact and a schema-2 document are accepted and validated first.
fn load_model(path: &Path) -> Result<PortableModel, Box<dyn Error>> {
    let contents = std::fs::read_to_string(path)
        .map_err(|error| format!("could not read model {}: {error}", path.display()))?;
    PortableModel::from_json(&contents).map_err(|error| {
        format!(
            "{} is not a usable StericX model: {error}. Run `stericx model validate` for the \
             full list of problems.",
            path.display()
        )
        .into()
    })
}

/// Column aliases accepted for each input, so both a descriptors CSV and a raw
/// reaction CSV work as a library without conversion.
const COLUMN_ALIASES: &[(&str, &[&str])] = &[
    (
        "label",
        &[
            "ligand",
            "Ligand",
            "Reaction_ID",
            "reaction_id",
            "Source_ID",
            "source_id",
            "id",
            "file",
            "Ligand_XYZ_Path",
        ],
    ),
    (
        "name",
        &[
            "name",
            "ligand_name",
            "Ligand_Name",
            "Ligand_SMILES",
            "ligand_smiles",
            "smiles",
            "SMILES",
        ],
    ),
    (
        "geometry",
        &[
            "Ligand_XYZ_Path",
            "ligand_xyz_path",
            "file",
            "xyz_path",
            "path",
        ],
    ),
    ("l", &["sterimol_l", "L_boltz", "l"]),
    ("b1", &["sterimol_b1", "B1_boltz", "b1"]),
    ("b5", &["sterimol_b5", "B5_boltz", "b5"]),
    ("nbo_charge", &["nbo_charge", "NBO_Charge"]),
    (
        "ir_frequency",
        &["ir_frequency", "IR_Frequency", "ir_freq", "IR_Freq"],
    ),
];

fn header_index(headers: &csv::StringRecord, key: &str) -> Option<usize> {
    let aliases = COLUMN_ALIASES
        .iter()
        .find(|(name, _)| *name == key)
        .map(|(_, aliases)| *aliases)?;
    headers.iter().position(|header| {
        aliases
            .iter()
            .any(|alias| alias.eq_ignore_ascii_case(header.trim()))
    })
}

fn load_library(
    library: &Path,
    donor_element: &str,
    sterimol_axis: SterimolAxis,
    config: BuriedVolumeConfig,
) -> Result<LoadedLibrary, Box<dyn Error>> {
    let mut excluded = Vec::new();
    if library.is_dir() {
        let mut paths = Vec::new();
        collect_coordinate_files(library, &mut paths)?;
        paths.sort();
        if paths.is_empty() {
            return Err(format!(
                "no .xyz/.sdf/.mol geometries found below {}",
                library.display()
            )
            .into());
        }
        let featurized = paths
            .par_iter()
            .map(|path| {
                match descriptors_for_file(path, donor_element, None, sterimol_axis, config) {
                    Ok(result) => Ok(Candidate {
                        library_index: 0,
                        label: result.file.clone(),
                        l: Some(result.sterimol_l),
                        b1: Some(result.sterimol_b1),
                        b5: Some(result.sterimol_b5),
                        ..Candidate::default()
                    }),
                    Err(message) => Err(Exclusion {
                        ligand: path.display().to_string(),
                        reason: "featurization_failed".to_owned(),
                        missing_descriptors: Vec::new(),
                        detail: message,
                    }),
                }
            })
            .collect::<Vec<_>>();
        let mut candidates: Vec<Candidate> = Vec::new();
        for outcome in featurized {
            match outcome {
                Ok(mut candidate) => {
                    // `paths` is sorted, so this is a stable position.
                    candidate.library_index = candidates.len();
                    candidates.push(candidate);
                }
                Err(exclusion) => excluded.push(exclusion),
            }
        }
        if candidates.is_empty() {
            return Err("no library geometry could be featurized".into());
        }
        // A geometry directory can never supply donor electronics.
        return Ok((
            candidates,
            Available {
                l: true,
                b1: true,
                b5: true,
                nbo_charge: false,
                ir_frequency: false,
            },
            excluded,
        ));
    }

    if !library.is_file() {
        return Err(format!("library does not exist: {}", library.display()).into());
    }
    let mut reader = csv::ReaderBuilder::new()
        .trim(csv::Trim::All)
        .from_path(library)?;
    let headers = reader.headers()?.clone();
    let columns = [
        "label",
        "name",
        "geometry",
        "l",
        "b1",
        "b5",
        "nbo_charge",
        "ir_frequency",
    ]
    .iter()
    .map(|key| (*key, header_index(&headers, key)))
    .collect::<HashMap<_, _>>();
    let available = Available {
        l: columns["l"].is_some(),
        b1: columns["b1"].is_some(),
        b5: columns["b5"].is_some(),
        nbo_charge: columns["nbo_charge"].is_some(),
        ir_frequency: columns["ir_frequency"].is_some(),
    };

    let number = |record: &csv::StringRecord, index: Option<usize>| -> Option<f32> {
        index
            .and_then(|index| record.get(index))
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .and_then(|value| value.parse::<f32>().ok())
            .filter(|value| value.is_finite())
    };

    let mut candidates = Vec::new();
    for (row, record) in reader.records().enumerate() {
        let record = record.map_err(|error| {
            format!(
                "{} row {} could not be parsed: {error}",
                library.display(),
                row + 2
            )
        })?;
        let label = columns["label"]
            .and_then(|index| record.get(index))
            .map(|value| value.trim().to_owned())
            .filter(|value| !value.is_empty())
            .unwrap_or_else(|| format!("row_{}", row + 2));
        let name = columns["name"]
            .and_then(|index| record.get(index))
            .map(|value| value.trim().to_owned())
            .filter(|value| !value.is_empty() && *value != label);
        candidates.push(Candidate {
            library_index: row,
            label,
            name,
            geometry: columns["geometry"]
                .and_then(|index| record.get(index))
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(PathBuf::from),
            l: number(&record, columns["l"]),
            b1: number(&record, columns["b1"]),
            b5: number(&record, columns["b5"]),
            nbo_charge: number(&record, columns["nbo_charge"]),
            ir_frequency: number(&record, columns["ir_frequency"]),
        });
    }
    if candidates.is_empty() {
        return Err(format!("library CSV contains no rows: {}", library.display()).into());
    }

    // A reaction CSV carries the electronics but not the Sterimol terms. When
    // the rows point at geometries, featurize them so a model mixing steric and
    // electronic terms can be screened from the CSV the user already has.
    let mut available = available;
    if !(available.l && available.b1 && available.b5) && columns["geometry"].is_some() {
        let base = library.parent().unwrap_or(Path::new("."));
        let featurized = candidates
            .par_iter()
            .map(|candidate| {
                let path = candidate.geometry.as_ref()?;
                let resolved = if path.is_file() {
                    path.clone()
                } else if base.join(path).is_file() {
                    base.join(path)
                } else {
                    // Left unfeaturized; the screening loop reports the ligand
                    // as missing whichever descriptors the model needs.
                    return None;
                };
                match descriptors_for_file(&resolved, donor_element, None, sterimol_axis, config) {
                    Ok(result) => Some((result.sterimol_l, result.sterimol_b1, result.sterimol_b5)),
                    Err(_) => None,
                }
            })
            .collect::<Vec<_>>();
        if featurized.iter().any(Option::is_some) {
            for (candidate, sterimol) in candidates.iter_mut().zip(featurized) {
                if let Some((l, b1, b5)) = sterimol {
                    candidate.l = candidate.l.or(Some(l));
                    candidate.b1 = candidate.b1.or(Some(b1));
                    candidate.b5 = candidate.b5.or(Some(b5));
                }
            }
            available.l = true;
            available.b1 = true;
            available.b5 = true;
        }
    }
    Ok((candidates, available, excluded))
}

fn collect_coordinate_files(
    directory: &Path,
    found: &mut Vec<PathBuf>,
) -> Result<(), Box<dyn Error>> {
    for entry in std::fs::read_dir(directory)? {
        let path = entry?.path();
        if path.is_dir() {
            collect_coordinate_files(&path, found)?;
        } else if path
            .extension()
            .and_then(|extension| extension.to_str())
            .is_some_and(|extension| {
                matches!(
                    extension.to_ascii_lowercase().as_str(),
                    "xyz" | "sdf" | "mol"
                )
            })
        {
            found.push(path);
        }
    }
    Ok(())
}

/// Propagate the bootstrap coefficient intervals through one feature vector.
///
/// This is deliberately interval arithmetic over each coefficient's 95 % band,
/// which ignores the correlation between coefficients and so is *conservative*
/// (wider than a joint region). It is a parameter-uncertainty band, not a
/// prediction interval: `model.json` does not carry the training design matrix
/// an OLS prediction interval would need.
fn coefficient_band(report: &ScientificFitReport, features: &[f32; 8]) -> Option<(f64, f64)> {
    if report.coefficient_intervals.is_empty() {
        return None;
    }
    let by_name = report
        .coefficient_intervals
        .iter()
        .map(|interval| (interval.feature.as_str(), interval))
        .collect::<HashMap<_, _>>();
    let mut low = 0.0_f64;
    let mut high = 0.0_f64;
    for (index, name) in MODEL_FEATURE_NAMES.iter().enumerate() {
        if report.weights[index] == 0.0 {
            continue;
        }
        let value = f64::from(features[index]);
        let interval = by_name.get(name)?;
        let a = interval.lower_95 * value;
        let b = interval.upper_95 * value;
        low += a.min(b);
        high += a.max(b);
    }
    Some((low, high))
}

fn domain_exceedances(report: &ScientificFitReport, features: &[f32; 8]) -> Vec<DomainExceedance> {
    report
        .selected_feature_indices
        .iter()
        .zip(&report.applicability_domain)
        .filter_map(|(column, domain)| {
            let value = f64::from(features[*column]);
            if value >= domain.minimum && value <= domain.maximum {
                return None;
            }
            let width = domain.maximum - domain.minimum;
            let distance = if value < domain.minimum {
                domain.minimum - value
            } else {
                value - domain.maximum
            };
            Some(DomainExceedance {
                feature: domain.feature.clone(),
                value,
                training_minimum: domain.minimum,
                training_maximum: domain.maximum,
                exceedance_fraction: if width > f64::EPSILON {
                    distance / width
                } else {
                    f64::INFINITY
                },
            })
        })
        .collect()
}

/// Grade how much a prediction deserves trust, from two independent signals:
/// whether every selected feature lies inside the training range (a 1-D box
/// check) and whether the ligand's leverage stays under `h* = 3p/n` (a check in
/// the correlated geometry the fit actually defines). They can disagree — a
/// ligand can sit inside every 1-D range yet still be far from the training
/// cloud — so both are reported and the worse one governs.
fn trust_grade(inside_range: bool, leverage_ratio: Option<f64>) -> String {
    match (inside_range, leverage_ratio) {
        (_, None) => {
            if inside_range {
                "range_only:inside".to_owned()
            } else {
                "range_only:outside".to_owned()
            }
        }
        (true, Some(ratio)) if ratio <= 1.0 => "reliable".to_owned(),
        (true, Some(_)) => "caution:high_leverage".to_owned(),
        (false, Some(ratio)) if ratio <= 1.0 => "caution:outside_range".to_owned(),
        (false, Some(_)) => "do_not_trust:extrapolation".to_owned(),
    }
}

pub(crate) fn screen_command(args: ScreenArgs<'_>) -> Result<(), Box<dyn Error>> {
    if !args.temperature.is_finite() || args.temperature <= 0.0 {
        return Err("--temperature must be a positive finite temperature".into());
    }
    let model = load_model(args.model)?;
    let report = &model.fit;
    // A portable model states how its features are built; use that mapping in
    // preference to inferring one from the weight vector.
    let required = match model.inference.as_ref() {
        Some(spec) => RequiredInputs::from_inference(spec)?,
        None => RequiredInputs::from_weights(&report.weights),
    };

    let (candidates, available, mut excluded) = load_library(
        args.library,
        args.donor_element,
        args.sterimol_axis,
        args.config,
    )?;
    let library_size = candidates.len();

    let missing = required.missing_from(&available);
    if !missing.is_empty() {
        return Err(format!(
            "model `{}` uses {} but the library does not provide {}. \
             StericX will not guess a missing input. Sterimol terms come from a \
             geometry directory or a CSV with a geometry-path column; donor \
             electronics (nbo_charge, ir_frequency) must be supplied as CSV \
             columns — a reaction CSV with NBO_Charge / IR_Frequency alongside \
             Ligand_XYZ_Path provides both.",
            report.model,
            report.selected_features.join(", "),
            missing.join(" and ")
        )
        .into());
    }

    let model_optimization = model
        .inference
        .as_ref()
        .map_or(Optimization::Unspecified, |spec| spec.response.optimization);
    let (order, overridden) =
        resolve_order(Some(model_optimization), args.ascending, args.descending)?;

    let mut hits = Vec::new();
    for candidate in &candidates {
        if !candidate.has(required) {
            let missing = candidate.missing_values(required);
            excluded.push(Exclusion {
                ligand: candidate.label.clone(),
                reason: "missing_descriptors".to_owned(),
                detail: format!(
                    "the model needs {} but this ligand does not supply {}",
                    required_input_names(required).join(", "),
                    missing.join(", ")
                ),
                missing_descriptors: missing,
            });
            continue;
        }
        let record = candidate.to_record();
        let features = expand_features(&record);
        let predicted = report
            .weights
            .iter()
            .zip(&features)
            .map(|(weight, feature)| f64::from(*weight) * f64::from(*feature))
            .sum::<f64>();
        if !predicted.is_finite() {
            excluded.push(Exclusion {
                ligand: candidate.label.clone(),
                reason: "non_finite_prediction".to_owned(),
                missing_descriptors: Vec::new(),
                detail: "the model produced a non-finite prediction for this ligand".to_owned(),
            });
            continue;
        }
        let applicability = assess_applicability(
            report.training_geometry.as_ref(),
            &report.applicability_domain,
            &report.selected_feature_indices,
            &features,
        );
        let outside = domain_exceedances(report, &features);
        let band = coefficient_band(report, &features);
        let geometry = report.training_geometry.as_ref();
        let leverage = geometry.and_then(|geometry| geometry.leverage(&features));
        let warning = geometry.map(|geometry| geometry.warning_leverage);
        let leverage_ratio = match (leverage, warning) {
            (Some(leverage), Some(warning)) if warning > 0.0 => Some(leverage / warning),
            _ => None,
        };
        let interval =
            geometry.and_then(|geometry| geometry.prediction_interval(predicted, &features));
        let trust = trust_grade(outside.is_empty(), leverage_ratio);
        hits.push(ScreenHit {
            rank: 0,
            library_index: candidate.library_index,
            ligand: candidate.label.clone(),
            ligand_name: candidate.name.clone(),
            predicted_ddg_kcal_mol: predicted,
            descriptors: candidate.descriptor_values(required),
            coefficient_band_low: band.map(|(low, _)| low),
            coefficient_band_high: band.map(|(_, high)| high),
            prediction_interval_low: interval.map(|(low, _)| low),
            prediction_interval_high: interval.map(|(_, high)| high),
            leverage,
            leverage_ratio,
            trust,
            // `calculate_enantiomeric_excess` is an absolute value, so carry the
            // ΔΔG‡ sign through: under this project's convention a positive
            // ΔΔG‡ favours R, and a negative prediction means the same excess
            // of the *opposite* enantiomer.
            predicted_ee_percent: f64::from(EyringKineticLink::calculate_enantiomeric_excess(
                predicted as f32,
                args.temperature,
            )) * if predicted < 0.0 { -1.0 } else { 1.0 },
            domain_verdict: applicability.verdict.label().to_owned(),
            nearest_training_distance: applicability.nearest_training_distance,
            nearest_training_threshold: applicability.nearest_training_threshold,
            nearest_training_ratio: applicability.nearest_training_ratio,
            mahalanobis_distance: applicability.mahalanobis_distance,
            maximum_extrapolation: applicability.maximum_extrapolation,
            applicability: if outside.is_empty() {
                "inside_training_range".to_owned()
            } else {
                format!(
                    "outside:{}",
                    outside
                        .iter()
                        .map(|exceedance| exceedance.feature.as_str())
                        .collect::<Vec<_>>()
                        .join("|")
                )
            },
            outside_domain: outside,
        });
    }

    let inside_domain = hits
        .iter()
        .filter(|hit| hit.outside_domain.is_empty())
        .count();
    if args.inside_domain_only {
        hits.retain(|hit| hit.outside_domain.is_empty());
    }
    if hits.is_empty() {
        return Err("no library member could be screened with this model".into());
    }

    // Ties resolve by ligand identifier, then by the order the library was
    // read in, so a repeated run reproduces the same table even when several
    // ligands share a prediction or an identifier.
    hits.sort_by(|left, right| {
        order
            .compare(left.predicted_ddg_kcal_mol, right.predicted_ddg_kcal_mol)
            .then_with(|| left.ligand.cmp(&right.ligand))
            .then_with(|| left.library_index.cmp(&right.library_index))
    });
    // Every candidate is predicted before anything is dropped; `--top` selects
    // from the finished ranking rather than truncating the work.
    let screened = hits.len();
    if let Some(top) = args.top {
        hits.truncate(top);
    }
    let returned = hits.len();
    for (index, hit) in hits.iter_mut().enumerate() {
        hit.rank = index + 1;
    }

    let report = ScreenReport {
        model_path: args.model.display().to_string(),
        model: report.model.clone(),
        selected_features: report.selected_features.clone(),
        required_inputs: required_input_names(required),
        training_count: report.training_count,
        training_r2: report.training.r2,
        training_rmse_kcal_mol: report.training.rmse,
        temperature_k: args.temperature,
        library_size,
        screened,
        returned,
        ranking_order: order.label().to_owned(),
        model_optimization: model_optimization.label().to_owned(),
        ranking_overridden: overridden,
        skipped: excluded.len(),
        exclusion_summary: {
            let mut counts: Vec<(String, usize)> = Vec::new();
            for exclusion in &excluded {
                match counts
                    .iter_mut()
                    .find(|(reason, _)| *reason == exclusion.reason)
                {
                    Some((_, count)) => *count += 1,
                    None => counts.push((exclusion.reason.clone(), 1)),
                }
            }
            counts.sort_by(|left, right| right.1.cmp(&left.1).then_with(|| left.0.cmp(&right.0)));
            counts
        },
        excluded,
        inside_domain,
        warning_leverage: report
            .training_geometry
            .as_ref()
            .map(|g| g.warning_leverage),
        high_leverage: hits
            .iter()
            .filter(|hit| hit.leverage_ratio.is_some_and(|ratio| ratio > 1.0))
            .count(),
        trust_summary: {
            let mut counts: Vec<(String, usize)> = Vec::new();
            for hit in &hits {
                match counts.iter_mut().find(|(grade, _)| *grade == hit.trust) {
                    Some((_, count)) => *count += 1,
                    None => counts.push((hit.trust.clone(), 1)),
                }
            }
            counts.sort_by(|left, right| right.1.cmp(&left.1).then_with(|| left.0.cmp(&right.0)));
            counts
        },
        domain_summary: {
            let mut counts: Vec<(String, usize)> = Vec::new();
            for hit in &hits {
                match counts
                    .iter_mut()
                    .find(|(verdict, _)| *verdict == hit.domain_verdict)
                {
                    Some((_, count)) => *count += 1,
                    None => counts.push((hit.domain_verdict.clone(), 1)),
                }
            }
            counts.sort_by(|left, right| right.1.cmp(&left.1).then_with(|| left.0.cmp(&right.0)));
            counts
        },
        neighbor_rule: report
            .training_geometry
            .as_ref()
            .and_then(|geometry| geometry.neighbor_calibration.as_ref())
            .map(|calibration| {
                format!(
                    "{} (threshold {:.4}, mean {:.4}, sd {:.4})",
                    calibration.rule,
                    calibration.threshold,
                    calibration.mean,
                    calibration.standard_deviation
                )
            }),
        uncertainty_note: if report.training_geometry.is_some() {
            "prediction_interval is a 95% Student-t interval y ± t(0.975, n−p)·s·√(1+h): it \
             widens with leverage, so a ligand far from the training set is reported with a \
             correspondingly wider band. Leverage h is compared against the conventional \
             warning leverage h* = 3p/n. coefficient_band remains the conservative bootstrap \
             propagation and is the weaker of the two signals."
                .to_owned()
        } else {
            "this model predates the recorded training geometry, so no leverage or prediction \
             interval can be computed — the domain check falls back to a per-feature range test \
             only. Refit with `stericx fit` to enable leverage-based extrapolation detection. \
             coefficient_band propagates the bootstrap 95% coefficient intervals by interval \
             arithmetic; it ignores coefficient correlation (so it is conservative) and is NOT \
             an OLS prediction interval."
                .to_owned()
        },
        hits,
    };

    match args.format {
        DescriptorFormat::Text => print_screen_text(&report),
        DescriptorFormat::Json => println!("{}", serde_json::to_string_pretty(&report)?),
        DescriptorFormat::Csv => print_screen_csv(&report),
    }
    Ok(())
}

fn required_input_names(required: RequiredInputs) -> Vec<String> {
    let mut names = Vec::new();
    for (needed, name) in [
        (required.l, "sterimol_l"),
        (required.b1, "sterimol_b1"),
        (required.b5, "sterimol_b5"),
        (required.nbo_charge, "nbo_charge"),
        (required.ir_frequency, "ir_frequency"),
    ] {
        if needed {
            names.push(name.to_owned());
        }
    }
    names
}

fn print_screen_text(report: &ScreenReport) {
    println!("model          {} ({})", report.model, report.model_path);
    println!("selected       {}", report.selected_features.join(", "));
    println!("requires       {}", report.required_inputs.join(", "));
    println!(
        "trained on     {} ligands   training R² {}   RMSE {:.3} kcal/mol",
        report.training_count,
        report
            .training_r2
            .map_or_else(|| "unavailable".to_owned(), |r2| format!("{r2:.4}")),
        report.training_rmse_kcal_mol
    );
    println!(
        "library        {} members, {} screened, {} returned, {} skipped, {} inside the \
         training domain",
        report.library_size, report.screened, report.returned, report.skipped, report.inside_domain
    );
    println!(
        "ranking        {} (model records: {})",
        report.ranking_order, report.model_optimization
    );
    if report.ranking_overridden {
        println!(
            "               ⚠ this reverses the direction the model records as better; \
             the top of this table is the model\'s worst"
        );
    }
    println!("temperature    {:.2} K", report.temperature_k);
    if let Some(warning) = report.warning_leverage {
        println!(
            "leverage       warning h* = {warning:.3} (3p/n); {} ligand(s) above it",
            report.high_leverage
        );
    } else {
        println!("leverage       unavailable — model carries no training geometry");
    }
    if !report.trust_summary.is_empty() {
        println!(
            "trust          {}",
            report
                .trust_summary
                .iter()
                .map(|(grade, count)| format!("{grade} {count}"))
                .collect::<Vec<_>>()
                .join(", ")
        );
    }
    println!(
        "\n{:>4}  {:>9}  {:>19}  {:>7}  {:>7}  {:<26}  ligand",
        "rank", "pred ddG", "95% pred interval", "leverage", "pred ee", "trust"
    );
    for hit in &report.hits {
        let interval = match (hit.prediction_interval_low, hit.prediction_interval_high) {
            (Some(low), Some(high)) => format!("[{low:>7.2}, {high:>7.2}]"),
            _ => match (hit.coefficient_band_low, hit.coefficient_band_high) {
                // Fall back to the weaker signal, marked so it is never mistaken
                // for a real prediction interval.
                (Some(low), Some(high)) => format!("~[{low:>6.2}, {high:>6.2}]"),
                _ => "unavailable".to_owned(),
            },
        };
        let leverage = match (hit.leverage, hit.leverage_ratio) {
            (Some(leverage), Some(ratio)) => {
                format!("{leverage:.3}{}", if ratio > 1.0 { "!" } else { " " })
            }
            (Some(leverage), None) => format!("{leverage:.3} "),
            _ => "     — ".to_owned(),
        };
        println!(
            "{:>4}  {:>9.3}  {:>19}  {:>7}  {:>6.1}%  {:<26}  {}",
            hit.rank,
            hit.predicted_ddg_kcal_mol,
            interval,
            leverage,
            hit.predicted_ee_percent,
            hit.trust,
            match hit.ligand_name.as_deref() {
                Some(name) => format!("{} ({name})", hit.ligand),
                None => hit.ligand.clone(),
            }
        );
        if !hit.descriptors.is_empty() {
            println!("{:>60}{}", "", describe_descriptors(&hit.descriptors));
        }
        println!(
            "{:>60}domain: {}{}{}",
            "",
            hit.domain_verdict,
            match hit.nearest_training_distance {
                Some(distance) => format!("  nearest training point {distance:.3}"),
                None => String::new(),
            },
            match hit.mahalanobis_distance {
                Some(distance) => format!("  mahalanobis {distance:.3}"),
                None => String::new(),
            }
        );
    }
    let untrusted = report
        .hits
        .iter()
        .filter(|hit| hit.trust.starts_with("do_not_trust"))
        .collect::<Vec<_>>();
    if !untrusted.is_empty() {
        println!(
            "\n{} ligand(s) are outside the training range AND above the warning leverage.",
            untrusted.len()
        );
        println!("For these the model is extrapolating and its prediction should not be trusted:");
        for hit in untrusted {
            println!(
                "  {} — leverage {:.3} = {:.1}x h*",
                hit.ligand,
                hit.leverage.unwrap_or(f64::NAN),
                hit.leverage_ratio.unwrap_or(f64::NAN)
            );
        }
    }
    let flagged = report
        .hits
        .iter()
        .filter(|hit| !hit.outside_domain.is_empty())
        .collect::<Vec<_>>();
    if !flagged.is_empty() {
        println!("\napplicability-domain warnings (extrapolation — treat with care):");
        for hit in flagged {
            for exceedance in &hit.outside_domain {
                println!(
                    "  {}: {} = {:.4} is outside the training range [{:.4}, {:.4}] by {:.0}% of its width",
                    hit.ligand,
                    exceedance.feature,
                    exceedance.value,
                    exceedance.training_minimum,
                    exceedance.training_maximum,
                    exceedance.exceedance_fraction * 100.0
                );
            }
        }
    }
    if !report.domain_summary.is_empty() {
        println!(
            "\napplicability: {}",
            report
                .domain_summary
                .iter()
                .map(|(verdict, count)| format!("{verdict} {count}"))
                .collect::<Vec<_>>()
                .join(", ")
        );
        if let Some(rule) = report.neighbor_rule.as_deref() {
            println!("  neighbour boundary: {rule}");
        }
    }
    if !report.excluded.is_empty() {
        println!(
            "\n{} candidate(s) were excluded and not screened:",
            report.excluded.len()
        );
        for (reason, count) in &report.exclusion_summary {
            println!("  {reason}: {count}");
        }
        // Name a bounded sample so a large library stays readable while the
        // machine-readable formats keep the complete list.
        for exclusion in report.excluded.iter().take(EXCLUSION_PREVIEW) {
            println!("    {} — {}", exclusion.ligand, exclusion.detail);
        }
        if report.excluded.len() > EXCLUSION_PREVIEW {
            println!(
                "    ... and {} more; use --format json or csv for the full list",
                report.excluded.len() - EXCLUSION_PREVIEW
            );
        }
    }
    println!("\nnote: {}", report.uncertainty_note);
}

/// How many excluded candidates the terminal report names individually.
const EXCLUSION_PREVIEW: usize = 10;

/// Renders the descriptors a prediction consumed as `name=value` pairs.
fn describe_descriptors(descriptors: &[DescriptorValue]) -> String {
    descriptors
        .iter()
        .map(|descriptor| format!("{}={:.4}", descriptor.name, descriptor.value))
        .collect::<Vec<_>>()
        .join(" ")
}

/// Descriptors at full precision, for the machine-readable formats.
fn describe_descriptors_exact(descriptors: &[DescriptorValue]) -> String {
    descriptors
        .iter()
        .map(|descriptor| format!("{}={}", descriptor.name, descriptor.value))
        .collect::<Vec<_>>()
        .join(" ")
}

fn print_screen_csv(report: &ScreenReport) {
    println!(
        "rank,ligand,ligand_name,predicted_ddg_kcal_mol,predicted_ee_percent,descriptors,\
         domain_verdict,nearest_training_distance,nearest_training_threshold,\
         nearest_training_ratio,mahalanobis_distance,maximum_extrapolation,\
         prediction_interval_low,prediction_interval_high,\
         leverage,leverage_ratio,trust,coefficient_band_low,coefficient_band_high,\
         applicability,outside_features"
    );
    for hit in &report.hits {
        // Full round-trip precision: a screened prediction must reproduce the
        // engine value exactly when read back, not a rounded copy of it.
        let format_bound =
            |value: Option<f64>| value.map_or_else(String::new, |value| format!("{value}"));
        println!(
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            hit.rank,
            crate::descriptors::csv_field(&hit.ligand),
            crate::descriptors::csv_field(hit.ligand_name.as_deref().unwrap_or_default()),
            hit.predicted_ddg_kcal_mol,
            hit.predicted_ee_percent,
            crate::descriptors::csv_field(&describe_descriptors_exact(&hit.descriptors)),
            hit.domain_verdict,
            format_bound(hit.nearest_training_distance),
            format_bound(hit.nearest_training_threshold),
            format_bound(hit.nearest_training_ratio),
            format_bound(hit.mahalanobis_distance),
            hit.maximum_extrapolation,
            format_bound(hit.prediction_interval_low),
            format_bound(hit.prediction_interval_high),
            format_bound(hit.leverage),
            format_bound(hit.leverage_ratio),
            hit.trust,
            format_bound(hit.coefficient_band_low),
            format_bound(hit.coefficient_band_high),
            hit.applicability,
            crate::descriptors::csv_field(
                &hit.outside_domain
                    .iter()
                    .map(|exceedance| exceedance.feature.as_str())
                    .collect::<Vec<_>>()
                    .join("|")
            )
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn weights(pairs: &[(usize, f32)]) -> [f32; 8] {
        let mut weights = [0.0_f32; 8];
        for (index, value) in pairs {
            weights[*index] = *value;
        }
        weights
    }

    #[test]
    fn geometry_only_models_do_not_require_electronics() {
        let required = RequiredInputs::from_weights(&weights(&[(0, 1.0), (F_B5, 0.3)])); // 0 = intercept
        assert!(required.b5);
        assert!(!required.nbo_charge);
        assert!(!required.ir_frequency);
        assert!(!required.l);
    }

    #[test]
    fn interaction_terms_pull_in_both_of_their_inputs() {
        let required = RequiredInputs::from_weights(&weights(&[(F_B5_NBO, -0.11)]));
        assert!(required.b5, "B5_x_nbo_charge needs B5");
        assert!(required.nbo_charge, "B5_x_nbo_charge needs the charge");
        assert!(!required.b1);
    }

    #[test]
    fn missing_inputs_are_named_precisely() {
        let required = RequiredInputs::from_weights(&weights(&[(F_B5_NBO, -0.11)]));
        let geometry_only = Available {
            l: true,
            b1: true,
            b5: true,
            nbo_charge: false,
            ir_frequency: false,
        };
        assert_eq!(required.missing_from(&geometry_only), vec!["nbo_charge"]);
        let complete = Available {
            nbo_charge: true,
            ir_frequency: true,
            ..geometry_only
        };
        assert!(required.missing_from(&complete).is_empty());
    }

    #[test]
    fn candidate_reports_whether_it_satisfies_the_model() {
        let required = RequiredInputs::from_weights(&weights(&[(F_B5_NBO, -0.11)]));
        let without_charge = Candidate {
            label: "a".into(),
            b5: Some(8.0),
            ..Candidate::default()
        };
        assert!(!without_charge.has(required));
        let with_charge = Candidate {
            nbo_charge: Some(-0.3),
            ..without_charge.clone()
        };
        assert!(with_charge.has(required));
    }

    #[test]
    fn header_lookup_accepts_descriptor_and_reaction_column_names() {
        let descriptors = csv::StringRecord::from(vec!["file", "sterimol_b5", "pyr_p"]);
        assert_eq!(header_index(&descriptors, "label"), Some(0));
        assert_eq!(header_index(&descriptors, "b5"), Some(1));
        assert_eq!(header_index(&descriptors, "nbo_charge"), None);

        let reaction = csv::StringRecord::from(vec!["Reaction_ID", "NBO_Charge", "IR_Frequency"]);
        assert_eq!(header_index(&reaction, "label"), Some(0));
        assert_eq!(header_index(&reaction, "nbo_charge"), Some(1));
        assert_eq!(header_index(&reaction, "ir_frequency"), Some(2));
    }
}
