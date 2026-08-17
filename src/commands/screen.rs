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
use steric_x::model::{MODEL_FEATURE_NAMES, expand_features};
use steric_x::{BuriedVolumeConfig, EyringKineticLink, PackedReactionRecord, ScientificFitReport};

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
    label: String,
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

#[derive(Clone, Debug, Serialize)]
struct ScreenHit {
    rank: usize,
    ligand: String,
    predicted_ddg_kcal_mol: f64,
    /// Conservative bounds from the bootstrap coefficient intervals. NOT an
    /// OLS prediction interval — see the note emitted with the report.
    coefficient_band_low: Option<f64>,
    coefficient_band_high: Option<f64>,
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
    screened: usize,
    skipped: usize,
    inside_domain: usize,
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
    pub(crate) donor_element: &'a str,
    pub(crate) sterimol_axis: SterimolAxis,
    pub(crate) format: DescriptorFormat,
    pub(crate) config: BuriedVolumeConfig,
}

fn load_model(path: &Path) -> Result<ScientificFitReport, Box<dyn Error>> {
    let contents = std::fs::read_to_string(path)
        .map_err(|error| format!("could not read model {}: {error}", path.display()))?;
    serde_json::from_str::<ScientificFitReport>(&contents).map_err(|error| {
        format!(
            "{} is not a StericX fit report produced by `stericx fit`: {error}",
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
            "file",
            "Reaction_ID",
            "reaction_id",
            "id",
            "Ligand_XYZ_Path",
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
) -> Result<(Vec<Candidate>, Available), Box<dyn Error>> {
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
        let candidates = paths
            .par_iter()
            .filter_map(|path| {
                match descriptors_for_file(path, donor_element, None, sterimol_axis, config) {
                    Ok(result) => Some(Candidate {
                        label: result.file.clone(),
                        l: Some(result.sterimol_l),
                        b1: Some(result.sterimol_b1),
                        b5: Some(result.sterimol_b5),
                        ..Candidate::default()
                    }),
                    Err(message) => {
                        eprintln!("skipped {}: {message}", path.display());
                        None
                    }
                }
            })
            .collect::<Vec<_>>();
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
        candidates.push(Candidate {
            label,
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
                    eprintln!(
                        "skipped {}: geometry not found: {}",
                        candidate.label,
                        path.display()
                    );
                    return None;
                };
                match descriptors_for_file(&resolved, donor_element, None, sterimol_axis, config) {
                    Ok(result) => Some((result.sterimol_l, result.sterimol_b1, result.sterimol_b5)),
                    Err(message) => {
                        eprintln!("skipped {}: {message}", candidate.label);
                        None
                    }
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
    Ok((candidates, available))
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

pub(crate) fn screen_command(args: ScreenArgs<'_>) -> Result<(), Box<dyn Error>> {
    if !args.temperature.is_finite() || args.temperature <= 0.0 {
        return Err("--temperature must be a positive finite temperature".into());
    }
    let report = load_model(args.model)?;
    let required = RequiredInputs::from_weights(&report.weights);

    let (candidates, available) = load_library(
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

    let mut skipped = 0_usize;
    let mut hits = Vec::new();
    for candidate in &candidates {
        if !candidate.has(required) {
            skipped += 1;
            eprintln!(
                "skipped {}: missing a required input value",
                candidate.label
            );
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
            skipped += 1;
            eprintln!("skipped {}: prediction is not finite", candidate.label);
            continue;
        }
        let outside = domain_exceedances(&report, &features);
        let band = coefficient_band(&report, &features);
        hits.push(ScreenHit {
            rank: 0,
            ligand: candidate.label.clone(),
            predicted_ddg_kcal_mol: predicted,
            coefficient_band_low: band.map(|(low, _)| low),
            coefficient_band_high: band.map(|(_, high)| high),
            // `calculate_enantiomeric_excess` is an absolute value, so carry the
            // ΔΔG‡ sign through: under this project's convention a positive
            // ΔΔG‡ favours R, and a negative prediction means the same excess
            // of the *opposite* enantiomer.
            predicted_ee_percent: f64::from(EyringKineticLink::calculate_enantiomeric_excess(
                predicted as f32,
                args.temperature,
            )) * if predicted < 0.0 { -1.0 } else { 1.0 },
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

    hits.sort_by(|left, right| {
        let ordering = left
            .predicted_ddg_kcal_mol
            .total_cmp(&right.predicted_ddg_kcal_mol);
        let ordering = if args.ascending {
            ordering
        } else {
            ordering.reverse()
        };
        ordering.then_with(|| left.ligand.cmp(&right.ligand))
    });
    let screened = hits.len();
    if let Some(top) = args.top {
        hits.truncate(top);
    }
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
        skipped,
        inside_domain,
        uncertainty_note:
            "coefficient_band propagates the bootstrap 95% coefficient intervals by interval \
             arithmetic; it ignores coefficient correlation (so it is conservative) and is NOT \
             an OLS prediction interval. Residual scatter is reported separately as \
             training_rmse_kcal_mol."
                .to_owned(),
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
        "library        {} members, {} screened, {} skipped, {} inside the training domain",
        report.library_size, report.screened, report.skipped, report.inside_domain
    );
    println!("temperature    {:.2} K", report.temperature_k);
    println!(
        "\n{:>4}  {:>9}  {:>19}  {:>7}  {:<28}  ligand",
        "rank", "pred ddG", "coef. band (95%)", "pred ee", "applicability"
    );
    for hit in &report.hits {
        let band = match (hit.coefficient_band_low, hit.coefficient_band_high) {
            (Some(low), Some(high)) => format!("[{low:>7.2}, {high:>7.2}]"),
            _ => "unavailable".to_owned(),
        };
        println!(
            "{:>4}  {:>9.3}  {:>19}  {:>6.1}%  {:<28}  {}",
            hit.rank,
            hit.predicted_ddg_kcal_mol,
            band,
            hit.predicted_ee_percent,
            hit.applicability,
            hit.ligand
        );
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
    println!("\nnote: {}", report.uncertainty_note);
}

fn print_screen_csv(report: &ScreenReport) {
    println!(
        "rank,ligand,predicted_ddg_kcal_mol,coefficient_band_low,coefficient_band_high,\
         predicted_ee_percent,applicability,outside_features"
    );
    for hit in &report.hits {
        let format_bound =
            |value: Option<f64>| value.map_or_else(String::new, |value| format!("{value:.6}"));
        println!(
            "{},{},{:.6},{},{},{:.4},{},{}",
            hit.rank,
            crate::descriptors::csv_field(&hit.ligand),
            hit.predicted_ddg_kcal_mol,
            format_bound(hit.coefficient_band_low),
            format_bound(hit.coefficient_band_high),
            hit.predicted_ee_percent,
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
