use clap::{Parser, Subcommand, ValueEnum};
use glam::Vec3;
use serde::{Deserialize, Serialize};
use std::error::Error;
use std::fs;
use std::mem::size_of;
use std::path::{Component, Path, PathBuf};
use std::time::{Duration, Instant};
use steric_x::model::{MODEL_FEATURE_COUNT, expand_features};
use steric_x::{
    BuriedVolumeCalculator, BuriedVolumeConfig, BuriedVolumeParams, EyringKineticLink, FitOptions,
    Molecule, PackedBuriedVolumeRecord, PackedReactionRecord, PackedReactionRecordV2,
    RegressXPredictor, ScientificFitReport, SigPackReader, SigPackV2Writer, SigPackWriter,
    SterimolCalculator, SterimolParams, coordination_center, covalent_radius,
    fit_scientific_model_grouped, parse_coordinate_file,
};

#[derive(Debug, Parser)]
#[command(
    name = "stericx",
    version,
    about = "SIMD physical-organic featurization and selectivity prediction"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Parse reaction metadata and XYZ files into a flat binary matrix.
    Parse {
        /// Reaction CSV produced by the data-preparation pipeline.
        #[arg(long)]
        csv: PathBuf,
        /// Root directory containing ligand XYZ files.
        #[arg(long)]
        xyz_dir: PathBuf,
        /// Destination `.sigpack` matrix.
        #[arg(long)]
        output: PathBuf,
    },
    /// Calculate coordination-aware buried volumes and write `.sigpack` v2.
    BuriedVolume {
        /// Reaction CSV produced by the data-preparation pipeline.
        #[arg(long)]
        csv: PathBuf,
        /// Root directory containing ligand conformer XYZ files.
        #[arg(long)]
        xyz_dir: PathBuf,
        /// Destination version-two `.sigpack` matrix.
        #[arg(long)]
        output: PathBuf,
        /// Optional per-conformer descriptor CSV for validation and audit.
        #[arg(long)]
        per_conformer_output: Option<PathBuf>,
        /// Metal-centred integration sphere radius in ångströms.
        #[arg(long, default_value_t = 3.5)]
        sphere_radius: f32,
        /// Grid density in Å³ per point (Kraken uses 0.01).
        #[arg(long, default_value_t = 0.01)]
        density: f32,
        /// Donor-to-virtual-metal distance in ångströms.
        #[arg(long, default_value_t = 2.1)]
        center_distance: f32,
        /// Bondi radius scale factor used by Morfeus.
        #[arg(long, default_value_t = 1.17)]
        radii_scale: f32,
        /// Fail if the CSV lacks one xTB-derived center per conformer.
        #[arg(long)]
        require_explicit_centers: bool,
    },
    /// Run parallel SIMD regression over a mapped `.sigpack` matrix.
    Predict {
        /// Input `.sigpack` matrix.
        #[arg(long)]
        data: PathBuf,
        /// JSON array of eight weights or object containing a `weights` array.
        #[arg(long)]
        weights: PathBuf,
    },
    /// Fit and freeze an interpretable physical-organic regression model.
    Fit {
        /// Input `.sigpack` matrix.
        #[arg(long)]
        data: PathBuf,
        /// Row-aligned CSV containing Reaction_ID, Dataset_Split, and Ligand_Group.
        #[arg(long)]
        metadata: PathBuf,
        /// Destination model and scientific diagnostics JSON.
        #[arg(long)]
        output: PathBuf,
        /// Frozen predictions for every non-training row.
        #[arg(long)]
        predictions: PathBuf,
        /// Maximum number of non-intercept model terms.
        #[arg(long, default_value_t = 3)]
        max_terms: usize,
        /// Bootstrap coefficient-interval replicates.
        #[arg(long, default_value_t = 1_000)]
        bootstrap: usize,
        /// Response-permutation null replicates.
        #[arg(long, default_value_t = 500)]
        permutations: usize,
        /// Deterministic resampling seed.
        #[arg(long, default_value_t = 20_260_725)]
        seed: u64,
    },
    /// Reveal and score previously frozen non-training predictions.
    Evaluate {
        /// Input `.sigpack` matrix containing experimental targets.
        #[arg(long)]
        data: PathBuf,
        /// Row-aligned split/provenance CSV used during fitting.
        #[arg(long)]
        metadata: PathBuf,
        /// Frozen model JSON produced by `stericx fit`.
        #[arg(long)]
        model: PathBuf,
        /// Frozen prediction CSV produced by `stericx fit`.
        #[arg(long)]
        predictions: PathBuf,
        /// Destination scored evaluation JSON.
        #[arg(long)]
        output: PathBuf,
    },
    /// Calculate an Eyring rate and enantiomeric product distribution.
    Simulate {
        /// ΔΔG‡ in kcal/mol.
        #[arg(long, allow_hyphen_values = true)]
        ddg: f32,
        /// Temperature in kelvin.
        #[arg(long)]
        temp: f32,
    },
    /// Compute Sterimol and buried-volume descriptors for one or more ligand files.
    ///
    /// Auto-detects the phosphine donor and its substituents straight from the
    /// geometry — no reaction CSV and no manual atom indices. Accepts `.xyz`
    /// (one geometry) and `.sdf`/`.mol` (one or more conformers); a multi-model
    /// file is treated as a conformer ensemble and reported as Kraken's tables
    /// do, including the `max_delta_qvbur_min` headline descriptor.
    Descriptors {
        /// Ligand coordinate files (`.xyz`, `.sdf`, `.mol`). Globs are expanded by the shell.
        #[arg(required = true, value_name = "FILE")]
        inputs: Vec<PathBuf>,
        /// Donor element to locate (defaults to phosphorus).
        #[arg(long, default_value = "P")]
        donor_element: String,
        /// Explicit zero-based donor atom index, overriding auto-detection.
        #[arg(long)]
        donor_index: Option<usize>,
        /// Sterimol axis: `bond` (donor→substituent) or `coordination`
        /// (metal→donor along the lone pair, the Kraken/ligand convention).
        #[arg(long, value_enum, default_value_t = SterimolAxis::Bond)]
        sterimol_axis: SterimolAxis,
        /// Output format.
        #[arg(long, value_enum, default_value_t = DescriptorFormat::Text)]
        format: DescriptorFormat,
        /// Metal-centred integration sphere radius in ångströms.
        #[arg(long, default_value_t = 3.5)]
        sphere_radius: f32,
        /// Grid density in Å³ per point.
        #[arg(long, default_value_t = 0.01)]
        density: f32,
        /// Donor-to-virtual-metal distance in ångströms (Kraken convention = 2.28).
        #[arg(long, default_value_t = 2.28)]
        center_distance: f32,
        /// Bondi radius scale factor used by Morfeus.
        #[arg(long, default_value_t = 1.17)]
        radii_scale: f32,
    },
}

/// Output format for the `descriptors` subcommand.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum DescriptorFormat {
    /// Human-readable summary, one block per file.
    Text,
    /// One JSON array of per-file descriptor records.
    Json,
    /// Comma-separated table, one row per file (spreadsheet-ready).
    Csv,
}

/// Which axis defines the Sterimol frame.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum SterimolAxis {
    /// Donor → nearest bonded substituent (a P–C bond direction).
    Bond,
    /// Virtual metal → donor, along the lone pair. This is the convention
    /// Kraken and the ligand-descriptor literature use, and it applies the
    /// historical +0.40 Å Verloop correction to `L`.
    Coordination,
}

/// Historical Verloop/Morfeus correction added to the raw geometric Sterimol
/// `L` under the coordination-axis convention.
const STERIMOL_L_CORRECTION: f32 = 0.40;

/// One row of the raw reaction CSV.
#[derive(Clone, Debug, Deserialize)]
struct ReactionCsvRow {
    #[serde(rename = "Reaction_ID", alias = "reaction_id")]
    reaction_id: String,
    #[serde(rename = "Ligand_XYZ_Path", alias = "ligand_xyz_path", alias = "file")]
    ligand_xyz_path: PathBuf,
    #[serde(rename = "Attach_Atom_Idx", alias = "attach_idx")]
    attach_atom_idx: usize,
    #[serde(
        rename = "Primary_Bond_Vector_Idx",
        alias = "primary_bond_vector_idx",
        alias = "axis_atom_idx",
        alias = "neighbor_idx"
    )]
    primary_bond_vector_idx: usize,
    #[serde(rename = "NBO_Charge", alias = "nbo_charge")]
    nbo_charge: f32,
    #[serde(rename = "IR_Frequency", alias = "ir_freq")]
    ir_freq: f32,
    #[serde(rename = "Temp_K", alias = "temp_k")]
    temp_k: f32,
    #[serde(rename = "Exp_ddG_kcal_mol", alias = "exp_ddg")]
    exp_ddg: f32,
    #[serde(rename = "Conformer_XYZ_Paths", alias = "conformer_xyz_paths", default)]
    conformer_xyz_paths: Option<String>,
    #[serde(
        rename = "Conformer_Relative_Energies_kcal_mol",
        alias = "conformer_relative_energies_kcal_mol",
        default
    )]
    conformer_relative_energies: Option<String>,
    #[serde(
        rename = "Conformer_Boltzmann_Weights",
        alias = "conformer_boltzmann_weights",
        default
    )]
    conformer_boltzmann_weights: Option<String>,
    #[serde(
        rename = "Conformer_Coordination_Centers_Angstrom",
        alias = "conformer_coordination_centers_angstrom",
        default
    )]
    conformer_coordination_centers: Option<String>,
    #[serde(
        rename = "Coordination_Center_Method",
        alias = "coordination_center_method",
        default
    )]
    coordination_center_method: Option<String>,
}

#[derive(Clone, Debug, Deserialize)]
struct ReactionMetadataRow {
    #[serde(rename = "Reaction_ID", alias = "reaction_id")]
    reaction_id: String,
    #[serde(rename = "Dataset_Split", alias = "dataset_split")]
    dataset_split: String,
    #[serde(rename = "Ligand_Group", alias = "ligand_group", default)]
    ligand_group: String,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct FrozenPredictionRow {
    #[serde(rename = "Reaction_ID")]
    reaction_id: String,
    #[serde(rename = "Ligand_Group")]
    ligand_group: String,
    #[serde(rename = "Dataset_Split")]
    dataset_split: String,
    #[serde(rename = "Predicted_ddG_kcal_mol")]
    predicted_ddg: f32,
    #[serde(rename = "Applicability_Domain")]
    applicability_domain: String,
}

#[derive(Debug, Serialize)]
struct EvaluationReport {
    schema_version: u32,
    model_path: String,
    frozen_predictions_path: String,
    frozen_predictions_fnv1a64: String,
    evaluated_records: usize,
    mae_kcal_mol: f64,
    rmse_kcal_mol: f64,
    r2: Option<f64>,
    applicability_warnings: usize,
    scored_predictions: Vec<ScoredPrediction>,
}

#[derive(Debug, Serialize)]
struct ScoredPrediction {
    reaction_id: String,
    ligand_group: String,
    dataset_split: String,
    predicted_ddg_kcal_mol: f32,
    experimental_ddg_kcal_mol: f32,
    residual_kcal_mol: f32,
    applicability_domain: String,
}

#[derive(Debug, Deserialize)]
#[serde(untagged)]
enum WeightDocument {
    Array([f32; MODEL_FEATURE_COUNT]),
    Object { weights: [f32; MODEL_FEATURE_COUNT] },
}

fn main() {
    if let Err(error) = run() {
        eprintln!("error: {error}");
        std::process::exit(2);
    }
}

fn run() -> Result<(), Box<dyn Error>> {
    match Cli::parse().command {
        Command::Parse {
            csv,
            xyz_dir,
            output,
        } => parse_command(&csv, &xyz_dir, &output),
        Command::BuriedVolume {
            csv,
            xyz_dir,
            output,
            per_conformer_output,
            sphere_radius,
            density,
            center_distance,
            radii_scale,
            require_explicit_centers,
        } => buried_volume_command(
            &csv,
            &xyz_dir,
            &output,
            per_conformer_output.as_deref(),
            BuriedVolumeConfig {
                sphere_radius,
                density,
                center_distance,
                radii_scale,
                include_hydrogens: false,
            },
            require_explicit_centers,
        ),
        Command::Predict { data, weights } => predict_command(&data, &weights),
        Command::Fit {
            data,
            metadata,
            output,
            predictions,
            max_terms,
            bootstrap,
            permutations,
            seed,
        } => fit_command(
            &data,
            &metadata,
            &output,
            &predictions,
            FitOptions {
                max_terms,
                bootstrap_samples: bootstrap,
                permutation_samples: permutations,
                seed,
            },
        ),
        Command::Evaluate {
            data,
            metadata,
            model,
            predictions,
            output,
        } => evaluate_command(&data, &metadata, &model, &predictions, &output),
        Command::Simulate { ddg, temp } => simulate_command(ddg, temp),
        Command::Descriptors {
            inputs,
            donor_element,
            donor_index,
            sterimol_axis,
            format,
            sphere_radius,
            density,
            center_distance,
            radii_scale,
        } => descriptors_command(
            &inputs,
            &donor_element,
            donor_index,
            sterimol_axis,
            format,
            BuriedVolumeConfig {
                sphere_radius,
                density,
                center_distance,
                radii_scale,
                include_hydrogens: false,
            },
        ),
    }
}

#[derive(Clone, Debug, Serialize)]
struct PerConformerBuriedVolumeRow {
    #[serde(rename = "Reaction_ID")]
    reaction_id: String,
    #[serde(rename = "Conformer_Index")]
    conformer_index: usize,
    #[serde(rename = "Conformer_XYZ_Path")]
    conformer_xyz_path: String,
    #[serde(rename = "Boltzmann_Weight")]
    boltzmann_weight: f32,
    #[serde(rename = "vbur")]
    vbur: f32,
    #[serde(rename = "percent_vbur")]
    percent_vbur: f32,
    #[serde(rename = "qvbur_min")]
    qvbur_min: f32,
    #[serde(rename = "qvbur_max")]
    qvbur_max: f32,
    #[serde(rename = "max_delta_qvbur")]
    max_delta_qvbur: f32,
    #[serde(rename = "ovbur_min")]
    ovbur_min: f32,
    #[serde(rename = "ovbur_max")]
    ovbur_max: f32,
    #[serde(rename = "near_vbur")]
    near_vbur: f32,
    #[serde(rename = "far_vbur")]
    far_vbur: f32,
    #[serde(rename = "Coordination_Center_Method")]
    coordination_center_method: String,
    #[serde(rename = "Coordination_Center_X")]
    coordination_center_x: Option<f32>,
    #[serde(rename = "Coordination_Center_Y")]
    coordination_center_y: Option<f32>,
    #[serde(rename = "Coordination_Center_Z")]
    coordination_center_z: Option<f32>,
}

fn parse_command(
    reactions_csv: &Path,
    xyz_dir: &Path,
    output: &Path,
) -> Result<(), Box<dyn Error>> {
    let total_started = Instant::now();
    let rss_start = resident_memory_bytes();
    if !reactions_csv.is_file() {
        return Err(format!("reaction CSV does not exist: {}", reactions_csv.display()).into());
    }
    if !xyz_dir.is_dir() {
        return Err(format!("XYZ directory does not exist: {}", xyz_dir.display()).into());
    }

    let input_bytes = fs::metadata(reactions_csv)?.len();
    let ingest_started = Instant::now();
    let mut csv_reader = csv::ReaderBuilder::new()
        .trim(csv::Trim::All)
        .from_path(reactions_csv)?;
    let mut records = Vec::new();
    let mut xyz_bytes = 0_u64;
    let mut geometry_time = Duration::ZERO;

    for (row_index, result) in csv_reader.deserialize::<ReactionCsvRow>().enumerate() {
        let row = result.map_err(|error| {
            format!(
                "{} row {} could not be parsed: {error}",
                reactions_csv.display(),
                row_index + 2
            )
        })?;
        validate_reaction_row(&row, row_index + 2)?;
        let geometry_started = Instant::now();
        let coordinate_paths = conformer_paths(&row)?;
        let weights = conformer_weights(&row, coordinate_paths.len())?;
        let energy_span = conformer_energy_span(&row, coordinate_paths.len())?;
        let mut conformer_params = Vec::with_capacity(coordinate_paths.len());
        for csv_coordinate_path in &coordinate_paths {
            let xyz_path = resolve_xyz_path(xyz_dir, csv_coordinate_path).map_err(|message| {
                format!(
                    "{} row {} ({}): {message}",
                    reactions_csv.display(),
                    row_index + 2,
                    row.reaction_id
                )
            })?;
            xyz_bytes = xyz_bytes.saturating_add(fs::metadata(&xyz_path)?.len());
            let molecule = Molecule::from_xyz_file(&xyz_path).map_err(|error| {
                format!(
                    "reaction {} failed to parse {}: {error}",
                    row.reaction_id,
                    xyz_path.display()
                )
            })?;
            conformer_params.push(sterimol_from_molecule(&molecule, &row).map_err(|message| {
                format!(
                    "reaction {} using {}: {message}",
                    row.reaction_id,
                    xyz_path.display()
                )
            })?);
        }
        let record = record_from_ensemble(&conformer_params, &weights, energy_span, &row)?;
        geometry_time += geometry_started.elapsed();
        records.push(record);
    }
    let ingest_time = ingest_started.elapsed();
    if records.is_empty() {
        return Err("reaction CSV contains no data rows".into());
    }

    if let Some(parent) = output
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)?;
    }
    let export_started = Instant::now();
    SigPackWriter::export(&records, output)?;
    let export_time = export_started.elapsed();
    let output_bytes = fs::metadata(output)?.len();
    let total_time = total_started.elapsed();
    let throughput = records.len() as f64 / total_time.as_secs_f64().max(f64::EPSILON);

    println!("command=parse");
    println!("records_processed={}", records.len());
    println!("csv_input={}", reactions_csv.display());
    println!("xyz_directory={}", xyz_dir.display());
    println!("sigpack_output={}", output.display());
    println!("csv_bytes={input_bytes}");
    println!("xyz_bytes_read={xyz_bytes}");
    println!("sigpack_bytes={output_bytes}");
    println!(
        "record_buffer_bytes={}",
        records.len() * size_of::<PackedReactionRecord>()
    );
    println!("csv_and_geometry_ms={:.3}", millis(ingest_time));
    println!("geometry_compute_ms={:.3}", millis(geometry_time));
    println!("binary_export_ms={:.3}", millis(export_time));
    println!("total_ms={:.3}", millis(total_time));
    println!("throughput_records_per_second={throughput:.1}");
    print_memory_metrics(rss_start, resident_memory_bytes());
    Ok(())
}

fn buried_volume_command(
    reactions_csv: &Path,
    xyz_dir: &Path,
    output: &Path,
    per_conformer_output: Option<&Path>,
    config: BuriedVolumeConfig,
    require_explicit_centers: bool,
) -> Result<(), Box<dyn Error>> {
    let total_started = Instant::now();
    let rss_start = resident_memory_bytes();
    if !reactions_csv.is_file() {
        return Err(format!("reaction CSV does not exist: {}", reactions_csv.display()).into());
    }
    if !xyz_dir.is_dir() {
        return Err(format!("XYZ directory does not exist: {}", xyz_dir.display()).into());
    }

    let mut csv_reader = csv::ReaderBuilder::new()
        .trim(csv::Trim::All)
        .from_path(reactions_csv)?;
    let mut packed_records = Vec::new();
    let mut audit_rows = Vec::new();
    let mut conformers_processed = 0_usize;
    let mut xyz_bytes = 0_u64;
    for (row_index, result) in csv_reader.deserialize::<ReactionCsvRow>().enumerate() {
        let row = result.map_err(|error| {
            format!(
                "{} row {} could not be parsed: {error}",
                reactions_csv.display(),
                row_index + 2
            )
        })?;
        validate_reaction_row(&row, row_index + 2)?;
        let coordinate_paths = conformer_paths(&row)?;
        let weights = conformer_weights(&row, coordinate_paths.len())?;
        let energy_span = conformer_energy_span(&row, coordinate_paths.len())?;
        let coordination_centers = conformer_coordination_centers(&row, coordinate_paths.len())?;
        if require_explicit_centers && coordination_centers.is_none() {
            return Err(format!(
                "reaction {} has no explicit conformer coordination centers",
                row.reaction_id
            )
            .into());
        }
        let center_method = row
            .coordination_center_method
            .as_deref()
            .filter(|value| !value.trim().is_empty())
            .unwrap_or("geometric_opposed_substituent_vector")
            .to_owned();
        let mut sterimol_params = Vec::with_capacity(coordinate_paths.len());
        let mut buried_volume_params = Vec::with_capacity(coordinate_paths.len());

        for (conformer_index, csv_coordinate_path) in coordinate_paths.iter().enumerate() {
            let xyz_path = resolve_xyz_path(xyz_dir, csv_coordinate_path).map_err(|message| {
                format!(
                    "{} row {} ({}): {message}",
                    reactions_csv.display(),
                    row_index + 2,
                    row.reaction_id
                )
            })?;
            xyz_bytes = xyz_bytes.saturating_add(fs::metadata(&xyz_path)?.len());
            let molecule = Molecule::from_xyz_file(&xyz_path).map_err(|error| {
                format!(
                    "reaction {} failed to parse {}: {error}",
                    row.reaction_id,
                    xyz_path.display()
                )
            })?;
            sterimol_params.push(sterimol_from_molecule(&molecule, &row)?);
            let explicit_center = coordination_centers
                .as_ref()
                .map(|centers| centers[conformer_index]);
            let buried = explicit_center
                .map_or_else(
                    || {
                        BuriedVolumeCalculator::compute(
                            &molecule,
                            row.attach_atom_idx,
                            row.primary_bond_vector_idx,
                            config,
                        )
                    },
                    |center| {
                        BuriedVolumeCalculator::compute_with_center(
                            &molecule,
                            row.attach_atom_idx,
                            row.primary_bond_vector_idx,
                            center,
                            config,
                        )
                    },
                )
                .map_err(|error| {
                    format!(
                        "reaction {} using {}: {error}",
                        row.reaction_id,
                        xyz_path.display()
                    )
                })?;
            audit_rows.push(PerConformerBuriedVolumeRow {
                reaction_id: row.reaction_id.clone(),
                conformer_index,
                conformer_xyz_path: csv_coordinate_path.display().to_string(),
                boltzmann_weight: weights[conformer_index],
                vbur: buried.buried_volume,
                percent_vbur: buried.percent_buried_volume,
                qvbur_min: buried.qvbur_min,
                qvbur_max: buried.qvbur_max,
                max_delta_qvbur: buried.max_delta_qvbur,
                ovbur_min: buried.ovbur_min,
                ovbur_max: buried.ovbur_max,
                near_vbur: buried.near_vbur,
                far_vbur: buried.far_vbur,
                coordination_center_method: center_method.clone(),
                coordination_center_x: explicit_center.map(|center| center.x),
                coordination_center_y: explicit_center.map(|center| center.y),
                coordination_center_z: explicit_center.map(|center| center.z),
            });
            buried_volume_params.push(buried);
            conformers_processed += 1;
        }

        let reaction = record_from_ensemble(&sterimol_params, &weights, energy_span, &row)?;
        let ensemble = BuriedVolumeCalculator::aggregate(&buried_volume_params, &weights)?;
        let buried_volume = PackedBuriedVolumeRecord {
            vbur_boltz: ensemble.vbur_boltz,
            vbur_min: ensemble.vbur_min,
            vbur_max: ensemble.vbur_max,
            vbur_delta: ensemble.vbur_delta,
            qvbur_min_boltz: ensemble.qvbur_min_boltz,
            qvbur_max_boltz: ensemble.qvbur_max_boltz,
            max_delta_qvbur_boltz: ensemble.max_delta_qvbur_boltz,
            max_delta_qvbur_min: ensemble.max_delta_qvbur_min,
            max_delta_qvbur_max: ensemble.max_delta_qvbur_max,
            max_delta_qvbur_delta: ensemble.max_delta_qvbur_delta,
            max_delta_qvbur_vburminconf: ensemble.max_delta_qvbur_vburminconf,
            near_vbur_boltz: ensemble.near_vbur_boltz,
            far_vbur_boltz: ensemble.far_vbur_boltz,
            conformer_count: ensemble.conformer_count as f32,
            sphere_radius: config.sphere_radius,
            grid_density: config.density,
        };
        packed_records.push(PackedReactionRecordV2 {
            reaction,
            buried_volume,
        });
        println!(
            "buried_volume_progress_rows={},reaction_id={},conformers={}",
            packed_records.len(),
            row.reaction_id,
            coordinate_paths.len()
        );
    }
    if packed_records.is_empty() {
        return Err("reaction CSV contains no data rows".into());
    }

    ensure_parent(output)?;
    let export_started = Instant::now();
    SigPackV2Writer::export(&packed_records, output)?;
    let export_time = export_started.elapsed();
    if let Some(path) = per_conformer_output {
        atomic_write_csv_rows(&audit_rows, path)?;
    }

    let total_time = total_started.elapsed();
    println!("command=buried-volume");
    println!("schema_version=2");
    println!("records_processed={}", packed_records.len());
    println!("conformers_processed={conformers_processed}");
    println!("sphere_radius_angstrom={:.4}", config.sphere_radius);
    println!("grid_density_angstrom3={:.6}", config.density);
    println!(
        "virtual_center_distance_angstrom={:.4}",
        config.center_distance
    );
    println!("radii_scale={:.4}", config.radii_scale);
    println!("explicit_coordination_centers_required={require_explicit_centers}");
    println!("official_kraken_center_method=xtb_lmo_center");
    println!("xyz_bytes_read={xyz_bytes}");
    println!("sigpack_v2_output={}", output.display());
    println!("sigpack_v2_bytes={}", fs::metadata(output)?.len());
    if let Some(path) = per_conformer_output {
        println!("per_conformer_output={}", path.display());
    }
    println!("binary_export_ms={:.3}", millis(export_time));
    println!("total_ms={:.3}", millis(total_time));
    println!(
        "throughput_conformers_per_second={:.1}",
        conformers_processed as f64 / total_time.as_secs_f64().max(f64::EPSILON)
    );
    print_memory_metrics(rss_start, resident_memory_bytes());
    Ok(())
}

fn predict_command(data: &Path, weights_path: &Path) -> Result<(), Box<dyn Error>> {
    let total_started = Instant::now();
    let rss_start = resident_memory_bytes();
    if !data.is_file() {
        return Err(format!("sigpack file does not exist: {}", data.display()).into());
    }
    if !weights_path.is_file() {
        return Err(format!("weight file does not exist: {}", weights_path.display()).into());
    }

    let weights_started = Instant::now();
    let weights = load_weights_json(weights_path)?;
    let weights_time = weights_started.elapsed();
    let predictor = RegressXPredictor::new(weights);

    let mapping_started = Instant::now();
    let reader = SigPackReader::open(data)?;
    let mapping_time = mapping_started.elapsed();
    let records = reader.records();
    if records.is_empty() {
        return Err("sigpack matrix contains no records".into());
    }

    let inference_started = Instant::now();
    let predictions = predictor.predict_batch(records);
    let inference_time = inference_started.elapsed();
    let mse = predictions
        .iter()
        .zip(records)
        .map(|(prediction, record)| {
            let residual = f64::from(*prediction) - f64::from(record.exp_ddg);
            residual * residual
        })
        .sum::<f64>()
        / records.len() as f64;
    let total_time = total_started.elapsed();
    let inference_throughput =
        records.len() as f64 / inference_time.as_secs_f64().max(f64::EPSILON);

    println!("command=predict");
    println!("records_predicted={}", records.len());
    println!("rayon_threads={}", rayon::current_num_threads());
    println!("sigpack_input={}", data.display());
    println!("weights_input={}", weights_path.display());
    println!("mapped_bytes={}", fs::metadata(data)?.len());
    println!(
        "prediction_buffer_bytes={}",
        predictions.len() * size_of::<f32>()
    );
    println!("weights_load_ms={:.3}", millis(weights_time));
    println!("memory_map_ms={:.3}", millis(mapping_time));
    println!("prediction_latency_ms={:.3}", millis(inference_time));
    println!("total_ms={:.3}", millis(total_time));
    println!("throughput_records_per_second={inference_throughput:.1}");
    println!("mse_kcal2_per_mol2={mse:.8}");
    for (index, prediction) in predictions.iter().take(5).enumerate() {
        println!(
            "prediction_preview[{index}]={prediction:.7},experimental={:.7}",
            records[index].exp_ddg
        );
    }
    print_memory_metrics(rss_start, resident_memory_bytes());
    Ok(())
}

fn fit_command(
    data: &Path,
    metadata_path: &Path,
    output: &Path,
    predictions_path: &Path,
    options: FitOptions,
) -> Result<(), Box<dyn Error>> {
    let total_started = Instant::now();
    let rss_start = resident_memory_bytes();
    let reader = SigPackReader::open(data)?;
    let records = reader.records();
    if records.is_empty() {
        return Err("sigpack matrix contains no records".into());
    }
    let metadata = load_reaction_metadata(metadata_path)?;
    if metadata.len() != records.len() {
        return Err(format!(
            "metadata has {} rows but sigpack contains {} records",
            metadata.len(),
            records.len()
        )
        .into());
    }
    let training_indices = metadata
        .iter()
        .enumerate()
        .filter(|(_, row)| row.dataset_split.eq_ignore_ascii_case("train"))
        .map(|(index, _)| index)
        .collect::<Vec<_>>();
    let frozen_indices = metadata
        .iter()
        .enumerate()
        .filter(|(_, row)| !row.dataset_split.eq_ignore_ascii_case("train"))
        .map(|(index, _)| index)
        .collect::<Vec<_>>();
    if frozen_indices.is_empty() {
        return Err("metadata contains no non-training rows to freeze".into());
    }
    let training_groups = training_indices
        .iter()
        .map(|&index| {
            let group = metadata[index].ligand_group.trim();
            if group.is_empty() {
                metadata[index].reaction_id.clone()
            } else {
                group.to_owned()
            }
        })
        .collect::<Vec<_>>();

    let fit_started = Instant::now();
    let report =
        fit_scientific_model_grouped(records, &training_indices, &training_groups, options)?;
    let fit_time = fit_started.elapsed();
    atomic_write_json(&report, output)?;

    let predictor = RegressXPredictor::new(report.weights);
    let frozen_records = frozen_indices
        .iter()
        .map(|&index| records[index])
        .collect::<Vec<_>>();
    let predicted = predictor.predict_batch(&frozen_records);
    let frozen_rows = frozen_indices
        .iter()
        .zip(predicted)
        .map(|(&index, predicted_ddg)| FrozenPredictionRow {
            reaction_id: metadata[index].reaction_id.clone(),
            ligand_group: metadata[index].ligand_group.clone(),
            dataset_split: metadata[index].dataset_split.clone(),
            predicted_ddg,
            applicability_domain: applicability_status(&records[index], &report),
        })
        .collect::<Vec<_>>();
    atomic_write_csv_rows(&frozen_rows, predictions_path)?;

    let total_time = total_started.elapsed();
    println!("command=fit");
    println!("records_total={}", records.len());
    println!("training_records={}", training_indices.len());
    println!("training_groups={}", report.training_group_count);
    println!("frozen_prediction_records={}", frozen_rows.len());
    println!("selected_features={}", report.selected_features.join(","));
    println!(
        "training_r2={}",
        display_optional_metric(report.training.r2)
    );
    println!(
        "fixed_feature_loo_r2={}",
        display_optional_metric(report.fixed_feature_loo.r2)
    );
    println!(
        "fixed_feature_group_loo_r2={}",
        display_optional_metric(report.fixed_feature_group_loo.r2)
    );
    println!("training_rmse={:.8}", report.training.rmse);
    println!(
        "fixed_feature_loo_rmse={:.8}",
        report.fixed_feature_loo.rmse
    );
    println!(
        "fixed_feature_group_loo_rmse={:.8}",
        report.fixed_feature_group_loo.rmse
    );
    println!(
        "response_permutation_p_value={:.8}",
        report.response_permutation_p_value
    );
    println!("model_output={}", output.display());
    println!("frozen_predictions_output={}", predictions_path.display());
    println!("fit_ms={:.3}", millis(fit_time));
    println!("total_ms={:.3}", millis(total_time));
    print_memory_metrics(rss_start, resident_memory_bytes());
    Ok(())
}

fn evaluate_command(
    data: &Path,
    metadata_path: &Path,
    model_path: &Path,
    predictions_path: &Path,
    output: &Path,
) -> Result<(), Box<dyn Error>> {
    let total_started = Instant::now();
    let rss_start = resident_memory_bytes();
    let reader = SigPackReader::open(data)?;
    let records = reader.records();
    let metadata = load_reaction_metadata(metadata_path)?;
    if metadata.len() != records.len() {
        return Err(format!(
            "metadata has {} rows but sigpack contains {} records",
            metadata.len(),
            records.len()
        )
        .into());
    }
    let model_contents = fs::read_to_string(model_path)?;
    let model: ScientificFitReport = serde_json::from_str(&model_contents)?;
    let frozen_bytes = fs::read(predictions_path)?;
    let mut frozen_reader = csv::ReaderBuilder::new()
        .trim(csv::Trim::All)
        .from_reader(frozen_bytes.as_slice());
    let frozen = frozen_reader
        .deserialize::<FrozenPredictionRow>()
        .collect::<Result<Vec<_>, _>>()?;
    if frozen.is_empty() {
        return Err("frozen prediction file contains no rows".into());
    }

    let predictor = RegressXPredictor::new(model.weights);
    let mut scored = Vec::with_capacity(frozen.len());
    let mut actual = Vec::with_capacity(frozen.len());
    let mut predicted = Vec::with_capacity(frozen.len());
    for frozen_row in &frozen {
        let index = metadata
            .iter()
            .position(|row| row.reaction_id == frozen_row.reaction_id)
            .ok_or_else(|| {
                format!(
                    "frozen reaction {} is absent from metadata",
                    frozen_row.reaction_id
                )
            })?;
        if metadata[index].dataset_split.eq_ignore_ascii_case("train") {
            return Err(format!(
                "frozen reaction {} is marked as training data",
                frozen_row.reaction_id
            )
            .into());
        }
        let recomputed = predictor.predict(&records[index]);
        if (recomputed - frozen_row.predicted_ddg).abs() > 1.0e-4 {
            return Err(format!(
                "frozen prediction for {} does not match the supplied model",
                frozen_row.reaction_id
            )
            .into());
        }
        let experimental = records[index].exp_ddg;
        if !experimental.is_finite() {
            return Err(format!(
                "reaction {} has no finite revealed target",
                frozen_row.reaction_id
            )
            .into());
        }
        actual.push(f64::from(experimental));
        predicted.push(f64::from(frozen_row.predicted_ddg));
        scored.push(ScoredPrediction {
            reaction_id: frozen_row.reaction_id.clone(),
            ligand_group: frozen_row.ligand_group.clone(),
            dataset_split: frozen_row.dataset_split.clone(),
            predicted_ddg_kcal_mol: frozen_row.predicted_ddg,
            experimental_ddg_kcal_mol: experimental,
            residual_kcal_mol: frozen_row.predicted_ddg - experimental,
            applicability_domain: frozen_row.applicability_domain.clone(),
        });
    }

    let residual_sum = actual
        .iter()
        .zip(&predicted)
        .map(|(actual, predicted)| (actual - predicted).powi(2))
        .sum::<f64>();
    let mean = actual.iter().sum::<f64>() / actual.len() as f64;
    let total_sum = actual
        .iter()
        .map(|actual| (actual - mean).powi(2))
        .sum::<f64>();
    let evaluation = EvaluationReport {
        schema_version: 1,
        model_path: model_path.display().to_string(),
        frozen_predictions_path: predictions_path.display().to_string(),
        frozen_predictions_fnv1a64: format!("{:016x}", fnv1a64(&frozen_bytes)),
        evaluated_records: actual.len(),
        mae_kcal_mol: actual
            .iter()
            .zip(&predicted)
            .map(|(actual, predicted)| (actual - predicted).abs())
            .sum::<f64>()
            / actual.len() as f64,
        rmse_kcal_mol: (residual_sum / actual.len() as f64).sqrt(),
        r2: (total_sum > f64::EPSILON).then_some(1.0 - residual_sum / total_sum),
        applicability_warnings: frozen
            .iter()
            .filter(|row| row.applicability_domain != "inside_training_range")
            .count(),
        scored_predictions: scored,
    };
    atomic_write_json(&evaluation, output)?;

    println!("command=evaluate");
    println!("evaluated_records={}", evaluation.evaluated_records);
    println!("mae_kcal_mol={:.8}", evaluation.mae_kcal_mol);
    println!("rmse_kcal_mol={:.8}", evaluation.rmse_kcal_mol);
    println!("r2={}", display_optional_metric(evaluation.r2));
    println!(
        "applicability_warnings={}",
        evaluation.applicability_warnings
    );
    println!(
        "frozen_predictions_fnv1a64={}",
        evaluation.frozen_predictions_fnv1a64
    );
    println!("evaluation_output={}", output.display());
    println!("total_ms={:.3}", millis(total_started.elapsed()));
    print_memory_metrics(rss_start, resident_memory_bytes());
    Ok(())
}

fn simulate_command(ddg_kcal: f32, temp_k: f32) -> Result<(), Box<dyn Error>> {
    if !ddg_kcal.is_finite() {
        return Err("--ddg must be finite".into());
    }
    if !temp_k.is_finite() || temp_k <= 0.0 {
        return Err("--temp must be a positive finite temperature".into());
    }

    let total_started = Instant::now();
    let rss_start = resident_memory_bytes();
    let rate = EyringKineticLink::calculate_rate_constant(ddg_kcal, temp_k);
    let (major_percent, minor_percent) =
        EyringKineticLink::calculate_enantiomeric_ratio(ddg_kcal, temp_k);
    let ee_percent = EyringKineticLink::calculate_enantiomeric_excess(ddg_kcal, temp_k);
    let distribution = EyringKineticLink::product_ratio(ddg_kcal, temp_k);
    let total_time = total_started.elapsed();

    println!("command=simulate");
    println!("ddg_kcal_mol={ddg_kcal:.7}");
    println!("temperature_k={temp_k:.2}");
    println!("rate_constant_s^-1={rate:.7e}");
    println!("major_enantiomer_percent={major_percent:.4}");
    println!("minor_enantiomer_percent={minor_percent:.4}");
    println!("percent_r={:.4}", distribution.percent_r);
    println!("percent_s={:.4}", distribution.percent_s);
    println!("ee_percent={ee_percent:.4}");
    println!("calculations_performed=3");
    println!("input_bytes={}", 2 * size_of::<f32>());
    println!("output_bytes={}", 6 * size_of::<f32>());
    println!("total_microseconds={:.3}", total_time.as_secs_f64() * 1e6);
    print_memory_metrics(rss_start, resident_memory_bytes());
    Ok(())
}

/// Covalent-radius tolerance for bond detection, matching the buried-volume kernel.
const DESCRIPTOR_BOND_TOLERANCE: f32 = 1.3;

/// The donor atom, its bonded substituents, and the substituent chosen as the
/// primary reference axis — everything the descriptor kernels need, recovered
/// from geometry alone.
#[derive(Clone, Copy, Debug)]
struct DonorTopology {
    donor_idx: usize,
    reference_idx: usize,
    substituents: [usize; 3],
}

/// Auto-detect the trivalent donor and its three bonded substituents.
///
/// The donor is the sole atom of `donor_element` unless `explicit_index` names
/// one. Substituents are the covalently bonded atoms (hydrogens included),
/// identical to the kernel's own frame construction, and the primary axis is the
/// nearest bonded heavy atom.
fn detect_donor(
    molecule: &Molecule,
    donor_element: &str,
    explicit_index: Option<usize>,
) -> Result<DonorTopology, String> {
    let donor_idx = match explicit_index {
        Some(index) => {
            if index >= molecule.atoms.len() {
                return Err(format!(
                    "--donor-index {index} is out of bounds for {} atoms",
                    molecule.atoms.len()
                ));
            }
            index
        }
        None => {
            let matches = molecule
                .atoms
                .iter()
                .enumerate()
                .filter(|(_, atom)| atom.element.eq_ignore_ascii_case(donor_element))
                .map(|(index, _)| index)
                .collect::<Vec<_>>();
            match matches.as_slice() {
                [single] => *single,
                [] => {
                    return Err(format!(
                        "no {donor_element} donor atom found; pass --donor-element or --donor-index"
                    ));
                }
                many => {
                    return Err(format!(
                        "found {} {donor_element} atoms; pass --donor-index to choose the donor",
                        many.len()
                    ));
                }
            }
        }
    };

    let donor = &molecule.atoms[donor_idx];
    let donor_covalent = covalent_radius(&donor.element);
    let mut bonded = molecule
        .atoms
        .iter()
        .enumerate()
        .filter(|(index, _)| *index != donor_idx)
        .filter_map(|(index, atom)| {
            let distance_squared = atom.position.distance_squared(donor.position);
            let cutoff =
                DESCRIPTOR_BOND_TOLERANCE * (donor_covalent + covalent_radius(&atom.element));
            (distance_squared <= cutoff * cutoff).then_some((distance_squared, index))
        })
        .collect::<Vec<_>>();
    if bonded.len() != 3 {
        return Err(format!(
            "{} donor (atom {donor_idx}) is not three-coordinate — found {} bonded atoms; \
             the descriptor is defined for trivalent (e.g. phosphine) donors",
            donor.element,
            bonded.len()
        ));
    }
    bonded.sort_by(|left, right| {
        left.0
            .total_cmp(&right.0)
            .then_with(|| left.1.cmp(&right.1))
    });
    let reference_idx = bonded
        .iter()
        .find(|(_, index)| !molecule.atoms[*index].element.eq_ignore_ascii_case("H"))
        .map(|(_, index)| *index)
        .ok_or_else(|| {
            format!(
                "{} donor (atom {donor_idx}) has no bonded heavy atom to define an axis",
                donor.element
            )
        })?;
    Ok(DonorTopology {
        donor_idx,
        reference_idx,
        substituents: [bonded[0].1, bonded[1].1, bonded[2].1],
    })
}

/// One file's descriptors, ready for text, JSON, or CSV emission.
#[derive(Debug, Serialize)]
struct DescriptorResult {
    file: String,
    conformers: usize,
    donor_element: String,
    donor_index: usize,
    substituents: Vec<String>,
    sterimol_l: f32,
    sterimol_b1: f32,
    sterimol_b5: f32,
    percent_buried_volume: f32,
    buried_volume: f32,
    qvbur_min: f32,
    qvbur_max: f32,
    max_delta_qvbur: f32,
    /// Kraken's headline `vbur_max_delta_qvbur_min`: the minimum over conformers.
    max_delta_qvbur_min: f32,
}

/// Compute ensemble-averaged descriptors for a single ligand file.
fn descriptors_for_file(
    path: &Path,
    donor_element: &str,
    donor_index: Option<usize>,
    sterimol_axis: SterimolAxis,
    config: BuriedVolumeConfig,
) -> Result<DescriptorResult, String> {
    let conformers = parse_coordinate_file(path)
        .map_err(|error| format!("failed to read {}: {error}", path.display()))?;
    if conformers.is_empty() {
        return Err(format!("{} contains no geometries", path.display()));
    }

    let topology = detect_donor(&conformers[0], donor_element, donor_index)?;
    let substituents = topology
        .substituents
        .iter()
        .map(|index| conformers[0].atoms[*index].element.clone())
        .collect::<Vec<_>>();

    let mut sterimol = Vec::with_capacity(conformers.len());
    let mut buried = Vec::with_capacity(conformers.len());
    for (conformer_index, molecule) in conformers.iter().enumerate() {
        // Re-detect per conformer so this tolerates files whose models differ.
        let topology = detect_donor(molecule, donor_element, donor_index).map_err(|message| {
            format!(
                "{} (conformer {conformer_index}): {message}",
                path.display()
            )
        })?;
        sterimol.push(match sterimol_axis {
            SterimolAxis::Bond => {
                SterimolCalculator::compute(molecule, topology.donor_idx, topology.reference_idx)
            }
            SterimolAxis::Coordination => {
                let center = coordination_center(
                    molecule,
                    topology.donor_idx,
                    topology.reference_idx,
                    config,
                )
                .map_err(|error| {
                    format!("{} (conformer {conformer_index}): {error}", path.display())
                })?;
                let mut params =
                    SterimolCalculator::compute_with_dummy(molecule, topology.donor_idx, center);
                params.l += STERIMOL_L_CORRECTION;
                params
            }
        });
        buried.push(
            BuriedVolumeCalculator::compute(
                molecule,
                topology.donor_idx,
                topology.reference_idx,
                config,
            )
            .map_err(|error| {
                format!("{} (conformer {conformer_index}): {error}", path.display())
            })?,
        );
    }

    let count = buried.len() as f32;
    let sterimol_mean =
        |select: fn(&SterimolParams) -> f32| sterimol.iter().map(select).sum::<f32>() / count;
    let buried_mean =
        |select: fn(&BuriedVolumeParams) -> f32| buried.iter().map(select).sum::<f32>() / count;
    let max_delta_qvbur_min = buried
        .iter()
        .map(|params| params.max_delta_qvbur)
        .fold(f32::INFINITY, f32::min);

    Ok(DescriptorResult {
        file: path.display().to_string(),
        conformers: conformers.len(),
        donor_element: conformers[0].atoms[topology.donor_idx].element.clone(),
        donor_index: topology.donor_idx,
        substituents,
        sterimol_l: sterimol_mean(|params| params.l),
        sterimol_b1: sterimol_mean(|params| params.b1),
        sterimol_b5: sterimol_mean(|params| params.b5),
        percent_buried_volume: buried_mean(|params| params.percent_buried_volume),
        buried_volume: buried_mean(|params| params.buried_volume),
        qvbur_min: buried_mean(|params| params.qvbur_min),
        qvbur_max: buried_mean(|params| params.qvbur_max),
        max_delta_qvbur: buried_mean(|params| params.max_delta_qvbur),
        max_delta_qvbur_min,
    })
}

fn print_descriptor_text(result: &DescriptorResult) {
    let ensemble = result.conformers > 1;
    let qualifier = if ensemble { " (conformer mean)" } else { "" };
    println!("{}", result.file);
    println!(
        "  donor          {} (atom {})",
        result.donor_element, result.donor_index
    );
    println!("  substituents   {}", result.substituents.join(", "));
    println!("  conformers     {}", result.conformers);
    println!(
        "  Sterimol{}      L {:.2}   B1 {:.2}   B5 {:.2}   Å",
        qualifier, result.sterimol_l, result.sterimol_b1, result.sterimol_b5
    );
    println!(
        "  buried volume{}  Vbur {:.1}%   ({:.1} Å³)",
        qualifier, result.percent_buried_volume, result.buried_volume
    );
    println!(
        "                 qvbur_min {:.2}   qvbur_max {:.2}   max_delta_qvbur {:.2}   Å³",
        result.qvbur_min, result.qvbur_max, result.max_delta_qvbur
    );
    if ensemble {
        println!(
            "                 max_delta_qvbur_min {:.2} Å³   [Kraken vbur_max_delta_qvbur_min]",
            result.max_delta_qvbur_min
        );
    }
}

fn print_descriptor_csv(results: &[DescriptorResult]) {
    println!(
        "file,conformers,donor_element,donor_index,substituents,sterimol_l,sterimol_b1,\
         sterimol_b5,percent_buried_volume,buried_volume,qvbur_min,qvbur_max,max_delta_qvbur,\
         max_delta_qvbur_min"
    );
    for result in results {
        println!(
            "{},{},{},{},{},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4}",
            csv_field(&result.file),
            result.conformers,
            result.donor_element,
            result.donor_index,
            csv_field(&result.substituents.join(" ")),
            result.sterimol_l,
            result.sterimol_b1,
            result.sterimol_b5,
            result.percent_buried_volume,
            result.buried_volume,
            result.qvbur_min,
            result.qvbur_max,
            result.max_delta_qvbur,
            result.max_delta_qvbur_min,
        );
    }
}

/// Quote a CSV field only when it contains a comma, quote, or newline.
fn csv_field(value: &str) -> String {
    if value.contains([',', '"', '\n']) {
        format!("\"{}\"", value.replace('"', "\"\""))
    } else {
        value.to_owned()
    }
}

fn descriptors_command(
    inputs: &[PathBuf],
    donor_element: &str,
    donor_index: Option<usize>,
    sterimol_axis: SterimolAxis,
    format: DescriptorFormat,
    config: BuriedVolumeConfig,
) -> Result<(), Box<dyn Error>> {
    if inputs.len() > 1 && donor_index.is_some() {
        return Err("--donor-index applies to a single file; omit it for batch runs".into());
    }
    let mut results = Vec::with_capacity(inputs.len());
    let mut failures = 0_usize;
    for path in inputs {
        match descriptors_for_file(path, donor_element, donor_index, sterimol_axis, config) {
            Ok(result) => results.push(result),
            Err(message) => {
                eprintln!("skipped {}: {message}", path.display());
                failures += 1;
            }
        }
    }
    if results.is_empty() {
        return Err("no ligand files could be featurized".into());
    }

    match format {
        DescriptorFormat::Text => {
            for (index, result) in results.iter().enumerate() {
                if index > 0 {
                    println!();
                }
                print_descriptor_text(result);
            }
        }
        DescriptorFormat::Json => {
            println!("{}", serde_json::to_string_pretty(&results)?);
        }
        DescriptorFormat::Csv => print_descriptor_csv(&results),
    }
    if failures > 0 {
        eprintln!(
            "featurized {} of {} files ({failures} skipped)",
            results.len(),
            inputs.len()
        );
    }
    Ok(())
}

fn sterimol_from_molecule(
    molecule: &Molecule,
    row: &ReactionCsvRow,
) -> Result<SterimolParams, String> {
    if molecule.atoms.len() < 2 {
        return Err("at least two atoms are required to establish the primary axis".into());
    }
    if row.attach_atom_idx == row.primary_bond_vector_idx {
        return Err("attachment and primary-vector indices must differ".into());
    }
    let attachment = molecule.atoms.get(row.attach_atom_idx).ok_or_else(|| {
        format!(
            "attachment atom index {} is out of bounds for {} atoms",
            row.attach_atom_idx,
            molecule.atoms.len()
        )
    })?;
    let neighbor = molecule
        .atoms
        .get(row.primary_bond_vector_idx)
        .ok_or_else(|| {
            format!(
                "primary-vector atom index {} is out of bounds for {} atoms",
                row.primary_bond_vector_idx,
                molecule.atoms.len()
            )
        })?;
    if (neighbor.position - attachment.position).length_squared() <= f32::EPSILON {
        return Err("attachment and primary-vector atoms occupy the same position".into());
    }

    Ok(SterimolCalculator::compute(
        molecule,
        row.attach_atom_idx,
        row.primary_bond_vector_idx,
    ))
}

fn record_from_ensemble(
    conformers: &[SterimolParams],
    weights: &[f32],
    energy_span: f32,
    row: &ReactionCsvRow,
) -> Result<PackedReactionRecord, String> {
    if conformers.is_empty() || conformers.len() != weights.len() {
        return Err("conformer parameters and weights must have equal non-zero lengths".into());
    }
    let weighted = conformers.iter().zip(weights).fold(
        SterimolParams::default(),
        |mut total, (params, weight)| {
            total.l += params.l * weight;
            total.b1 += params.b1 * weight;
            total.b5 += params.b5 * weight;
            total
        },
    );
    let minimum = conformers.iter().fold(
        SterimolParams {
            l: f32::INFINITY,
            b1: f32::INFINITY,
            b5: f32::INFINITY,
        },
        |minimum, params| SterimolParams {
            l: minimum.l.min(params.l),
            b1: minimum.b1.min(params.b1),
            b5: minimum.b5.min(params.b5),
        },
    );
    let maximum = conformers
        .iter()
        .fold(SterimolParams::default(), |maximum, params| {
            SterimolParams {
                l: maximum.l.max(params.l),
                b1: maximum.b1.max(params.b1),
                b5: maximum.b5.max(params.b5),
            }
        });
    if [
        weighted.l,
        weighted.b1,
        weighted.b5,
        minimum.l,
        minimum.b1,
        minimum.b5,
        maximum.l,
        maximum.b1,
        maximum.b5,
        energy_span,
    ]
    .iter()
    .any(|value| !value.is_finite())
    {
        return Err("ensemble aggregation produced a non-finite descriptor".into());
    }

    let mut record = PackedReactionRecord::from_ensemble(
        weighted.l,
        weighted.b1,
        weighted.b5,
        minimum.l,
        maximum.l,
        minimum.b1,
        maximum.b1,
        minimum.b5,
        maximum.b5,
        conformers.len(),
        energy_span,
    );
    record.nbo_charge = row.nbo_charge;
    record.ir_freq = row.ir_freq;
    record.temp_k = row.temp_k;
    record.exp_ddg = row.exp_ddg;
    Ok(record)
}

fn conformer_paths(row: &ReactionCsvRow) -> Result<Vec<PathBuf>, String> {
    let paths = row
        .conformer_xyz_paths
        .as_deref()
        .map(|value| {
            value
                .split(';')
                .map(str::trim)
                .filter(|path| !path.is_empty())
                .map(PathBuf::from)
                .collect::<Vec<_>>()
        })
        .filter(|paths| !paths.is_empty())
        .unwrap_or_else(|| vec![row.ligand_xyz_path.clone()]);
    for path in &paths {
        if path.as_os_str().is_empty()
            || path
                .components()
                .any(|component| component == Component::ParentDir)
        {
            return Err(format!(
                "invalid conformer coordinate path `{}`",
                path.display()
            ));
        }
    }
    Ok(paths)
}

fn parse_semicolon_floats(value: &str, label: &str) -> Result<Vec<f32>, String> {
    value
        .split(';')
        .map(str::trim)
        .filter(|field| !field.is_empty())
        .map(|field| {
            field
                .parse::<f32>()
                .map_err(|_| format!("{label} contains invalid float `{field}`"))
        })
        .collect()
}

fn conformer_weights(row: &ReactionCsvRow, conformer_count: usize) -> Result<Vec<f32>, String> {
    let mut weights = match row.conformer_boltzmann_weights.as_deref() {
        Some(value) if !value.trim().is_empty() => {
            parse_semicolon_floats(value, "Conformer_Boltzmann_Weights")?
        }
        _ => vec![1.0 / conformer_count as f32; conformer_count],
    };
    if weights.len() != conformer_count {
        return Err(format!(
            "found {} conformer weights for {conformer_count} coordinate paths",
            weights.len()
        ));
    }
    if weights
        .iter()
        .any(|weight| !weight.is_finite() || *weight < 0.0)
    {
        return Err("conformer weights must be finite and non-negative".into());
    }
    let total = weights.iter().sum::<f32>();
    if !total.is_finite() || total <= f32::EPSILON {
        return Err("conformer weights must have a positive sum".into());
    }
    weights.iter_mut().for_each(|weight| *weight /= total);
    Ok(weights)
}

fn conformer_coordination_centers(
    row: &ReactionCsvRow,
    conformer_count: usize,
) -> Result<Option<Vec<Vec3>>, String> {
    let Some(value) = row.conformer_coordination_centers.as_deref() else {
        return Ok(None);
    };
    if value.trim().is_empty() {
        return Ok(None);
    }
    let centers = value
        .split(';')
        .map(str::trim)
        .filter(|field| !field.is_empty())
        .map(|field| {
            let components = field
                .split(',')
                .map(str::trim)
                .map(|component| {
                    component.parse::<f32>().map_err(|_| {
                        format!(
                            "Conformer_Coordination_Centers_Angstrom contains invalid float `{component}`"
                        )
                    })
                })
                .collect::<Result<Vec<_>, _>>()?;
            if components.len() != 3 {
                return Err(format!(
                    "coordination center `{field}` must contain exactly x,y,z"
                ));
            }
            let center = Vec3::new(components[0], components[1], components[2]);
            if !center.is_finite() {
                return Err("coordination centers must be finite".to_owned());
            }
            Ok(center)
        })
        .collect::<Result<Vec<_>, String>>()?;
    if centers.len() != conformer_count {
        return Err(format!(
            "found {} coordination centers for {conformer_count} coordinate paths",
            centers.len()
        ));
    }
    Ok(Some(centers))
}

fn conformer_energy_span(row: &ReactionCsvRow, conformer_count: usize) -> Result<f32, String> {
    let Some(value) = row.conformer_relative_energies.as_deref() else {
        return Ok(0.0);
    };
    if value.trim().is_empty() {
        return Ok(0.0);
    }
    let energies = parse_semicolon_floats(value, "Conformer_Relative_Energies_kcal_mol")?;
    if energies.len() != conformer_count {
        return Err(format!(
            "found {} conformer energies for {conformer_count} coordinate paths",
            energies.len()
        ));
    }
    if energies
        .iter()
        .any(|energy| !energy.is_finite() || *energy < 0.0)
    {
        return Err("relative conformer energies must be finite and non-negative".into());
    }
    Ok(energies.into_iter().fold(0.0_f32, f32::max))
}

fn validate_reaction_row(row: &ReactionCsvRow, csv_line: usize) -> Result<(), Box<dyn Error>> {
    if row.reaction_id.trim().is_empty() {
        return Err(format!("CSV line {csv_line} has an empty Reaction_ID").into());
    }
    if row.ligand_xyz_path.as_os_str().is_empty() {
        return Err(format!("CSV line {csv_line} has an empty Ligand_XYZ_Path").into());
    }
    if row
        .ligand_xyz_path
        .components()
        .any(|component| component == Component::ParentDir)
    {
        return Err(format!("CSV line {csv_line} contains a parent-directory XYZ path").into());
    }
    if !row.nbo_charge.is_finite()
        || !row.ir_freq.is_finite()
        || !row.temp_k.is_finite()
        || !row.exp_ddg.is_finite()
    {
        return Err(format!("CSV line {csv_line} contains a non-finite numeric value").into());
    }
    if row.temp_k <= 0.0 {
        return Err(format!("CSV line {csv_line} has a non-positive temperature").into());
    }
    Ok(())
}

fn resolve_xyz_path(xyz_dir: &Path, csv_path: &Path) -> Result<PathBuf, String> {
    if csv_path.is_absolute() {
        return csv_path
            .is_file()
            .then(|| csv_path.to_owned())
            .ok_or_else(|| format!("XYZ file does not exist: {}", csv_path.display()));
    }

    let direct = xyz_dir.join(csv_path);
    if direct.is_file() {
        return Ok(direct);
    }

    if let (Some(directory_name), Some(first_component)) = (
        xyz_dir.file_name(),
        csv_path
            .components()
            .next()
            .and_then(|component| match component {
                Component::Normal(value) => Some(value),
                _ => None,
            }),
    ) && directory_name == first_component
    {
        let stripped: PathBuf = csv_path.components().skip(1).collect();
        let candidate = xyz_dir.join(stripped);
        if candidate.is_file() {
            return Ok(candidate);
        }
    }

    if let Some(filename) = csv_path.file_name() {
        let candidate = xyz_dir.join(filename);
        if candidate.is_file() {
            return Ok(candidate);
        }
    }
    Err(format!(
        "could not resolve XYZ path `{}` below {}",
        csv_path.display(),
        xyz_dir.display()
    ))
}

fn load_weights_json(path: &Path) -> Result<[f32; MODEL_FEATURE_COUNT], Box<dyn Error>> {
    let contents = fs::read_to_string(path)?;
    let document: WeightDocument = serde_json::from_str(&contents).map_err(|error| {
        format!(
            "weight file {} must be an eight-number JSON array or an object with an eight-number `weights` array: {error}",
            path.display()
        )
    })?;
    let weights = match document {
        WeightDocument::Array(weights) | WeightDocument::Object { weights } => weights,
    };
    if weights.iter().any(|weight| !weight.is_finite()) {
        return Err(format!("weight file {} contains NaN or infinity", path.display()).into());
    }
    Ok(weights)
}

fn load_reaction_metadata(path: &Path) -> Result<Vec<ReactionMetadataRow>, Box<dyn Error>> {
    let mut reader = csv::ReaderBuilder::new()
        .trim(csv::Trim::All)
        .from_path(path)?;
    let rows = reader
        .deserialize::<ReactionMetadataRow>()
        .enumerate()
        .map(|(index, row)| {
            let row = row.map_err(|error| {
                format!(
                    "{} row {} could not be parsed: {error}",
                    path.display(),
                    index + 2
                )
            })?;
            if row.reaction_id.trim().is_empty() {
                return Err(format!(
                    "{} row {} has no Reaction_ID",
                    path.display(),
                    index + 2
                ));
            }
            if !matches!(
                row.dataset_split.to_ascii_lowercase().as_str(),
                "train" | "external" | "blind" | "test"
            ) {
                return Err(format!(
                    "{} row {} has unsupported Dataset_Split `{}`",
                    path.display(),
                    index + 2,
                    row.dataset_split
                ));
            }
            Ok(row)
        })
        .collect::<Result<Vec<_>, String>>()?;
    if rows.is_empty() {
        return Err(format!("metadata CSV contains no rows: {}", path.display()).into());
    }
    Ok(rows)
}

fn applicability_status(record: &PackedReactionRecord, report: &ScientificFitReport) -> String {
    let features = expand_features(record);
    let outside = report
        .selected_feature_indices
        .iter()
        .zip(&report.applicability_domain)
        .filter(|(column, domain)| {
            let value = f64::from(features[**column]);
            value < domain.minimum || value > domain.maximum
        })
        .map(|(_, domain)| domain.feature.as_str())
        .collect::<Vec<_>>();
    if outside.is_empty() {
        "inside_training_range".into()
    } else {
        format!("outside:{}", outside.join("|"))
    }
}

fn temporary_sibling(path: &Path) -> PathBuf {
    let mut value = path.as_os_str().to_os_string();
    value.push(".tmp");
    PathBuf::from(value)
}

fn ensure_parent(path: &Path) -> std::io::Result<()> {
    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)?;
    }
    Ok(())
}

fn atomic_write_json<T: Serialize>(value: &T, path: &Path) -> Result<(), Box<dyn Error>> {
    ensure_parent(path)?;
    let temporary = temporary_sibling(path);
    let json = serde_json::to_string_pretty(value)?;
    fs::write(&temporary, format!("{json}\n"))?;
    fs::rename(temporary, path)?;
    Ok(())
}

fn atomic_write_csv_rows<T: Serialize>(rows: &[T], path: &Path) -> Result<(), Box<dyn Error>> {
    ensure_parent(path)?;
    let temporary = temporary_sibling(path);
    {
        let mut writer = csv::Writer::from_path(&temporary)?;
        for row in rows {
            writer.serialize(row)?;
        }
        writer.flush()?;
    }
    fs::rename(temporary, path)?;
    Ok(())
}

fn display_optional_metric(value: Option<f64>) -> String {
    value.map_or_else(|| "unavailable".into(), |value| format!("{value:.8}"))
}

fn fnv1a64(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325, |hash, byte| {
        (hash ^ u64::from(*byte)).wrapping_mul(0x0000_0100_0000_01b3)
    })
}

fn resident_memory_bytes() -> Option<u64> {
    let status = fs::read_to_string("/proc/self/status").ok()?;
    let line = status.lines().find(|line| line.starts_with("VmRSS:"))?;
    let kibibytes = line.split_whitespace().nth(1)?.parse::<u64>().ok()?;
    kibibytes.checked_mul(1_024)
}

fn print_memory_metrics(start: Option<u64>, end: Option<u64>) {
    match start {
        Some(bytes) => println!("rss_start_bytes={bytes}"),
        None => println!("rss_start_bytes=unavailable"),
    }
    match end {
        Some(bytes) => println!("rss_end_bytes={bytes}"),
        None => println!("rss_end_bytes=unavailable"),
    }
    match (start, end) {
        (Some(start), Some(end)) => {
            println!("rss_delta_bytes={}", i128::from(end) - i128::from(start));
        }
        _ => println!("rss_delta_bytes=unavailable"),
    }
}

fn millis(duration: Duration) -> f64 {
    duration.as_secs_f64() * 1_000.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::Vec3;
    use std::time::{SystemTime, UNIX_EPOCH};
    use steric_x::Atom;

    fn temporary_directory(test_name: &str) -> PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!(
            "steric_x_cli_{test_name}_{}_{nonce}",
            std::process::id()
        ))
    }

    fn sample_row() -> ReactionCsvRow {
        ReactionCsvRow {
            reaction_id: "RXN-001".into(),
            ligand_xyz_path: "ligand.xyz".into(),
            attach_atom_idx: 0,
            primary_bond_vector_idx: 1,
            nbo_charge: -0.3,
            ir_freq: 1_650.0,
            temp_k: 298.15,
            exp_ddg: 1.2,
            conformer_xyz_paths: None,
            conformer_relative_energies: None,
            conformer_boltzmann_weights: None,
            conformer_coordination_centers: None,
            coordination_center_method: None,
        }
    }

    #[test]
    fn creates_record_from_csv_row() {
        let molecule = Molecule {
            atoms: vec![
                Atom::new("H", Vec3::ZERO),
                Atom::new("C", Vec3::new(2.0, 0.0, 0.0)),
            ],
        };
        let row = sample_row();
        let params = sterimol_from_molecule(&molecule, &row).unwrap();
        let record = record_from_ensemble(&[params], &[1.0], 0.0, &row).unwrap();
        assert!((record.l - 3.7).abs() < 1.0e-5);
        assert_eq!(record.nbo_charge, -0.3);
        assert_eq!(record.exp_ddg, 1.2);
        assert_eq!(record.conformer_count(), 1);
    }

    #[test]
    fn parse_command_packs_csv_and_xyz_geometry() {
        let directory = temporary_directory("parse");
        let xyz_dir = directory.join("xyz");
        fs::create_dir_all(&xyz_dir).unwrap();
        fs::write(xyz_dir.join("ligand.xyz"), "2\nligand\nH 0 0 0\nC 2 0 0\n").unwrap();
        let csv_path = directory.join("reactions.csv");
        fs::write(
            &csv_path,
            "Reaction_ID,Ligand_XYZ_Path,Attach_Atom_Idx,Primary_Bond_Vector_Idx,NBO_Charge,IR_Frequency,Temp_K,Exp_ddG_kcal_mol\n\
             RXN-001,xyz/ligand.xyz,0,1,-0.25,1650,298.15,1.5\n",
        )
        .unwrap();
        let output = directory.join("dataset.sigpack");

        parse_command(&csv_path, &xyz_dir, &output).unwrap();
        let reader = SigPackReader::open(&output).unwrap();
        assert_eq!(reader.records().len(), 1);
        assert_eq!(reader.records()[0].nbo_charge, -0.25);
        assert_eq!(reader.records()[0].exp_ddg, 1.5);
        drop(reader);
        fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn parse_command_boltzmann_averages_conformer_envelopes() {
        let directory = temporary_directory("ensemble");
        let xyz_dir = directory.join("xyz");
        fs::create_dir_all(&xyz_dir).unwrap();
        fs::write(
            xyz_dir.join("compact.xyz"),
            "3\ncompact\nH 0 0 0\nH 0 0 1\nC 0 0 1\n",
        )
        .unwrap();
        fs::write(
            xyz_dir.join("wide.xyz"),
            "3\nwide\nH 0 0 0\nH 0 0 1\nC 2 0 1\n",
        )
        .unwrap();
        let csv_path = directory.join("reactions.csv");
        fs::write(
            &csv_path,
            "Reaction_ID,Ligand_XYZ_Path,Attach_Atom_Idx,Primary_Bond_Vector_Idx,NBO_Charge,IR_Frequency,Temp_K,Exp_ddG_kcal_mol,Conformer_XYZ_Paths,Conformer_Relative_Energies_kcal_mol,Conformer_Boltzmann_Weights\n\
             RXN-E,xyz/compact.xyz,0,1,-0.25,1650,298.15,1.5,xyz/compact.xyz;xyz/wide.xyz,0;1,0.25;0.75\n",
        )
        .unwrap();
        let output = directory.join("ensemble.sigpack");

        parse_command(&csv_path, &xyz_dir, &output).unwrap();
        let reader = SigPackReader::open(&output).unwrap();
        let record = reader.records()[0];
        assert_eq!(record.conformer_count(), 2);
        assert!((record.b5 - 3.2).abs() < 1.0e-5);
        assert!((record.b5_min() - 1.7).abs() < 1.0e-5);
        assert!((record.b5_max() - 3.7).abs() < 1.0e-5);
        assert_eq!(record.ensemble_energy_span(), 1.0);
        drop(reader);
        fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn buried_volume_command_writes_version_two_and_audit_csv() {
        let directory = temporary_directory("buried_volume");
        let xyz_dir = directory.join("xyz");
        fs::create_dir_all(&xyz_dir).unwrap();
        fs::write(
            xyz_dir.join("phosphine.xyz"),
            "5\nphosphine\n\
             P 0 0 0\n\
             C 1.4 0 0.45\n\
             C -0.7 1.212 0.45\n\
             C -0.7 -1.212 0.45\n\
             C 2.8 0 0.7\n",
        )
        .unwrap();
        let csv_path = directory.join("reactions.csv");
        fs::write(
            &csv_path,
            "Reaction_ID,Ligand_XYZ_Path,Attach_Atom_Idx,Primary_Bond_Vector_Idx,NBO_Charge,IR_Frequency,Temp_K,Exp_ddG_kcal_mol,Conformer_Coordination_Centers_Angstrom,Coordination_Center_Method\n\
             RXN-BV,xyz/phosphine.xyz,0,1,0.8,1650,298.15,1.2,\"0,0,-2.1\",xtb_lmo_kraken\n",
        )
        .unwrap();
        let output = directory.join("dataset_v2.sigpack");
        let audit = directory.join("conformers.csv");

        buried_volume_command(
            &csv_path,
            &xyz_dir,
            &output,
            Some(&audit),
            BuriedVolumeConfig::default(),
            true,
        )
        .unwrap();
        let reader = steric_x::SigPackV2Reader::open(&output).unwrap();
        assert_eq!(reader.len(), 1);
        assert!(reader.records()[0].buried_volume.max_delta_qvbur_boltz > 0.0);
        assert_eq!(reader.records()[0].buried_volume.conformer_count, 1.0);
        assert!(
            fs::read_to_string(&audit)
                .unwrap()
                .contains("max_delta_qvbur")
        );
        assert!(
            fs::read_to_string(&audit)
                .unwrap()
                .contains("xtb_lmo_kraken")
        );
        drop(reader);
        fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn weight_loader_accepts_array_and_object_json() {
        let directory = temporary_directory("weights");
        fs::create_dir_all(&directory).unwrap();
        let array_path = directory.join("array.json");
        let object_path = directory.join("object.json");
        fs::write(&array_path, "[0,1,2,3,4,5,6,7]").unwrap();
        fs::write(&object_path, r#"{"weights":[7,6,5,4,3,2,1,0]}"#).unwrap();

        assert_eq!(load_weights_json(&array_path).unwrap()[7], 7.0);
        assert_eq!(load_weights_json(&object_path).unwrap()[0], 7.0);
        fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn clap_accepts_all_command_contracts() {
        assert!(
            Cli::try_parse_from([
                "stericx",
                "parse",
                "--csv",
                "reactions.csv",
                "--xyz-dir",
                "xyz",
                "--output",
                "data.sigpack",
            ])
            .is_ok()
        );
        assert!(
            Cli::try_parse_from([
                "stericx",
                "buried-volume",
                "--csv",
                "reactions.csv",
                "--xyz-dir",
                "xyz",
                "--output",
                "data_v2.sigpack",
                "--per-conformer-output",
                "conformers.csv",
            ])
            .is_ok()
        );
        assert!(
            Cli::try_parse_from([
                "stericx",
                "fit",
                "--data",
                "data.sigpack",
                "--metadata",
                "reactions.csv",
                "--output",
                "model.json",
                "--predictions",
                "frozen.csv",
            ])
            .is_ok()
        );
        assert!(
            Cli::try_parse_from([
                "stericx",
                "evaluate",
                "--data",
                "data.sigpack",
                "--metadata",
                "reactions.csv",
                "--model",
                "model.json",
                "--predictions",
                "frozen.csv",
                "--output",
                "evaluation.json",
            ])
            .is_ok()
        );
        assert!(
            Cli::try_parse_from([
                "stericx",
                "predict",
                "--data",
                "data.sigpack",
                "--weights",
                "weights.json",
            ])
            .is_ok()
        );
        assert!(
            Cli::try_parse_from(["stericx", "simulate", "--ddg", "1.82", "--temp", "298.15",])
                .is_ok()
        );
        assert!(
            Cli::try_parse_from([
                "stericx",
                "descriptors",
                "ligand_a.xyz",
                "ligand_b.sdf",
                "--sterimol-axis",
                "coordination",
                "--format",
                "csv",
            ])
            .is_ok()
        );
        // `descriptors` requires at least one input file.
        assert!(Cli::try_parse_from(["stericx", "descriptors"]).is_err());
    }

    fn molecule_from(atoms: &[(&str, f32, f32, f32)]) -> Molecule {
        Molecule {
            atoms: atoms
                .iter()
                .map(|(element, x, y, z)| Atom::new(*element, Vec3::new(*x, *y, *z)))
                .collect(),
        }
    }

    fn tertiary_phosphine() -> Molecule {
        molecule_from(&[
            ("P", 0.0, 0.0, 0.0),
            ("C", 1.4, 0.0, 0.45),
            ("C", -0.7, 1.212, 0.45),
            ("C", -0.7, -1.212, 0.45),
            // A distal, non-bonded carbon that must never be taken as a substituent.
            ("C", 2.8, 0.0, 0.7),
        ])
    }

    #[test]
    fn detect_donor_finds_the_single_phosphine() {
        let topology = detect_donor(&tertiary_phosphine(), "P", None).unwrap();
        assert_eq!(topology.donor_idx, 0);
        let mut substituents = topology.substituents;
        substituents.sort_unstable();
        assert_eq!(substituents, [1, 2, 3]);
        // The primary axis is a bonded heavy atom, never the distal carbon (4).
        assert!([1, 2, 3].contains(&topology.reference_idx));
    }

    #[test]
    fn detect_donor_counts_bonded_hydrogens_as_substituents() {
        // Primary phosphine R-PH2: the two bonded hydrogens are substituents,
        // and the reference axis falls on the lone bonded carbon.
        let molecule = molecule_from(&[
            ("P", 0.0, 0.0, 0.0),
            ("C", 1.5, 0.0, 0.6),
            ("H", -0.6, 1.1, 0.55),
            ("H", -0.6, -1.1, 0.55),
        ]);
        let topology = detect_donor(&molecule, "P", None).unwrap();
        let mut substituents = topology.substituents;
        substituents.sort_unstable();
        assert_eq!(substituents, [1, 2, 3]);
        assert_eq!(topology.reference_idx, 1);
    }

    #[test]
    fn detect_donor_reports_missing_and_ambiguous_donors() {
        let no_phosphorus = molecule_from(&[("C", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 1.1)]);
        let error = detect_donor(&no_phosphorus, "P", None).unwrap_err();
        assert!(error.contains("no P donor"));

        let two_phosphorus = molecule_from(&[
            ("P", 0.0, 0.0, 0.0),
            ("C", 1.5, 0.0, 0.0),
            ("C", -1.5, 0.0, 0.0),
            ("C", 0.0, 1.5, 0.0),
            ("P", 6.0, 0.0, 0.0),
            ("C", 7.5, 0.0, 0.0),
            ("C", 4.5, 0.0, 0.0),
            ("C", 6.0, 1.5, 0.0),
        ]);
        let error = detect_donor(&two_phosphorus, "P", None).unwrap_err();
        assert!(error.contains("--donor-index"));
        // Naming one of them resolves the ambiguity.
        assert_eq!(
            detect_donor(&two_phosphorus, "P", Some(4))
                .unwrap()
                .donor_idx,
            4
        );
    }

    #[test]
    fn detect_donor_rejects_a_non_trivalent_donor() {
        let two_coordinate = molecule_from(&[
            ("P", 0.0, 0.0, 0.0),
            ("C", 1.5, 0.0, 0.0),
            ("H", -0.9, 0.9, 0.0),
        ]);
        let error = detect_donor(&two_coordinate, "P", None).unwrap_err();
        assert!(error.contains("three-coordinate"));
    }

    #[test]
    fn detect_donor_supports_a_non_phosphorus_element() {
        // A trivalent amine donor located by --donor-element N.
        let amine = molecule_from(&[
            ("N", 0.0, 0.0, 0.0),
            ("C", 1.47, 0.0, 0.3),
            ("C", -0.7, 1.2, 0.3),
            ("C", -0.7, -1.2, 0.3),
        ]);
        let topology = detect_donor(&amine, "N", None).unwrap();
        assert_eq!(topology.donor_idx, 0);
    }

    #[test]
    fn descriptors_for_file_featurizes_an_xyz_geometry() {
        let directory = temporary_directory("descriptors_xyz");
        std::fs::create_dir_all(&directory).unwrap();
        let path = directory.join("phosphine.xyz");
        std::fs::write(
            &path,
            "5\nphosphine\nP 0 0 0\nC 1.4 0 0.45\nC -0.7 1.212 0.45\n\
             C -0.7 -1.212 0.45\nC 2.8 0 0.7\n",
        )
        .unwrap();
        let result = descriptors_for_file(
            &path,
            "P",
            None,
            SterimolAxis::Bond,
            BuriedVolumeConfig::default(),
        )
        .unwrap();
        assert_eq!(result.donor_element, "P");
        assert_eq!(result.conformers, 1);
        assert_eq!(result.substituents.len(), 3);
        assert!(result.buried_volume.is_finite() && result.buried_volume > 0.0);
        assert!(result.max_delta_qvbur_min.is_finite());
        // With one conformer the ensemble minimum equals the single value.
        assert_eq!(result.max_delta_qvbur, result.max_delta_qvbur_min);

        // The coordination axis is a different frame, so it yields a different
        // (also valid) Sterimol L, and applies the +0.40 Å correction.
        let coordination = descriptors_for_file(
            &path,
            "P",
            None,
            SterimolAxis::Coordination,
            BuriedVolumeConfig::default(),
        )
        .unwrap();
        assert!(coordination.sterimol_l.is_finite() && coordination.sterimol_l > 0.0);
        assert!((coordination.sterimol_l - result.sterimol_l).abs() > 1.0e-4);
        std::fs::remove_dir_all(&directory).ok();
    }

    #[test]
    fn csv_field_quotes_only_when_required() {
        assert_eq!(csv_field("plain"), "plain");
        assert_eq!(csv_field("a,b"), "\"a,b\"");
        assert_eq!(csv_field("say \"hi\""), "\"say \"\"hi\"\"\"");
    }
}
