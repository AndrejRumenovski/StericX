//! Command-line surface: the `clap` argument types for every subcommand.

use clap::{Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(
    name = "stericx",
    version,
    about = "SIMD physical-organic featurization and selectivity prediction"
)]
pub(crate) struct Cli {
    #[command(subcommand)]
    pub(crate) command: Command,
}

#[derive(Debug, Subcommand)]
pub(crate) enum Command {
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
    /// Rank a ligand library by steric similarity to a query ligand.
    ///
    /// Featurizes the query, then ranks every library member by Euclidean
    /// distance in *standardized* descriptor space (each descriptor z-scored
    /// against the library, so Å-scale and percent-scale features count
    /// comparably). Constraints narrow the field before ranking, so questions
    /// like "a less bulky ligand of similar shape" or "%Vbur between 30 and 35
    /// with B5 under 7 Å" are one command. The ranking is steric similarity,
    /// not a prediction of reactivity.
    Search {
        /// Query ligand geometry (`.xyz`, `.sdf`, `.mol`).
        #[arg(long, value_name = "FILE")]
        ligand: PathBuf,
        /// Library to search: a directory of geometries, or a CSV previously
        /// written by `stericx descriptors --format csv`.
        #[arg(long, value_name = "DIR|CSV")]
        library: PathBuf,
        /// Number of hits to report.
        #[arg(long, default_value_t = 10)]
        top: usize,
        /// Comma-separated descriptors defining similarity (default: shape
        /// envelope, %Vbur, quadrant asymmetry, and pyramidalization).
        #[arg(long)]
        features: Option<String>,
        /// Constraint such as `vbur=30..35`, `b5<7`, or `l>=8`. Repeatable.
        #[arg(long = "filter", value_name = "EXPR")]
        filters: Vec<String>,
        /// Keep only candidates less buried than the query.
        #[arg(long)]
        less_bulky: bool,
        /// Keep only candidates more buried than the query.
        #[arg(long)]
        more_bulky: bool,
        /// Donor element to locate (defaults to phosphorus).
        #[arg(long, default_value = "P")]
        donor_element: String,
        /// Sterimol axis used for both query and library.
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
        /// Donor-to-virtual-metal distance in ångströms.
        #[arg(long, default_value_t = 2.28)]
        center_distance: f32,
        /// Bondi radius scale factor used by Morfeus.
        #[arg(long, default_value_t = 1.17)]
        radii_scale: f32,
    },
}

/// Output format for the `descriptors` subcommand.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub(crate) enum DescriptorFormat {
    /// Human-readable summary, one block per file.
    Text,
    /// One JSON array of per-file descriptor records.
    Json,
    /// Comma-separated table, one row per file (spreadsheet-ready).
    Csv,
}

/// Which axis defines the Sterimol frame.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub(crate) enum SterimolAxis {
    /// Donor → nearest bonded substituent (a P–C bond direction).
    Bond,
    /// Virtual metal → donor, along the lone pair. This is the convention
    /// Kraken and the ligand-descriptor literature use, and it applies the
    /// historical +0.40 Å Verloop correction to `L`.
    Coordination,
}

/// Historical Verloop/Morfeus correction added to the raw geometric Sterimol
/// `L` under the coordination-axis convention.
pub(crate) const STERIMOL_L_CORRECTION: f32 = 0.40;

#[cfg(test)]
mod tests {
    use super::*;

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
}
