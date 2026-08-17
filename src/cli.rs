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

// clap derives one variant per subcommand with its flags as flat fields, so the
// variants are inherently uneven in size and cannot be boxed without giving up
// the derive. The enum is constructed once per process, so the size spread costs
// nothing measurable.
#[allow(clippy::large_enum_variant)]
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
        /// Query ligand geometry (`.xyz`, `.sdf`, `.mol`). Omit it to run a
        /// pure constraint query, where the filters alone select the ligands.
        #[arg(long = "similar-to", visible_alias = "ligand", value_name = "FILE")]
        ligand: Option<PathBuf>,
        /// Database to search: a `stericx db build` table, a descriptors CSV, or
        /// a directory of geometries to featurize on the fly.
        #[arg(long = "database", visible_alias = "library", value_name = "DB|DIR")]
        library: PathBuf,
        /// Number of hits to report.
        #[arg(long, default_value_t = 10)]
        top: usize,
        /// Descriptor to order a constraint-only query by (default: the first
        /// constrained descriptor). Ignored when --similar-to is given.
        #[arg(long, value_name = "DESCRIPTOR")]
        sort_by: Option<String>,
        /// Sort a constraint-only query descending.
        #[arg(long)]
        descending: bool,
        /// Comma-separated descriptors defining similarity (default: shape
        /// envelope, %Vbur, quadrant asymmetry, and pyramidalization).
        #[arg(long)]
        features: Option<String>,
        /// Constraint such as `vbur=30..35`, `b5<7`, or `l>=8`. Repeatable.
        #[arg(long = "filter", value_name = "EXPR")]
        filters: Vec<String>,
        /// Keep ligands whose %Vbur falls in this range, e.g. `30:35`.
        #[arg(long, value_name = "LOW:HIGH")]
        vbur: Option<String>,
        /// Keep ligands whose Sterimol L falls in this range.
        #[arg(long, value_name = "LOW:HIGH")]
        l: Option<String>,
        /// Keep ligands whose Sterimol B1 falls in this range.
        #[arg(long, value_name = "LOW:HIGH")]
        b1: Option<String>,
        /// Keep ligands whose Sterimol B5 falls in this range.
        #[arg(long, value_name = "LOW:HIGH")]
        b5: Option<String>,
        /// Minimum %Vbur (inclusive).
        #[arg(long)]
        vbur_min: Option<f32>,
        /// Maximum %Vbur (inclusive).
        #[arg(long)]
        vbur_max: Option<f32>,
        /// Minimum Sterimol L (inclusive).
        #[arg(long)]
        l_min: Option<f32>,
        /// Maximum Sterimol L (inclusive).
        #[arg(long)]
        l_max: Option<f32>,
        /// Minimum Sterimol B1 (inclusive).
        #[arg(long)]
        b1_min: Option<f32>,
        /// Maximum Sterimol B1 (inclusive).
        #[arg(long)]
        b1_max: Option<f32>,
        /// Minimum Sterimol B5 (inclusive).
        #[arg(long)]
        b5_min: Option<f32>,
        /// Maximum Sterimol B5 (inclusive).
        #[arg(long)]
        b5_max: Option<f32>,
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
    /// Rank a ligand library with a fitted reaction model.
    ///
    /// Reports predicted performance, a conservative uncertainty band from the
    /// model's bootstrap coefficient intervals, and an applicability-domain
    /// warning for every ligand that falls outside the range the model was
    /// trained on. The model decides what the library must supply: a model that
    /// selected an electronic term (`nbo_charge`, `ir_frequency`, or an
    /// interaction) cannot be screened from geometry alone, and `screen` says so
    /// rather than guessing the missing quantity.
    Screen {
        /// Fitted model JSON produced by `stericx fit`.
        #[arg(value_name = "MODEL")]
        model: PathBuf,
        /// Library: a directory of geometries, a descriptors CSV, or a reaction
        /// CSV carrying `NBO_Charge` / `IR_Frequency`.
        #[arg(value_name = "LIBRARY")]
        library: PathBuf,
        /// Report only the best N ligands (default: all).
        #[arg(long)]
        top: Option<usize>,
        /// Temperature used to convert predicted ΔΔG‡ into an ee.
        #[arg(long, default_value_t = 298.15)]
        temperature: f32,
        /// Drop ligands that fall outside the training domain.
        #[arg(long)]
        inside_domain_only: bool,
        /// Rank ascending (smallest predicted value first).
        #[arg(long)]
        ascending: bool,
        /// Donor element to locate when featurizing a geometry library.
        #[arg(long, default_value = "P")]
        donor_element: String,
        /// Sterimol axis used when featurizing a geometry library.
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
    /// Compare two or more ligands side by side.
    ///
    /// Prints every descriptor for each ligand, the spread across them, and —
    /// with `--database` — that spread in library standard deviations plus a
    /// standardized pairwise distance. The σ view is what makes a raw difference
    /// interpretable: 0.5 Å of B5 is small in a library spanning 4 Å and large in
    /// one spanning 0.6 Å.
    Compare {
        /// Ligand geometries to compare (two or more).
        #[arg(required = true, value_name = "FILE")]
        inputs: Vec<PathBuf>,
        /// Database supplying the scale for σ columns and distances.
        #[arg(long = "database", visible_alias = "library", value_name = "DB|DIR")]
        database: Option<PathBuf>,
        /// Donor element to locate (defaults to phosphorus).
        #[arg(long, default_value = "P")]
        donor_element: String,
        /// Sterimol axis.
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
    /// Build and inspect precomputed ligand descriptor databases.
    Db {
        #[command(subcommand)]
        action: DbCommand,
    },
}

/// Which path component becomes a database row's ligand label.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub(crate) enum DbLabel {
    /// The file stem, e.g. `49973.sdf` becomes `49973`.
    Stem,
    /// The parent directory, e.g. `1088/49973.sdf` becomes `1088` — the Kraken
    /// cache layout, where the directory is the ligand identifier.
    Parent,
    /// The path relative to the source root.
    Path,
}

#[derive(Debug, Subcommand)]
pub(crate) enum DbCommand {
    /// Featurize a directory of geometries into a reusable database.
    Build {
        /// Directory of ligand geometries to featurize (searched recursively).
        #[arg(long, value_name = "DIR")]
        source: PathBuf,
        /// Destination CSV; a `.manifest.json` sibling records provenance.
        #[arg(long, value_name = "CSV")]
        output: PathBuf,
        /// Treat each parent directory as one ligand and aggregate its
        /// conformers (`max_delta_qvbur_min` takes the minimum, per Kraken).
        #[arg(long)]
        group_by_parent: bool,
        /// Restrict to these geometry extensions (repeatable). Useful when a
        /// source tree mirrors the same conformers in more than one format.
        #[arg(long = "extension", value_name = "EXT")]
        extensions: Vec<String>,
        /// Which path component becomes the ligand label.
        #[arg(long, value_enum, default_value_t = DbLabel::Stem)]
        label_from: DbLabel,
        /// Donor element to locate (defaults to phosphorus).
        #[arg(long, default_value = "P")]
        donor_element: String,
        /// Sterimol axis recorded in the database.
        #[arg(long, value_enum, default_value_t = SterimolAxis::Bond)]
        sterimol_axis: SterimolAxis,
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
