//! Export ordered descriptors and a conformer sensitivity report using Rust.
use clap::Parser;
use std::{error::Error, fs, io::Write, path::PathBuf};
use steric_x::ensemble::{
    DescriptorOptions, EnergyMetadata, Ensemble, MissingConformerPolicy, StatisticsOptions,
};

#[derive(Parser)]
#[command(about = "Export conformer descriptors and sensitivity as JSON and CSV")]
struct Args {
    /// SDF, MOL or multi-frame XYZ input; original order is retained.
    input: PathBuf,
    /// Output stem; creates STEM.json and STEM.csv without overwriting files.
    output: PathBuf,
    #[arg(long)]
    ligand_id: Option<String>,
    /// JSON object containing DescriptorOptions; omitted fields use CLI defaults.
    #[arg(long)]
    options: Option<PathBuf>,
    /// Ordered JSON array of explicit EnergyMetadata objects or nulls.
    #[arg(long)]
    energies: Option<PathBuf>,
    /// Opt into energy weighting at this temperature in kelvin.
    #[arg(long, requires = "population_assumptions")]
    temperature: Option<f64>,
    /// State membership, degeneracy and thermochemical assumptions.
    #[arg(long, requires = "temperature")]
    population_assumptions: Option<String>,
    /// Explicitly condition statistics on successful conformers.
    #[arg(long)]
    available_cases: bool,
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    let mut ensemble = Ensemble::from_file(&args.input)?;
    let options: DescriptorOptions = match args.options {
        Some(path) => serde_json::from_str(&fs::read_to_string(path)?)?,
        None => DescriptorOptions::default(),
    };
    if let Some(path) = args.energies {
        let energies: Vec<Option<EnergyMetadata>> =
            serde_json::from_str(&fs::read_to_string(path)?)?;
        if energies.len() != ensemble.conformers.len() {
            return Err("energy metadata must match every ordered source record".into());
        }
        for (index, energy) in energies.into_iter().enumerate() {
            if let Some(energy) = energy {
                ensemble.set_energy(index, energy)?;
            }
        }
    }
    let population = args
        .temperature
        .map(|temperature| {
            ensemble.boltzmann_population(
                temperature,
                args.population_assumptions.as_deref().unwrap_or_default(),
            )
        })
        .transpose()?;
    let descriptors = ensemble.descriptors(&options);
    let statistics = descriptors.statistics(
        &StatisticsOptions {
            missing: if args.available_cases {
                MissingConformerPolicy::AvailableCases
            } else {
                MissingConformerPolicy::RequireComplete
            },
        },
        population.as_ref(),
    )?;
    let report = statistics.sensitivity(
        args.ligand_id
            .as_deref()
            .unwrap_or(&args.input.to_string_lossy()),
    );
    let json = serde_json::to_string_pretty(
        &serde_json::json!({"descriptors":descriptors,"sensitivity":report}),
    )?;
    let csv = report.to_csv()?;
    let paths = [
        args.output.with_extension("json"),
        args.output.with_extension("csv"),
    ];
    if paths.iter().any(|path| path.exists()) {
        return Err("output already exists; choose a new stem to preserve prior evidence".into());
    }
    for (path, data) in paths.iter().zip([json, csv]) {
        fs::File::create_new(path)?.write_all(data.as_bytes())?;
    }
    Ok(())
}
