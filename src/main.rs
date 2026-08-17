mod cli;
mod commands;
mod descriptors;
mod output;
mod reaction;
#[cfg(test)]
mod test_support;

use clap::Parser;
use cli::{Cli, Command};
use std::error::Error;
use steric_x::{BuriedVolumeConfig, FitOptions};

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
        } => commands::parse::parse_command(&csv, &xyz_dir, &output),
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
        } => commands::buried_volume::buried_volume_command(
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
        Command::Predict { data, weights } => commands::predict::predict_command(&data, &weights),
        Command::Fit {
            data,
            metadata,
            output,
            predictions,
            max_terms,
            bootstrap,
            permutations,
            seed,
        } => commands::fit::fit_command(
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
        } => commands::evaluate::evaluate_command(&data, &metadata, &model, &predictions, &output),
        Command::Simulate { ddg, temp } => commands::simulate::simulate_command(ddg, temp),
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
        } => descriptors::descriptors_command(
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
        Command::Search {
            ligand,
            library,
            top,
            features,
            filters,
            less_bulky,
            more_bulky,
            donor_element,
            sterimol_axis,
            format,
            sphere_radius,
            density,
            center_distance,
            radii_scale,
        } => commands::search::search_command(
            &ligand,
            &library,
            top,
            features.as_deref(),
            &filters,
            less_bulky,
            more_bulky,
            &donor_element,
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
        Command::Screen {
            model,
            library,
            top,
            temperature,
            inside_domain_only,
            ascending,
            donor_element,
            sterimol_axis,
            format,
            sphere_radius,
            density,
            center_distance,
            radii_scale,
        } => commands::screen::screen_command(commands::screen::ScreenArgs {
            model: &model,
            library: &library,
            top,
            temperature,
            inside_domain_only,
            ascending,
            donor_element: &donor_element,
            sterimol_axis,
            format,
            config: BuriedVolumeConfig {
                sphere_radius,
                density,
                center_distance,
                radii_scale,
                include_hydrogens: false,
            },
        }),
    }
}
