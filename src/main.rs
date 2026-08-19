mod cli;
mod commands;
mod descriptors;
mod digest;
mod output;
mod reaction;
#[cfg(test)]
mod test_support;

use clap::Parser;
use cli::{Cli, Command};
use std::error::Error;
use steric_x::model::ReactionProvenance;
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
            portable_model,
            model_id,
            reaction_family,
            catalyst_metal,
            ligand_class,
            source_url,
            response_temp_k,
            omit_bootstrap_ensemble,
            optimize,
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
            commands::fit::PortableModelRequest {
                path: portable_model,
                model_id,
                reaction: ReactionProvenance {
                    reaction_family,
                    catalyst_metal,
                    ligand_class,
                    source_url,
                    notes: Vec::new(),
                },
                response_temp_k,
                omit_bootstrap_ensemble,
                optimization: optimize.into(),
            },
        ),
        Command::Evaluate {
            data,
            metadata,
            model,
            predictions,
            output,
        } => commands::evaluate::evaluate_command(&data, &metadata, &model, &predictions, &output),
        Command::Model { action } => match action {
            cli::ModelCommand::Inspect { model, format } => {
                commands::model::inspect_command(&model, format)
            }
            cli::ModelCommand::Validate {
                model,
                format,
                strict,
            } => commands::model::validate_command(&model, format, strict),
        },
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
            sort_by,
            descending,
            features,
            filters,
            vbur,
            l,
            b1,
            b5,
            vbur_min,
            vbur_max,
            l_min,
            l_max,
            b1_min,
            b1_max,
            b5_min,
            b5_max,
            less_bulky,
            more_bulky,
            donor_element,
            sterimol_axis,
            format,
            sphere_radius,
            density,
            center_distance,
            radii_scale,
        } => {
            let shorthand = commands::search::RangeFlags {
                vbur,
                l,
                b1,
                b5,
                vbur_min,
                vbur_max,
                l_min,
                l_max,
                b1_min,
                b1_max,
                b5_min,
                b5_max,
            }
            .to_filters()?;
            commands::search::search_command(commands::search::SearchArgs {
                ligand: ligand.as_deref(),
                library: &library,
                top,
                sort_by: sort_by.as_deref(),
                descending,
                feature_names: features.as_deref(),
                filter_expressions: &filters,
                shorthand_filters: shorthand,
                less_bulky,
                more_bulky,
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
            })
        }
        Command::Screen {
            model,
            library,
            library_flag,
            top,
            temperature,
            in_domain_only,
            ascending,
            descending,
            domain_rule,
            exclude_tested,
            diverse,
            diversity_weight,
            donor_element,
            sterimol_axis,
            format,
            sphere_radius,
            density,
            center_distance,
            radii_scale,
        } => {
            let library = library
                .or(library_flag)
                .ok_or("a library is required: pass it positionally or as --library <PATH>")?;
            commands::screen::screen_command(commands::screen::ScreenArgs {
                model: &model,
                library: &library,
                top,
                temperature,
                in_domain_only,
                ascending,
                descending,
                domain_rule: domain_rule.into(),
                exclude_tested: exclude_tested.as_deref(),
                diverse,
                diversity_weight,
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
            })
        }
        Command::Compare {
            inputs,
            database,
            donor_element,
            sterimol_axis,
            format,
            sphere_radius,
            density,
            center_distance,
            radii_scale,
        } => commands::compare::compare_command(commands::compare::CompareArgs {
            inputs: &inputs,
            database: database.as_deref(),
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
        Command::Db { action } => match action {
            cli::DbCommand::Build {
                source,
                output,
                group_by_parent,
                extensions,
                label_from,
                donor_element,
                sterimol_axis,
                sphere_radius,
                density,
                center_distance,
                radii_scale,
            } => commands::db::db_build_command(commands::db::DbBuildArgs {
                source: &source,
                output: &output,
                group_by_parent,
                extensions,
                label_from,
                donor_element: &donor_element,
                sterimol_axis,
                config: BuriedVolumeConfig {
                    sphere_radius,
                    density,
                    center_distance,
                    radii_scale,
                    include_hydrogens: false,
                },
            }),
        },
    }
}
