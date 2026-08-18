//! `stericx compare`: put two or more ligands side by side.
//!
//! Reports every descriptor for each ligand, the pairwise difference, and — when
//! a database is supplied — that difference expressed in standard deviations of
//! the library. The σ view is what makes a raw number interpretable: a 0.5 Å
//! difference in `B5` is small if the library spans 4 Å and large if it spans
//! 0.6 Å, and only the library can say which.

use crate::cli::{DescriptorFormat, SterimolAxis};
use crate::commands::search::{
    FEATURES, LibraryEntry, default_feature_names, library_statistics, load_library_entries,
};
use crate::descriptors::descriptors_for_file;
use serde::Serialize;
use std::error::Error;
use std::path::{Path, PathBuf};
use steric_x::BuriedVolumeConfig;

#[derive(Debug, Serialize)]
struct DescriptorRow {
    descriptor: String,
    values: Vec<f32>,
    /// `max − min` across the compared ligands.
    spread: f32,
    /// `spread` in library standard deviations, when a database is supplied.
    spread_sigma: Option<f64>,
}

#[derive(Debug, Serialize)]
struct PairDistance {
    left: String,
    right: String,
    /// Euclidean distance in library-standardized descriptor space.
    standardized_distance: f64,
}

#[derive(Debug, Serialize)]
struct CompareReport {
    ligands: Vec<LibraryEntry>,
    descriptors: Vec<DescriptorRow>,
    /// Present only when a database provided the standardization.
    database: Option<String>,
    database_size: Option<usize>,
    similarity_features: Vec<String>,
    pairwise_distances: Vec<PairDistance>,
    note: String,
}

pub(crate) struct CompareArgs<'a> {
    pub(crate) inputs: &'a [PathBuf],
    pub(crate) database: Option<&'a Path>,
    pub(crate) donor_element: &'a str,
    pub(crate) sterimol_axis: SterimolAxis,
    pub(crate) format: DescriptorFormat,
    pub(crate) config: BuriedVolumeConfig,
}

pub(crate) fn compare_command(args: CompareArgs<'_>) -> Result<(), Box<dyn Error>> {
    if args.inputs.len() < 2 {
        return Err("compare needs at least two ligand files".into());
    }
    let mut ligands = Vec::with_capacity(args.inputs.len());
    for path in args.inputs {
        let result = descriptors_for_file(
            path,
            args.donor_element,
            None,
            args.sterimol_axis,
            args.config,
        )
        .map_err(|message| format!("{}: {message}", path.display()))?;
        ligands.push(LibraryEntry::from(&result));
    }

    // A database is optional: without it the raw differences still stand, they
    // just cannot be placed on a scale.
    let statistics = match args.database {
        Some(path) => {
            let entries =
                load_library_entries(path, args.donor_element, args.sterimol_axis, args.config)?;
            if entries.len() < 2 {
                return Err("a database needs at least two members to standardize".into());
            }
            Some((
                path.display().to_string(),
                entries.len(),
                library_statistics(&entries),
            ))
        }
        None => None,
    };

    let descriptors = FEATURES
        .iter()
        .map(|feature| {
            let (name, read) = (feature.name, feature.get);
            let values = ligands.iter().map(read).collect::<Vec<_>>();
            let maximum = values.iter().copied().fold(f32::NEG_INFINITY, f32::max);
            let minimum = values.iter().copied().fold(f32::INFINITY, f32::min);
            let spread = maximum - minimum;
            let spread_sigma = statistics.as_ref().and_then(|(_, _, stats)| {
                stats
                    .get(name)
                    .filter(|(_, deviation)| *deviation > f64::EPSILON)
                    .map(|(_, deviation)| f64::from(spread) / deviation)
            });
            DescriptorRow {
                descriptor: name.to_owned(),
                values,
                spread,
                spread_sigma,
            }
        })
        .collect::<Vec<_>>();

    // Pairwise distances use the same default similarity space as `search`, so a
    // `compare` distance and a `search` distance mean the same thing.
    let similarity_features = default_feature_names();
    let mut pairwise_distances = Vec::new();
    if let Some((_, _, stats)) = statistics.as_ref() {
        for left in 0..ligands.len() {
            for right in (left + 1)..ligands.len() {
                let mut sum = 0.0_f64;
                for name in &similarity_features {
                    let Some((mean, deviation)) = stats.get(name.as_str()) else {
                        continue;
                    };
                    if *deviation <= f64::EPSILON {
                        continue;
                    }
                    let Some(read) = FEATURES
                        .iter()
                        .find(|feature| feature.name == name.as_str())
                        .map(|feature| feature.get)
                    else {
                        continue;
                    };
                    let a = (f64::from(read(&ligands[left])) - mean) / deviation;
                    let b = (f64::from(read(&ligands[right])) - mean) / deviation;
                    sum += (a - b).powi(2);
                }
                pairwise_distances.push(PairDistance {
                    left: ligands[left].file.clone(),
                    right: ligands[right].file.clone(),
                    standardized_distance: sum.sqrt(),
                });
            }
        }
    }

    let report = CompareReport {
        ligands: ligands.clone(),
        descriptors,
        database: statistics.as_ref().map(|(path, _, _)| path.clone()),
        database_size: statistics.as_ref().map(|(_, size, _)| *size),
        similarity_features,
        pairwise_distances,
        note: if statistics.is_some() {
            "sigma columns and distances are scaled by the supplied database's spread; a \
             difference is only large or small relative to the library it sits in."
                .to_owned()
        } else {
            "no --database supplied, so differences are reported in raw units only. Pass a \
             database to express them in library standard deviations and to compute a \
             standardized distance."
                .to_owned()
        },
    };

    match args.format {
        DescriptorFormat::Text => print_compare_text(&report),
        DescriptorFormat::Json => println!("{}", serde_json::to_string_pretty(&report)?),
        DescriptorFormat::Csv => print_compare_csv(&report),
    }
    Ok(())
}

fn short_label(path: &str) -> String {
    Path::new(path)
        .file_stem()
        .and_then(|stem| stem.to_str())
        .unwrap_or(path)
        .to_owned()
}

fn print_compare_text(report: &CompareReport) {
    let labels = report
        .ligands
        .iter()
        .map(|entry| short_label(&entry.file))
        .collect::<Vec<_>>();
    let width = labels.iter().map(String::len).max().unwrap_or(8).max(10);

    for (entry, label) in report.ligands.iter().zip(&labels) {
        println!(
            "{label:<width$}  {}  donor {} (atom {}), substituents {}, {} conformer(s)",
            entry.file,
            entry.donor_element,
            entry.donor_index,
            entry.substituents,
            entry.conformers
        );
    }
    match (&report.database, report.database_size) {
        (Some(path), Some(size)) => println!("\nscale: {path} ({size} members)"),
        _ => println!("\nscale: none (raw units only — pass --database for σ)"),
    }

    print!("\n{:<24}", "descriptor");
    for label in &labels {
        print!("{label:>width$}");
    }
    print!("{:>10}", "spread");
    if report.database.is_some() {
        print!("{:>9}", "σ");
    }
    println!();

    for row in &report.descriptors {
        print!("{:<24}", row.descriptor);
        for value in &row.values {
            print!("{value:>width$.3}");
        }
        print!("{:>10.3}", row.spread);
        if report.database.is_some() {
            match row.spread_sigma {
                Some(sigma) => print!("{sigma:>9.2}"),
                None => print!("{:>9}", "—"),
            }
        }
        println!();
    }

    if !report.pairwise_distances.is_empty() {
        println!(
            "\nstandardized distance over {}:",
            report.similarity_features.join(", ")
        );
        for pair in &report.pairwise_distances {
            println!(
                "  {:<width$} vs {:<width$}  {:.3}",
                short_label(&pair.left),
                short_label(&pair.right),
                pair.standardized_distance
            );
        }
    }
    println!("\nnote: {}", report.note);
}

fn print_compare_csv(report: &CompareReport) {
    let labels = report
        .ligands
        .iter()
        .map(|entry| crate::descriptors::csv_field(&short_label(&entry.file)))
        .collect::<Vec<_>>();
    println!("descriptor,{},spread,spread_sigma", labels.join(","));
    for row in &report.descriptors {
        let values = row
            .values
            .iter()
            .map(|value| format!("{value:.4}"))
            .collect::<Vec<_>>()
            .join(",");
        let sigma = row
            .spread_sigma
            .map_or_else(String::new, |sigma| format!("{sigma:.4}"));
        println!("{},{},{:.4},{}", row.descriptor, values, row.spread, sigma);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shortens_paths_to_readable_labels() {
        assert_eq!(short_label("ligands/xphos_1.sdf"), "xphos_1");
        assert_eq!(short_label("plain"), "plain");
    }
}
