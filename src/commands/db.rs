//! `stericx db build`: precompute a reusable ligand descriptor database.
//!
//! Featurizing is the expensive step — grid integration over every conformer —
//! while reading a table back is trivially cheap. A database therefore stores
//! computed descriptors as a plain CSV next to a JSON manifest recording exactly
//! which settings produced it, how many ligands it covers, and the SHA-256 of
//! the table. Text was chosen deliberately over a packed binary: it costs
//! nothing measurable to parse at this scale, stays diff-able and inspectable,
//! and avoids the endianness caveat the `.sigpack` format carries.

use crate::cli::{DbLabel, SterimolAxis};
use crate::commands::search::{FEATURES, LibraryEntry};
use crate::descriptors::{csv_field, descriptors_for_file};
use rayon::prelude::*;
use serde::Serialize;
use std::collections::BTreeMap;
use std::error::Error;
use std::path::{Path, PathBuf};
use std::time::Instant;
use steric_x::BuriedVolumeConfig;

#[derive(Debug, Serialize)]
struct DatabaseManifest {
    schema_version: u32,
    generator: String,
    source: String,
    grouping: String,
    extensions: Vec<String>,
    label_from: String,
    donor_element: String,
    sterimol_axis: String,
    sphere_radius: f32,
    grid_density: f32,
    center_distance: f32,
    radii_scale: f32,
    ligands: usize,
    geometries_featurized: usize,
    geometries_skipped: usize,
    descriptor_columns: Vec<String>,
    table: String,
    table_sha256: String,
    build_seconds: f64,
}

pub(crate) struct DbBuildArgs<'a> {
    pub(crate) source: &'a Path,
    pub(crate) output: &'a Path,
    pub(crate) group_by_parent: bool,
    /// Geometry extensions to include, lower-case and without the dot.
    pub(crate) extensions: Vec<String>,
    pub(crate) label_from: DbLabel,
    pub(crate) donor_element: &'a str,
    pub(crate) sterimol_axis: SterimolAxis,
    pub(crate) config: BuriedVolumeConfig,
}

/// One ligand's row: a stable label plus the descriptor columns.
#[derive(Clone, Debug)]
struct DatabaseRow {
    label: String,
    entry: LibraryEntry,
    /// Geometries that contributed, > 1 when conformers were aggregated.
    geometries: usize,
}

fn collect_coordinate_files(
    directory: &Path,
    allowed: &[String],
    found: &mut Vec<PathBuf>,
) -> std::io::Result<()> {
    let mut entries = std::fs::read_dir(directory)?
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .map(|entry| entry.path())
        .collect::<Vec<_>>();
    entries.sort();
    for path in entries {
        if path.is_dir() {
            collect_coordinate_files(&path, allowed, found)?;
        } else if path
            .extension()
            .and_then(|extension| extension.to_str())
            .is_some_and(|extension| {
                allowed
                    .iter()
                    .any(|allow| allow.eq_ignore_ascii_case(extension))
            })
        {
            found.push(path);
        }
    }
    Ok(())
}

fn label_for(path: &Path, source: &Path, label_from: DbLabel) -> String {
    match label_from {
        DbLabel::Stem => path
            .file_stem()
            .and_then(|stem| stem.to_str())
            .unwrap_or("unknown")
            .to_owned(),
        DbLabel::Parent => path
            .parent()
            .and_then(Path::file_name)
            .and_then(|name| name.to_str())
            .unwrap_or("unknown")
            .to_owned(),
        DbLabel::Path => path
            .strip_prefix(source)
            .unwrap_or(path)
            .display()
            .to_string(),
    }
}

/// Aggregate several conformers of one ligand into a single row.
///
/// Most descriptors average, but `max_delta_qvbur_min` takes the minimum across
/// conformers — that is Kraken's own `*_min` convention for the headline
/// quadrant-asymmetry descriptor, so aggregating any other way would silently
/// redefine it.
fn aggregate(label: String, group: &[LibraryEntry]) -> DatabaseRow {
    let count = group.len() as f32;
    let mean = |read: fn(&LibraryEntry) -> f32| group.iter().map(read).sum::<f32>() / count;
    let representative = &group[0];
    DatabaseRow {
        entry: LibraryEntry {
            ligand: label.clone(),
            file: representative.file.clone(),
            conformers: group.iter().map(|entry| entry.conformers).sum(),
            donor_element: representative.donor_element.clone(),
            donor_index: representative.donor_index,
            substituents: representative.substituents.clone(),
            sterimol_l: mean(|entry| entry.sterimol_l),
            sterimol_b1: mean(|entry| entry.sterimol_b1),
            sterimol_b5: mean(|entry| entry.sterimol_b5),
            percent_buried_volume: mean(|entry| entry.percent_buried_volume),
            buried_volume: mean(|entry| entry.buried_volume),
            qvbur_min: mean(|entry| entry.qvbur_min),
            qvbur_max: mean(|entry| entry.qvbur_max),
            max_delta_qvbur: mean(|entry| entry.max_delta_qvbur),
            max_delta_qvbur_min: group
                .iter()
                .map(|entry| entry.max_delta_qvbur_min)
                .fold(f32::INFINITY, f32::min),
            pyr_p: mean(|entry| entry.pyr_p),
            pyr_alpha: mean(|entry| entry.pyr_alpha),
        },
        geometries: group.len(),
        label,
    }
}

pub(crate) fn db_build_command(args: DbBuildArgs<'_>) -> Result<(), Box<dyn Error>> {
    let started = Instant::now();
    if !args.source.is_dir() {
        return Err(format!("source directory does not exist: {}", args.source.display()).into());
    }
    let extensions = if args.extensions.is_empty() {
        vec!["xyz".to_owned(), "sdf".to_owned(), "mol".to_owned()]
    } else {
        args.extensions
            .iter()
            .map(|e| e.trim_start_matches('.').to_ascii_lowercase())
            .collect()
    };
    let mut paths = Vec::new();
    collect_coordinate_files(args.source, &extensions, &mut paths)?;
    if paths.is_empty() {
        return Err(format!(
            "no .xyz/.sdf/.mol geometries found below {}",
            args.source.display()
        )
        .into());
    }

    let featurized = paths
        .par_iter()
        .map(|path| {
            let label = label_for(path, args.source, args.label_from);
            match descriptors_for_file(
                path,
                args.donor_element,
                None,
                args.sterimol_axis,
                args.config,
            ) {
                Ok(result) => Some((label, LibraryEntry::from(&result))),
                Err(message) => {
                    eprintln!("skipped {}: {message}", path.display());
                    None
                }
            }
        })
        .collect::<Vec<_>>();
    let skipped = featurized.iter().filter(|slot| slot.is_none()).count();

    // BTreeMap keeps the emitted order deterministic regardless of how rayon
    // scheduled the work, so rebuilding an unchanged source is reproducible.
    let mut grouped: BTreeMap<String, Vec<LibraryEntry>> = BTreeMap::new();
    for (label, entry) in featurized.into_iter().flatten() {
        grouped.entry(label).or_default().push(entry);
    }
    if grouped.is_empty() {
        return Err("no geometry could be featurized into the database".into());
    }

    let rows = if args.group_by_parent {
        grouped
            .into_iter()
            .map(|(label, group)| aggregate(label, &group))
            .collect::<Vec<_>>()
    } else {
        grouped
            .into_iter()
            .flat_map(|(label, group)| {
                // Without grouping each geometry stands alone; a label collision
                // gets an index suffix so every row stays addressable.
                let multiple = group.len() > 1;
                group.into_iter().enumerate().map(move |(index, entry)| {
                    let label = if multiple {
                        format!("{label}#{index}")
                    } else {
                        label.clone()
                    };
                    DatabaseRow {
                        label,
                        entry,
                        geometries: 1,
                    }
                })
            })
            .collect::<Vec<_>>()
    };

    let descriptor_columns = FEATURES
        .iter()
        .map(|feature| feature.name.to_owned())
        .collect::<Vec<_>>();
    let mut table = String::new();
    table.push_str(
        "ligand,file,geometries,conformers,donor_element,donor_index,substituents,\
         sterimol_l,sterimol_b1,sterimol_b5,percent_buried_volume,buried_volume,\
         qvbur_min,qvbur_max,max_delta_qvbur,max_delta_qvbur_min,pyr_p,pyr_alpha\n",
    );
    for row in &rows {
        let entry = &row.entry;
        table.push_str(&format!(
            "{},{},{},{},{},{},{},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4}\n",
            csv_field(&row.label),
            csv_field(&entry.file),
            row.geometries,
            entry.conformers,
            entry.donor_element,
            entry.donor_index,
            csv_field(&entry.substituents),
            entry.sterimol_l,
            entry.sterimol_b1,
            entry.sterimol_b5,
            entry.percent_buried_volume,
            entry.buried_volume,
            entry.qvbur_min,
            entry.qvbur_max,
            entry.max_delta_qvbur,
            entry.max_delta_qvbur_min,
            entry.pyr_p,
            entry.pyr_alpha
        ));
    }

    crate::output::write_atomic_text(args.output, &table)?;
    let manifest_path = args.output.with_extension("manifest.json");
    let manifest = DatabaseManifest {
        schema_version: 1,
        generator: format!("stericx {}", env!("CARGO_PKG_VERSION")),
        source: args.source.display().to_string(),
        grouping: if args.group_by_parent {
            "parent_directory_is_one_ligand".to_owned()
        } else {
            "one_row_per_geometry".to_owned()
        },
        extensions: extensions.clone(),
        label_from: format!("{:?}", args.label_from).to_lowercase(),
        donor_element: args.donor_element.to_owned(),
        sterimol_axis: format!("{:?}", args.sterimol_axis).to_lowercase(),
        sphere_radius: args.config.sphere_radius,
        grid_density: args.config.density,
        center_distance: args.config.center_distance,
        radii_scale: args.config.radii_scale,
        ligands: rows.len(),
        geometries_featurized: paths.len() - skipped,
        geometries_skipped: skipped,
        descriptor_columns,
        table: args
            .output
            .file_name()
            .map(|name| name.to_string_lossy().into_owned())
            .unwrap_or_default(),
        table_sha256: crate::digest::sha256_hex(table.as_bytes()),
        build_seconds: started.elapsed().as_secs_f64(),
    };
    crate::output::atomic_write_json(&manifest, &manifest_path)?;

    println!("command=db build");
    println!("ligands={}", manifest.ligands);
    println!("geometries_featurized={}", manifest.geometries_featurized);
    println!("geometries_skipped={}", manifest.geometries_skipped);
    println!("database={}", args.output.display());
    println!("manifest={}", manifest_path.display());
    println!("table_sha256={}", manifest.table_sha256);
    println!("build_seconds={:.3}", manifest.build_seconds);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_published_sha256_vectors() {
        assert_eq!(
            crate::digest::sha256_hex(b""),
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        );
        assert_eq!(
            crate::digest::sha256_hex(b"abc"),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
        // Spans a multi-block message so the padding path is exercised.
        assert_eq!(
            crate::digest::sha256_hex(b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"),
            "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"
        );
    }

    fn entry(l: f32, b5: f32, vbur: f32, dq: f32) -> LibraryEntry {
        LibraryEntry {
            ligand: String::new(),
            file: "x.sdf".into(),
            conformers: 1,
            donor_element: "P".into(),
            donor_index: 0,
            substituents: "C C C".into(),
            sterimol_l: l,
            sterimol_b1: 1.7,
            sterimol_b5: b5,
            percent_buried_volume: vbur,
            buried_volume: vbur * 1.7959,
            qvbur_min: 9.0,
            qvbur_max: 13.0,
            max_delta_qvbur: dq,
            max_delta_qvbur_min: dq,
            pyr_p: 0.93,
            pyr_alpha: 17.0,
        }
    }

    #[test]
    fn aggregation_averages_but_takes_the_minimum_of_the_kraken_min_descriptor() {
        let group = vec![entry(8.0, 7.0, 30.0, 5.0), entry(10.0, 9.0, 40.0, 3.0)];
        let row = aggregate("lig".into(), &group);
        assert!((row.entry.sterimol_l - 9.0).abs() < 1e-6);
        assert!((row.entry.percent_buried_volume - 35.0).abs() < 1e-6);
        // Kraken's `*_min` convention: the minimum over conformers, not the mean.
        assert!((row.entry.max_delta_qvbur_min - 3.0).abs() < 1e-6);
        assert_eq!(row.geometries, 2);
        assert_eq!(row.entry.conformers, 2);
    }

    #[test]
    fn labels_derive_from_the_requested_component() {
        let path = Path::new("/lib/1088/49973.sdf");
        let source = Path::new("/lib");
        assert_eq!(label_for(path, source, DbLabel::Stem), "49973");
        assert_eq!(label_for(path, source, DbLabel::Parent), "1088");
        assert_eq!(label_for(path, source, DbLabel::Path), "1088/49973.sdf");
    }
}
