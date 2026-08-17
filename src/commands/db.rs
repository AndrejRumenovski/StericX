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

fn sha256_hex(bytes: &[u8]) -> String {
    // Compact SHA-256 so a database carries the same kind of integrity anchor
    // the project's frozen artifacts already use, without a new dependency.
    const K: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4,
        0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe,
        0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f,
        0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
        0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc,
        0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
        0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116,
        0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7,
        0xc67178f2,
    ];
    let mut hash: [u32; 8] = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab,
        0x5be0cd19,
    ];
    let mut message = bytes.to_vec();
    let bit_length = (bytes.len() as u64) * 8;
    message.push(0x80);
    while message.len() % 64 != 56 {
        message.push(0);
    }
    message.extend_from_slice(&bit_length.to_be_bytes());

    for chunk in message.chunks(64) {
        let mut w = [0_u32; 64];
        for (index, word) in chunk.chunks(4).enumerate() {
            w[index] = u32::from_be_bytes([word[0], word[1], word[2], word[3]]);
        }
        for index in 16..64 {
            let s0 = w[index - 15].rotate_right(7)
                ^ w[index - 15].rotate_right(18)
                ^ (w[index - 15] >> 3);
            let s1 = w[index - 2].rotate_right(17)
                ^ w[index - 2].rotate_right(19)
                ^ (w[index - 2] >> 10);
            w[index] = w[index - 16]
                .wrapping_add(s0)
                .wrapping_add(w[index - 7])
                .wrapping_add(s1);
        }
        let mut v = hash;
        for index in 0..64 {
            let s1 = v[4].rotate_right(6) ^ v[4].rotate_right(11) ^ v[4].rotate_right(25);
            let choice = (v[4] & v[5]) ^ (!v[4] & v[6]);
            let temp1 = v[7]
                .wrapping_add(s1)
                .wrapping_add(choice)
                .wrapping_add(K[index])
                .wrapping_add(w[index]);
            let s0 = v[0].rotate_right(2) ^ v[0].rotate_right(13) ^ v[0].rotate_right(22);
            let majority = (v[0] & v[1]) ^ (v[0] & v[2]) ^ (v[1] & v[2]);
            let temp2 = s0.wrapping_add(majority);
            v[7] = v[6];
            v[6] = v[5];
            v[5] = v[4];
            v[4] = v[3].wrapping_add(temp1);
            v[3] = v[2];
            v[2] = v[1];
            v[1] = v[0];
            v[0] = temp1.wrapping_add(temp2);
        }
        for (slot, value) in hash.iter_mut().zip(v) {
            *slot = slot.wrapping_add(value);
        }
    }
    hash.iter().map(|word| format!("{word:08x}")).collect()
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
        table_sha256: sha256_hex(table.as_bytes()),
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
            sha256_hex(b""),
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        );
        assert_eq!(
            sha256_hex(b"abc"),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
        // Spans a multi-block message so the padding path is exercised.
        assert_eq!(
            sha256_hex(b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"),
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
