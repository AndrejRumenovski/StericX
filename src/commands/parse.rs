//! `stericx parse`: reaction metadata + XYZ geometries → flat `.sigpack` matrix.

use crate::output::{millis, print_memory_metrics, resident_memory_bytes};
use crate::reaction::{
    ReactionCsvRow, conformer_energy_span, conformer_paths, conformer_weights,
    record_from_ensemble, resolve_xyz_path, sterimol_from_molecule, validate_reaction_row,
};
use std::error::Error;
use std::fs;
use std::mem::size_of;
use std::path::Path;
use std::time::{Duration, Instant};
use steric_x::{Molecule, PackedReactionRecord, SigPackWriter};

pub(crate) fn parse_command(
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::temporary_directory;
    use steric_x::SigPackReader;

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
}
