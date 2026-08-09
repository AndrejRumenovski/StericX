//! `stericx buried-volume`: coordination-aware buried volumes → `.sigpack` v2.

use crate::output::{
    atomic_write_csv_rows, ensure_parent, millis, print_memory_metrics, resident_memory_bytes,
};
use crate::reaction::{
    ReactionCsvRow, conformer_coordination_centers, conformer_energy_span, conformer_paths,
    conformer_weights, record_from_ensemble, resolve_xyz_path, sterimol_from_molecule,
    validate_reaction_row,
};
use serde::Serialize;
use std::error::Error;
use std::fs;
use std::path::Path;
use std::time::Instant;
use steric_x::{
    BuriedVolumeCalculator, BuriedVolumeConfig, Molecule, PackedBuriedVolumeRecord,
    PackedReactionRecordV2, SigPackV2Writer,
};

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

pub(crate) fn buried_volume_command(
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::temporary_directory;

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
}
