//! Reaction-CSV ingestion: the row schemas plus the helpers that resolve
//! conformer geometries, aggregate Sterimol ensembles, and validate input.

use glam::Vec3;
use serde::Deserialize;
use std::error::Error;
use std::path::{Component, Path, PathBuf};
use steric_x::ReactionLabel;
use steric_x::model::is_supported_split;
use steric_x::{Molecule, PackedReactionRecord, SterimolCalculator, SterimolParams};

/// One row of the raw reaction CSV.
#[derive(Clone, Debug, Deserialize)]
pub(crate) struct ReactionCsvRow {
    #[serde(rename = "Reaction_ID", alias = "reaction_id")]
    pub(crate) reaction_id: String,
    #[serde(rename = "Ligand_XYZ_Path", alias = "ligand_xyz_path", alias = "file")]
    pub(crate) ligand_xyz_path: PathBuf,
    #[serde(rename = "Attach_Atom_Idx", alias = "attach_idx")]
    pub(crate) attach_atom_idx: usize,
    #[serde(
        rename = "Primary_Bond_Vector_Idx",
        alias = "primary_bond_vector_idx",
        alias = "axis_atom_idx",
        alias = "neighbor_idx"
    )]
    pub(crate) primary_bond_vector_idx: usize,
    #[serde(rename = "NBO_Charge", alias = "nbo_charge")]
    pub(crate) nbo_charge: f32,
    #[serde(rename = "IR_Frequency", alias = "ir_freq")]
    pub(crate) ir_freq: f32,
    #[serde(rename = "Temp_K", alias = "temp_k")]
    pub(crate) temp_k: f32,
    #[serde(rename = "Exp_ddG_kcal_mol", alias = "exp_ddg")]
    pub(crate) exp_ddg: f32,
    #[serde(rename = "Conformer_XYZ_Paths", alias = "conformer_xyz_paths", default)]
    pub(crate) conformer_xyz_paths: Option<String>,
    #[serde(
        rename = "Conformer_Relative_Energies_kcal_mol",
        alias = "conformer_relative_energies_kcal_mol",
        default
    )]
    pub(crate) conformer_relative_energies: Option<String>,
    #[serde(
        rename = "Conformer_Boltzmann_Weights",
        alias = "conformer_boltzmann_weights",
        default
    )]
    pub(crate) conformer_boltzmann_weights: Option<String>,
    #[serde(
        rename = "Conformer_Coordination_Centers_Angstrom",
        alias = "conformer_coordination_centers_angstrom",
        default
    )]
    pub(crate) conformer_coordination_centers: Option<String>,
    #[serde(
        rename = "Coordination_Center_Method",
        alias = "coordination_center_method",
        default
    )]
    pub(crate) coordination_center_method: Option<String>,
}

/// Row-aligned split metadata. The schema lives in the model layer so the
/// library can derive a training partition without the CLI.
pub(crate) type ReactionMetadataRow = ReactionLabel;

pub(crate) fn sterimol_from_molecule(
    molecule: &Molecule,
    row: &ReactionCsvRow,
) -> Result<SterimolParams, String> {
    steric_x::profile_scope!("sterimol", "reaction::sterimol_from_molecule");
    if molecule.atoms.len() < 2 {
        return Err("at least two atoms are required to establish the primary axis".into());
    }
    if row.attach_atom_idx == row.primary_bond_vector_idx {
        return Err("attachment and primary-vector indices must differ".into());
    }
    let attachment = molecule.atoms.get(row.attach_atom_idx).ok_or_else(|| {
        format!(
            "attachment atom index {} is out of bounds for {} atoms",
            row.attach_atom_idx,
            molecule.atoms.len()
        )
    })?;
    let neighbor = molecule
        .atoms
        .get(row.primary_bond_vector_idx)
        .ok_or_else(|| {
            format!(
                "primary-vector atom index {} is out of bounds for {} atoms",
                row.primary_bond_vector_idx,
                molecule.atoms.len()
            )
        })?;
    if neighbor.position == attachment.position {
        return Err("attachment and primary-vector atoms occupy the same position".into());
    }

    SterimolCalculator::compute(molecule, row.attach_atom_idx, row.primary_bond_vector_idx)
        .map_err(|error| error.to_string())
}

pub(crate) fn record_from_ensemble(
    conformers: &[SterimolParams],
    weights: &[f32],
    energy_span: f32,
    row: &ReactionCsvRow,
) -> Result<PackedReactionRecord, String> {
    steric_x::profile_scope!("conformer_processing", "reaction::record_from_ensemble");
    if !energy_span.is_finite() || energy_span < 0.0 {
        return Err("ensemble energy span must be finite and non-negative".into());
    }
    let weighted = weighted_sterimol(conformers, weights)?;
    let minimum = conformers.iter().fold(
        SterimolParams {
            l: f32::INFINITY,
            b1: f32::INFINITY,
            b5: f32::INFINITY,
        },
        |minimum, params| SterimolParams {
            l: minimum.l.min(params.l),
            b1: minimum.b1.min(params.b1),
            b5: minimum.b5.min(params.b5),
        },
    );
    let maximum = conformers
        .iter()
        .fold(conformers[0], |maximum, params| SterimolParams {
            l: maximum.l.max(params.l),
            b1: maximum.b1.max(params.b1),
            b5: maximum.b5.max(params.b5),
        });
    if [
        weighted.l,
        weighted.b1,
        weighted.b5,
        minimum.l,
        minimum.b1,
        minimum.b5,
        maximum.l,
        maximum.b1,
        maximum.b5,
        energy_span,
    ]
    .iter()
    .any(|value| !value.is_finite())
    {
        return Err("ensemble aggregation produced a non-finite descriptor".into());
    }

    let mut record = PackedReactionRecord::from_ensemble(
        weighted.l,
        weighted.b1,
        weighted.b5,
        minimum.l,
        maximum.l,
        minimum.b1,
        maximum.b1,
        minimum.b5,
        maximum.b5,
        conformers.len(),
        energy_span,
    );
    record.nbo_charge = row.nbo_charge;
    record.ir_freq = row.ir_freq;
    record.temp_k = row.temp_k;
    record.exp_ddg = row.exp_ddg;
    Ok(record)
}

/// Average supplied finite descriptors with explicit nonnegative weights.
/// This operation does not infer thermodynamic provenance from the weights.
pub(crate) fn weighted_sterimol(
    conformers: &[SterimolParams],
    weights: &[f32],
) -> Result<SterimolParams, String> {
    if conformers.is_empty() || conformers.len() != weights.len() {
        return Err("conformer parameters and weights must have equal non-zero lengths".into());
    }
    if conformers
        .iter()
        .any(|p| !p.l.is_finite() || !p.b1.is_finite() || !p.b5.is_finite())
    {
        return Err("conformer descriptors must be finite".into());
    }
    let normalized = normalized_weights(weights)?;
    let mean = |select: fn(&SterimolParams) -> f32| {
        conformers
            .iter()
            .zip(&normalized)
            .map(|(params, weight)| f64::from(select(params)) * weight)
            .sum::<f64>() as f32
    };
    let result = SterimolParams {
        l: mean(|p| p.l),
        b1: mean(|p| p.b1),
        b5: mean(|p| p.b5),
    };
    if !result.l.is_finite() || !result.b1.is_finite() || !result.b5.is_finite() {
        return Err("ensemble aggregation exceeds the finite f32 descriptor range".into());
    }
    Ok(result)
}

pub(crate) fn conformer_paths(row: &ReactionCsvRow) -> Result<Vec<PathBuf>, String> {
    steric_x::profile_scope!("conformer_processing", "reaction::conformer_paths");
    let paths = row
        .conformer_xyz_paths
        .as_deref()
        .map(|value| {
            value
                .split(';')
                .map(str::trim)
                .filter(|path| !path.is_empty())
                .map(PathBuf::from)
                .collect::<Vec<_>>()
        })
        .filter(|paths| !paths.is_empty())
        .unwrap_or_else(|| vec![row.ligand_xyz_path.clone()]);
    for path in &paths {
        if path.as_os_str().is_empty()
            || path
                .components()
                .any(|component| component == Component::ParentDir)
        {
            return Err(format!(
                "invalid conformer coordinate path `{}`",
                path.display()
            ));
        }
    }
    Ok(paths)
}

fn parse_semicolon_floats(value: &str, label: &str) -> Result<Vec<f32>, String> {
    value
        .split(';')
        .map(str::trim)
        .filter(|field| !field.is_empty())
        .map(|field| {
            field
                .parse::<f32>()
                .map_err(|_| format!("{label} contains invalid float `{field}`"))
        })
        .collect()
}

/// Normalization is invariant to a positive common scale, including subnormal
/// inputs. Convert before dividing so small ratios and large sums use f64.
pub(crate) fn normalized_weights(weights: &[f32]) -> Result<Vec<f64>, String> {
    if weights.is_empty()
        || weights
            .iter()
            .any(|weight| !weight.is_finite() || *weight < 0.0)
    {
        return Err("conformer weights must be finite and non-negative".into());
    }
    let scale = weights.iter().copied().fold(0.0_f32, f32::max);
    if scale == 0.0 {
        return Err("conformer weights must have a positive sum".into());
    }
    let scaled = weights
        .iter()
        .map(|weight| f64::from(*weight) / f64::from(scale))
        .collect::<Vec<_>>();
    let total = scaled.iter().sum::<f64>();
    Ok(scaled.into_iter().map(|weight| weight / total).collect())
}

pub(crate) fn conformer_weights(
    row: &ReactionCsvRow,
    conformer_count: usize,
) -> Result<Vec<f32>, String> {
    steric_x::profile_scope!("conformer_processing", "reaction::conformer_weights");
    if conformer_count == 0 {
        return Err("conformer weights require at least one conformer".into());
    }
    let weights = match row.conformer_boltzmann_weights.as_deref() {
        Some(value) if !value.trim().is_empty() => {
            parse_semicolon_floats(value, "Conformer_Boltzmann_Weights")?
        }
        // Missing weights explicitly mean a uniform average, even when energy
        // metadata is present. No thermodynamic provenance is inferred here.
        _ => vec![1.0; conformer_count],
    };
    if weights.len() != conformer_count {
        return Err(format!(
            "found {} conformer weights for {conformer_count} coordinate paths",
            weights.len()
        ));
    }
    // Preserve raw values until multiplication in f64. Serializing normalized
    // probabilities as f32 first can erase small but meaningful contributions.
    normalized_weights(&weights)?;
    Ok(weights)
}

pub(crate) fn conformer_coordination_centers(
    row: &ReactionCsvRow,
    conformer_count: usize,
) -> Result<Option<Vec<Vec3>>, String> {
    steric_x::profile_scope!(
        "conformer_processing",
        "reaction::conformer_coordination_centers"
    );
    let Some(value) = row.conformer_coordination_centers.as_deref() else {
        return Ok(None);
    };
    if value.trim().is_empty() {
        return Ok(None);
    }
    let centers = value
        .split(';')
        .map(str::trim)
        .filter(|field| !field.is_empty())
        .map(|field| {
            let components = field
                .split(',')
                .map(str::trim)
                .map(|component| {
                    component.parse::<f32>().map_err(|_| {
                        format!(
                            "Conformer_Coordination_Centers_Angstrom contains invalid float `{component}`"
                        )
                    })
                })
                .collect::<Result<Vec<_>, _>>()?;
            if components.len() != 3 {
                return Err(format!(
                    "coordination center `{field}` must contain exactly x,y,z"
                ));
            }
            let center = Vec3::new(components[0], components[1], components[2]);
            if !center.is_finite() {
                return Err("coordination centers must be finite".to_owned());
            }
            Ok(center)
        })
        .collect::<Result<Vec<_>, String>>()?;
    if centers.len() != conformer_count {
        return Err(format!(
            "found {} coordination centers for {conformer_count} coordinate paths",
            centers.len()
        ));
    }
    Ok(Some(centers))
}

pub(crate) fn conformer_energy_span(
    row: &ReactionCsvRow,
    conformer_count: usize,
) -> Result<f32, String> {
    steric_x::profile_scope!("conformer_processing", "reaction::conformer_energy_span");
    if conformer_count == 0 {
        return Err("ensemble energy span requires at least one conformer".into());
    }
    let Some(value) = row.conformer_relative_energies.as_deref() else {
        return Ok(0.0);
    };
    if value.trim().is_empty() {
        return Ok(0.0);
    }
    let energies = parse_semicolon_floats(value, "Conformer_Relative_Energies_kcal_mol")?;
    if energies.len() != conformer_count {
        return Err(format!(
            "found {} conformer energies for {conformer_count} coordinate paths",
            energies.len()
        ));
    }
    if energies
        .iter()
        .any(|energy| !energy.is_finite() || *energy < 0.0)
    {
        return Err("relative conformer energies must be finite and non-negative".into());
    }
    let minimum = energies.iter().copied().fold(f32::INFINITY, f32::min);
    let maximum = energies.iter().copied().fold(f32::NEG_INFINITY, f32::max);
    Ok(maximum - minimum)
}

pub(crate) fn validate_reaction_row(
    row: &ReactionCsvRow,
    csv_line: usize,
) -> Result<(), Box<dyn Error>> {
    if row.reaction_id.trim().is_empty() {
        return Err(format!("CSV line {csv_line} has an empty Reaction_ID").into());
    }
    if row.ligand_xyz_path.as_os_str().is_empty() {
        return Err(format!("CSV line {csv_line} has an empty Ligand_XYZ_Path").into());
    }
    if row
        .ligand_xyz_path
        .components()
        .any(|component| component == Component::ParentDir)
    {
        return Err(format!("CSV line {csv_line} contains a parent-directory XYZ path").into());
    }
    if !row.nbo_charge.is_finite()
        || !row.ir_freq.is_finite()
        || !row.temp_k.is_finite()
        || !row.exp_ddg.is_finite()
    {
        return Err(format!("CSV line {csv_line} contains a non-finite numeric value").into());
    }
    if row.temp_k <= 0.0 {
        return Err(format!("CSV line {csv_line} has a non-positive temperature").into());
    }
    Ok(())
}

pub(crate) fn resolve_xyz_path(xyz_dir: &Path, csv_path: &Path) -> Result<PathBuf, String> {
    steric_x::profile_scope!("file_parsing", "reaction::resolve_xyz_path");
    if csv_path.is_absolute() {
        return csv_path
            .is_file()
            .then(|| csv_path.to_owned())
            .ok_or_else(|| format!("XYZ file does not exist: {}", csv_path.display()));
    }

    let direct = xyz_dir.join(csv_path);
    if direct.is_file() {
        return Ok(direct);
    }

    if let (Some(directory_name), Some(first_component)) = (
        xyz_dir.file_name(),
        csv_path
            .components()
            .next()
            .and_then(|component| match component {
                Component::Normal(value) => Some(value),
                _ => None,
            }),
    ) && directory_name == first_component
    {
        let stripped: PathBuf = csv_path.components().skip(1).collect();
        let candidate = xyz_dir.join(stripped);
        if candidate.is_file() {
            return Ok(candidate);
        }
    }

    if let Some(filename) = csv_path.file_name() {
        let candidate = xyz_dir.join(filename);
        if candidate.is_file() {
            return Ok(candidate);
        }
    }
    Err(format!(
        "could not resolve XYZ path `{}` below {}",
        csv_path.display(),
        xyz_dir.display()
    ))
}

pub(crate) fn load_reaction_metadata(
    path: &Path,
) -> Result<Vec<ReactionMetadataRow>, Box<dyn Error>> {
    steric_x::profile_scope!("file_parsing", "reaction::load_reaction_metadata");
    let mut reader = csv::ReaderBuilder::new()
        .trim(csv::Trim::All)
        .from_path(path)?;
    let rows = reader
        .deserialize::<ReactionMetadataRow>()
        .enumerate()
        .map(|(index, row)| {
            let row = row.map_err(|error| {
                format!(
                    "{} row {} could not be parsed: {error}",
                    path.display(),
                    index + 2
                )
            })?;
            if row.reaction_id.trim().is_empty() {
                return Err(format!(
                    "{} row {} has no Reaction_ID",
                    path.display(),
                    index + 2
                ));
            }
            if !is_supported_split(&row.dataset_split) {
                return Err(format!(
                    "{} row {} has unsupported Dataset_Split `{}`",
                    path.display(),
                    index + 2,
                    row.dataset_split
                ));
            }
            Ok(row)
        })
        .collect::<Result<Vec<_>, String>>()?;
    if rows.is_empty() {
        return Err(format!("metadata CSV contains no rows: {}", path.display()).into());
    }
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::Vec3;
    use steric_x::{Atom, Molecule};

    fn sample_row() -> ReactionCsvRow {
        ReactionCsvRow {
            reaction_id: "RXN-001".into(),
            ligand_xyz_path: "ligand.xyz".into(),
            attach_atom_idx: 0,
            primary_bond_vector_idx: 1,
            nbo_charge: -0.3,
            ir_freq: 1_650.0,
            temp_k: 298.15,
            exp_ddg: 1.2,
            conformer_xyz_paths: None,
            conformer_relative_energies: None,
            conformer_boltzmann_weights: None,
            conformer_coordination_centers: None,
            coordination_center_method: None,
        }
    }

    #[test]
    fn creates_record_from_csv_row() {
        let molecule = Molecule {
            atoms: vec![
                Atom::new("H", Vec3::ZERO),
                Atom::new("C", Vec3::new(2.0, 0.0, 0.0)),
            ],
        };
        let row = sample_row();
        let params = sterimol_from_molecule(&molecule, &row).unwrap();
        let record = record_from_ensemble(&[params], &[1.0], 0.0, &row).unwrap();
        assert!((record.l - 3.7).abs() < 1.0e-5);
        assert_eq!(record.nbo_charge, -0.3);
        assert_eq!(record.exp_ddg, 1.2);
        assert_eq!(record.conformer_count(), 1);
    }
    #[test]
    fn supplied_weights_keep_scale_and_small_contributions() {
        let params = [
            SterimolParams {
                l: 2.0,
                b1: 4.0,
                b5: 6.0,
            },
            SterimolParams {
                l: 6.0,
                b1: 8.0,
                b5: 10.0,
            },
        ];
        for text in ["1;3", "1e-10;3e-10", "1e38;3e38"] {
            let mut row = sample_row();
            row.conformer_boltzmann_weights = Some(text.into());
            let weights = conformer_weights(&row, 2).unwrap();
            let mean = weighted_sterimol(&params, &weights).unwrap();
            assert!((mean.l - 5.0).abs() < 1e-6);
            assert!((mean.b1 - 7.0).abs() < 1e-6);
            assert!((mean.b5 - 9.0).abs() < 1e-6);
        }
        let mut row = sample_row();
        row.conformer_boltzmann_weights = Some("1.40129846e-45;2".into());
        let weights = conformer_weights(&row, 2).unwrap();
        let mean = weighted_sterimol(
            &[
                SterimolParams {
                    l: f32::MAX,
                    b1: 0.0,
                    b5: 0.0,
                },
                SterimolParams::default(),
            ],
            &weights,
        )
        .unwrap();
        assert_eq!(
            mean.l,
            (f64::from(f32::MAX) * f64::from(f32::from_bits(1)) / 2.0) as f32
        );
        assert!(mean.l > 0.0);
    }

    #[test]
    fn invalid_weights_and_descriptors_are_not_silently_accepted() {
        for weights in [
            [-0.5, 1.5],
            [0.0, 0.0],
            [f32::NAN, 1.0],
            [f32::INFINITY, 1.0],
        ] {
            assert!(weighted_sterimol(&[SterimolParams::default(); 2], &weights).is_err());
        }
        for bad in [f32::NAN, f32::INFINITY, f32::NEG_INFINITY] {
            for p in [
                SterimolParams {
                    l: bad,
                    ..SterimolParams::default()
                },
                SterimolParams {
                    b1: bad,
                    ..SterimolParams::default()
                },
                SterimolParams {
                    b5: bad,
                    ..SterimolParams::default()
                },
            ] {
                assert!(weighted_sterimol(&[p], &[1.0]).is_err());
            }
        }
    }

    #[test]
    fn energy_span_is_offset_invariant_and_missing_weights_remain_uniform() {
        let mut row = sample_row();
        for text in ["0;1", "1;2", "10;11"] {
            row.conformer_relative_energies = Some(text.into());
            assert_eq!(conformer_energy_span(&row, 2).unwrap(), 1.0);
            assert_eq!(
                normalized_weights(&conformer_weights(&row, 2).unwrap()).unwrap(),
                vec![0.5, 0.5]
            );
        }
        assert!(conformer_energy_span(&row, 0).is_err());
        row.conformer_relative_energies = Some("-1;0".into());
        assert!(conformer_energy_span(&row, 2).is_err());
    }
    #[test]
    fn reaction_axis_validation_accepts_nonzero_short_and_subnormal_axes() {
        let row = sample_row();
        for length in [0.000_345_266_95_f32, f32::from_bits(1)] {
            let molecule = Molecule {
                atoms: vec![
                    Atom::new("H", Vec3::ZERO),
                    Atom::new("C", Vec3::new(0.0, 0.0, length)),
                ],
            };
            let result = sterimol_from_molecule(&molecule, &row).unwrap();
            // A sphere centered on the explicit +Z axis has L=z+r, B1=B5=r.
            assert_eq!(result.l, length + 1.7_f32);
            assert_eq!(result.b1, 1.7);
            assert_eq!(result.b5, 1.7);
        }
        let coincident = Molecule {
            atoms: vec![Atom::new("H", Vec3::ZERO), Atom::new("C", Vec3::ZERO)],
        };
        assert!(sterimol_from_molecule(&coincident, &row).is_err());
    }
}
