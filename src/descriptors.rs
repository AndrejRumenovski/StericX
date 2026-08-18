//! The `descriptors` subcommand: geometry-only donor detection, ensemble-averaged
//! Sterimol / buried-volume / pyramidalization descriptors, and their text, JSON,
//! and CSV emitters.

use crate::cli::{DescriptorFormat, STERIMOL_L_CORRECTION, SterimolAxis};
use serde::Serialize;
use std::error::Error;
use std::path::{Path, PathBuf};
use steric_x::{
    BuriedVolumeCalculator, BuriedVolumeConfig, Molecule, PyramidalizationCalculator,
    SterimolCalculator, bonded_neighbors, coordination_center, parse_coordinate_file,
};

/// The donor atom, its bonded substituents, and the substituent chosen as the
/// primary reference axis — everything the descriptor kernels need, recovered
/// from geometry alone.
#[derive(Clone, Copy, Debug)]
struct DonorTopology {
    donor_idx: usize,
    reference_idx: usize,
    substituents: [usize; 3],
}

/// Auto-detect the trivalent donor and its three bonded substituents.
///
/// The donor is the sole atom of `donor_element` unless `explicit_index` names
/// one. Substituents are the covalently bonded atoms (hydrogens included),
/// identical to the kernel's own frame construction, and the primary axis is the
/// nearest bonded heavy atom.
fn detect_donor(
    molecule: &Molecule,
    donor_element: &str,
    explicit_index: Option<usize>,
) -> Result<DonorTopology, String> {
    let donor_idx = match explicit_index {
        Some(index) => {
            if index >= molecule.atoms.len() {
                return Err(format!(
                    "--donor-index {index} is out of bounds for {} atoms",
                    molecule.atoms.len()
                ));
            }
            index
        }
        None => {
            let matches = molecule
                .atoms
                .iter()
                .enumerate()
                .filter(|(_, atom)| atom.element.eq_ignore_ascii_case(donor_element))
                .map(|(index, _)| index)
                .collect::<Vec<_>>();
            match matches.as_slice() {
                [single] => *single,
                [] => {
                    return Err(format!(
                        "no {donor_element} donor atom found; pass --donor-element or --donor-index"
                    ));
                }
                many => {
                    return Err(format!(
                        "found {} {donor_element} atoms; pass --donor-index to choose the donor",
                        many.len()
                    ));
                }
            }
        }
    };

    let donor = &molecule.atoms[donor_idx];
    // Shared covalent-radius frame (already sorted by ascending distance/index),
    // identical to the buried-volume kernel's own substituent detection.
    let bonded = bonded_neighbors(molecule, donor_idx);
    if bonded.len() != 3 {
        return Err(format!(
            "{} donor (atom {donor_idx}) is not three-coordinate — found {} bonded atoms; \
             the descriptor is defined for trivalent (e.g. phosphine) donors",
            donor.element,
            bonded.len()
        ));
    }
    let reference_idx = bonded
        .iter()
        .find(|(_, index)| !molecule.atoms[*index].element.eq_ignore_ascii_case("H"))
        .map(|(_, index)| *index)
        .ok_or_else(|| {
            format!(
                "{} donor (atom {donor_idx}) has no bonded heavy atom to define an axis",
                donor.element
            )
        })?;
    Ok(DonorTopology {
        donor_idx,
        reference_idx,
        substituents: [bonded[0].1, bonded[1].1, bonded[2].1],
    })
}

/// One file's descriptors, ready for text, JSON, or CSV emission.
#[derive(Debug, Serialize)]
pub(crate) struct DescriptorResult {
    pub(crate) file: String,
    pub(crate) conformers: usize,
    pub(crate) donor_element: String,
    pub(crate) donor_index: usize,
    pub(crate) substituents: Vec<String>,
    pub(crate) sterimol_l: f32,
    pub(crate) sterimol_b1: f32,
    pub(crate) sterimol_b5: f32,
    pub(crate) percent_buried_volume: f32,
    pub(crate) buried_volume: f32,
    pub(crate) qvbur_min: f32,
    pub(crate) qvbur_max: f32,
    pub(crate) max_delta_qvbur: f32,
    /// Kraken's headline `vbur_max_delta_qvbur_min`: the minimum over conformers.
    pub(crate) max_delta_qvbur_min: f32,
    /// Radhakrishnan pyramidalization `P` (Kraken `pyr_P`), conformer mean.
    pub(crate) pyr_p: f32,
    /// Mean pyramidalization angle in degrees (Kraken `pyr_alpha`), conformer mean.
    pub(crate) pyr_alpha: f32,
}

/// Compute ensemble-averaged descriptors for a single ligand file.
pub(crate) fn descriptors_for_file(
    path: &Path,
    donor_element: &str,
    donor_index: Option<usize>,
    sterimol_axis: SterimolAxis,
    config: BuriedVolumeConfig,
) -> Result<DescriptorResult, String> {
    let conformers = parse_coordinate_file(path)
        .map_err(|error| format!("failed to read {}: {error}", path.display()))?;
    if conformers.is_empty() {
        return Err(format!("{} contains no geometries", path.display()));
    }

    let topology = detect_donor(&conformers[0], donor_element, donor_index)?;
    let substituents = topology
        .substituents
        .iter()
        .map(|index| conformers[0].atoms[*index].element.clone())
        .collect::<Vec<_>>();

    let mut sterimol = Vec::with_capacity(conformers.len());
    let mut buried = Vec::with_capacity(conformers.len());
    let mut pyramidalization = Vec::with_capacity(conformers.len());
    for (conformer_index, molecule) in conformers.iter().enumerate() {
        // Re-detect per conformer so this tolerates files whose models differ.
        let topology = detect_donor(molecule, donor_element, donor_index).map_err(|message| {
            format!(
                "{} (conformer {conformer_index}): {message}",
                path.display()
            )
        })?;
        pyramidalization.push(PyramidalizationCalculator::compute(
            molecule,
            topology.donor_idx,
            topology.substituents,
        ));
        sterimol.push(match sterimol_axis {
            SterimolAxis::Bond => {
                SterimolCalculator::compute(molecule, topology.donor_idx, topology.reference_idx)
            }
            SterimolAxis::Coordination => {
                let center = coordination_center(
                    molecule,
                    topology.donor_idx,
                    topology.reference_idx,
                    config,
                )
                .map_err(|error| {
                    format!("{} (conformer {conformer_index}): {error}", path.display())
                })?;
                let mut params =
                    SterimolCalculator::compute_with_dummy(molecule, topology.donor_idx, center);
                params.l += STERIMOL_L_CORRECTION;
                params
            }
        });
        buried.push(
            BuriedVolumeCalculator::compute(
                molecule,
                topology.donor_idx,
                topology.reference_idx,
                config,
            )
            .map_err(|error| {
                format!("{} (conformer {conformer_index}): {error}", path.display())
            })?,
        );
    }

    let max_delta_qvbur_min = buried
        .iter()
        .map(|params| params.max_delta_qvbur)
        .fold(f32::INFINITY, f32::min);

    Ok(DescriptorResult {
        file: path.display().to_string(),
        conformers: conformers.len(),
        donor_element: conformers[0].atoms[topology.donor_idx].element.clone(),
        donor_index: topology.donor_idx,
        substituents,
        sterimol_l: conformer_mean(&sterimol, |params| params.l),
        sterimol_b1: conformer_mean(&sterimol, |params| params.b1),
        sterimol_b5: conformer_mean(&sterimol, |params| params.b5),
        percent_buried_volume: conformer_mean(&buried, |params| params.percent_buried_volume),
        buried_volume: conformer_mean(&buried, |params| params.buried_volume),
        qvbur_min: conformer_mean(&buried, |params| params.qvbur_min),
        qvbur_max: conformer_mean(&buried, |params| params.qvbur_max),
        max_delta_qvbur: conformer_mean(&buried, |params| params.max_delta_qvbur),
        max_delta_qvbur_min,
        pyr_p: conformer_mean(&pyramidalization, |params| params.pyr_p),
        pyr_alpha: conformer_mean(&pyramidalization, |params| params.pyr_alpha),
    })
}

/// Arithmetic mean of one descriptor field over a conformer ensemble.
fn conformer_mean<T>(items: &[T], select: impl Fn(&T) -> f32) -> f32 {
    items.iter().map(select).sum::<f32>() / items.len() as f32
}

fn print_descriptor_text(result: &DescriptorResult) {
    let ensemble = result.conformers > 1;
    let qualifier = if ensemble { " (conformer mean)" } else { "" };
    println!("{}", result.file);
    println!(
        "  donor          {} (atom {})",
        result.donor_element, result.donor_index
    );
    println!("  substituents   {}", result.substituents.join(", "));
    println!("  conformers     {}", result.conformers);
    println!(
        "  Sterimol{}      L {:.2}   B1 {:.2}   B5 {:.2}   Å",
        qualifier, result.sterimol_l, result.sterimol_b1, result.sterimol_b5
    );
    println!(
        "  buried volume{}  Vbur {:.1}%   ({:.1} Å³)",
        qualifier, result.percent_buried_volume, result.buried_volume
    );
    println!(
        "                 qvbur_min {:.2}   qvbur_max {:.2}   max_delta_qvbur {:.2}   Å³",
        result.qvbur_min, result.qvbur_max, result.max_delta_qvbur
    );
    if ensemble {
        println!(
            "                 max_delta_qvbur_min {:.2} Å³   [Kraken vbur_max_delta_qvbur_min]",
            result.max_delta_qvbur_min
        );
    }
    println!(
        "  pyramidalization{} pyr_P {:.3}   pyr_alpha {:.2}°",
        qualifier, result.pyr_p, result.pyr_alpha
    );
}

/// A CSV column: its header and a formatter for one result's value.
type CsvColumn = (&'static str, fn(&DescriptorResult) -> String);

fn print_descriptor_csv(results: &[DescriptorResult]) {
    // Single source of truth: each column pairs its header with its value
    // formatter, so the header row and every data row derive from one list and
    // adding a descriptor cannot desynchronise column names, order, or count.
    let columns: [CsvColumn; 16] = [
        ("file", |r| csv_field(&r.file)),
        ("conformers", |r| r.conformers.to_string()),
        ("donor_element", |r| r.donor_element.clone()),
        ("donor_index", |r| r.donor_index.to_string()),
        ("substituents", |r| csv_field(&r.substituents.join(" "))),
        ("sterimol_l", |r| format!("{:.4}", r.sterimol_l)),
        ("sterimol_b1", |r| format!("{:.4}", r.sterimol_b1)),
        ("sterimol_b5", |r| format!("{:.4}", r.sterimol_b5)),
        ("percent_buried_volume", |r| {
            format!("{:.4}", r.percent_buried_volume)
        }),
        ("buried_volume", |r| format!("{:.4}", r.buried_volume)),
        ("qvbur_min", |r| format!("{:.4}", r.qvbur_min)),
        ("qvbur_max", |r| format!("{:.4}", r.qvbur_max)),
        ("max_delta_qvbur", |r| format!("{:.4}", r.max_delta_qvbur)),
        ("max_delta_qvbur_min", |r| {
            format!("{:.4}", r.max_delta_qvbur_min)
        }),
        ("pyr_p", |r| format!("{:.4}", r.pyr_p)),
        ("pyr_alpha", |r| format!("{:.4}", r.pyr_alpha)),
    ];
    println!(
        "{}",
        columns
            .iter()
            .map(|(name, _)| *name)
            .collect::<Vec<_>>()
            .join(",")
    );
    for result in results {
        println!(
            "{}",
            columns
                .iter()
                .map(|(_, format)| format(result))
                .collect::<Vec<_>>()
                .join(",")
        );
    }
}

/// Quote a CSV field only when it contains a comma, quote, or newline.
pub(crate) fn csv_field(value: &str) -> String {
    if value.contains([',', '"', '\n']) {
        format!("\"{}\"", value.replace('"', "\"\""))
    } else {
        value.to_owned()
    }
}

pub(crate) fn descriptors_command(
    inputs: &[PathBuf],
    donor_element: &str,
    donor_index: Option<usize>,
    sterimol_axis: SterimolAxis,
    format: DescriptorFormat,
    config: BuriedVolumeConfig,
) -> Result<(), Box<dyn Error>> {
    if inputs.len() > 1 && donor_index.is_some() {
        return Err("--donor-index applies to a single file; omit it for batch runs".into());
    }
    let mut results = Vec::with_capacity(inputs.len());
    let mut failures = 0_usize;
    for path in inputs {
        match descriptors_for_file(path, donor_element, donor_index, sterimol_axis, config) {
            Ok(result) => results.push(result),
            Err(message) => {
                eprintln!("skipped {}: {message}", path.display());
                failures += 1;
            }
        }
    }
    if results.is_empty() {
        return Err("no ligand files could be featurized".into());
    }

    match format {
        DescriptorFormat::Text => {
            for (index, result) in results.iter().enumerate() {
                if index > 0 {
                    println!();
                }
                print_descriptor_text(result);
            }
        }
        DescriptorFormat::Json => {
            println!("{}", serde_json::to_string_pretty(&results)?);
        }
        DescriptorFormat::Csv => print_descriptor_csv(&results),
    }
    if failures > 0 {
        eprintln!(
            "featurized {} of {} files ({failures} skipped)",
            results.len(),
            inputs.len()
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::temporary_directory;
    use glam::Vec3;
    use steric_x::{Atom, BuriedVolumeConfig, Molecule};

    fn molecule_from(atoms: &[(&str, f32, f32, f32)]) -> Molecule {
        Molecule {
            atoms: atoms
                .iter()
                .map(|(element, x, y, z)| Atom::new(*element, Vec3::new(*x, *y, *z)))
                .collect(),
        }
    }

    fn tertiary_phosphine() -> Molecule {
        molecule_from(&[
            ("P", 0.0, 0.0, 0.0),
            ("C", 1.4, 0.0, 0.45),
            ("C", -0.7, 1.212, 0.45),
            ("C", -0.7, -1.212, 0.45),
            // A distal, non-bonded carbon that must never be taken as a substituent.
            ("C", 2.8, 0.0, 0.7),
        ])
    }

    #[test]
    fn detect_donor_finds_the_single_phosphine() {
        let topology = detect_donor(&tertiary_phosphine(), "P", None).unwrap();
        assert_eq!(topology.donor_idx, 0);
        let mut substituents = topology.substituents;
        substituents.sort_unstable();
        assert_eq!(substituents, [1, 2, 3]);
        // The primary axis is a bonded heavy atom, never the distal carbon (4).
        assert!([1, 2, 3].contains(&topology.reference_idx));
    }

    #[test]
    fn detect_donor_counts_bonded_hydrogens_as_substituents() {
        // Primary phosphine R-PH2: the two bonded hydrogens are substituents,
        // and the reference axis falls on the lone bonded carbon.
        let molecule = molecule_from(&[
            ("P", 0.0, 0.0, 0.0),
            ("C", 1.5, 0.0, 0.6),
            ("H", -0.6, 1.1, 0.55),
            ("H", -0.6, -1.1, 0.55),
        ]);
        let topology = detect_donor(&molecule, "P", None).unwrap();
        let mut substituents = topology.substituents;
        substituents.sort_unstable();
        assert_eq!(substituents, [1, 2, 3]);
        assert_eq!(topology.reference_idx, 1);
    }

    #[test]
    fn detect_donor_reports_missing_and_ambiguous_donors() {
        let no_phosphorus = molecule_from(&[("C", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 1.1)]);
        let error = detect_donor(&no_phosphorus, "P", None).unwrap_err();
        assert!(error.contains("no P donor"));

        let two_phosphorus = molecule_from(&[
            ("P", 0.0, 0.0, 0.0),
            ("C", 1.5, 0.0, 0.0),
            ("C", -1.5, 0.0, 0.0),
            ("C", 0.0, 1.5, 0.0),
            ("P", 6.0, 0.0, 0.0),
            ("C", 7.5, 0.0, 0.0),
            ("C", 4.5, 0.0, 0.0),
            ("C", 6.0, 1.5, 0.0),
        ]);
        let error = detect_donor(&two_phosphorus, "P", None).unwrap_err();
        assert!(error.contains("--donor-index"));
        // Naming one of them resolves the ambiguity.
        assert_eq!(
            detect_donor(&two_phosphorus, "P", Some(4))
                .unwrap()
                .donor_idx,
            4
        );
    }

    #[test]
    fn detect_donor_rejects_a_non_trivalent_donor() {
        let two_coordinate = molecule_from(&[
            ("P", 0.0, 0.0, 0.0),
            ("C", 1.5, 0.0, 0.0),
            ("H", -0.9, 0.9, 0.0),
        ]);
        let error = detect_donor(&two_coordinate, "P", None).unwrap_err();
        assert!(error.contains("three-coordinate"));
    }

    #[test]
    fn detect_donor_supports_a_non_phosphorus_element() {
        // A trivalent amine donor located by --donor-element N.
        let amine = molecule_from(&[
            ("N", 0.0, 0.0, 0.0),
            ("C", 1.47, 0.0, 0.3),
            ("C", -0.7, 1.2, 0.3),
            ("C", -0.7, -1.2, 0.3),
        ]);
        let topology = detect_donor(&amine, "N", None).unwrap();
        assert_eq!(topology.donor_idx, 0);
    }

    #[test]
    fn descriptors_for_file_featurizes_an_xyz_geometry() {
        let directory = temporary_directory("descriptors_xyz");
        std::fs::create_dir_all(&directory).unwrap();
        let path = directory.join("phosphine.xyz");
        std::fs::write(
            &path,
            "5\nphosphine\nP 0 0 0\nC 1.4 0 0.45\nC -0.7 1.212 0.45\n\
             C -0.7 -1.212 0.45\nC 2.8 0 0.7\n",
        )
        .unwrap();
        let result = descriptors_for_file(
            &path,
            "P",
            None,
            SterimolAxis::Bond,
            BuriedVolumeConfig::default(),
        )
        .unwrap();
        assert_eq!(result.donor_element, "P");
        assert_eq!(result.conformers, 1);
        assert_eq!(result.substituents.len(), 3);
        assert!(result.buried_volume.is_finite() && result.buried_volume > 0.0);
        assert!(result.max_delta_qvbur_min.is_finite());
        // With one conformer the ensemble minimum equals the single value.
        assert_eq!(result.max_delta_qvbur, result.max_delta_qvbur_min);

        // The coordination axis is a different frame, so it yields a different
        // (also valid) Sterimol L, and applies the +0.40 Å correction.
        let coordination = descriptors_for_file(
            &path,
            "P",
            None,
            SterimolAxis::Coordination,
            BuriedVolumeConfig::default(),
        )
        .unwrap();
        assert!(coordination.sterimol_l.is_finite() && coordination.sterimol_l > 0.0);
        assert!((coordination.sterimol_l - result.sterimol_l).abs() > 1.0e-4);
        std::fs::remove_dir_all(&directory).ok();
    }

    #[test]
    fn csv_field_quotes_only_when_required() {
        assert_eq!(csv_field("plain"), "plain");
        assert_eq!(csv_field("a,b"), "\"a,b\"");
        assert_eq!(csv_field("say \"hi\""), "\"say \"\"hi\"\"\"");
    }
}
