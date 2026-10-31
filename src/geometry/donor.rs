//! Shared donor selection for the legacy CLI and conformer API.

use super::{Molecule, bonded_neighbors};

/// The donor atom, its bonded substituents, and the substituent chosen as the
/// primary reference axis — everything the descriptor kernels need, recovered
/// from geometry alone.
#[derive(Clone, Copy, Debug)]
pub struct DonorTopology {
    pub donor_idx: usize,
    pub reference_idx: usize,
    pub substituents: [usize; 3],
}

/// Auto-detect the trivalent donor and its three bonded substituents.
///
/// The donor is the sole atom of `donor_element` unless `explicit_index` names
/// one. Substituents are the covalently bonded atoms (hydrogens included),
/// identical to the kernel's own frame construction, and the primary axis is the
/// nearest bonded heavy atom, unless tied (explicit reference required).
pub fn detect_trivalent_donor(
    molecule: &Molecule,
    bonds: Option<&[[usize; 2]]>,
    donor_element: &str,
    explicit_index: Option<usize>,
    explicit_reference: Option<usize>,
    coordination_axis: bool,
) -> Result<DonorTopology, String> {
    crate::profile_scope!("donor_bond_detection", "descriptors::detect_donor");
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

    if let Some(bonds) = bonds {
        if bonds
            .iter()
            .any(|&[a, b]| a >= molecule.atoms.len() || b >= molecule.atoms.len() || a == b)
        {
            return Err("invalid explicit bond indices".into());
        }
    }
    let donor = &molecule.atoms[donor_idx];
    // Shared covalent-radius frame (already sorted by ascending distance/index),
    // identical to the buried-volume kernel's own substituent detection.
    let mut bonded = if let Some(bonds) = bonds {
        bonds
            .iter()
            .filter_map(|&[a, b]| {
                let index = if a == donor_idx {
                    b
                } else if b == donor_idx {
                    a
                } else {
                    return None;
                };
                Some((
                    molecule.atoms[index]
                        .position
                        .distance_squared(donor.position),
                    index,
                ))
            })
            .collect::<Vec<_>>()
    } else {
        bonded_neighbors(molecule, donor_idx)
    };
    bonded.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
    if bonded.len() != 3 {
        return Err(format!(
            "{} donor (atom {donor_idx}) is not three-coordinate — found {} bonded atoms; \
             the descriptor is defined for trivalent (e.g. phosphine) donors",
            donor.element,
            bonded.len()
        ));
    }
    let reference_idx = if let Some(index) = explicit_reference {
        if !bonded.iter().any(|(_, neighbor)| *neighbor == index) {
            return Err("--reference-index must name a bonded donor substituent".into());
        }
        index
    } else if coordination_axis {
        // All three planes are reduced symmetrically; no unique heavy bond axis
        // is needed for this physically distinct coordination-axis descriptor.
        bonded[0].1
    } else {
        let heavy = bonded
            .iter()
            .filter(|(_, index)| !molecule.atoms[*index].element.eq_ignore_ascii_case("H"))
            .collect::<Vec<_>>();
        let Some(first) = heavy.first() else {
            return Err(format!(
                "{} donor (atom {donor_idx}) has no bonded heavy atom to define an axis; pass --reference-index or use coordination mode",
                donor.element
            ));
        };
        if heavy.get(1).is_some_and(|next| next.0 == first.0) {
            return Err(
                "nearest bonded heavy-atom axis is ambiguous; pass --reference-index".into(),
            );
        }
        first.1
    };
    Ok(DonorTopology {
        donor_idx,
        reference_idx,
        substituents: [bonded[0].1, bonded[1].1, bonded[2].1],
    })
}
