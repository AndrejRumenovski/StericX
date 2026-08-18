use super::{Molecule, covalent_radius};
use glam::Vec3;
use std::error::Error;
use std::f32::consts::PI;
use std::fmt::{Display, Formatter};

/// Kraken-compatible geometric settings for a metal-centred buried-volume scan.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BuriedVolumeConfig {
    /// Radius of the integration sphere in ångströms.
    pub sphere_radius: f32,
    /// Volume represented by one grid point in Å³.
    pub density: f32,
    /// Distance from the donor atom to the virtual metal centre in ångströms.
    pub center_distance: f32,
    /// Scale applied to Bondi-style atomic radii.
    pub radii_scale: f32,
    /// Whether hydrogen atoms contribute to occupied volume.
    pub include_hydrogens: bool,
}

impl Default for BuriedVolumeConfig {
    fn default() -> Self {
        Self {
            sphere_radius: 3.5,
            density: 0.01,
            center_distance: 2.1,
            radii_scale: 1.17,
            include_hydrogens: false,
        }
    }
}

/// Buried-volume descriptors for one conformer.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BuriedVolumeParams {
    /// Ligand-occupied volume inside the integration sphere in Å³.
    pub buried_volume: f32,
    /// Occupied fraction of the integration sphere in percent.
    pub percent_buried_volume: f32,
    /// Smallest quadrant occupied volume across all donor-substituent orientations.
    pub qvbur_min: f32,
    /// Largest quadrant occupied volume across all donor-substituent orientations.
    pub qvbur_max: f32,
    /// Largest difference between adjacent quadrant occupied volumes.
    pub max_delta_qvbur: f32,
    /// Smallest octant occupied volume across all orientations.
    pub ovbur_min: f32,
    /// Largest octant occupied volume across all orientations.
    pub ovbur_max: f32,
    /// Occupied volume in the donor-facing hemisphere.
    pub near_vbur: f32,
    /// Occupied volume in the distal hemisphere.
    pub far_vbur: f32,
}

/// Conformer-ensemble aggregation used by Kraken-style descriptor tables.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BuriedVolumeEnsembleParams {
    pub vbur_boltz: f32,
    pub vbur_min: f32,
    pub vbur_max: f32,
    pub vbur_delta: f32,
    pub qvbur_min_boltz: f32,
    pub qvbur_max_boltz: f32,
    pub max_delta_qvbur_boltz: f32,
    pub max_delta_qvbur_min: f32,
    pub max_delta_qvbur_max: f32,
    pub max_delta_qvbur_delta: f32,
    pub max_delta_qvbur_vburminconf: f32,
    pub near_vbur_boltz: f32,
    pub far_vbur_boltz: f32,
    pub conformer_count: usize,
}

/// Input or geometry failure in a buried-volume calculation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BuriedVolumeError(String);

impl Display for BuriedVolumeError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        formatter.write_str(&self.0)
    }
}

impl Error for BuriedVolumeError {}

/// Deterministic voxel implementation of Kraken's Morfeus buried-volume protocol.
#[derive(Clone, Copy, Debug, Default)]
pub struct BuriedVolumeCalculator;

impl BuriedVolumeCalculator {
    /// Calculate one conformer's coordination-aware buried-volume descriptors.
    ///
    /// The virtual centre is placed 2.1 Å from the donor along a geometrically
    /// inferred lone-pair direction. This is a documented approximation to
    /// Kraken's xTB localized-molecular-orbital centre. The three covalently
    /// bonded donor substituents (see [`donor_neighbor_indices`], hydrogens
    /// included) are each used to define the XZ plane, matching Kraken's
    /// three-orientation quadrant scan.
    pub fn compute(
        molecule: &Molecule,
        donor_idx: usize,
        reference_neighbor_idx: usize,
        config: BuriedVolumeConfig,
    ) -> Result<BuriedVolumeParams, BuriedVolumeError> {
        validate_config(config)?;
        let donor = molecule.atoms.get(donor_idx).ok_or_else(|| {
            BuriedVolumeError(format!("donor index {donor_idx} is out of bounds"))
        })?;
        if !donor.position.is_finite() {
            return Err(BuriedVolumeError(
                "donor coordinate is not finite".to_owned(),
            ));
        }

        let neighbor_indices = donor_neighbor_indices(molecule, donor_idx, reference_neighbor_idx)?;
        let lone_pair_direction =
            infer_lone_pair_direction(molecule, donor_idx, &neighbor_indices)?;
        let center = donor.position + config.center_distance * lone_pair_direction;
        Self::compute_from_center(molecule, donor_idx, neighbor_indices, center, config)
    }

    /// Calculate buried volume around an explicitly supplied coordination center.
    ///
    /// This path is used for xTB localized-molecular-orbital centers and avoids
    /// the geometric lone-pair approximation in [`Self::compute`].
    pub fn compute_with_center(
        molecule: &Molecule,
        donor_idx: usize,
        reference_neighbor_idx: usize,
        center: Vec3,
        config: BuriedVolumeConfig,
    ) -> Result<BuriedVolumeParams, BuriedVolumeError> {
        validate_config(config)?;
        let donor = molecule.atoms.get(donor_idx).ok_or_else(|| {
            BuriedVolumeError(format!("donor index {donor_idx} is out of bounds"))
        })?;
        if !donor.position.is_finite() || !center.is_finite() {
            return Err(BuriedVolumeError(
                "donor or coordination-center coordinate is not finite".to_owned(),
            ));
        }
        if center.distance_squared(donor.position) <= f32::EPSILON {
            return Err(BuriedVolumeError(
                "coordination center coincides with donor atom".to_owned(),
            ));
        }
        let neighbor_indices = donor_neighbor_indices(molecule, donor_idx, reference_neighbor_idx)?;
        Self::compute_from_center(molecule, donor_idx, neighbor_indices, center, config)
    }

    fn compute_from_center(
        molecule: &Molecule,
        donor_idx: usize,
        neighbor_indices: [usize; 3],
        center: Vec3,
        config: BuriedVolumeConfig,
    ) -> Result<BuriedVolumeParams, BuriedVolumeError> {
        let sphere = integration_grid(config);
        if sphere.is_empty() {
            return Err(BuriedVolumeError(
                "integration grid contains no points".to_owned(),
            ));
        }

        let mut result = BuriedVolumeParams {
            qvbur_min: f32::INFINITY,
            ovbur_min: f32::INFINITY,
            ..BuriedVolumeParams::default()
        };
        for (orientation_index, &plane_idx) in neighbor_indices.iter().enumerate() {
            let basis = coordinate_basis(molecule, donor_idx, plane_idx, center)?;
            let atoms = aligned_atoms(molecule, center, basis, config)?;
            let volumes = occupied_volumes(&sphere, &atoms, config.sphere_radius);
            if orientation_index == 0 {
                result.buried_volume = volumes.buried_volume;
                result.percent_buried_volume =
                    100.0 * volumes.buried_volume / sphere_volume(config.sphere_radius);
                result.near_vbur = volumes.near_vbur;
                result.far_vbur = volumes.far_vbur;
            }
            result.qvbur_min = result
                .qvbur_min
                .min(volumes.quadrants.into_iter().fold(f32::INFINITY, f32::min));
            result.qvbur_max = result
                .qvbur_max
                .max(volumes.quadrants.into_iter().fold(0.0, f32::max));
            result.ovbur_min = result
                .ovbur_min
                .min(volumes.octants.into_iter().fold(f32::INFINITY, f32::min));
            result.ovbur_max = result
                .ovbur_max
                .max(volumes.octants.into_iter().fold(0.0, f32::max));
            let q = volumes.quadrants;
            let adjacent_delta = (0..4)
                .map(|index| (q[index] - q[(index + 3) % 4]).abs())
                .fold(0.0_f32, f32::max);
            result.max_delta_qvbur = result.max_delta_qvbur.max(adjacent_delta);
        }
        // A donor that occupies volume but yields perfectly equal quadrants in
        // every orientation is unphysical: it signals a collapsed frame (e.g.
        // only the donor atom fell inside the sphere) rather than a genuine
        // symmetric ligand. Refuse to emit that spurious zero so the ensemble
        // minimum is never poisoned by it.
        if result.buried_volume > 0.0 && result.max_delta_qvbur == 0.0 {
            return Err(BuriedVolumeError(
                "degenerate coordination frame produced a symmetric zero max_delta_qvbur"
                    .to_owned(),
            ));
        }
        Ok(result)
    }

    /// Aggregate conformer descriptors using normalized non-negative weights.
    pub fn aggregate(
        conformers: &[BuriedVolumeParams],
        weights: &[f32],
    ) -> Result<BuriedVolumeEnsembleParams, BuriedVolumeError> {
        if conformers.is_empty() || conformers.len() != weights.len() {
            return Err(BuriedVolumeError(
                "conformers and weights must have equal non-zero lengths".to_owned(),
            ));
        }
        if weights
            .iter()
            .any(|weight| !weight.is_finite() || *weight < 0.0)
        {
            return Err(BuriedVolumeError(
                "weights must be finite and non-negative".to_owned(),
            ));
        }
        let weight_sum = weights.iter().sum::<f32>();
        if !weight_sum.is_finite() || weight_sum <= f32::EPSILON {
            return Err(BuriedVolumeError(
                "conformer weights have zero total".to_owned(),
            ));
        }
        let normalized = weights.iter().map(|weight| weight / weight_sum);
        let mut ensemble = BuriedVolumeEnsembleParams {
            vbur_min: f32::INFINITY,
            max_delta_qvbur_min: f32::INFINITY,
            conformer_count: conformers.len(),
            ..BuriedVolumeEnsembleParams::default()
        };
        let mut minimum_vbur_index = 0;
        for (index, (params, weight)) in conformers.iter().zip(normalized).enumerate() {
            if params.buried_volume < conformers[minimum_vbur_index].buried_volume {
                minimum_vbur_index = index;
            }
            ensemble.vbur_boltz += params.buried_volume * weight;
            ensemble.qvbur_min_boltz += params.qvbur_min * weight;
            ensemble.qvbur_max_boltz += params.qvbur_max * weight;
            ensemble.max_delta_qvbur_boltz += params.max_delta_qvbur * weight;
            ensemble.near_vbur_boltz += params.near_vbur * weight;
            ensemble.far_vbur_boltz += params.far_vbur * weight;
            ensemble.vbur_min = ensemble.vbur_min.min(params.buried_volume);
            ensemble.vbur_max = ensemble.vbur_max.max(params.buried_volume);
            ensemble.max_delta_qvbur_min = ensemble.max_delta_qvbur_min.min(params.max_delta_qvbur);
            ensemble.max_delta_qvbur_max = ensemble.max_delta_qvbur_max.max(params.max_delta_qvbur);
        }
        ensemble.vbur_delta = ensemble.vbur_max - ensemble.vbur_min;
        ensemble.max_delta_qvbur_delta =
            ensemble.max_delta_qvbur_max - ensemble.max_delta_qvbur_min;
        ensemble.max_delta_qvbur_vburminconf = conformers[minimum_vbur_index].max_delta_qvbur;
        Ok(ensemble)
    }
}

/// The virtual coordination centre for a donor: the point a fixed distance from
/// the donor along its geometrically inferred lone-pair direction, where a
/// coordinating metal (or a Sterimol dummy) sits.
///
/// This is the same centre [`BuriedVolumeCalculator::compute`] integrates
/// around, exposed so callers can place a metal probe or Sterimol dummy there.
pub fn coordination_center(
    molecule: &Molecule,
    donor_idx: usize,
    reference_neighbor_idx: usize,
    config: BuriedVolumeConfig,
) -> Result<Vec3, BuriedVolumeError> {
    validate_config(config)?;
    let donor = molecule
        .atoms
        .get(donor_idx)
        .ok_or_else(|| BuriedVolumeError(format!("donor index {donor_idx} is out of bounds")))?;
    if !donor.position.is_finite() {
        return Err(BuriedVolumeError(
            "donor coordinate is not finite".to_owned(),
        ));
    }
    let neighbor_indices = donor_neighbor_indices(molecule, donor_idx, reference_neighbor_idx)?;
    let lone_pair_direction = infer_lone_pair_direction(molecule, donor_idx, &neighbor_indices)?;
    Ok(donor.position + config.center_distance * lone_pair_direction)
}

#[derive(Clone, Copy)]
struct Basis {
    x: Vec3,
    y: Vec3,
    z: Vec3,
}

#[derive(Clone, Copy)]
struct AlignedAtom {
    position: Vec3,
    radius_squared: f32,
}

struct OccupiedVolumes {
    buried_volume: f32,
    quadrants: [f32; 4],
    octants: [f32; 8],
    near_vbur: f32,
    far_vbur: f32,
}

fn validate_config(config: BuriedVolumeConfig) -> Result<(), BuriedVolumeError> {
    if !config.sphere_radius.is_finite() || config.sphere_radius <= 0.0 {
        return Err(BuriedVolumeError(
            "sphere radius must be positive and finite".to_owned(),
        ));
    }
    if !config.density.is_finite() || config.density <= 0.0 {
        return Err(BuriedVolumeError(
            "density must be positive and finite".to_owned(),
        ));
    }
    if !config.center_distance.is_finite() || config.center_distance <= 0.0 {
        return Err(BuriedVolumeError(
            "centre distance must be positive and finite".to_owned(),
        ));
    }
    if !config.radii_scale.is_finite() || config.radii_scale <= 0.0 {
        return Err(BuriedVolumeError(
            "radii scale must be positive and finite".to_owned(),
        ));
    }
    Ok(())
}

/// Tolerance applied to summed covalent radii when inferring bonds from
/// Cartesian coordinates. Real P–X bonds sit near 1.0× the summed radii while
/// the nearest non-bonded contact is ~1.5×, so 1.3 separates them with margin.
const BOND_TOLERANCE_FACTOR: f32 = 1.3;

/// Covalently bonded neighbours of `donor_idx` as `(distance_squared, index)`,
/// sorted by ascending distance then index.
///
/// The bond cutoff is [`BOND_TOLERANCE_FACTOR`] × the summed covalent radii.
/// `donor_idx` must be a valid atom index. This is the single source of the
/// covalent-radius frame used by both the buried-volume quadrant frame and the
/// `descriptors` command's donor detection, so the two never drift apart.
pub fn bonded_neighbors(molecule: &Molecule, donor_idx: usize) -> Vec<(f32, usize)> {
    let donor = &molecule.atoms[donor_idx];
    let donor_covalent = covalent_radius(&donor.element);
    let mut bonded = molecule
        .atoms
        .iter()
        .enumerate()
        .filter(|(index, _)| *index != donor_idx)
        .filter_map(|(index, atom)| {
            let distance_squared = atom.position.distance_squared(donor.position);
            let bond_cutoff =
                BOND_TOLERANCE_FACTOR * (donor_covalent + covalent_radius(&atom.element));
            (distance_squared <= bond_cutoff * bond_cutoff).then_some((distance_squared, index))
        })
        .collect::<Vec<_>>();
    bonded.sort_by(|left, right| {
        left.0
            .total_cmp(&right.0)
            .then_with(|| left.1.cmp(&right.1))
    });
    bonded
}

/// Identify the three atoms covalently bonded to a trivalent donor.
///
/// Substituents are selected by covalent-radius bond detection, **including any
/// bonded hydrogens**, rather than by taking the three nearest heavy atoms.
/// The nearest-heavy heuristic silently misidentifies the frame of primary and
/// secondary phosphines (R–PH2, R2P–H): it discards the bonded hydrogens and
/// reaches for distant non-bonded heavy atoms, which places the lone-pair
/// centre in empty space (spurious `max_delta_qvbur = 0`) or skews it (gross
/// overestimates). For a genuinely trisubstituted donor the two rules coincide,
/// because no non-bonded heavy atom can lie closer than a real P–X bond.
///
/// Hydrogens participate only in defining the geometric frame here; whether
/// they contribute occupied volume remains governed by [`BuriedVolumeConfig`].
fn donor_neighbor_indices(
    molecule: &Molecule,
    donor_idx: usize,
    reference_neighbor_idx: usize,
) -> Result<[usize; 3], BuriedVolumeError> {
    if reference_neighbor_idx == donor_idx || reference_neighbor_idx >= molecule.atoms.len() {
        return Err(BuriedVolumeError(
            "reference neighbor index is invalid".to_owned(),
        ));
    }
    // `bonded_neighbors` returns the covalent-radius frame already sorted by
    // ascending distance then index.
    let bonded = bonded_neighbors(molecule, donor_idx);
    if bonded.len() != 3 {
        return Err(BuriedVolumeError(format!(
            "donor must be trivalent for the quadrant frame, found {} bonded substituents",
            bonded.len()
        )));
    }
    if !bonded
        .iter()
        .any(|(_, index)| *index == reference_neighbor_idx)
    {
        return Err(BuriedVolumeError(
            "reference neighbor is not bonded to the donor".to_owned(),
        ));
    }
    // Deterministic order: reference substituent first, then the remaining two
    // by ascending distance and index. Order does not affect the min/max
    // quadrant descriptors, but keeps the audit output reproducible.
    let mut neighbors = [reference_neighbor_idx; 3];
    let mut cursor = 1;
    for (_, index) in bonded {
        if index != reference_neighbor_idx {
            neighbors[cursor] = index;
            cursor += 1;
        }
    }
    Ok(neighbors)
}

fn infer_lone_pair_direction(
    molecule: &Molecule,
    donor_idx: usize,
    neighbors: &[usize; 3],
) -> Result<Vec3, BuriedVolumeError> {
    let donor = molecule.atoms[donor_idx].position;
    let vectors = neighbors.map(|index| (molecule.atoms[index].position - donor).normalize());
    if vectors.iter().any(|vector| !vector.is_finite()) {
        return Err(BuriedVolumeError(
            "donor and substituent coordinates are coincident".to_owned(),
        ));
    }
    let opposing_sum = -(vectors[0] + vectors[1] + vectors[2]);
    if opposing_sum.length_squared() > 1.0e-4 {
        return Ok(opposing_sum.normalize());
    }

    let normal = vectors[0].cross(vectors[1]);
    if normal.length_squared() <= 1.0e-6 {
        return Err(BuriedVolumeError(
            "donor substituents do not define a lone-pair direction".to_owned(),
        ));
    }
    let candidate = normal.normalize();
    let positive_clearance = center_clearance(molecule, donor_idx, donor + candidate);
    let negative_clearance = center_clearance(molecule, donor_idx, donor - candidate);
    Ok(if positive_clearance >= negative_clearance {
        candidate
    } else {
        -candidate
    })
}

fn center_clearance(molecule: &Molecule, donor_idx: usize, point: Vec3) -> f32 {
    molecule
        .atoms
        .iter()
        .enumerate()
        .filter(|(index, atom)| *index != donor_idx && !atom.element.eq_ignore_ascii_case("H"))
        .map(|(_, atom)| point.distance_squared(atom.position))
        .fold(f32::INFINITY, f32::min)
}

fn coordinate_basis(
    molecule: &Molecule,
    donor_idx: usize,
    plane_idx: usize,
    center: Vec3,
) -> Result<Basis, BuriedVolumeError> {
    // Morfeus maps the centre->donor vector onto negative Z.
    let z = -(molecule.atoms[donor_idx].position - center).normalize();
    let plane_vector = molecule.atoms[plane_idx].position - center;
    let x_projection = plane_vector - z * plane_vector.dot(z);
    if !z.is_finite() || x_projection.length_squared() <= 1.0e-8 {
        return Err(BuriedVolumeError(
            "donor plane atom is collinear with the coordination axis".to_owned(),
        ));
    }
    let x = x_projection.normalize();
    let y = z.cross(x).normalize();
    Ok(Basis { x, y, z })
}

fn aligned_atoms(
    molecule: &Molecule,
    center: Vec3,
    basis: Basis,
    config: BuriedVolumeConfig,
) -> Result<Vec<AlignedAtom>, BuriedVolumeError> {
    molecule
        .atoms
        .iter()
        .filter(|atom| config.include_hydrogens || !atom.element.eq_ignore_ascii_case("H"))
        .map(|atom| {
            if !atom.position.is_finite() || !atom.vdw_radius.is_finite() || atom.vdw_radius <= 0.0
            {
                return Err(BuriedVolumeError(
                    "atom has invalid coordinates or radius".to_owned(),
                ));
            }
            let relative = atom.position - center;
            let position = Vec3::new(
                relative.dot(basis.x),
                relative.dot(basis.y),
                relative.dot(basis.z),
            );
            let radius = atom.vdw_radius * config.radii_scale;
            Ok(AlignedAtom {
                position,
                radius_squared: radius * radius,
            })
        })
        .collect()
}

fn integration_grid(config: BuriedVolumeConfig) -> Vec<Vec3> {
    let volume = sphere_volume(config.sphere_radius);
    let side_points = ((volume / config.density * 6.0 / PI).cbrt().round() as usize).max(2);
    let step = 2.0 * config.sphere_radius / (side_points - 1) as f32;
    let radius_squared = config.sphere_radius * config.sphere_radius;
    let mut points = Vec::with_capacity(side_points.pow(3));
    for ix in 0..side_points {
        let x = -config.sphere_radius + ix as f32 * step;
        for iy in 0..side_points {
            let y = -config.sphere_radius + iy as f32 * step;
            for iz in 0..side_points {
                let z = -config.sphere_radius + iz as f32 * step;
                let point = Vec3::new(x, y, z);
                if point.length_squared() <= radius_squared {
                    points.push(point);
                }
            }
        }
    }
    points
}

fn occupied_volumes(sphere: &[Vec3], atoms: &[AlignedAtom], sphere_radius: f32) -> OccupiedVolumes {
    let mut occupied_total = 0_usize;
    let mut quadrant_total = [0_usize; 4];
    let mut quadrant_occupied = [0_usize; 4];
    let mut octant_total = [0_usize; 8];
    let mut octant_occupied = [0_usize; 8];
    for point in sphere {
        let quadrant = quadrant_index(*point);
        let octant = octant_index(*point);
        quadrant_total[quadrant] += 1;
        octant_total[octant] += 1;
        let occupied = atoms
            .iter()
            .any(|atom| point.distance_squared(atom.position) <= atom.radius_squared);
        if occupied {
            occupied_total += 1;
            quadrant_occupied[quadrant] += 1;
            octant_occupied[octant] += 1;
        }
    }
    let volume = sphere_volume(sphere_radius);
    let quadrants = std::array::from_fn(|index| {
        occupied_fraction(quadrant_occupied[index], quadrant_total[index]) * volume / 4.0
    });
    let octants = std::array::from_fn(|index| {
        occupied_fraction(octant_occupied[index], octant_total[index]) * volume / 8.0
    });
    let near_vbur = octants[4..].iter().sum();
    let far_vbur = octants[..4].iter().sum();
    OccupiedVolumes {
        buried_volume: occupied_fraction(occupied_total, sphere.len()) * volume,
        quadrants,
        octants,
        near_vbur,
        far_vbur,
    }
}

fn quadrant_index(point: Vec3) -> usize {
    match (point.x >= 0.0, point.y >= 0.0) {
        (true, true) => 0,
        (false, true) => 1,
        (false, false) => 2,
        (true, false) => 3,
    }
}

fn octant_index(point: Vec3) -> usize {
    let quadrant = quadrant_index(point);
    if point.z >= 0.0 {
        quadrant
    } else {
        4 + quadrant
    }
}

fn occupied_fraction(occupied: usize, total: usize) -> f32 {
    if total == 0 {
        0.0
    } else {
        occupied as f32 / total as f32
    }
}

fn sphere_volume(radius: f32) -> f32 {
    4.0 * PI * radius.powi(3) / 3.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Atom;

    fn atom(element: &str, position: Vec3) -> Atom {
        Atom::new(element, position)
    }

    fn phosphine() -> Molecule {
        Molecule {
            atoms: vec![
                atom("P", Vec3::ZERO),
                atom("C", Vec3::new(1.4, 0.0, 0.45)),
                atom("C", Vec3::new(-0.7, 1.212, 0.45)),
                atom("C", Vec3::new(-0.7, -1.212, 0.45)),
                atom("C", Vec3::new(2.8, 0.0, 0.7)),
            ],
        }
    }

    #[test]
    fn descriptor_is_finite_and_volume_is_partitioned() {
        let params =
            BuriedVolumeCalculator::compute(&phosphine(), 0, 1, BuriedVolumeConfig::default())
                .unwrap();
        assert!(params.buried_volume.is_finite() && params.buried_volume > 0.0);
        assert!(params.percent_buried_volume > 0.0 && params.percent_buried_volume < 100.0);
        assert!((params.near_vbur + params.far_vbur - params.buried_volume).abs() < 0.5);
        assert!(params.qvbur_min <= params.qvbur_max);
        assert!(params.ovbur_min <= params.ovbur_max);
    }

    #[test]
    fn descriptor_is_rotation_and_translation_invariant() {
        let source = phosphine();
        let rotation = glam::Quat::from_rotation_y(0.73) * glam::Quat::from_rotation_z(-1.21);
        let translation = Vec3::new(5.0, -3.0, 1.5);
        let transformed = Molecule {
            atoms: source
                .atoms
                .iter()
                .map(|source_atom| {
                    atom(
                        source_atom.element.as_str(),
                        rotation * source_atom.position + translation,
                    )
                })
                .collect(),
        };
        let first =
            BuriedVolumeCalculator::compute(&source, 0, 1, BuriedVolumeConfig::default()).unwrap();
        let second =
            BuriedVolumeCalculator::compute(&transformed, 0, 1, BuriedVolumeConfig::default())
                .unwrap();
        assert!((first.buried_volume - second.buried_volume).abs() < 1.0e-4);
        assert!((first.max_delta_qvbur - second.max_delta_qvbur).abs() < 1.0e-4);
    }

    #[test]
    fn accepts_an_explicit_quantum_coordination_center() {
        let molecule = phosphine();
        let inferred =
            BuriedVolumeCalculator::compute(&molecule, 0, 1, BuriedVolumeConfig::default())
                .unwrap();
        let explicit = BuriedVolumeCalculator::compute_with_center(
            &molecule,
            0,
            1,
            Vec3::new(0.0, 0.0, -2.1),
            BuriedVolumeConfig::default(),
        )
        .unwrap();
        assert!((inferred.buried_volume - explicit.buried_volume).abs() < 1.0e-4);
        assert!((inferred.max_delta_qvbur - explicit.max_delta_qvbur).abs() < 1.0e-4);
    }

    #[test]
    fn ensemble_aggregation_tracks_vbur_minimum_conformer() {
        let conformers = [
            BuriedVolumeParams {
                buried_volume: 20.0,
                max_delta_qvbur: 4.0,
                ..BuriedVolumeParams::default()
            },
            BuriedVolumeParams {
                buried_volume: 15.0,
                max_delta_qvbur: 7.0,
                ..BuriedVolumeParams::default()
            },
        ];
        let ensemble = BuriedVolumeCalculator::aggregate(&conformers, &[0.25, 0.75]).unwrap();
        assert_eq!(ensemble.vbur_boltz, 16.25);
        assert_eq!(ensemble.vbur_delta, 5.0);
        assert_eq!(ensemble.max_delta_qvbur_vburminconf, 7.0);
        assert_eq!(ensemble.conformer_count, 2);
    }

    fn primary_phosphine() -> Molecule {
        // R–PH2: one heavy (C) substituent plus two bonded hydrogens. The two
        // nearest *heavy* atoms after the bonded carbon are non-bonded ring
        // carbons, so the legacy nearest-heavy rule mis-framed this donor.
        Molecule {
            atoms: vec![
                atom("P", Vec3::ZERO),
                atom("C", Vec3::new(1.5, 0.0, 0.6)),
                atom("H", Vec3::new(-0.6, 1.1, 0.55)),
                atom("H", Vec3::new(-0.6, -1.1, 0.55)),
                // Distal ring carbons, closer than nothing but not bonded to P.
                atom("C", Vec3::new(2.6, 0.9, 0.9)),
                atom("C", Vec3::new(2.6, -0.9, 0.9)),
                atom("C", Vec3::new(3.9, 0.0, 1.2)),
            ],
        }
    }

    #[test]
    fn primary_phosphine_uses_bonded_hydrogens_not_distal_heavy_atoms() {
        // The bonded set is {C, H, H}; the two distal ring carbons must be
        // excluded. A correct frame yields a finite, strictly positive
        // asymmetry instead of the spurious zero the nearest-heavy rule gave.
        let neighbors = donor_neighbor_indices(&primary_phosphine(), 0, 1).unwrap();
        let mut sorted = neighbors;
        sorted.sort_unstable();
        assert_eq!(sorted, [1, 2, 3]);
        let params = BuriedVolumeCalculator::compute(
            &primary_phosphine(),
            0,
            1,
            BuriedVolumeConfig::default(),
        )
        .unwrap();
        assert!(params.max_delta_qvbur > 0.0);
        assert!(params.buried_volume.is_finite() && params.buried_volume > 0.0);
    }

    #[test]
    fn non_trivalent_donor_is_rejected() {
        // Only the bonded carbon and one hydrogen: a two-coordinate donor that
        // cannot define the three-orientation frame.
        let molecule = Molecule {
            atoms: vec![
                atom("P", Vec3::ZERO),
                atom("C", Vec3::new(1.5, 0.0, 0.6)),
                atom("H", Vec3::new(-0.7, 1.0, 0.5)),
                atom("C", Vec3::new(3.9, 0.0, 1.2)),
            ],
        };
        let error = donor_neighbor_indices(&molecule, 0, 1).unwrap_err();
        assert!(error.to_string().contains("trivalent"));
    }

    #[test]
    fn rejects_invalid_inputs() {
        let error =
            BuriedVolumeCalculator::compute(&phosphine(), 8, 1, BuriedVolumeConfig::default())
                .unwrap_err();
        assert!(error.to_string().contains("out of bounds"));
    }
}
