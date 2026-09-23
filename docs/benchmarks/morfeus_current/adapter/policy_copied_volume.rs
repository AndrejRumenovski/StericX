use super::{Molecule, covalent_radius};
use glam::Vec3;
use std::error::Error;
use std::f32::consts::PI;
use std::fmt::{Display, Formatter};

/// Geometric settings for a metal-centred buried-volume scan.
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
    /// Mean ligand-occupied volume over the three donor-plane grids in Å³.
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
    /// Mean occupied volume in the donor-facing hemisphere over three planes.
    pub near_vbur: f32,
    /// Mean occupied volume in the distal hemisphere over three planes.
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

/// Deterministic f32 voxel descriptors with an explicit geometric-center convention.
#[derive(Clone, Copy, Debug, Default)]
pub struct BuriedVolumeCalculator;

impl BuriedVolumeCalculator {
    /// Calculate one conformer's coordination-aware buried-volume descriptors.
    ///
    /// The virtual centre is placed 2.1 Å from the donor along a geometrically
    /// inferred lone-pair direction (negative sum of unit bond directions).
    /// This differs from both the raw-bond-vector DFT construction and an
    /// electronic localized-orbital center. The three covalently
    /// bonded donor substituents (see [`bonded_neighbors`], hydrogens
    /// included) are each used to define the XZ plane, matching Kraken's
    /// three-orientation quadrant scan.
    pub fn compute(
        molecule: &Molecule,
        donor_idx: usize,
        reference_neighbor_idx: usize,
        config: BuriedVolumeConfig,
    ) -> Result<BuriedVolumeParams, BuriedVolumeError> {
        crate::profile_scope!("buried_volume", "BuriedVolumeCalculator::compute");
        let (neighbor_indices, center) =
            donor_geometry(molecule, donor_idx, reference_neighbor_idx, config)?;
        Self::compute_from_center(molecule, donor_idx, neighbor_indices, center, config)
    }

    /// Check the same frame, input, and grid preconditions without calculating occupancy.
    pub fn validate_for_sterimol(
        molecule: &Molecule,
        donor_idx: usize,
        reference_neighbor_idx: usize,
        config: BuriedVolumeConfig,
    ) -> Result<(), BuriedVolumeError> {
        let (neighbors, center) =
            donor_geometry(molecule, donor_idx, reference_neighbor_idx, config)?;
        Self::validate_with_neighbors(molecule, donor_idx, neighbors, center, config)
    }

    /// Validate an explicitly supplied trivalent topology and coordination center.
    pub fn validate_with_neighbors(
        molecule: &Molecule,
        donor_idx: usize,
        neighbors: [usize; 3],
        center: Vec3,
        config: BuriedVolumeConfig,
    ) -> Result<(), BuriedVolumeError> {
        validate_config(config)?;
        validate_neighbors(molecule, donor_idx, neighbors)?;
        if integration_grid(config).is_empty() {
            return Err(BuriedVolumeError(
                "integration grid contains no points".into(),
            ));
        }
        for plane in neighbors {
            let basis = coordinate_basis(molecule, donor_idx, plane, center)?;
            aligned_atoms(molecule, center, basis, config)?;
        }
        Ok(())
    }

    /// Calculate using three explicit bonded neighbors, without distance inference.
    /// Supplied topology is authoritative; neighbors must be distinct and valid.
    pub fn compute_with_neighbors(
        molecule: &Molecule,
        donor_idx: usize,
        neighbors: [usize; 3],
        config: BuriedVolumeConfig,
    ) -> Result<BuriedVolumeParams, BuriedVolumeError> {
        let center = coordination_center_with_neighbors(molecule, donor_idx, neighbors, config)?;
        Self::compute_from_center(molecule, donor_idx, neighbors, center, config)
    }

    /// Calculate using explicit bonded neighbors and an explicit coordination center.
    /// This also resolves an otherwise ambiguous planar-center sign.
    pub fn compute_with_center_and_neighbors(
        molecule: &Molecule,
        donor_idx: usize,
        neighbors: [usize; 3],
        center: Vec3,
        config: BuriedVolumeConfig,
    ) -> Result<BuriedVolumeParams, BuriedVolumeError> {
        validate_config(config)?;
        validate_neighbors(molecule, donor_idx, neighbors)?;
        Self::compute_from_center(molecule, donor_idx, neighbors, center, config)
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
        crate::profile_scope!(
            "buried_volume",
            "BuriedVolumeCalculator::compute_with_center"
        );
        validate_config(config)?;
        let donor = molecule.atoms.get(donor_idx).ok_or_else(|| {
            BuriedVolumeError(format!("donor index {donor_idx} is out of bounds"))
        })?;
        if !donor.position.is_finite() || !center.is_finite() {
            return Err(BuriedVolumeError(
                "donor or coordination-center coordinate is not finite".to_owned(),
            ));
        }
        if center == donor.position {
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
        crate::profile_scope!(
            "buried_volume",
            "BuriedVolumeCalculator::compute_from_center"
        );
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
        let mut totals = [0.0_f32; 3];
        let mut near = [0.0_f32; 3];
        let mut far = [0.0_f32; 3];
        for (orientation_index, &plane_idx) in neighbor_indices.iter().enumerate() {
            let basis = coordinate_basis(molecule, donor_idx, plane_idx, center)?;
            let atoms = aligned_atoms(molecule, center, basis, config)?;
            let volumes = occupied_volumes(&sphere, &atoms, config.sphere_radius);
            totals[orientation_index] = volumes.buried_volume;
            near[orientation_index] = volumes.near_vbur;
            far[orientation_index] = volumes.far_vbur;
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
        // The previous first-plane selection changed under atom permutations.
        // Average the same three lattice estimates symmetrically. Sorting before
        // reduction fixes f32 accumulation order without choosing a chemical axis.
        fn mean(mut values: [f32; 3]) -> f32 {
            values.sort_by(f32::total_cmp);
            values.into_iter().map(|value| value / 3.0).sum()
        }
        result.buried_volume = mean(totals);
        result.percent_buried_volume =
            100.0 * (result.buried_volume / sphere_volume(config.sphere_radius));
        result.near_vbur = mean(near);
        result.far_vbur = mean(far);
        Ok(result)
    }

    /// Aggregate conformer descriptors using normalized non-negative weights.
    pub fn aggregate(
        conformers: &[BuriedVolumeParams],
        weights: &[f32],
    ) -> Result<BuriedVolumeEnsembleParams, BuriedVolumeError> {
        crate::profile_scope!("conformer_processing", "BuriedVolumeCalculator::aggregate");
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
        if conformers.iter().any(|params| {
            [
                params.buried_volume,
                params.percent_buried_volume,
                params.qvbur_min,
                params.qvbur_max,
                params.max_delta_qvbur,
                params.ovbur_min,
                params.ovbur_max,
                params.near_vbur,
                params.far_vbur,
            ]
            .iter()
            .any(|value| !value.is_finite())
        }) {
            return Err(BuriedVolumeError(
                "conformer descriptors must be finite".to_owned(),
            ));
        }
        let scale = weights.iter().copied().fold(0.0_f32, f32::max);
        if scale == 0.0 {
            return Err(BuriedVolumeError(
                "conformer weights have zero total".to_owned(),
            ));
        }
        // Normalize before the sum can overflow, without an absolute epsilon
        // cutoff that would reject a valid common scaling of the weights.
        let scaled = weights
            .iter()
            .map(|weight| f64::from(*weight) / f64::from(scale))
            .collect::<Vec<_>>();
        let weight_sum = scaled.iter().sum::<f64>();
        let mean = |select: fn(&BuriedVolumeParams) -> f32| {
            (conformers
                .iter()
                .zip(&scaled)
                .map(|(params, weight)| f64::from(select(params)) * (weight / weight_sum))
                .sum::<f64>()) as f32
        };
        let mut ensemble = BuriedVolumeEnsembleParams {
            vbur_boltz: mean(|params| params.buried_volume),
            qvbur_min_boltz: mean(|params| params.qvbur_min),
            qvbur_max_boltz: mean(|params| params.qvbur_max),
            max_delta_qvbur_boltz: mean(|params| params.max_delta_qvbur),
            near_vbur_boltz: mean(|params| params.near_vbur),
            far_vbur_boltz: mean(|params| params.far_vbur),
            vbur_min: conformers[0].buried_volume,
            vbur_max: conformers[0].buried_volume,
            max_delta_qvbur_min: conformers[0].max_delta_qvbur,
            max_delta_qvbur_max: conformers[0].max_delta_qvbur,
            conformer_count: conformers.len(),
            ..BuriedVolumeEnsembleParams::default()
        };
        let mut minimum_vbur_index = 0;
        for (index, params) in conformers.iter().enumerate() {
            if params.buried_volume < conformers[minimum_vbur_index].buried_volume {
                minimum_vbur_index = index;
            }
            ensemble.vbur_min = ensemble.vbur_min.min(params.buried_volume);
            ensemble.vbur_max = ensemble.vbur_max.max(params.buried_volume);
            ensemble.max_delta_qvbur_min = ensemble.max_delta_qvbur_min.min(params.max_delta_qvbur);
            ensemble.max_delta_qvbur_max = ensemble.max_delta_qvbur_max.max(params.max_delta_qvbur);
        }
        ensemble.vbur_delta = ensemble.vbur_max - ensemble.vbur_min;
        ensemble.max_delta_qvbur_delta =
            ensemble.max_delta_qvbur_max - ensemble.max_delta_qvbur_min;
        ensemble.max_delta_qvbur_vburminconf = conformers[minimum_vbur_index].max_delta_qvbur;
        if [
            ensemble.vbur_boltz,
            ensemble.qvbur_min_boltz,
            ensemble.qvbur_max_boltz,
            ensemble.max_delta_qvbur_boltz,
            ensemble.near_vbur_boltz,
            ensemble.far_vbur_boltz,
            ensemble.vbur_delta,
            ensemble.max_delta_qvbur_delta,
        ]
        .iter()
        .any(|value| !value.is_finite())
        {
            return Err(BuriedVolumeError(
                "ensemble aggregation exceeds the finite f32 descriptor range".to_owned(),
            ));
        }
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
    crate::profile_scope!("donor_bond_detection", "coordination_center");
    donor_geometry(molecule, donor_idx, reference_neighbor_idx, config).map(|(_, center)| center)
}

/// Geometric coordination center using authoritative trivalent connectivity.
/// The negative sum of normalized bonds defines the direction. A planar
/// clearance tie has no unique sign and requires an explicit center instead.
pub fn coordination_center_with_neighbors(
    molecule: &Molecule,
    donor_idx: usize,
    neighbors: [usize; 3],
    config: BuriedVolumeConfig,
) -> Result<Vec3, BuriedVolumeError> {
    validate_config(config)?;
    validate_neighbors(molecule, donor_idx, neighbors)?;
    let donor = molecule.atoms[donor_idx].position;
    let direction = infer_lone_pair_direction(molecule, donor_idx, &neighbors)?;
    let center = donor + config.center_distance * direction;
    if !center.is_finite() || center == donor {
        return Err(BuriedVolumeError(
            "coordination center is non-finite or indistinguishable from donor at f32 precision"
                .into(),
        ));
    }
    Ok(center)
}

fn validate_neighbors(
    molecule: &Molecule,
    donor_idx: usize,
    neighbors: [usize; 3],
) -> Result<(), BuriedVolumeError> {
    let donor = molecule
        .atoms
        .get(donor_idx)
        .ok_or_else(|| BuriedVolumeError("donor index is out of bounds".into()))?;
    if !donor.position.is_finite() {
        return Err(BuriedVolumeError("donor coordinate is not finite".into()));
    }
    if neighbors[0] == neighbors[1] || neighbors[0] == neighbors[2] || neighbors[1] == neighbors[2]
    {
        return Err(BuriedVolumeError("donor neighbors must be distinct".into()));
    }
    for index in neighbors {
        let atom = molecule
            .atoms
            .get(index)
            .ok_or_else(|| BuriedVolumeError("neighbor index is out of bounds".into()))?;
        if index == donor_idx || !atom.position.is_finite() || atom.position == donor.position {
            return Err(BuriedVolumeError(
                "invalid or coincident donor neighbor".into(),
            ));
        }
    }
    Ok(())
}

/// Shared preconditions in the same order as the full buried-volume calculation.
fn donor_geometry(
    molecule: &Molecule,
    donor_idx: usize,
    reference_neighbor_idx: usize,
    config: BuriedVolumeConfig,
) -> Result<([usize; 3], Vec3), BuriedVolumeError> {
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
    let center = coordination_center_with_neighbors(molecule, donor_idx, neighbor_indices, config)?;
    Ok((neighbor_indices, center))
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
    let volume = sphere_volume(config.sphere_radius);
    let side = (volume / config.density * 6.0 / PI).cbrt().round().max(2.0);
    let max_points = (isize::MAX as usize) / std::mem::size_of::<Vec3>();
    if !volume.is_finite()
        || volume <= 0.0
        || !side.is_finite()
        || side > (max_points as f64).cbrt() as f32
        || !(config.sphere_radius * config.sphere_radius).is_finite()
    {
        return Err(BuriedVolumeError(
            "integration grid exceeds finite representable range".into(),
        ));
    }

    Ok(())
}

/// Tolerance applied to summed covalent radii when inferring bonds from
/// Cartesian coordinates. This heuristic can include non-bonded close contacts;
/// supplied connectivity should be used when available.
const BOND_TOLERANCE_FACTOR: f32 = 1.3;

/// Covalently bonded neighbours of `donor_idx` as `(distance_squared, index)`,
/// sorted by ascending distance then index.
///
/// The bond cutoff is `BOND_TOLERANCE_FACTOR` × the summed covalent radii.
/// `donor_idx` must be a valid atom index. This is the single source of the
/// covalent-radius frame used by both the buried-volume quadrant frame and the
/// `descriptors` command's donor detection, so the two never drift apart.
pub fn bonded_neighbors(molecule: &Molecule, donor_idx: usize) -> Vec<(f32, usize)> {
    crate::profile_scope!("donor_bond_detection", "bonded_neighbors");
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
/// overestimates). Distance inference itself can still include close non-bonded
/// contacts; explicit-neighbor APIs avoid that ambiguity.
///
/// Hydrogens participate only in defining the geometric frame here; whether
/// they contribute occupied volume remains governed by [`BuriedVolumeConfig`].
fn donor_neighbor_indices(
    molecule: &Molecule,
    donor_idx: usize,
    reference_neighbor_idx: usize,
) -> Result<[usize; 3], BuriedVolumeError> {
    crate::profile_scope!("donor_bond_detection", "donor_neighbor_indices");
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
    crate::profile_scope!("donor_bond_detection", "infer_lone_pair_direction");
    let donor = molecule.atoms[donor_idx].position;
    let mut vectors = neighbors.map(|index| {
        let displacement = molecule.atoms[index].position - donor;
        let squared = displacement.length_squared();
        if squared.is_finite() && squared > f32::EPSILON {
            displacement.normalize()
        } else {
            // Only the demonstrated small/overflow-prone normalization path
            // needs wider intermediates; ordinary donor arithmetic stays f32.
            (molecule.atoms[index].position.as_dvec3() - donor.as_dvec3())
                .normalize()
                .as_vec3()
        }
    });
    // Reproducible sum and planar reference under neighbor/atom permutations.
    vectors.sort_by(|a, b| {
        a.x.total_cmp(&b.x)
            .then(a.y.total_cmp(&b.y))
            .then(a.z.total_cmp(&b.z))
    });
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
    // In the planar fallback, an exact f32 equality comparison can select a
    // sign solely because a rigid transform rounded the input coordinates or
    // because adding a unit probe to a translated donor lost precision.
    // Resolve the fixed represented trial normal only when coordinate-rounding
    // clearance intervals are disjoint. This does not propagate uncertainty in
    // the inferred normal itself or establish an experimental coordinate error.
    let positive = center_clearance_interval(molecule, donor_idx, candidate)?;
    let negative = center_clearance_interval(molecule, donor_idx, -candidate)?;
    if positive.0 > negative.1 {
        Ok(candidate)
    } else if negative.0 > positive.1 {
        Ok(-candidate)
    } else {
        Err(BuriedVolumeError(
            "planar coordination direction is ambiguous at input coordinate precision; supply an explicit center".into(),
        ))
    }
}

/// Half the larger adjacent f32 spacing: a conservative symmetric rounding cell.
fn coordinate_halfwidth(value: f32) -> f64 {
    let magnitude = value.abs();
    let bits = magnitude.to_bits();
    let gap = if bits == f32::MAX.to_bits() {
        f64::from(magnitude) - f64::from(f32::from_bits(bits - 1))
    } else {
        f64::from(f32::from_bits(bits + 1)) - f64::from(magnitude)
    };
    gap * 0.5
}

// Explicit adjacent-value operations preserve the crate's Rust 1.85 support.
fn outward_down(value: f64) -> f64 {
    if value == 0.0 {
        return -f64::from_bits(1);
    }
    f64::from_bits(if value > 0.0 {
        value.to_bits() - 1
    } else {
        value.to_bits() + 1
    })
}

fn outward_up(value: f64) -> f64 {
    if value == 0.0 {
        return f64::from_bits(1);
    }
    f64::from_bits(if value > 0.0 {
        value.to_bits() + 1
    } else {
        value.to_bits() - 1
    })
}

/// Squared nearest-heavy-atom distance intervals for a fixed one-Å trial normal.
/// Only coordinate rounding is propagated; the represented probe stays fixed.
fn center_clearance_interval(
    molecule: &Molecule,
    donor_idx: usize,
    probe: Vec3,
) -> Result<(f64, f64), BuriedVolumeError> {
    crate::profile_scope!("donor_bond_detection", "center_clearance");
    let donor = molecule.atoms[donor_idx].position.to_array();
    let mut nearest = (f64::INFINITY, f64::INFINITY);
    for (index, atom) in molecule.atoms.iter().enumerate() {
        if index == donor_idx || atom.element.eq_ignore_ascii_case("H") {
            continue;
        }
        if !atom.position.is_finite() {
            return Err(BuriedVolumeError("atom coordinate is not finite".into()));
        }
        let mut distance = (0.0, 0.0);
        for ((a, d), p) in atom
            .position
            .to_array()
            .into_iter()
            .zip(donor)
            .zip(probe.to_array())
        {
            let da = coordinate_halfwidth(a);
            let dd = coordinate_halfwidth(d);
            let a = f64::from(a);
            let d = f64::from(d);
            let p = f64::from(p);
            let low = outward_down(outward_down(outward_down(a - da) - outward_up(d + dd)) - p);
            let high = outward_up(outward_up(outward_up(a + da) - outward_down(d - dd)) - p);
            let closest = if low <= 0.0 && high >= 0.0 {
                0.0
            } else {
                low.abs().min(high.abs())
            };
            let farthest = low.abs().max(high.abs());
            distance.0 =
                outward_down(distance.0 + outward_down(closest * closest).max(0.0)).max(0.0);
            distance.1 = outward_up(distance.1 + outward_up(farthest * farthest));
        }
        nearest.0 = nearest.0.min(distance.0);
        nearest.1 = nearest.1.min(distance.1);
    }
    Ok(nearest)
}

fn coordinate_basis(
    molecule: &Molecule,
    donor_idx: usize,
    plane_idx: usize,
    center: Vec3,
) -> Result<Basis, BuriedVolumeError> {
    crate::profile_scope!("buried_volume", "buried_volume::coordinate_basis");
    // Morfeus maps the centre->donor vector onto negative Z.
    if !center.is_finite() || center == molecule.atoms[donor_idx].position {
        return Err(BuriedVolumeError(
            "coordination center is invalid or coincident with donor".into(),
        ));
    }
    let z = -(molecule.atoms[donor_idx].position - center).normalize();
    let plane_vector = molecule.atoms[plane_idx].position - center;
    let x_projection = plane_vector - z * plane_vector.dot(z);
    if !z.is_finite() || !x_projection.is_finite() || x_projection.length_squared() <= 1.0e-8 {
        // Resolve a genuinely nonzero small plane rather than imposing the
        // old absolute threshold. A cross-product construction avoids subtracting
        // nearly equal axial projections in the audited nearly-collinear case.
        let axis = center.as_dvec3() - molecule.atoms[donor_idx].position.as_dvec3();
        let plane = molecule.atoms[plane_idx].position.as_dvec3() - center.as_dvec3();
        let normal = plane.cross(axis);
        if normal.length_squared() == 0.0 || !normal.is_finite() {
            return Err(BuriedVolumeError(
                "donor plane atom is collinear with the coordination axis".into(),
            ));
        }
        let z = axis.normalize();
        let y = normal.normalize();
        let x = z.cross(y);
        return Ok(Basis {
            x: x.as_vec3(),
            y: -y.as_vec3(),
            z: z.as_vec3(),
        });
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
    crate::profile_scope!("buried_volume", "buried_volume::aligned_atoms");
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
            let radius_squared = radius * radius;
            if !position.is_finite() || !radius_squared.is_finite() || radius_squared == 0.0 {
                return Err(BuriedVolumeError(
                    "aligned coordinate or squared radius exceeds finite f32 range".into(),
                ));
            }
            Ok(AlignedAtom {
                position,
                radius_squared,
            })
        })
        .collect()
}

fn integration_grid(config: BuriedVolumeConfig) -> Vec<Vec3> {
    crate::profile_scope!("buried_volume", "buried_volume::integration_grid");
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
    crate::profile_scope!("buried_volume", "buried_volume::occupied_volumes");
    struct RowCandidate {
        z: f32,
        radius_squared: f32,
        xy_squared: f32,
    }

    let mut occupied_total = 0_usize;
    let mut quadrant_total = [0_usize; 4];
    let mut quadrant_occupied = [0_usize; 4];
    let mut octant_total = [0_usize; 8];
    let mut octant_occupied = [0_usize; 8];
    let mut row_xy = None;
    let mut next_atom = 0;
    let mut candidates: Vec<RowCandidate> = Vec::with_capacity(atoms.len().min(64));
    for point in sphere {
        let xy = (point.x.to_bits(), point.y.to_bits());
        if row_xy != Some(xy) {
            row_xy = Some(xy);
            next_atom = 0;
            candidates.clear();
        }
        let quadrant = quadrant_index(*point);
        let octant = octant_index(*point);
        quadrant_total[quadrant] += 1;
        octant_total[octant] += 1;
        // glam::Vec3 uses (dx*dx + dy*dy) + dz*dz. A row reuses the
        // bit-identical first sum without changing that f32 operation order.
        let mut occupied = candidates.iter().any(|candidate| {
            let dz = point.z - candidate.z;
            candidate.xy_squared + dz * dz <= candidate.radius_squared
        });
        // The cache covers retained atoms from [0, next_atom), in original
        // order. Extend it only after all cached candidates miss this point;
        // a dense first-atom hit therefore never prepares the unused suffix.
        while !occupied && next_atom < atoms.len() {
            let atom = &atoms[next_atom];
            next_atom += 1;
            let dx = point.x - atom.position.x;
            let dy = point.y - atom.position.y;
            let xy_squared = dx * dx + dy * dy;
            // Adding a non-negative z square cannot turn this miss into a hit.
            // A NaN comparison is false and retains the candidate; inf <= inf
            // behavior is preserved by the final test.
            if xy_squared > atom.radius_squared {
                continue;
            }
            candidates.push(RowCandidate {
                z: atom.position.z,
                radius_squared: atom.radius_squared,
                xy_squared,
            });
            let dz = point.z - atom.position.z;
            occupied = xy_squared + dz * dz <= atom.radius_squared;
        }
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

    fn assert_sterimol_validation_matches_full(
        molecule: &Molecule,
        donor: usize,
        reference: usize,
        config: BuriedVolumeConfig,
    ) -> Result<(), BuriedVolumeError> {
        let expected =
            BuriedVolumeCalculator::compute(molecule, donor, reference, config).map(|_| ());
        let actual =
            BuriedVolumeCalculator::validate_for_sterimol(molecule, donor, reference, config);
        assert_eq!(actual, expected, "configuration: {config:?}");
        actual
    }

    #[test]
    fn sterimol_validation_accepts_valid_symmetric_and_empty_occupancy() {
        let molecule = phosphine();
        // The donor covers every voxel, so every orientation has exactly equal
        // quadrants. Occupancy symmetry is not a degenerate geometric frame.
        let symmetric = BuriedVolumeConfig {
            sphere_radius: 0.5,
            density: 0.001,
            center_distance: 0.1,
            ..BuriedVolumeConfig::default()
        };
        assert_sterimol_validation_matches_full(&molecule, 0, 1, symmetric).unwrap();
        let result = BuriedVolumeCalculator::compute(&molecule, 0, 1, symmetric).unwrap();
        assert!(result.buried_volume > 0.0);
        assert_eq!(result.max_delta_qvbur, 0.0);

        // Zero occupied volume is accepted even though quadrant differences
        // are also zero. It must not be mistaken for the positive-volume case.
        let empty = BuriedVolumeConfig {
            center_distance: 1000.0,
            ..BuriedVolumeConfig::default()
        };
        let full = BuriedVolumeCalculator::compute(&molecule, 0, 1, empty).unwrap();
        assert_eq!(full.buried_volume, 0.0);
        assert_eq!(full.max_delta_qvbur, 0.0);
        assert_sterimol_validation_matches_full(&molecule, 0, 1, empty).unwrap();
    }

    #[test]
    fn sterimol_validation_checks_later_frames_after_an_occupancy_witness() {
        let molecule = Molecule {
            atoms: vec![
                atom("P", Vec3::ZERO),
                atom("C", Vec3::new(1.5, 0.0, 0.0)),
                atom("C", Vec3::new(-1.5, 0.0, 0.0)),
                atom("C", Vec3::new(0.0, 0.0, 1.5)),
                atom("C", Vec3::new(2.7, 0.8, 0.4)),
            ],
        };
        let config = BuriedVolumeConfig::default();
        let (neighbors, center) = donor_geometry(&molecule, 0, 1, config).unwrap();
        assert_eq!(neighbors, [1, 2, 3]);
        let basis = coordinate_basis(&molecule, 0, neighbors[0], center).unwrap();
        let atoms = aligned_atoms(&molecule, center, basis, config).unwrap();
        let sphere = integration_grid(config);
        let volumes = occupied_volumes(&sphere, &atoms, config.sphere_radius);
        let q = volumes.quadrants;
        let first_delta = (0..4)
            .map(|index| (q[index] - q[(index + 3) % 4]).abs())
            .fold(0.0_f32, f32::max);
        assert!(first_delta > 0.0, "fixture must prove a first-scan witness");
        let error = assert_sterimol_validation_matches_full(&molecule, 0, 1, config).unwrap_err();
        assert!(error.to_string().contains("collinear"));
    }

    #[test]
    fn sterimol_validation_preserves_configuration_and_atom_error_order() {
        let molecule = phosphine();
        let default = BuriedVolumeConfig::default();
        for invalid in [0.0, -1.0, f32::NAN, f32::INFINITY] {
            for config in [
                BuriedVolumeConfig {
                    sphere_radius: invalid,
                    ..default
                },
                BuriedVolumeConfig {
                    density: invalid,
                    ..default
                },
                BuriedVolumeConfig {
                    center_distance: invalid,
                    ..default
                },
                BuriedVolumeConfig {
                    radii_scale: invalid,
                    ..default
                },
            ] {
                assert!(
                    assert_sterimol_validation_matches_full(&molecule, 99, 99, config).is_err()
                );
            }
        }
        for (donor, reference) in [(99, 1), (0, 99), (0, 0), (0, 4)] {
            assert!(
                assert_sterimol_validation_matches_full(&molecule, donor, reference, default)
                    .is_err()
            );
        }
        let mut invalid_donor = molecule.clone();
        invalid_donor.atoms[0].position.x = f32::NAN;
        assert!(
            assert_sterimol_validation_matches_full(&invalid_donor, 0, 1, default)
                .unwrap_err()
                .to_string()
                .contains("donor coordinate")
        );
        let mut invalid_radius = molecule.clone();
        invalid_radius.atoms[4].vdw_radius = 0.0;
        assert!(
            assert_sterimol_validation_matches_full(&invalid_radius, 0, 1, default)
                .unwrap_err()
                .to_string()
                .contains("coordinates or radius")
        );

        let mut invalid_hydrogen = molecule;
        invalid_hydrogen
            .atoms
            .push(atom("H", Vec3::splat(f32::NAN)));
        assert_sterimol_validation_matches_full(&invalid_hydrogen, 0, 1, default).unwrap();
        assert!(
            assert_sterimol_validation_matches_full(
                &invalid_hydrogen,
                0,
                1,
                BuriedVolumeConfig {
                    include_hydrogens: true,
                    ..default
                },
            )
            .unwrap_err()
            .to_string()
            .contains("coordinates or radius")
        );
    }

    #[test]
    fn sterimol_validation_matches_extreme_but_bounded_grid_configurations() {
        let molecule = phosphine();
        let default = BuriedVolumeConfig::default();
        // Keep grid allocations bounded while exercising underflow, overflow
        // in scaled radii/coordinates, an empty grid, and adjacent f32 settings.
        for config in [
            BuriedVolumeConfig {
                sphere_radius: f32::MIN_POSITIVE,
                ..default
            },
            BuriedVolumeConfig {
                sphere_radius: f32::from_bits(1),
                ..default
            },
            BuriedVolumeConfig {
                radii_scale: f32::MAX,
                ..default
            },
            BuriedVolumeConfig {
                center_distance: f32::MAX,
                ..default
            },
            BuriedVolumeConfig {
                density: f32::MAX,
                ..default
            },
            BuriedVolumeConfig {
                radii_scale: f32::from_bits(default.radii_scale.to_bits() - 1),
                ..default
            },
            BuriedVolumeConfig {
                radii_scale: f32::from_bits(default.radii_scale.to_bits() + 1),
                ..default
            },
        ] {
            let _ = assert_sterimol_validation_matches_full(&molecule, 0, 1, config);
        }
    }

    #[test]
    fn sterimol_validation_matches_deterministic_perturbed_geometries() {
        let mut state = 0x4d595df4d0f33173_u64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((state >> 40) as f32 / (1_u32 << 24) as f32 - 0.5) * 0.3
        };
        for index in 0..32 {
            let mut molecule = if index % 2 == 0 {
                phosphine()
            } else {
                primary_phosphine()
            };
            for atom in &mut molecule.atoms {
                atom.position += Vec3::new(next(), next(), next());
            }
            let config = BuriedVolumeConfig {
                density: 0.04,
                center_distance: 2.1 + next(),
                radii_scale: 1.17 + next(),
                include_hydrogens: index % 3 == 0,
                ..BuriedVolumeConfig::default()
            };
            let _ = assert_sterimol_validation_matches_full(&molecule, 0, 1, config);
        }
    }

    // Frozen pre-optimization loop: keep the original glam distance predicate
    // independent of the row cache and its scalar prefix implementation.
    fn occupied_volumes_reference(
        sphere: &[Vec3],
        atoms: &[AlignedAtom],
        sphere_radius: f32,
    ) -> OccupiedVolumes {
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

    fn assert_occupancy_bits_match(sphere: &[Vec3], atoms: &[AlignedAtom], radius: f32) {
        let actual = occupied_volumes(sphere, atoms, radius);
        let expected = occupied_volumes_reference(sphere, atoms, radius);
        let all_bits = |result: OccupiedVolumes| {
            let mut fields = vec![result.buried_volume, result.near_vbur, result.far_vbur];
            fields.extend(result.quadrants);
            fields.extend(result.octants);
            fields.into_iter().map(f32::to_bits).collect::<Vec<_>>()
        };
        assert_eq!(all_bits(actual), all_bits(expected));
    }

    #[test]
    fn row_occupancy_preserves_adjacent_float_boundary_decisions() {
        let points = [
            Vec3::new(0.5, 0.75, -0.25),
            Vec3::new(0.5, 0.75, 0.0),
            Vec3::new(0.5, 0.75, 0.25),
            Vec3::new(0.5, 0.75, f32::from_bits(0.25_f32.to_bits() + 1)),
            Vec3::new(-0.5, -0.75, -0.25),
        ];
        for point in points {
            let distance_squared = point.distance_squared(Vec3::ZERO);
            for radius_squared in [
                f32::from_bits(distance_squared.to_bits() - 1),
                distance_squared,
                f32::from_bits(distance_squared.to_bits() + 1),
            ] {
                let atoms = [AlignedAtom {
                    position: Vec3::ZERO,
                    radius_squared,
                }];
                // Compare the individual predicate, and the same atom across
                // repeated XY rows with boundaries on either side of equality.
                assert_occupancy_bits_match(&[point], &atoms, 3.5);
                assert_occupancy_bits_match(&points, &atoms, 3.5);
            }
        }
    }

    #[test]
    fn row_occupancy_preserves_signed_zero_subnormal_and_nonfinite_arithmetic() {
        let values = [
            0.0,
            -0.0,
            f32::from_bits(1),
            -f32::from_bits(1),
            f32::MIN_POSITIVE,
            -f32::MIN_POSITIVE,
            1.0,
            -1.0,
            1.0e20,
            -1.0e20,
            f32::MAX,
            -f32::MAX,
            f32::INFINITY,
            f32::NEG_INFINITY,
            f32::NAN,
            f32::from_bits(0x7fc00001),
        ];
        let mut points = Vec::new();
        for (index, &x) in values.iter().enumerate() {
            for &z in &values {
                points.push(Vec3::new(x, values[(index + 1) % values.len()], z));
            }
        }
        for (index, &x) in values.iter().enumerate() {
            for radius_squared in [
                0.0,
                -0.0,
                -1.0,
                f32::from_bits(1),
                f32::MIN_POSITIVE,
                1.0,
                f32::MAX,
                f32::INFINITY,
                f32::NAN,
            ] {
                let atom = AlignedAtom {
                    position: Vec3::new(x, values[(index + 2) % values.len()], 0.0),
                    radius_squared,
                };
                assert_occupancy_bits_match(&points, &[atom], 3.5);
            }
        }
        // Cover an underflowed XY prefix followed by a nonzero subnormal sum.
        let atom = AlignedAtom {
            position: Vec3::ZERO,
            radius_squared: 1.0e-40,
        };
        let tiny = [Vec3::new(1.0e-22, -1.0e-22, 1.0e-20), Vec3::ZERO];
        assert_occupancy_bits_match(&tiny, &[atom], 3.5);
    }

    #[test]
    fn row_occupancy_handles_lazy_extension_row_resets_and_dense_first_hits() {
        let atoms = [
            AlignedAtom {
                position: Vec3::new(0.0, 0.0, 1.0),
                radius_squared: 0.01,
            },
            AlignedAtom {
                position: Vec3::new(20.0, 20.0, 0.0),
                radius_squared: 0.01,
            },
            AlignedAtom {
                position: Vec3::new(0.0, 0.0, -1.0),
                radius_squared: 0.01,
            },
            AlignedAtom {
                position: Vec3::new(0.0, 0.0, 2.0),
                radius_squared: 0.01,
            },
        ];
        let points = [
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 0.0, -1.0),
            Vec3::new(0.0, 0.0, 2.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 0.0, 10.0),
            Vec3::new(20.0, 20.0, 0.0),
            Vec3::new(-0.0, 0.0, -1.0),
            Vec3::new(0.0, 0.0, 2.0),
        ];
        assert_occupancy_bits_match(&points, &atoms, 3.5);
        let reversed: Vec<_> = points.into_iter().rev().collect();
        assert_occupancy_bits_match(&reversed, &atoms, 3.5);
        assert_occupancy_bits_match(&points, &[], 3.5);
        assert_occupancy_bits_match(&[], &atoms, 3.5);
        let mut dense = vec![AlignedAtom {
            position: Vec3::ZERO,
            radius_squared: 1000.0,
        }];
        dense.extend((0..2000).map(|index| AlignedAtom {
            position: Vec3::splat(index as f32),
            radius_squared: 1.0,
        }));
        assert_occupancy_bits_match(&points, &dense, 3.5);
    }

    #[test]
    fn row_occupancy_matches_random_atoms_and_non_grid_point_order() {
        let mut state = 0xfedcba9876543210_u64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (state >> 40) as f32 / (1_u32 << 24) as f32
        };
        for case in 0..48 {
            let atoms = (0..case % 31 + 1)
                .map(|_| AlignedAtom {
                    position: Vec3::new(next() * 6.0 - 3.0, next() * 6.0 - 3.0, next() * 6.0 - 3.0),
                    radius_squared: next() * 5.0,
                })
                .collect::<Vec<_>>();
            let mut points = Vec::new();
            for _ in 0..12 {
                let x = next() * 7.0 - 3.5;
                let y = next() * 7.0 - 3.5;
                for _ in 0..16 {
                    points.push(Vec3::new(x, y, next() * 7.0 - 3.5));
                }
            }
            assert_occupancy_bits_match(&points, &atoms, 3.5);
            // Rotate and swap arbitrary rows/points without restoring grid order.
            points.rotate_left(case + 1);
            points.swap(0, 100);
            points.swap(55, 140);
            assert_occupancy_bits_match(&points, &atoms, 3.5);
        }
    }

    #[test]
    fn row_occupancy_matches_all_real_conformer_orientations_bitwise() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("data/conformers");
        let config = BuriedVolumeConfig::default();
        let sphere = integration_grid(config);
        let mut geometries = 0;
        for directory in std::fs::read_dir(root).unwrap() {
            let directory = directory.unwrap().path();
            if !directory.is_dir() {
                continue;
            }
            for entry in std::fs::read_dir(directory).unwrap() {
                let path = entry.unwrap().path();
                if path.extension().and_then(|value| value.to_str()) != Some("xyz") {
                    continue;
                }
                let molecule = super::super::parse_coordinate_file(&path)
                    .unwrap()
                    .remove(0);
                let donor = molecule
                    .atoms
                    .iter()
                    .position(|atom| atom.element == "P")
                    .unwrap();
                let reference = bonded_neighbors(&molecule, donor)
                    .into_iter()
                    .find(|(_, index)| molecule.atoms[*index].element != "H")
                    .unwrap()
                    .1;
                let (neighbors, center) =
                    donor_geometry(&molecule, donor, reference, config).unwrap();
                for plane in neighbors {
                    let basis = coordinate_basis(&molecule, donor, plane, center).unwrap();
                    let atoms = aligned_atoms(&molecule, center, basis, config).unwrap();
                    assert_occupancy_bits_match(&sphere, &atoms, config.sphere_radius);
                }
                geometries += 1;
            }
        }
        assert!(geometries >= 56);
    }
}

#[cfg(test)]
#[path = "buried_volume/aggregation_remediation_tests.rs"]
mod aggregation_remediation_tests;

pub fn policy_probe(m: &Molecule,d:usize,ns:[usize;3],cfg:BuriedVolumeConfig)->serde_json::Value {
    let center=coordination_center_with_neighbors(m,d,ns,cfg).unwrap();
    let grid=integration_grid(cfg);
    let points:Vec<_>=grid.iter().map(|p| {
        let index=p.to_array().map(|x| ((f64::from(x)+3.5)*31.0/7.0).round() as usize);
        serde_json::json!({"index":index,"coordinates":p.to_array(),"region":octant_index(*p)})
    }).collect();
    let frames:Vec<_>=ns.iter().map(|&plane| {
        let basis=coordinate_basis(m,d,plane,center).unwrap();
        let atoms=aligned_atoms(m,center,basis,cfg).unwrap();
        let occupied:Vec<_>=grid.iter().map(|p| atoms.iter().any(|a| p.distance_squared(a.position)<=a.radius_squared)).collect();
        let v=occupied_volumes(&grid,&atoms,cfg.sphere_radius);
        serde_json::json!({"plane":plane,"occupied":occupied,"buried_volume":v.buried_volume,"quadrants":v.quadrants,"octants":v.octants,"near_vbur":v.near_vbur,"far_vbur":v.far_vbur,
            "aligned_atoms":atoms.iter().map(|a|serde_json::json!({"position":a.position.to_array(),"radius_squared":a.radius_squared})).collect::<Vec<_>>()})
    }).collect();
    serde_json::json!({"donor":d,"neighbors":ns,"center":center.to_array(),"points":points,"frames":frames,"sphere_volume":sphere_volume(cfg.sphere_radius)})
}
