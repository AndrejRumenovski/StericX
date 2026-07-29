use super::Molecule;
use glam::Vec3;

/// Pyramidalization descriptors for a trivalent donor atom.
///
/// Reproduces Kraken's published `pyr_P` and `pyr_alpha` (Gensch et al., *J. Am.
/// Chem. Soc.* **2022**, *144*, 1205), which are computed by `morfeus`'
/// `Pyramidalization` class from the three donor→substituent unit vectors.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PyramidalizationParams {
    /// Radhakrishnan pyramidalization `P` (dimensionless).
    ///
    /// Equal to the absolute scalar triple product of the three donor→neighbor
    /// unit vectors: `1` for an ideal tetrahedral-like apex through `0` for a
    /// planar centre, with `morfeus`' `2 - P` correction for an inverted
    /// ("acute") pyramid.
    pub pyr_p: f32,
    /// Mean pyramidalization angle `alpha` (degrees).
    ///
    /// Average over the three substituents of the angle between a substituent
    /// vector and the plane-normal of the other two: `0°` for a fully
    /// pyramidal apex, `90°` for a planar centre.
    pub pyr_alpha: f32,
}

/// Computes pyramidalization descriptors from Cartesian coordinates.
#[derive(Clone, Copy, Debug, Default)]
pub struct PyramidalizationCalculator;

impl PyramidalizationCalculator {
    /// Computes `pyr_P` and `pyr_alpha` for a donor and its three substituents.
    ///
    /// `neighbor_indices` are the three atoms bonded to the donor (the same
    /// covalently-bonded substituents used for the buried-volume quadrant
    /// frame). The descriptors are invariant to the order of the three
    /// neighbours. Out-of-bounds indices, a neighbour coincident with the
    /// donor, non-finite coordinates, or two collinear substituent vectors
    /// yield the zero-valued default.
    #[must_use]
    pub fn compute(
        molecule: &Molecule,
        donor_idx: usize,
        neighbor_indices: [usize; 3],
    ) -> PyramidalizationParams {
        let Some(donor) = molecule.atoms.get(donor_idx) else {
            return PyramidalizationParams::default();
        };
        if !donor.position.is_finite() {
            return PyramidalizationParams::default();
        }

        // Unit vectors from the donor to each substituent.
        let mut unit = [Vec3::ZERO; 3];
        for (slot, &neighbor_idx) in unit.iter_mut().zip(neighbor_indices.iter()) {
            if neighbor_idx == donor_idx {
                return PyramidalizationParams::default();
            }
            let Some(neighbor) = molecule.atoms.get(neighbor_idx) else {
                return PyramidalizationParams::default();
            };
            let bond = neighbor.position - donor.position;
            if !bond.is_finite() || bond.length_squared() <= f32::EPSILON {
                return PyramidalizationParams::default();
            }
            *slot = bond.normalize();
        }
        let [a, b, c] = unit;

        // pyr_P is the absolute scalar triple product of the three unit
        // vectors — identical to morfeus' `sin(theta) * cos(alpha)` evaluated on
        // the permutation with a positive determinant, but order-invariant.
        let mut pyr_p = a.dot(b.cross(c)).abs();

        // pyr_alpha averages the signed out-of-plane angle for each substituent
        // taken in turn as the apex vector; the other two define the reference
        // plane, ordered so their cross product points toward the apex.
        let mut alpha_sum = 0.0_f32;
        for apex in 0..3 {
            let v3 = unit[apex];
            let mut v1 = unit[(apex + 1) % 3];
            let mut v2 = unit[(apex + 2) % 3];
            let mut normal = v1.cross(v2);
            if v3.dot(normal) < 0.0 {
                std::mem::swap(&mut v1, &mut v2);
                normal = -normal;
            }
            if normal.length_squared() <= f32::EPSILON {
                return PyramidalizationParams::default();
            }
            let cos_alpha = v3.dot(normal.normalize()).clamp(-1.0, 1.0);
            let mut alpha = cos_alpha.acos();
            // A pyramid whose apex leans onto the same side as the other two
            // vectors' bisector is "acute": the angle is signed negative. Only
            // the sign of the projection matters, so the bisector need not be
            // normalised (its length is positive by the guard).
            let bisector = v1 + v2;
            if bisector.length_squared() > f32::EPSILON && bisector.dot(v3) > 0.0 {
                alpha = -alpha;
            }
            alpha_sum += alpha;
        }
        let mean_alpha = alpha_sum / 3.0;

        // Correct P for an inverted pyramid, matching morfeus.
        if mean_alpha < 0.0 {
            pyr_p = 2.0 - pyr_p;
        }

        let params = PyramidalizationParams {
            pyr_p,
            pyr_alpha: mean_alpha.to_degrees(),
        };
        if params.pyr_p.is_finite() && params.pyr_alpha.is_finite() {
            params
        } else {
            PyramidalizationParams::default()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Atom;
    use std::f32::consts::PI;

    fn molecule_from(donor: Vec3, neighbors: [Vec3; 3]) -> Molecule {
        let mut atoms = vec![Atom::new("P", donor)];
        for position in neighbors {
            atoms.push(Atom::new("C", position));
        }
        Molecule { atoms }
    }

    #[test]
    fn planar_center_is_unpyramidalized() {
        let neighbors = [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new((2.0 * PI / 3.0).cos(), (2.0 * PI / 3.0).sin(), 0.0),
            Vec3::new((4.0 * PI / 3.0).cos(), (4.0 * PI / 3.0).sin(), 0.0),
        ];
        let molecule = molecule_from(Vec3::ZERO, neighbors);
        let params = PyramidalizationCalculator::compute(&molecule, 0, [1, 2, 3]);
        assert!(params.pyr_p.abs() < 1.0e-5, "P = {}", params.pyr_p);
        assert!(
            (params.pyr_alpha - 90.0).abs() < 1.0e-3,
            "alpha = {}",
            params.pyr_alpha
        );
    }

    #[test]
    fn tetrahedral_apex_matches_morfeus_reference() {
        // Three of the four tetrahedral vertices; morfeus reports
        // pyr_P = 0.769800, pyr_alpha = 35.264390.
        let neighbors = [
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(1.0, -1.0, -1.0),
            Vec3::new(-1.0, 1.0, -1.0),
        ];
        let molecule = molecule_from(Vec3::ZERO, neighbors);
        let params = PyramidalizationCalculator::compute(&molecule, 0, [1, 2, 3]);
        assert!(
            (params.pyr_p - 0.769_8).abs() < 1.0e-4,
            "P = {}",
            params.pyr_p
        );
        assert!(
            (params.pyr_alpha - 35.264_39).abs() < 1.0e-3,
            "alpha = {}",
            params.pyr_alpha
        );
    }

    #[test]
    fn descriptors_are_invariant_to_neighbor_order() {
        let neighbors = [
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(1.0, -1.0, -1.0),
            Vec3::new(-1.0, 1.0, -1.0),
        ];
        let molecule = molecule_from(Vec3::ZERO, neighbors);
        let base = PyramidalizationCalculator::compute(&molecule, 0, [1, 2, 3]);
        for order in [[2, 3, 1], [3, 1, 2], [3, 2, 1], [1, 3, 2], [2, 1, 3]] {
            let permuted = PyramidalizationCalculator::compute(&molecule, 0, order);
            assert!((permuted.pyr_p - base.pyr_p).abs() < 1.0e-6);
            assert!((permuted.pyr_alpha - base.pyr_alpha).abs() < 1.0e-4);
        }
    }

    #[test]
    fn descriptors_are_translation_invariant() {
        let neighbors = [
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(1.0, -1.0, -1.0),
            Vec3::new(-1.0, 1.0, -1.0),
        ];
        let origin = molecule_from(Vec3::ZERO, neighbors);
        let shift = Vec3::new(-3.0, 7.5, 2.0);
        let shifted = molecule_from(shift, neighbors.map(|position| position + shift));
        let first = PyramidalizationCalculator::compute(&origin, 0, [1, 2, 3]);
        let second = PyramidalizationCalculator::compute(&shifted, 0, [1, 2, 3]);
        assert!((first.pyr_p - second.pyr_p).abs() < 1.0e-6);
        assert!((first.pyr_alpha - second.pyr_alpha).abs() < 1.0e-4);
    }

    #[test]
    fn invalid_inputs_return_default() {
        let molecule = molecule_from(
            Vec3::ZERO,
            [
                Vec3::new(1.0, 1.0, 1.0),
                Vec3::new(1.0, -1.0, -1.0),
                Vec3::new(-1.0, 1.0, -1.0),
            ],
        );
        // Out-of-bounds neighbour.
        assert_eq!(
            PyramidalizationCalculator::compute(&molecule, 0, [1, 2, 9]),
            PyramidalizationParams::default()
        );
        // Neighbour coincident with the donor.
        assert_eq!(
            PyramidalizationCalculator::compute(&molecule, 0, [0, 2, 3]),
            PyramidalizationParams::default()
        );
    }
}
