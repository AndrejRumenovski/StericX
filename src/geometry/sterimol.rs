use super::Molecule;
use glam::{Quat, Vec2, Vec3};

/// Three-dimensional Sterimol dimensions in ångströms.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SterimolParams {
    /// Longitudinal extent along the attachment axis.
    pub l: f32,
    /// Minimum perpendicular radial half-width.
    pub b1: f32,
    /// Maximum perpendicular radial half-width.
    pub b5: f32,
}

/// Computes Sterimol descriptors from Cartesian coordinates and atomic radii.
#[derive(Clone, Copy, Debug, Default)]
pub struct SterimolCalculator;

impl SterimolCalculator {
    /// Computes `L`, `B1`, and `B5` for a molecule.
    ///
    /// The attachment atom is translated to the origin and the vector from
    /// `attach_idx` to `neighbor_idx` is rotated onto positive Z. Atomic van der
    /// Waals spheres are then projected onto the axial and perpendicular
    /// directions. The attachment atom is treated as the conventional dummy
    /// atom and excluded from the substituent envelope. `B1` is the minimum
    /// signed support radius found by a one-degree azimuthal scan; `B5` is the
    /// exact maximum radial support.
    ///
    /// Invalid indices, identical indices, coincident axis atoms, non-finite
    /// coordinates, or fewer than two atoms produce zero-valued parameters.
    #[must_use]
    pub fn compute(molecule: &Molecule, attach_idx: usize, neighbor_idx: usize) -> SterimolParams {
        if attach_idx == neighbor_idx || molecule.atoms.len() < 2 {
            return SterimolParams::default();
        }
        let Some(attachment) = molecule.atoms.get(attach_idx) else {
            return SterimolParams::default();
        };
        let Some(neighbor) = molecule.atoms.get(neighbor_idx) else {
            return SterimolParams::default();
        };

        let primary_axis = neighbor.position - attachment.position;
        if !primary_axis.is_finite() || primary_axis.length_squared() <= f32::EPSILON {
            return SterimolParams::default();
        }

        let rotation = Quat::from_rotation_arc(primary_axis.normalize(), Vec3::Z);
        let mut projected = Vec::with_capacity(molecule.atoms.len() - 1);
        for (atom_index, atom) in molecule.atoms.iter().enumerate() {
            if !atom.position.is_finite() || !atom.vdw_radius.is_finite() || atom.vdw_radius < 0.0 {
                return SterimolParams::default();
            }
            // The attachment atom is the conventional Sterimol dummy/origin,
            // not part of the substituent surface envelope.
            if atom_index == attach_idx {
                continue;
            }
            let aligned = rotation * (atom.position - attachment.position);
            projected.push((aligned.truncate(), aligned.z, atom.vdw_radius));
        }

        let l = projected
            .iter()
            .map(|(_, z, radius)| z + radius)
            .fold(f32::NEG_INFINITY, f32::max);

        let b5 = projected
            .iter()
            .map(|(xy, _, radius)| xy.length() + radius)
            .fold(0.0_f32, f32::max);

        let b1 = (0_u16..360)
            .map(|degrees| {
                let angle = f32::from(degrees).to_radians();
                let scan = Vec2::new(angle.cos(), angle.sin());
                projected
                    .iter()
                    .map(|(xy, _, radius)| xy.dot(scan) + radius)
                    .fold(0.0_f32, f32::max)
            })
            .fold(f32::INFINITY, f32::min);

        SterimolParams { l, b1, b5 }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Atom;

    fn atom(element: &str, x: f32, y: f32, z: f32) -> Atom {
        Atom::new(element, Vec3::new(x, y, z))
    }

    #[test]
    fn computes_linear_geometry_after_axis_alignment() {
        let molecule = Molecule {
            atoms: vec![atom("H", 0.0, 0.0, 0.0), atom("C", 2.0, 0.0, 0.0)],
        };

        let params = SterimolCalculator::compute(&molecule, 0, 1);

        assert!((params.l - 3.7).abs() < 1.0e-5);
        assert!((params.b1 - 1.7).abs() < 1.0e-5);
        assert!((params.b5 - 1.7).abs() < 1.0e-5);
    }

    #[test]
    fn computes_known_off_axis_widths() {
        let molecule = Molecule {
            atoms: vec![
                atom("H", 0.0, 0.0, 0.0),
                atom("H", 0.0, 0.0, 1.0),
                atom("C", 2.0, 0.0, 1.0),
            ],
        };

        let params = SterimolCalculator::compute(&molecule, 0, 1);

        assert!((params.l - 2.7).abs() < 1.0e-5);
        assert!((params.b1 - 1.2).abs() < 1.0e-5);
        assert!((params.b5 - 3.7).abs() < 1.0e-5);
    }

    #[test]
    fn attachment_dummy_radius_is_excluded() {
        let molecule = Molecule {
            atoms: vec![atom("I", 0.0, 0.0, 0.0), atom("H", 0.0, 0.0, 1.0)],
        };

        let params = SterimolCalculator::compute(&molecule, 0, 1);

        assert!((params.l - 2.2).abs() < 1.0e-5);
        assert!((params.b1 - 1.2).abs() < 1.0e-5);
        assert!((params.b5 - 1.2).abs() < 1.0e-5);
    }

    #[test]
    fn descriptors_are_translation_and_rotation_invariant() {
        let along_z = Molecule {
            atoms: vec![
                atom("H", 0.0, 0.0, 0.0),
                atom("H", 0.0, 0.0, 1.0),
                atom("O", 1.5, 0.0, 2.0),
            ],
        };
        let translated_along_x = Molecule {
            atoms: vec![
                atom("H", 4.0, -2.0, 3.0),
                atom("H", 5.0, -2.0, 3.0),
                atom("O", 6.0, -0.5, 3.0),
            ],
        };

        let first = SterimolCalculator::compute(&along_z, 0, 1);
        let second = SterimolCalculator::compute(&translated_along_x, 0, 1);

        assert!((first.l - second.l).abs() < 1.0e-5);
        assert!((first.b1 - second.b1).abs() < 1.0e-5);
        assert!((first.b5 - second.b5).abs() < 1.0e-5);
    }

    #[test]
    fn invalid_inputs_return_zeroes() {
        let one_atom = Molecule {
            atoms: vec![atom("C", 0.0, 0.0, 0.0)],
        };
        assert_eq!(
            SterimolCalculator::compute(&one_atom, 0, 0),
            SterimolParams::default()
        );

        let coincident = Molecule {
            atoms: vec![atom("C", 0.0, 0.0, 0.0), atom("H", 0.0, 0.0, 0.0)],
        };
        assert_eq!(
            SterimolCalculator::compute(&coincident, 0, 1),
            SterimolParams::default()
        );
        assert_eq!(
            SterimolCalculator::compute(&coincident, 2, 1),
            SterimolParams::default()
        );
    }
}
