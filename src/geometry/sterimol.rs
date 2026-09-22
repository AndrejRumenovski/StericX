use super::{DescriptorError, Molecule};
use glam::{DVec3, Vec2, Vec3};

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
    /// Computes raw geometric `L`, `B1`, and `B5` in ångströms.
    ///
    /// The axis points from attachment (excluded dummy) to neighbor. `B1`
    /// retains the one-degree signed-support scan in the shortest-arc transverse
    /// frame; its angular discretization is not exactly rotation invariant.
    /// `L` and `B5` use the supplied axis directly, without a near-axis rotation
    /// approximation. Invalid indices, non-finite coordinates/radii, negative
    /// radii, coincident axis atoms, or unrepresentable outputs return an error.
    pub fn compute(
        molecule: &Molecule,
        attach_idx: usize,
        neighbor_idx: usize,
    ) -> Result<SterimolParams, DescriptorError> {
        crate::profile_scope!("sterimol", "SterimolCalculator::compute");
        if attach_idx == neighbor_idx || molecule.atoms.len() < 2 {
            return Err(DescriptorError(
                "Sterimol requires two distinct axis atoms".into(),
            ));
        }
        let attachment = molecule
            .atoms
            .get(attach_idx)
            .ok_or_else(|| DescriptorError("attachment index is out of bounds".into()))?;
        let neighbor = molecule
            .atoms
            .get(neighbor_idx)
            .ok_or_else(|| DescriptorError("neighbor index is out of bounds".into()))?;
        project(
            molecule,
            attachment.position,
            neighbor.position,
            Some(attach_idx),
        )
    }

    /// Computes against an explicit virtual dummy, including every real atom.
    ///
    /// The axis points from `dummy` to `base_idx`. This returns raw geometric
    /// `L`; a caller using the historical Verloop convention adds 0.40 Å.
    /// Input and output validation matches [`Self::compute`].
    pub fn compute_with_dummy(
        molecule: &Molecule,
        base_idx: usize,
        dummy: Vec3,
    ) -> Result<SterimolParams, DescriptorError> {
        crate::profile_scope!("sterimol", "SterimolCalculator::compute_with_dummy");
        let base = molecule
            .atoms
            .get(base_idx)
            .ok_or_else(|| DescriptorError("base index is out of bounds".into()))?;
        project(molecule, dummy, base.position, None)
    }
}

fn checked_f32(value: f64) -> Result<f32, DescriptorError> {
    let result = value as f32;
    if result.is_finite() {
        Ok(result)
    } else {
        Err(DescriptorError(
            "Sterimol projection or descriptor exceeds finite f32 range".into(),
        ))
    }
}

fn project(
    molecule: &Molecule,
    origin: Vec3,
    base: Vec3,
    excluded: Option<usize>,
) -> Result<SterimolParams, DescriptorError> {
    if !origin.is_finite() || !base.is_finite() {
        return Err(DescriptorError(
            "Sterimol axis coordinate is not finite".into(),
        ));
    }
    // Wider intermediates specifically prevent the audited finite-coordinate
    // square overflow and short-axis underflow. Inputs and reported quantities
    // remain f32; the B1 scan below retains its f32 angles and support arithmetic.
    let axis = base.as_dvec3() - origin.as_dvec3();
    if axis.length_squared() == 0.0 {
        return Err(DescriptorError("Sterimol axis atoms coincide".into()));
    }
    let z = axis.normalize();
    // Stable exact shortest-arc basis. The two algebraic branches avoid
    // cancellation at either pole, unlike Quat::from_rotation_arc's deliberate
    // near-parallel approximation. An exactly antiparallel axis uses a Y half-turn.
    let transverse = z.x.hypot(z.y);
    let (x, y) = if transverse == 0.0 {
        (if z.z >= 0.0 { DVec3::X } else { -DVec3::X }, DVec3::Y)
    } else {
        let rotation_axis = DVec3::new(-z.y / transverse, z.x / transverse, 0.0);
        let (sin_half, cos_half) = if z.z >= 0.0 {
            let c = ((1.0 + z.z) * 0.5).sqrt();
            (transverse / (2.0 * c), c)
        } else {
            let s = ((1.0 - z.z) * 0.5).sqrt();
            (s, transverse / (2.0 * s))
        };
        // Rotate +Z to the supplied axis; these columns form the transverse
        // basis whose dot products rotate coordinates back onto +Z.
        let q = glam::DQuat::from_xyzw(
            rotation_axis.x * sin_half,
            rotation_axis.y * sin_half,
            0.0,
            cos_half,
        );
        (q * DVec3::X, q * DVec3::Y)
    };
    let mut projected = Vec::with_capacity(molecule.atoms.len());
    let mut b5 = 0.0_f64;
    for (index, atom) in molecule.atoms.iter().enumerate() {
        if !atom.position.is_finite() || !atom.vdw_radius.is_finite() || atom.vdw_radius < 0.0 {
            return Err(DescriptorError(
                "atom has non-finite coordinates/radius or negative radius".into(),
            ));
        }
        if excluded == Some(index) {
            continue;
        }
        let relative = atom.position.as_dvec3() - origin.as_dvec3();
        let axial = checked_f32(relative.dot(z))?;
        let xy = Vec2::new(checked_f32(relative.dot(x))?, checked_f32(relative.dot(y))?);
        // Direct cross-product radius is independent of transverse-frame phase.
        b5 = b5.max(relative.cross(z).length() + f64::from(atom.vdw_radius));
        projected.push((xy, axial, atom.vdw_radius));
    }
    params_from_projection(&projected, checked_f32(b5)?)
}

fn params_from_projection(
    projected: &[(Vec2, f32, f32)],
    b5: f32,
) -> Result<SterimolParams, DescriptorError> {
    crate::profile_scope!("sterimol", "sterimol::params_from_projection");
    if projected.is_empty() {
        return Err(DescriptorError(
            "Sterimol envelope contains no atoms".into(),
        ));
    }
    let l = projected
        .iter()
        .map(|(_, z, radius)| z + radius)
        .fold(f32::NEG_INFINITY, f32::max);
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
    if !l.is_finite() || !b1.is_finite() {
        return Err(DescriptorError(
            "Sterimol descriptor exceeds finite f32 range".into(),
        ));
    }
    Ok(SterimolParams { l, b1, b5 })
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

        let params = SterimolCalculator::compute(&molecule, 0, 1).unwrap();

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

        let params = SterimolCalculator::compute(&molecule, 0, 1).unwrap();

        assert!((params.l - 2.7).abs() < 1.0e-5);
        assert!((params.b1 - 1.2).abs() < 1.0e-5);
        assert!((params.b5 - 3.7).abs() < 1.0e-5);
    }

    #[test]
    fn attachment_dummy_radius_is_excluded() {
        let molecule = Molecule {
            atoms: vec![atom("I", 0.0, 0.0, 0.0), atom("H", 0.0, 0.0, 1.0)],
        };

        let params = SterimolCalculator::compute(&molecule, 0, 1).unwrap();

        assert!((params.l - 2.2).abs() < 1.0e-5);
        assert!((params.b1 - 1.2).abs() < 1.0e-5);
        assert!((params.b5 - 1.2).abs() < 1.0e-5);
    }

    #[test]
    fn dummy_axis_measures_length_from_the_virtual_origin() {
        // Donor P at the origin, a carbon 2 Å beyond it, and a virtual metal
        // dummy 2.28 Å below along the same line. The dummy is not a real atom,
        // so both P and C contribute; L is measured from the dummy origin.
        let molecule = Molecule {
            atoms: vec![atom("P", 0.0, 0.0, 0.0), atom("C", 0.0, 0.0, 2.0)],
        };

        let params =
            SterimolCalculator::compute_with_dummy(&molecule, 0, Vec3::new(0.0, 0.0, -2.28))
                .unwrap();

        // C at axial 4.28 + radius 1.70 = 5.98 is the farthest; P (2.28 + 1.80)
        // is nearer. Both lie on the axis, so B1 = B5 = the widest radius (P).
        assert!((params.l - 5.98).abs() < 1.0e-4);
        assert!((params.b1 - 1.80).abs() < 1.0e-4);
        assert!((params.b5 - 1.80).abs() < 1.0e-4);
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

        let first = SterimolCalculator::compute(&along_z, 0, 1).unwrap();
        let second = SterimolCalculator::compute(&translated_along_x, 0, 1).unwrap();

        assert!((first.l - second.l).abs() < 1.0e-5);
        assert!((first.b1 - second.b1).abs() < 1.0e-5);
        assert!((first.b5 - second.b5).abs() < 1.0e-5);
    }

    #[test]
    fn invalid_inputs_return_errors() {
        let molecule = Molecule {
            atoms: vec![atom("C", 0.0, 0.0, 0.0), atom("H", 0.0, 0.0, 0.0)],
        };
        assert!(SterimolCalculator::compute(&molecule, 0, 0).is_err());
        assert!(SterimolCalculator::compute(&molecule, 0, 1).is_err());
        assert!(SterimolCalculator::compute(&molecule, 2, 1).is_err());
    }
}
