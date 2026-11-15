//! Continuous angular support, explicitly separate from the default f32 scan.
use super::{DescriptorError, Molecule, SterimolCalculator};
use glam::{DVec2, DVec3, Vec3};

/// Continuous-angle Sterimol envelope, computed and returned in f64 (Å).
///
/// Atomic coordinates and radii still have the source molecule's f32 precision.
/// This removes angular discretization, not input or floating-point error.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize)]
pub struct ContinuousSterimolParams {
    /// Raw longitudinal support, without the historical 0.40 Å correction.
    pub l: f64,
    /// Minimum signed transverse support over the full circle, clipped at zero.
    pub b1: f64,
    /// Maximum transverse distance plus atomic radius.
    pub b5: f64,
    /// Number of directions whose complete support envelope was evaluated.
    pub support_directions: usize,
}

impl SterimolCalculator {
    /// Evaluate the continuous angular envelope, excluding the attachment atom.
    ///
    /// Enumerates support stationary points and pairwise equality directions.
    /// The candidate set contains a global minimum in real arithmetic. This is
    /// an O(N³) reference method; no certified floating-point bound is claimed.
    /// Existing `compute` calls retain their original angular scan and precision.
    pub fn compute_continuous(
        molecule: &Molecule,
        attach_idx: usize,
        neighbor_idx: usize,
    ) -> Result<ContinuousSterimolParams, DescriptorError> {
        if attach_idx == neighbor_idx || molecule.atoms.len() < 2 {
            return Err(DescriptorError(
                "Sterimol requires two distinct axis atoms".into(),
            ));
        }
        let origin = molecule
            .atoms
            .get(attach_idx)
            .ok_or_else(|| DescriptorError("attachment index is out of bounds".into()))?
            .position;
        let base = molecule
            .atoms
            .get(neighbor_idx)
            .ok_or_else(|| DescriptorError("neighbor index is out of bounds".into()))?
            .position;
        continuous_projection(molecule, origin, base, Some(attach_idx))
    }

    /// Evaluate the continuous envelope against a virtual dummy, including all
    /// real atoms. Uses the same quantity as [`Self::compute_continuous`].
    pub fn compute_continuous_with_dummy(
        molecule: &Molecule,
        base_idx: usize,
        dummy: Vec3,
    ) -> Result<ContinuousSterimolParams, DescriptorError> {
        let base = molecule
            .atoms
            .get(base_idx)
            .ok_or_else(|| DescriptorError("base index is out of bounds".into()))?
            .position;
        continuous_projection(molecule, dummy, base, None)
    }
}

fn continuous_projection(
    molecule: &Molecule,
    origin: Vec3,
    base: Vec3,
    excluded: Option<usize>,
) -> Result<ContinuousSterimolParams, DescriptorError> {
    crate::profile_scope!("sterimol", "sterimol::continuous_support");
    if !origin.is_finite() || !base.is_finite() {
        return Err(DescriptorError(
            "Sterimol axis coordinate is not finite".into(),
        ));
    }
    let axis = base.as_dvec3() - origin.as_dvec3();
    if axis.length_squared() == 0.0 {
        return Err(DescriptorError("Sterimol axis atoms coincide".into()));
    }
    let z = axis.normalize();
    // Cross with the least-aligned coordinate axis to avoid a small denominator.
    // A continuous support minimum does not depend on transverse-frame phase.
    let seed = if z.x.abs() <= z.y.abs() && z.x.abs() <= z.z.abs() {
        DVec3::X
    } else if z.y.abs() <= z.z.abs() {
        DVec3::Y
    } else {
        DVec3::Z
    };
    let x = z.cross(seed).normalize();
    let y = z.cross(x);
    let mut discs = Vec::with_capacity(molecule.atoms.len() + 1);
    // Treat the clipping at zero as another support branch so minima where it
    // becomes active also have a stationary/equality candidate.
    discs.push((DVec2::ZERO, 0.0));
    let mut l = f64::NEG_INFINITY;
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
        let radius = f64::from(atom.vdw_radius);
        l = l.max(relative.dot(z) + radius);
        b5 = b5.max(relative.cross(z).length() + radius);
        discs.push((DVec2::new(relative.dot(x), relative.dot(y)), radius));
    }
    if discs.len() == 1 {
        return Err(DescriptorError(
            "Sterimol envelope contains no atoms".into(),
        ));
    }
    let (b1, support_directions) = continuous_support(&discs);
    Ok(ContinuousSterimolParams {
        l,
        b1,
        b5,
        support_directions,
    })
}

fn continuous_support(discs: &[(DVec2, f64)]) -> (f64, usize) {
    let mut minimum = f64::INFINITY;
    let mut directions = 0;
    let mut evaluate = |direction: DVec2| {
        let support = discs
            .iter()
            .map(|(position, radius)| position.dot(direction) + radius)
            .fold(0.0_f64, f64::max);
        minimum = minimum.min(support);
        directions += 1;
    };
    evaluate(DVec2::X);
    for (i, &(position, radius)) in discs.iter().enumerate() {
        let distance = position.x.hypot(position.y);
        if distance > 0.0 {
            evaluate(-position / distance);
        }
        for &(other, other_radius) in &discs[..i] {
            let delta = position - other;
            let separation = delta.x.hypot(delta.y);
            let radius_delta = other_radius - radius;
            if separation == 0.0 || radius_delta.abs() > separation {
                continue;
            }
            let along = delta / separation;
            let cosine = (radius_delta / separation).clamp(-1.0, 1.0);
            let sine = ((1.0 - cosine) * (1.0 + cosine)).max(0.0).sqrt();
            let across = DVec2::new(-along.y, along.x);
            evaluate((cosine * along + sine * across).normalize());
            evaluate((cosine * along - sine * across).normalize());
        }
    }
    (minimum, directions)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Atom;

    fn molecule(discs: &[(f32, f32, f32)]) -> Molecule {
        let mut atoms = vec![Atom {
            element: "C".into(),
            position: Vec3::ZERO,
            vdw_radius: 0.0,
        }];
        atoms.extend(discs.iter().map(|&(x, y, r)| Atom {
            element: "C".into(),
            position: Vec3::new(x, y, 0.0),
            vdw_radius: r,
        }));
        Molecule { atoms }
    }
    fn compute(discs: &[(f32, f32, f32)]) -> ContinuousSterimolParams {
        SterimolCalculator::compute_continuous_with_dummy(&molecule(discs), 0, -Vec3::Z).unwrap()
    }

    #[test]
    fn analytical_continuous_envelopes() {
        let single = compute(&[(2.0, 0.0, 3.0)]);
        assert_eq!((single.l, single.b1, single.b5), (4.0, 1.0, 5.0));
        assert_eq!(compute(&[(2.0, 0.0, 0.5), (-2.0, 0.0, 0.5)]).b1, 0.5);
        assert_eq!(compute(&[(2.0, 0.0, 0.5)]).b1, 0.0);
        assert_eq!(compute(&[(0.0, 0.0, 3.0), (1.0, 0.0, 1.0)]).b1, 3.0);
        let square = compute(&[
            (2.0, 2.0, 0.5),
            (-2.0, 2.0, 0.5),
            (-2.0, -2.0, 0.5),
            (2.0, -2.0, 0.5),
        ]);
        assert_eq!(square.b1, 2.5);
        assert_eq!(square.b5, 8.0_f64.sqrt() + 0.5);
    }

    #[test]
    fn continuous_envelope_preserves_permutations_and_exact_rigid_transforms() {
        let mut m = molecule(&[(1.5, 0.5, 1.0), (-0.5, 2.0, 0.5), (-2.0, -1.0, 0.25)]);
        let original = SterimolCalculator::compute_continuous_with_dummy(&m, 0, -Vec3::Z).unwrap();
        m.atoms[1..].reverse();
        assert_eq!(
            SterimolCalculator::compute_continuous_with_dummy(&m, 0, -Vec3::Z).unwrap(),
            original
        );
        let transform = |p: Vec3| Vec3::new(p.z + 2.0, p.y + 4.0, -p.x + 1.0);
        for atom in &mut m.atoms {
            atom.position = transform(atom.position);
        }
        let transformed =
            SterimolCalculator::compute_continuous_with_dummy(&m, 0, transform(-Vec3::Z)).unwrap();
        for (a, b) in [
            (original.l, transformed.l),
            (original.b1, transformed.b1),
            (original.b5, transformed.b5),
        ] {
            assert!((a - b).abs() <= 32.0 * f64::EPSILON * a.abs().max(1.0));
        }
    }

    #[test]
    fn continuous_axis_and_input_validation() {
        let mut m = molecule(&[(0.0, 0.0, 1.0)]);
        assert!(SterimolCalculator::compute_continuous_with_dummy(&m, 0, Vec3::ZERO).is_err());
        assert!(SterimolCalculator::compute_continuous_with_dummy(&m, 3, -Vec3::Z).is_err());
        assert!(
            SterimolCalculator::compute_continuous_with_dummy(&m, 0, Vec3::splat(f32::NAN))
                .is_err()
        );
        assert!(
            SterimolCalculator::compute_continuous_with_dummy(&m, 0, -Vec3::Z * f32::from_bits(1))
                .is_ok()
        );
        m.atoms[1].vdw_radius = -1.0;
        assert!(SterimolCalculator::compute_continuous_with_dummy(&m, 0, -Vec3::Z).is_err());
    }

    #[test]
    fn continuous_real_attachment_is_excluded() {
        let molecule = Molecule {
            atoms: vec![
                Atom {
                    element: "C".into(),
                    position: Vec3::ZERO,
                    vdw_radius: 100.0,
                },
                Atom {
                    element: "C".into(),
                    position: Vec3::Z,
                    vdw_radius: 1.0,
                },
            ],
        };
        let value = SterimolCalculator::compute_continuous(&molecule, 0, 1).unwrap();
        assert_eq!((value.l, value.b1, value.b5), (2.0, 1.0, 1.0));
        assert!(SterimolCalculator::compute_continuous(&molecule, 0, 0).is_err());
    }
}
