//! Optional analytical union certificates; the original chord API stays unchanged.
use super::*;
use std::time::Instant;

/// Validate and canonically order the exact represented balls for new APIs.
pub(super) fn prepare(
    region: AtomicSphere,
    atoms: &[AtomicSphere],
    options: CertifiedVolumeOptions,
) -> Result<Vec<AtomicSphere>, DescriptorError> {
    if !options.absolute_tolerance.is_finite()
        || options.absolute_tolerance <= 0.0
        || options.max_cells == 0
    {
        return Err(error(
            "certified tolerance and cell budget must be positive and finite",
        ));
    }
    for s in std::iter::once(&region).chain(atoms) {
        if !s.radius.is_finite()
            || s.radius < 0.0
            || s.radius > 1e100
            || s.center.iter().any(|x| !x.is_finite() || x.abs() > 1e100)
        {
            return Err(error(
                "continuous sphere inputs must be finite, |coordinate| and radius <= 1e100, radius >= 0",
            ));
        }
    }
    if region.radius == 0.0 {
        return Err(error("integration radius must be positive"));
    }
    let mut atoms: Vec<_> = atoms.iter().copied().filter(|s| s.radius > 0.0).collect();
    atoms.sort_by(|a, b| {
        a.center[0]
            .total_cmp(&b.center[0])
            .then(a.center[1].total_cmp(&b.center[1]))
            .then(a.center[2].total_cmp(&b.center[2]))
            .then(a.radius.total_cmp(&b.radius))
    });
    atoms.dedup();
    Ok(atoms)
}

// Knuth TwoSum, with no overflow in the validated coordinate domain. Gradual
// underflow is required. The exact sum is s+e; the sign of e selects the only
// direction that needs widening. Unlike unconditional widening, exact dyadic
// contact distances remain exact. Also valid for subtraction via b -> -b.
pub(super) fn tight_add(a: f64, b: f64) -> Interval {
    let s = a + b;
    let v = s - a;
    let e = (a - (s - v)) + (b - v);
    Interval {
        lo: if e < 0.0 { s.down() } else { s },
        hi: if e > 0.0 { s.up() } else { s },
    }
}

fn disjoint(a: AtomicSphere, b: AtomicSphere) -> bool {
    let sum = tight_add(a.radius, b.radius);
    // A coordinate separation is a lower bound on the Euclidean separation.
    // Equality is sufficient: a tangency has zero Lebesgue volume.
    for i in 0..3 {
        let d = tight_add(a.center[i].max(b.center[i]), -a.center[i].min(b.center[i]));
        if d.lo >= sum.hi {
            return true;
        }
    }
    distance(a, b).lo >= sum.hi
}

fn triangle_free(
    region: AtomicSphere,
    atoms: &[AtomicSphere],
    evaluations: &mut usize,
) -> Option<Interval> {
    // Do not use pair-only inclusion-exclusion on clipped balls or when any
    // triple may have positive overlap. False negatives just use chord bounds.
    for &atom in atoms {
        *evaluations += 1;
        if !contained(atom, region) {
            return None;
        }
    }
    let mut neighbors: Vec<Vec<usize>> = vec![Vec::new(); atoms.len()];
    let mut overlaps = Interval::point(0.0);
    for (i, &a) in atoms.iter().enumerate() {
        for (j, &b) in atoms[..i].iter().enumerate() {
            *evaluations += 1;
            if disjoint(a, b) {
                continue;
            }
            if neighbors[i].iter().any(|k| neighbors[j].contains(k)) {
                return None;
            }
            neighbors[i].push(j);
            neighbors[j].push(i);
            overlaps = overlaps.add(intersection(a, b)?);
        }
    }
    let sum = atoms
        .iter()
        .fold(Interval::point(0.0), |v, a| v.add(ball_volume(a.radius)));
    let mut result = sum.sub(overlaps);
    result.lo = result.lo.max(0.0);
    result.hi = result.hi.min(ball_volume(region.radius).hi);
    (result.lo <= result.hi && result.hi.is_finite()).then_some(result)
}

pub(super) fn analytical_result(
    bounds: Interval,
    options: CertifiedVolumeOptions,
    evaluations: usize,
    elapsed_seconds: f64,
) -> CertifiedVolumeResult {
    let lo = bounds.lo;
    let hi = bounds.hi;
    let value = lo + (hi - lo) * 0.5;
    // Error-free subtraction avoids adding an unnecessary ULP when an endpoint
    // difference is exact, including [0, minimum-positive-subnormal].
    let absolute_error = tight_add(value, -lo).hi.max(tight_add(hi, -value).hi);
    let target_met = absolute_error <= options.absolute_tolerance;
    let relative = if lo > 0.0 {
        (absolute_error / lo).up()
    } else {
        f64::INFINITY
    };
    CertifiedVolumeResult {
        options,
        volume_angstrom3: value,
        lower_bound_angstrom3: lo,
        upper_bound_angstrom3: hi,
        absolute_error_bound_angstrom3: absolute_error,
        relative_error_bound: relative.is_finite().then_some(relative),
        target_met,
        status: if target_met {
            CertifiedVolumeStatus::ToleranceAchieved
        } else {
            CertifiedVolumeStatus::FloatingPointLimit
        },
        geometric_evaluations: evaluations,
        subdivisions: 0,
        cells_evaluated: 0,
        analytical: true,
        smooth_cells: 0,
        range_cells: 0,
        empty_cells: 0,
        peak_heap_cells: 0,
        elapsed_seconds,
    }
}

/// Explicit geometric extension of [`sphere_union_volume_certified`].
///
/// Certifies contained unions by pair inclusion-exclusion only if a conservative
/// overlap graph is triangle-free. Exact axis contacts count as disjoint in
/// volume. If this proof fails, uses the original certified chord integrator.
/// Integration spheres with entire volume below the minimum positive f64 use
/// a proven global volume bound. No rescaling, molecular tolerance change, or
/// inference of zero volume from underflow occurs.
///
/// Preprocessing has O(n²) pair queries; graph triangle checks can take O(n³)
/// time in the worst case. Its work is counted separately from the cell budget.
/// Original finite-grid, reference and certified-chord entry points are unchanged.
pub fn sphere_union_volume_certified_geometric(
    region: AtomicSphere,
    atoms: &[AtomicSphere],
    options: CertifiedVolumeOptions,
) -> Result<CertifiedVolumeResult, DescriptorError> {
    let start = Instant::now();
    let atoms = prepare(region, atoms, options)?;
    let mut evaluations = 0;
    let no_overlap = atoms.iter().all(|&atom| {
        evaluations += 1;
        disjoint(region, atom)
    });
    if no_overlap {
        return Ok(analytical_result(
            Interval::point(0.0),
            options,
            evaluations,
            start.elapsed().as_secs_f64(),
        ));
    }
    // pi < 6, hence 4*pi*R³/3 < 8R³. For R <= 2^-359 this is <= 2^-1074.
    // Construct the exact power of two by exponent bits (no pow rounding).
    let tiny_radius = f64::from_bits((1023_u64 - 359) << 52);
    if region.radius <= tiny_radius {
        return Ok(analytical_result(
            Interval {
                lo: 0.0,
                hi: f64::from_bits(1),
            },
            options,
            evaluations,
            start.elapsed().as_secs_f64(),
        ));
    }
    if let Some(bounds) = triangle_free(region, &atoms, &mut evaluations) {
        return Ok(analytical_result(
            bounds,
            options,
            evaluations,
            start.elapsed().as_secs_f64(),
        ));
    }
    let mut result = sphere_union_volume_certified(region, &atoms, options)?;
    result.geometric_evaluations += evaluations;
    result.elapsed_seconds = start.elapsed().as_secs_f64();
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn tight_add_matches_independent_exact_fraction_oracle() {
        let rows: Vec<[u64; 4]> = serde_json::from_str(include_str!(
            "../../../tests/fixtures/geometric_sum_oracle.json"
        ))
        .unwrap();
        for [a, b, lo, hi] in rows {
            let result = tight_add(f64::from_bits(a), f64::from_bits(b));
            assert_eq!(result.lo, f64::from_bits(lo), "{a:x} {b:x}");
            assert_eq!(result.hi, f64::from_bits(hi), "{a:x} {b:x}");
        }
    }
}
