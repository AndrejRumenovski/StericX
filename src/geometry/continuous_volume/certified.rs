//! Adaptive interval-enclosed chord integration with certified smooth-patch remainders.
use super::*;
use std::time::Instant;

#[path = "regional.rs"]
mod regional;
pub use regional::{
    CertifiedPercentage, CertifiedRadiusProfile, CertifiedRadiusProfileEntry,
    CertifiedRegionalMeasure, CertifiedRegionalVolumeResult,
    sphere_union_volume_certified_radius_profile, sphere_union_volume_certified_regions,
};

/// Work and absolute-error requirements for the explicit certified integrator.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct CertifiedVolumeOptions {
    /// Requested absolute error in Å³. Must be positive and finite.
    pub absolute_tolerance: f64,
    /// Maximum number of classified XY rectangles, including the root.
    pub max_cells: usize,
}
impl Default for CertifiedVolumeOptions {
    fn default() -> Self {
        Self {
            absolute_tolerance: 0.05,
            max_cells: 262143,
        }
    }
}

/// Why the certified calculation stopped.
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum CertifiedVolumeStatus {
    ToleranceAchieved,
    CellBudgetExhausted,
    FloatingPointLimit,
}

/// Certified continuous volume for exact supplied binary64 spheres.
///
/// Requires IEEE round-to-nearest and gradual underflow. Bounds include numerical
/// integration and arithmetic error, but not coordinate/radius uncertainty.
/// Runtime is diagnostic and is not part of deterministic scientific output.
#[derive(Clone, Debug, serde::Serialize)]
pub struct CertifiedVolumeResult {
    pub options: CertifiedVolumeOptions,
    pub volume_angstrom3: f64,
    pub lower_bound_angstrom3: f64,
    pub upper_bound_angstrom3: f64,
    pub absolute_error_bound_angstrom3: f64,
    /// Error divided by the certified positive lower volume; None at zero or
    /// when that relative bound is not representable as finite binary64.
    pub relative_error_bound: Option<f64>,
    pub target_met: bool,
    pub status: CertifiedVolumeStatus,
    /// Sphere/rectangle and sphere/point chord evaluations, including remainder
    /// evaluations, plus preprocessing sphere-pair queries.
    /// This counts geometric queries, not individual arithmetic operations.
    pub geometric_evaluations: usize,
    pub subdivisions: usize,
    pub cells_evaluated: usize,
    pub analytical: bool,
    pub smooth_cells: usize,
    pub range_cells: usize,
    pub empty_cells: usize,
    pub peak_heap_cells: usize,
    pub elapsed_seconds: f64,
}

fn relative_error(error: f64, lower: f64) -> Option<f64> {
    if lower <= 0.0 {
        return None;
    }
    let bound = (error / lower).up();
    bound.is_finite().then_some(bound)
}

fn analytical_bounds(
    region: AtomicSphere,
    atoms: &[AtomicSphere],
    evaluations: &mut usize,
) -> Option<Interval> {
    if atoms.is_empty() {
        return Some(Interval::point(0.0));
    }
    for &atom in atoms {
        *evaluations += 1;
        if contained(region, atom) {
            return Some(ball_volume(region.radius));
        }
    }
    if atoms.len() == 1 {
        *evaluations += 1;
        return intersection(region, atoms[0]);
    }
    if atoms.len() == 2 {
        for &atom in atoms {
            *evaluations += 1;
            if !contained(atom, region) {
                return None;
            }
        }
        *evaluations += 1;
        return intersection(atoms[0], atoms[1]).map(|overlap| {
            ball_volume(atoms[0].radius)
                .add(ball_volume(atoms[1].radius))
                .sub(overlap)
        });
    }
    None
}

/// Integrate the intersection of a ball with an atomic-sphere union.
///
/// Uses outward-rounded analytical cases when available. Otherwise, adaptively
/// integrates Z-chord unions over XY rectangles. A derivative remainder is used
/// only after interval predicates prove a fixed smooth endpoint topology;
/// unresolved rectangles retain conservative chord-range bounds. Refinement
/// returns nested global bounds for a fixed input as the work budget increases.
/// Failure to reach the tolerance is explicit in both status and target_met.
///
/// Input coordinates/radii have the same domain as [`sphere_union_volume`].
/// Existing finite-grid and continuous-reference methods remain unchanged.
pub fn sphere_union_volume_certified(
    region: AtomicSphere,
    atoms: &[AtomicSphere],
    options: CertifiedVolumeOptions,
) -> Result<CertifiedVolumeResult, DescriptorError> {
    let start = Instant::now();
    let tolerance = options.absolute_tolerance;
    let budget = options.max_cells;
    if !tolerance.is_finite() || tolerance <= 0.0 || budget == 0 {
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
    let mut preprocessing_evaluations = 0;
    if let Some(bounds) = analytical_bounds(region, &atoms, &mut preprocessing_evaluations)
        .filter(|b| b.lo.is_finite() && b.hi.is_finite() && b.lo <= b.hi)
    {
        let result = bounded_result(
            bounds,
            ContinuousVolumeMethod::Reference {
                absolute_tolerance: tolerance,
                max_boxes: 1,
            },
            1,
            true,
            tolerance,
        );
        let lower = result.lower_bound_angstrom3.unwrap();
        let error = result.absolute_error_bound_angstrom3.unwrap();
        let met = result.target_met.unwrap();
        return Ok(CertifiedVolumeResult {
            options,
            volume_angstrom3: result.volume_angstrom3,
            lower_bound_angstrom3: lower,
            upper_bound_angstrom3: result.upper_bound_angstrom3.unwrap(),
            absolute_error_bound_angstrom3: error,
            relative_error_bound: relative_error(error, lower),
            target_met: met,
            status: if met {
                CertifiedVolumeStatus::ToleranceAchieved
            } else {
                CertifiedVolumeStatus::FloatingPointLimit
            },
            geometric_evaluations: preprocessing_evaluations,
            subdivisions: 0,
            cells_evaluated: 0,
            analytical: true,
            smooth_cells: 0,
            range_cells: 0,
            empty_cells: 0,
            peak_heap_cells: 0,
            elapsed_seconds: start.elapsed().as_secs_f64(),
        });
    }
    let spheres: Vec<_> = std::iter::once(region).chain(atoms).collect();
    let bounds: Vec<_> = spheres
        .iter()
        .map(|&s| relative_sphere(s, region.center))
        .collect();
    let mut scratch = ChordScratch {
        counters: ChordCounters {
            evaluations: preprocessing_evaluations,
            ..Default::default()
        },
        ..Default::default()
    };
    let root = chord_cell(
        [-region.radius; 2],
        [region.radius; 2],
        &(1..spheres.len()).collect::<Vec<_>>(),
        &spheres,
        &bounds,
        0,
        &mut scratch,
    )?;
    let mut lower = root.bounds.lo;
    let mut upper = root.bounds.hi;
    let mut best_lower = lower;
    let mut best_upper = upper.min(ball_volume(region.radius).hi);
    let mut heap = BinaryHeap::new();
    heap.push(root);
    let mut cells = 1usize;
    let mut subdivisions = 0usize;
    let mut peak = 1usize;
    let (estimate, error, met, status) = loop {
        best_lower = best_lower.max(lower);
        best_upper = best_upper.min(upper);
        if best_lower > best_upper {
            return Err(error("contradictory certified chord bounds"));
        }
        let value = best_lower + (best_upper - best_lower) * 0.5;
        let error = if best_lower == best_upper {
            0.0
        } else {
            (value - best_lower).max(best_upper - value).up()
        };
        if error <= tolerance {
            break (value, error, true, CertifiedVolumeStatus::ToleranceAchieved);
        }
        if cells > budget.saturating_sub(2) {
            break (
                value,
                error,
                false,
                CertifiedVolumeStatus::CellBudgetExhausted,
            );
        }
        let Some(parent) = heap.pop() else {
            break (
                value,
                error,
                false,
                CertifiedVolumeStatus::FloatingPointLimit,
            );
        };
        let axis = usize::from(parent.hi[1] - parent.lo[1] >= parent.hi[0] - parent.lo[0]);
        let middle = parent.lo[axis] + (parent.hi[axis] - parent.lo[axis]) * 0.5;
        if middle <= parent.lo[axis] || middle >= parent.hi[axis] {
            continue;
        }
        let mut left = parent.hi;
        left[axis] = middle;
        let mut right = parent.lo;
        right[axis] = middle;
        lower = (lower - parent.bounds.lo).down();
        upper = (upper - parent.bounds.hi).up();
        for (lo, hi) in [(parent.lo, left), (right, parent.hi)] {
            let child = chord_cell(
                lo,
                hi,
                &parent.candidates,
                &spheres,
                &bounds,
                cells,
                &mut scratch,
            )?;
            cells += 1;
            lower = (lower + child.bounds.lo).down();
            upper = (upper + child.bounds.hi).up();
            if child.bounds.hi > child.bounds.lo {
                heap.push(child);
            }
        }
        subdivisions += 1;
        peak = peak.max(heap.len());
    };
    let counters = scratch.counters;
    Ok(CertifiedVolumeResult {
        options,
        volume_angstrom3: estimate,
        lower_bound_angstrom3: best_lower,
        upper_bound_angstrom3: best_upper,
        absolute_error_bound_angstrom3: error,
        relative_error_bound: relative_error(error, best_lower),
        target_met: met,
        status,
        geometric_evaluations: counters.evaluations,
        subdivisions,
        cells_evaluated: cells,
        analytical: false,
        smooth_cells: counters.smooth,
        range_cells: counters.range,
        empty_cells: counters.empty,
        peak_heap_cells: peak,
        elapsed_seconds: start.elapsed().as_secs_f64(),
    })
}

impl BuriedVolumeCalculator {
    /// Explicit certified total volume using the same f64 scaling and explicit
    /// center convention as compute_continuous_with_center. No frame averaging.
    pub fn compute_certified_with_center(
        molecule: &Molecule,
        center: [f64; 3],
        sphere_radius: f64,
        radii_scale: f64,
        include_hydrogens: bool,
        options: CertifiedVolumeOptions,
    ) -> Result<CertifiedVolumeResult, DescriptorError> {
        let start = Instant::now();
        if !radii_scale.is_finite() || radii_scale <= 0.0 {
            return Err(error("continuous radii scale must be finite and positive"));
        }
        let mut atoms = Vec::with_capacity(molecule.atoms.len());
        for atom in &molecule.atoms {
            if !include_hydrogens && atom.element.eq_ignore_ascii_case("H") {
                continue;
            }
            if !atom.position.is_finite() || !atom.vdw_radius.is_finite() || atom.vdw_radius <= 0.0
            {
                return Err(error("atom has invalid coordinates or radius"));
            }
            atoms.push(AtomicSphere {
                center: atom.position.as_dvec3().to_array(),
                radius: f64::from(atom.vdw_radius) * radii_scale,
            });
        }
        let mut result = sphere_union_volume_certified(
            AtomicSphere {
                center,
                radius: sphere_radius,
            },
            &atoms,
            options,
        )?;
        result.elapsed_seconds = start.elapsed().as_secs_f64();
        Ok(result)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct ChordEndpoint {
    sphere: usize,
    sign: i8,
}
#[derive(Clone, Copy)]
struct ChordSpan {
    low: Interval,
    high: Interval,
    lower: Option<ChordEndpoint>,
    upper: Option<ChordEndpoint>,
    regular: bool,
}
#[derive(Clone)]
struct ChordCell {
    lo: [f64; 2],
    hi: [f64; 2],
    bounds: Interval,
    candidates: Vec<usize>,
    ordinal: usize,
}
impl PartialEq for ChordCell {
    fn eq(&self, b: &Self) -> bool {
        self.ordinal == b.ordinal
    }
}
impl Eq for ChordCell {}
impl PartialOrd for ChordCell {
    fn partial_cmp(&self, b: &Self) -> Option<Ordering> {
        Some(self.cmp(b))
    }
}
impl Ord for ChordCell {
    fn cmp(&self, b: &Self) -> Ordering {
        (self.bounds.hi - self.bounds.lo)
            .total_cmp(&(b.bounds.hi - b.bounds.lo))
            .then_with(|| b.ordinal.cmp(&self.ordinal))
    }
}
#[derive(Default)]
struct ChordCounters {
    evaluations: usize,
    smooth: usize,
    range: usize,
    empty: usize,
}

fn chord_q(s: BoundedSphere, lo: [f64; 2], hi: [f64; 2]) -> Interval {
    let mut q = s.radius_squared;
    for i in 0..2 {
        q = q.sub(
            Interval {
                lo: lo[i],
                hi: hi[i],
            }
            .sub(s.center[i])
            .square(),
        );
    }
    q
}
fn chord_span(index: usize, s: BoundedSphere, lo: [f64; 2], hi: [f64; 2]) -> Option<ChordSpan> {
    let q = chord_q(s, lo, hi);
    if q.hi <= 0.0 {
        return None;
    }
    let extent = q.sqrt();
    Some(ChordSpan {
        low: s.center[2].sub(extent),
        high: s.center[2].add(extent),
        lower: Some(ChordEndpoint {
            sphere: index,
            sign: -1,
        }),
        upper: Some(ChordEndpoint {
            sphere: index,
            sign: 1,
        }),
        regular: q.lo > 0.0,
    })
}
fn chord_max(a: Interval, b: Interval) -> Interval {
    Interval {
        lo: a.lo.max(b.lo),
        hi: a.hi.max(b.hi),
    }
}
fn chord_min(a: Interval, b: Interval) -> Interval {
    Interval {
        lo: a.lo.min(b.lo),
        hi: a.hi.min(b.hi),
    }
}
fn chord_length(ranges: &mut Vec<(f64, f64)>, upper: bool) -> f64 {
    ranges.retain(|&(a, b)| a < b);
    ranges.sort_unstable_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
    let Some(&(mut a, mut b)) = ranges.first() else {
        return 0.0;
    };
    let mut sum = 0.0;
    let round = |x: f64| if upper { x.up() } else { x.down() };
    for &(c, d) in &ranges[1..] {
        if c > b {
            sum = round(sum + round(b - a));
            a = c;
            b = d;
        } else {
            b = b.max(d);
        }
    }
    round(sum + round(b - a)).max(0.0)
}
fn same_endpoint(
    a: Option<ChordEndpoint>,
    b: Option<ChordEndpoint>,
    spheres: &[AtomicSphere],
) -> bool {
    match (a, b) {
        (Some(a), Some(b)) => a.sign == b.sign && spheres[a.sphere] == spheres[b.sphere],
        _ => false,
    }
}
#[derive(Default)]
// Reused only within one integration. No geometry or decisions survive calls.
struct ChordScratch {
    counters: ChordCounters,
    spans: Vec<ChordSpan>,
    ranges: Vec<(f64, f64)>,
    order: Vec<usize>,
    components: Vec<(usize, usize)>,
    endpoints: Vec<(ChordEndpoint, ChordEndpoint)>,
}
fn chord_smooth(
    spans: &[ChordSpan],
    spheres: &[AtomicSphere],
    order: &mut Vec<usize>,
    components: &mut Vec<(usize, usize)>,
    endpoints: &mut Vec<(ChordEndpoint, ChordEndpoint)>,
) -> bool {
    if spans
        .iter()
        .any(|s| !s.regular || s.lower.is_none() || s.upper.is_none() || s.low.hi > s.high.lo)
    {
        return false;
    }
    order.clear();
    order.extend(0..spans.len());
    order.sort_unstable_by(|&a, &b| spans[a].low.hi.total_cmp(&spans[b].low.hi).then(a.cmp(&b)));
    components.clear();
    let mut inner_high = f64::NEG_INFINITY;
    for (position, &index) in order.iter().enumerate() {
        let s = spans[index];
        if s.low.hi > inner_high {
            components.push((position, position + 1));
            inner_high = s.high.lo;
        } else {
            inner_high = inner_high.max(s.high.lo);
            let Some(component) = components.last_mut() else {
                return false;
            };
            component.1 = position + 1;
        }
    }
    let mut previous_upper = f64::NEG_INFINITY;
    endpoints.clear();
    for &(start, end) in components.iter() {
        let group = &order[start..end];
        let Some(&lower) = group
            .iter()
            .min_by(|&&a, &&b| spans[a].low.lo.total_cmp(&spans[b].low.lo))
        else {
            return false;
        };
        let Some(&upper) = group
            .iter()
            .max_by(|&&a, &&b| spans[a].high.hi.total_cmp(&spans[b].high.hi))
        else {
            return false;
        };
        let low = spans[lower];
        let high = spans[upper];
        if low.low.lo < previous_upper {
            return false;
        }
        for &i in group {
            if i != lower
                && low.low.hi > spans[i].low.lo
                && !same_endpoint(low.lower, spans[i].lower, spheres)
            {
                return false;
            }
            if i != upper
                && high.high.lo < spans[i].high.hi
                && !same_endpoint(high.upper, spans[i].upper, spheres)
            {
                return false;
            }
        }
        previous_upper = high.high.hi;
        let (Some(a), Some(b)) = (low.lower, high.upper) else {
            return false;
        };
        endpoints.push((a, b));
    }
    true
}

fn chord_endpoint_value(
    endpoint: ChordEndpoint,
    spheres: &[BoundedSphere],
    mid: [Interval; 2],
) -> Interval {
    let s = spheres[endpoint.sphere];
    let mut q = s.radius_squared;
    for (i, m) in mid.into_iter().enumerate() {
        q = q.sub(m.sub(s.center[i]).square());
    }
    let extent = q.sqrt();
    if endpoint.sign < 0 {
        s.center[2].sub(extent)
    } else {
        s.center[2].add(extent)
    }
}
fn chord_cell(
    lo: [f64; 2],
    hi: [f64; 2],
    indices: &[usize],
    spheres: &[AtomicSphere],
    bounds: &[BoundedSphere],
    ordinal: usize,
    scratch: &mut ChordScratch,
) -> Result<ChordCell, DescriptorError> {
    let counters = &mut scratch.counters;
    let mut cell = ChordCell {
        lo,
        hi,
        bounds: Interval::point(0.0),
        candidates: Vec::new(),
        ordinal,
    };
    counters.evaluations += 1;
    let Some(region) = chord_span(0, bounds[0], lo, hi) else {
        counters.empty += 1;
        return Ok(cell);
    };
    let widths =
        std::array::from_fn::<_, 2, _>(|i| Interval::point(hi[i]).sub(Interval::point(lo[i])));
    let area = widths[0].chord_mul(widths[1]);
    let spans = &mut scratch.spans;
    spans.clear();
    for &i in indices {
        counters.evaluations += 1;
        let Some(mut s) = chord_span(i, bounds[i], lo, hi) else {
            continue;
        };
        let same = spheres[i] == spheres[0];
        s.lower = if same || s.low.lo >= region.low.hi {
            s.lower
        } else if region.low.lo >= s.low.hi {
            region.lower
        } else {
            None
        };
        s.upper = if same || s.high.hi <= region.high.lo {
            s.upper
        } else if region.high.hi <= s.high.lo {
            region.upper
        } else {
            None
        };
        s.low = chord_max(s.low, region.low);
        s.high = chord_min(s.high, region.high);
        s.regular &= region.regular;
        if s.low.lo >= s.high.hi {
            continue;
        }
        cell.candidates.push(i);
        spans.push(s);
    }
    if spans.is_empty() {
        counters.empty += 1;
        return Ok(cell);
    }
    scratch.ranges.clear();
    scratch.ranges.extend(
        spans
            .iter()
            .filter(|s| s.regular)
            .map(|s| (s.low.hi, s.high.lo)),
    );
    let lower = chord_length(&mut scratch.ranges, false);
    scratch.ranges.clear();
    scratch
        .ranges
        .extend(spans.iter().map(|s| (s.low.lo, s.high.hi)));
    let upper = chord_length(&mut scratch.ranges, true);
    cell.bounds = Interval {
        lo: lower,
        hi: upper,
    }
    .chord_mul(area);
    cell.bounds.lo = cell.bounds.lo.max(0.0);
    if !chord_smooth(
        spans,
        spheres,
        &mut scratch.order,
        &mut scratch.components,
        &mut scratch.endpoints,
    ) {
        counters.range += 1;
        return Ok(cell);
    };
    let mid = std::array::from_fn(|i| {
        Interval::point(lo[i])
            .add(Interval::point(hi[i]))
            .chord_mul(Interval::point(0.5))
    });
    let mut length = Interval::point(0.0);
    let mut curvature = Interval::point(0.0);
    for &(a, b) in &scratch.endpoints {
        counters.evaluations += 2;
        length = length
            .add(chord_endpoint_value(b, bounds, mid).sub(chord_endpoint_value(a, bounds, mid)));
        for endpoint in [a, b] {
            counters.evaluations += 1;
            let s = bounds[endpoint.sphere];
            let q = chord_q(s, lo, hi);
            let qmin = Interval::point(q.lo);
            let denominator = qmin.chord_mul(qmin.sqrt());
            if denominator.lo <= 0.0 || !denominator.hi.is_finite() {
                counters.range += 1;
                return Ok(cell);
            }
            curvature =
                curvature.add(Interval::point(s.radius_squared.hi).chord_div_positive(denominator));
        }
    }
    let remainder = area
        .chord_mul(widths[0].square().add(widths[1].square()))
        .chord_mul(curvature)
        .chord_div_positive(Interval::point(24.0));
    let integral = length.chord_mul(area).add(Interval {
        lo: -remainder.hi,
        hi: remainder.hi,
    });
    if integral.lo.is_finite() && integral.hi.is_finite() {
        cell.bounds.lo = cell.bounds.lo.max(integral.lo);
        cell.bounds.hi = cell.bounds.hi.min(integral.hi);
        if cell.bounds.lo > cell.bounds.hi {
            return Err(error("contradictory certified smooth-patch bounds"));
        }
        counters.smooth += 1;
    } else {
        counters.range += 1;
    }
    Ok(cell)
}

// For nonnegative finite endpoints, monotonic rounded multiplication selects
// the same extrema as the general four-product implementation. Keep its exact
// zero handling and fall back for signed or infinite intervals. Only C1 uses
// these methods; legacy continuous-reference arithmetic is unchanged.
trait PositiveChordArithmetic {
    fn chord_mul(self, other: Self) -> Self;
    fn chord_div_positive(self, other: Self) -> Self;
}
impl PositiveChordArithmetic for Interval {
    fn chord_mul(self, b: Self) -> Self {
        if self.lo >= 0.0 && b.lo >= 0.0 && self.hi.is_finite() && b.hi.is_finite() {
            if (self.lo == 0.0 && self.hi == 0.0) || (b.lo == 0.0 && b.hi == 0.0) {
                return Self::point(0.0);
            }
            Self {
                lo: (self.lo * b.lo).down(),
                hi: (self.hi * b.hi).up(),
            }
        } else {
            self.mul(b)
        }
    }
    fn chord_div_positive(self, b: Self) -> Self {
        debug_assert!(b.lo > 0.0);
        self.chord_mul(Self {
            lo: (1.0 / b.hi).down(),
            hi: (1.0 / b.lo).up(),
        })
    }
}

#[cfg(test)]
mod positive_interval_tests {
    use super::*;
    fn equal(a: Interval, b: Interval) {
        assert_eq!(a.lo.to_bits(), b.lo.to_bits(), "lower {a:?} {b:?}");
        assert_eq!(a.hi.to_bits(), b.hi.to_bits(), "upper {a:?} {b:?}");
    }
    #[test]
    fn specialized_bounds_match_four_product_reference_across_exponents() {
        let mut state = 0xa45e_2689_b213_54df_u64;
        let mut random = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            f64::from_bits(state & 0x7fef_ffff_ffff_ffff)
        };
        for _ in 0..100000 {
            let (a, b, c, d) = (random(), random(), random(), random());
            let x = Interval {
                lo: a.min(b),
                hi: a.max(b),
            };
            let y = Interval {
                lo: c.min(d),
                hi: c.max(d),
            };
            equal(x.chord_mul(y), x.mul(y));
            if y.lo > 0.0 {
                equal(x.chord_div_positive(y), x.div_positive(y));
            }
        }
        let cases = [
            Interval::point(0.0),
            Interval::point(-0.0),
            Interval::point(f64::from_bits(1)),
            Interval {
                lo: 0.0,
                hi: f64::MIN_POSITIVE,
            },
            Interval { lo: -1.0, hi: 1.0 },
            Interval { lo: -2.0, hi: -1.0 },
            Interval {
                lo: 1.0,
                hi: f64::MAX,
            },
            Interval {
                lo: 0.0,
                hi: f64::INFINITY,
            },
            Interval {
                lo: f64::NEG_INFINITY,
                hi: f64::INFINITY,
            },
        ];
        for x in cases {
            for y in cases {
                equal(x.chord_mul(y), x.mul(y));
                if y.lo > 0.0 {
                    equal(x.chord_div_positive(y), x.div_positive(y));
                }
            }
        }
    }
}
