//! Shared XY subdivisions for certified quadrant and octant integrals.
use super::super::geometric::{prepare, tight_add};
use super::*;

/// A certified percentage with the exact spherical-region volume as denominator.
#[derive(Clone, Copy, Debug, serde::Serialize)]
pub struct CertifiedPercentage {
    pub estimate: f64,
    pub lower_bound: f64,
    pub upper_bound: f64,
    pub absolute_error_bound: f64,
}

/// One member of a shared regional integration. Work and runtime belong to the
/// enclosing result because one rectangle contributes to two octants and a quadrant.
#[derive(Clone, Copy, Debug, serde::Serialize)]
pub struct CertifiedRegionalMeasure {
    pub volume_angstrom3: f64,
    pub lower_bound_angstrom3: f64,
    pub upper_bound_angstrom3: f64,
    pub absolute_error_bound_angstrom3: f64,
    pub relative_error_bound: Option<f64>,
    /// None if the positive regional denominator cannot be enclosed away from zero.
    pub percent_buried_volume: Option<CertifiedPercentage>,
    pub target_met: bool,
    pub status: CertifiedVolumeStatus,
}

/// Certified volumes in a specified Cartesian frame, with planes through the
/// integration sphere center. Quadrants are (++),(-+),(--),(+-); octants are these
/// four at positive Z followed by the four at negative Z. Boundary planes have
/// measure zero. Bounds include rounding but not source/frame/chemical uncertainty.
#[derive(Clone, Debug, serde::Serialize)]
pub struct CertifiedRegionalVolumeResult {
    pub options: CertifiedVolumeOptions,
    pub integration_sphere: AtomicSphere,
    pub total: CertifiedRegionalMeasure,
    pub quadrants: [CertifiedRegionalMeasure; 4],
    pub octants: [CertifiedRegionalMeasure; 8],
    pub target_met: bool,
    pub status: CertifiedVolumeStatus,
    /// Shared budget counts XY rectangles, including four quadrant roots.
    pub cells_evaluated: usize,
    pub subdivisions: usize,
    pub geometric_evaluations: usize,
    /// Piece counters count hemisphere integrals, at most two per rectangle.
    pub smooth_pieces: usize,
    pub range_pieces: usize,
    pub empty_pieces: usize,
    pub peak_heap_cells: usize,
    /// Whole shared calculation, including preprocessing. Not per-bin timing.
    pub elapsed_seconds: f64,
}

pub(super) fn measure(
    bounds: Interval,
    denominator: Interval,
    tolerance: f64,
    stopped: CertifiedVolumeStatus,
) -> CertifiedRegionalMeasure {
    let volume = bounds.lo + (bounds.hi - bounds.lo) * 0.5;
    let error = tight_add(volume, -bounds.lo)
        .hi
        .max(tight_add(bounds.hi, -volume).hi);
    let percent = (denominator.lo > 0.0).then(|| {
        let p = bounds
            .chord_div_positive(denominator)
            .chord_mul(Interval::point(100.0));
        let lo = p.lo.max(0.0);
        let hi = p.hi.min(100.0);
        let value = lo + (hi - lo) * 0.5;
        CertifiedPercentage {
            estimate: value,
            lower_bound: lo,
            upper_bound: hi,
            absolute_error_bound: tight_add(value, -lo).hi.max(tight_add(hi, -value).hi),
        }
    });
    let met = error <= tolerance;
    CertifiedRegionalMeasure {
        volume_angstrom3: volume,
        lower_bound_angstrom3: bounds.lo,
        upper_bound_angstrom3: bounds.hi,
        absolute_error_bound_angstrom3: error,
        relative_error_bound: relative_error(error, bounds.lo),
        percent_buried_volume: percent,
        target_met: met,
        status: if met {
            CertifiedVolumeStatus::ToleranceAchieved
        } else {
            stopped
        },
    }
}

#[derive(Clone)]
struct RegionalCell {
    lo: [f64; 2],
    hi: [f64; 2],
    bounds: [Interval; 2],
    candidates: Vec<usize>,
    quadrant: usize,
    ordinal: usize,
}
impl RegionalCell {
    fn width(&self) -> f64 {
        (self.bounds[0].hi - self.bounds[0].lo) + (self.bounds[1].hi - self.bounds[1].lo)
    }
}
impl PartialEq for RegionalCell {
    fn eq(&self, b: &Self) -> bool {
        self.ordinal == b.ordinal
    }
}
impl Eq for RegionalCell {}
impl PartialOrd for RegionalCell {
    fn partial_cmp(&self, b: &Self) -> Option<Ordering> {
        Some(self.cmp(b))
    }
}
impl Ord for RegionalCell {
    fn cmp(&self, b: &Self) -> Ordering {
        self.width()
            .total_cmp(&b.width())
            .then_with(|| b.ordinal.cmp(&self.ordinal))
    }
}
#[derive(Default)]
struct RegionalScratch {
    chord: ChordScratch,
    base: Vec<ChordSpan>,
}

fn clip_plane(mut s: ChordSpan, positive: bool) -> Option<ChordSpan> {
    // sign=0 identifies the exact constant endpoint z=0, with zero derivatives.
    let plane = Some(ChordEndpoint { sphere: 0, sign: 0 });
    if positive {
        s.lower = if s.low.lo >= 0.0 {
            s.lower
        } else if s.low.hi <= 0.0 {
            plane
        } else {
            None
        };
        s.low = chord_max(s.low, Interval::point(0.0));
    } else {
        s.upper = if s.high.hi <= 0.0 {
            s.upper
        } else if s.high.lo >= 0.0 {
            plane
        } else {
            None
        };
        s.high = chord_min(s.high, Interval::point(0.0));
    }
    (s.low.lo < s.high.hi).then_some(s)
}
fn endpoint_value(e: ChordEndpoint, bounds: &[BoundedSphere], mid: [Interval; 2]) -> Interval {
    if e.sign == 0 {
        Interval::point(0.0)
    } else {
        chord_endpoint_value(e, bounds, mid)
    }
}
fn bound_piece(
    lo: [f64; 2],
    hi: [f64; 2],
    spheres: &[AtomicSphere],
    bounds: &[BoundedSphere],
    scratch: &mut ChordScratch,
) -> Result<Interval, DescriptorError> {
    let spans = &scratch.spans;
    if spans.is_empty() {
        scratch.counters.empty += 1;
        return Ok(Interval::point(0.0));
    }
    let widths =
        std::array::from_fn::<_, 2, _>(|i| Interval::point(hi[i]).sub(Interval::point(lo[i])));
    let area = widths[0].chord_mul(widths[1]);
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
    let mut enclosure = Interval {
        lo: lower,
        hi: upper,
    }
    .chord_mul(area);
    enclosure.lo = enclosure.lo.max(0.0);
    if !chord_smooth(
        spans,
        spheres,
        &mut scratch.order,
        &mut scratch.components,
        &mut scratch.endpoints,
    ) {
        scratch.counters.range += 1;
        return Ok(enclosure);
    }
    let mid = std::array::from_fn(|i| {
        Interval::point(lo[i])
            .add(Interval::point(hi[i]))
            .chord_mul(Interval::point(0.5))
    });
    let mut length = Interval::point(0.0);
    let mut curvature = Interval::point(0.0);
    for &(a, b) in &scratch.endpoints {
        scratch.counters.evaluations += usize::from(a.sign != 0) + usize::from(b.sign != 0);
        length = length.add(endpoint_value(b, bounds, mid).sub(endpoint_value(a, bounds, mid)));
        for e in [a, b] {
            if e.sign == 0 {
                continue;
            }
            scratch.counters.evaluations += 1;
            let s = bounds[e.sphere];
            let q = chord_q(s, lo, hi);
            let qmin = Interval::point(q.lo);
            let denominator = qmin.chord_mul(qmin.sqrt());
            if denominator.lo <= 0.0 || !denominator.hi.is_finite() {
                scratch.counters.range += 1;
                return Ok(enclosure);
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
        enclosure.lo = enclosure.lo.max(integral.lo);
        enclosure.hi = enclosure.hi.min(integral.hi);
        if enclosure.lo > enclosure.hi {
            return Err(error("contradictory regional smooth-patch bounds"));
        }
        scratch.counters.smooth += 1;
    } else {
        scratch.counters.range += 1;
    }
    Ok(enclosure)
}

#[allow(clippy::too_many_arguments)]
fn regional_cell(
    lo: [f64; 2],
    hi: [f64; 2],
    indices: &[usize],
    spheres: &[AtomicSphere],
    bounds: &[BoundedSphere],
    quadrant: usize,
    ordinal: usize,
    scratch: &mut RegionalScratch,
) -> Result<RegionalCell, DescriptorError> {
    let mut cell = RegionalCell {
        lo,
        hi,
        bounds: [Interval::point(0.0); 2],
        candidates: Vec::new(),
        quadrant,
        ordinal,
    };
    scratch.chord.counters.evaluations += 1;
    let Some(region) = chord_span(0, bounds[0], lo, hi) else {
        scratch.chord.counters.empty += 2;
        return Ok(cell);
    };
    scratch.base.clear();
    for &i in indices {
        scratch.chord.counters.evaluations += 1;
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
        scratch.base.push(s);
    }
    for side in 0..2 {
        scratch.chord.spans.clear();
        scratch.chord.spans.extend(
            scratch
                .base
                .iter()
                .filter_map(|&s| clip_plane(s, side == 0)),
        );
        cell.bounds[side] = bound_piece(lo, hi, spheres, bounds, &mut scratch.chord)?;
    }
    Ok(cell)
}

fn sum_intervals(xs: impl IntoIterator<Item = Interval>) -> Interval {
    let mut sum = xs.into_iter().fold(Interval::point(0.0), |a, b| a.add(b));
    sum.lo = sum.lo.max(0.0);
    sum
}

/// Integrate all quadrants/octants with shared chord evaluations and subdivisions.
///
/// The tolerance must hold for the total and every region. One global cell budget
/// covers the four quadrant roots and both Z halves; max_cells must be >=4.
/// Positive/negative Z bounds share each cell's geometric chord queries, candidate
/// list and subdivision. Regional percentages use V(sphere)/4 or /8, not a lattice
/// count. Total and quadrant intervals are outward sums of the constituent octants
/// (with the known spherical cap on the total). Estimates need not sum bitwise.
pub fn sphere_union_volume_certified_regions(
    region: AtomicSphere,
    atoms: &[AtomicSphere],
    options: CertifiedVolumeOptions,
) -> Result<CertifiedRegionalVolumeResult, DescriptorError> {
    let start = Instant::now();
    let atoms = prepare(region, atoms, options)?;
    if options.max_cells < 4 {
        return Err(error(
            "regional cell budget must include four quadrant roots",
        ));
    }
    let spheres: Vec<_> = std::iter::once(region).chain(atoms).collect();
    let bounds: Vec<_> = spheres
        .iter()
        .map(|&s| relative_sphere(s, region.center))
        .collect();
    integrate_regions(region, &spheres, &bounds, options, start)
}

fn integrate_regions(
    region: AtomicSphere,
    spheres: &[AtomicSphere],
    bounds: &[BoundedSphere],
    options: CertifiedVolumeOptions,
    start: Instant,
) -> Result<CertifiedRegionalVolumeResult, DescriptorError> {
    let indices: Vec<_> = (1..spheres.len()).collect();
    let sphere_volume = ball_volume(region.radius);
    let quadrant_volume = sphere_volume.chord_mul(Interval::point(0.25));
    let octant_volume = sphere_volume.chord_mul(Interval::point(0.125));
    let mut scratch = RegionalScratch::default();
    let mut heap = BinaryHeap::new();
    let mut lower = [0.; 8];
    let mut upper = [0.; 8];
    let r = region.radius;
    let roots = [
        ([0., 0.], [r, r]),
        ([-r, 0.], [0., r]),
        ([-r, -r], [0., 0.]),
        ([0., -r], [r, 0.]),
    ];
    for (q, (lo, hi)) in roots.into_iter().enumerate() {
        let cell = regional_cell(lo, hi, &indices, spheres, bounds, q, q, &mut scratch)?;
        for side in 0..2 {
            lower[q + 4 * side] = cell.bounds[side].lo;
            upper[q + 4 * side] = cell.bounds[side].hi;
        }
        if cell.width() > 0.0 {
            heap.push(cell);
        }
    }
    let mut best_lo = [0.0_f64; 8];
    let mut best_hi = [octant_volume.hi; 8];
    let mut cells = 4;
    let mut subdivisions = 0;
    let mut peak = heap.len();
    let (total, quadrants, octants, status) = loop {
        for i in 0..8 {
            best_lo[i] = best_lo[i].max(lower[i]);
            best_hi[i] = best_hi[i].min(upper[i]);
            if best_lo[i] > best_hi[i] {
                return Err(error("contradictory certified regional bounds"));
            }
        }
        let intervals = std::array::from_fn::<_, 8, _>(|i| Interval {
            lo: best_lo[i],
            hi: best_hi[i],
        });
        let qb =
            std::array::from_fn::<_, 4, _>(|i| sum_intervals([intervals[i], intervals[i + 4]]));
        let mut tb = sum_intervals(qb);
        tb.lo = tb.lo.max(0.0);
        tb.hi = tb.hi.min(sphere_volume.hi);
        if tb.lo > tb.hi {
            return Err(error("contradictory regional total bounds"));
        }
        let stopped = if heap.is_empty() {
            CertifiedVolumeStatus::FloatingPointLimit
        } else {
            CertifiedVolumeStatus::CellBudgetExhausted
        };
        let total = measure(tb, sphere_volume, options.absolute_tolerance, stopped);
        let quadrants =
            qb.map(|b| measure(b, quadrant_volume, options.absolute_tolerance, stopped));
        let octants =
            intervals.map(|b| measure(b, octant_volume, options.absolute_tolerance, stopped));
        if total.target_met
            && quadrants.iter().all(|x| x.target_met)
            && octants.iter().all(|x| x.target_met)
        {
            break (
                total,
                quadrants,
                octants,
                CertifiedVolumeStatus::ToleranceAchieved,
            );
        }
        if cells > options.max_cells.saturating_sub(2) || heap.is_empty() {
            break (total, quadrants, octants, stopped);
        }
        let parent = heap.pop().expect("nonempty checked");
        let axis = usize::from(parent.hi[1] - parent.lo[1] >= parent.hi[0] - parent.lo[0]);
        let middle = parent.lo[axis] + (parent.hi[axis] - parent.lo[axis]) * 0.5;
        if middle <= parent.lo[axis] || middle >= parent.hi[axis] {
            continue;
        }
        let mut left = parent.hi;
        left[axis] = middle;
        let mut right = parent.lo;
        right[axis] = middle;
        for side in 0..2 {
            let i = parent.quadrant + 4 * side;
            lower[i] = (lower[i] - parent.bounds[side].lo).down();
            upper[i] = (upper[i] - parent.bounds[side].hi).up();
        }
        for (lo, hi) in [(parent.lo, left), (right, parent.hi)] {
            let child = regional_cell(
                lo,
                hi,
                &parent.candidates,
                spheres,
                bounds,
                parent.quadrant,
                cells,
                &mut scratch,
            )?;
            cells += 1;
            for side in 0..2 {
                let i = parent.quadrant + 4 * side;
                lower[i] = (lower[i] + child.bounds[side].lo).down();
                upper[i] = (upper[i] + child.bounds[side].hi).up();
            }
            if child.width() > 0.0 {
                heap.push(child);
            }
        }
        subdivisions += 1;
        peak = peak.max(heap.len());
    };
    let c = scratch.chord.counters;
    Ok(CertifiedRegionalVolumeResult {
        options,
        integration_sphere: region,
        total,
        quadrants,
        octants,
        target_met: status == CertifiedVolumeStatus::ToleranceAchieved,
        status,
        cells_evaluated: cells,
        subdivisions,
        geometric_evaluations: c.evaluations,
        smooth_pieces: c.smooth,
        range_pieces: c.range,
        empty_pieces: c.empty,
        peak_heap_cells: peak,
        elapsed_seconds: start.elapsed().as_secs_f64(),
    })
}
