//! Explicit continuous sphere-union integration, separate from the native grid.
use super::{BuriedVolumeCalculator, DescriptorError, Molecule};
use std::{cmp::Ordering, collections::BinaryHeap};

mod certified;
pub use certified::{
    CertifiedVolumeOptions, CertifiedVolumeResult, CertifiedVolumeStatus,
    sphere_union_volume_certified,
};

/// A mathematical closed ball. Coordinates/radius are in Å.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct AtomicSphere {
    pub center: [f64; 3],
    pub radius: f64,
}

/// Explicit numerical definitions; neither changes the existing finite grid.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(tag = "method", rename_all = "snake_case")]
pub enum ContinuousVolumeMethod {
    /// Midpoint XY quadrature of exact merged Z chords. No certified error bound.
    MidpointSlices { subdivisions: usize },
    /// Outward-rounded analytical cases, otherwise adaptive box enclosures.
    Reference {
        absolute_tolerance: f64,
        max_boxes: usize,
    },
}

/// Volume of the integration ball intersected with the union of atomic balls.
/// Bounds, when present, enclose this continuous quantity for the supplied f64
/// inputs under IEEE round-to-nearest with gradual underflow. They do not cover
/// coordinate/radius uncertainty or the physical adequacy of a hard-sphere model.
#[derive(Clone, Debug, PartialEq, serde::Serialize)]
pub struct ContinuousVolumeResult {
    pub method: ContinuousVolumeMethod,
    pub volume_angstrom3: f64,
    pub lower_bound_angstrom3: Option<f64>,
    pub upper_bound_angstrom3: Option<f64>,
    pub absolute_error_bound_angstrom3: Option<f64>,
    /// None for fixed-resolution quadrature; false when a reference target was
    /// not attained within the work or floating-point range/precision limit.
    pub target_met: Option<bool>,
    /// XY sample count or number of classified boxes (one for analytical cases).
    pub work_units: usize,
    pub analytical: bool,
}

impl BuriedVolumeCalculator {
    /// Continuous total volume at an explicit center. No donor-plane rotation
    /// or regional discretization is used. Radii are scaled in f64; this explicit
    /// method therefore has a distinct precision contract from the native grid.
    /// Bounds refer to the resulting rounded f64 radii, not to exact-real
    /// multiplication of an uncertain physical radius and scaling factor.
    pub fn compute_continuous_with_center(
        molecule: &Molecule,
        center: [f64; 3],
        sphere_radius: f64,
        radii_scale: f64,
        include_hydrogens: bool,
        method: ContinuousVolumeMethod,
    ) -> Result<ContinuousVolumeResult, DescriptorError> {
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
        sphere_union_volume(
            AtomicSphere {
                center,
                radius: sphere_radius,
            },
            &atoms,
            method,
        )
    }
}

/// Integrate an explicitly specified union of balls inside another ball.
///
/// Values must be finite and at most 1e100 in magnitude, to keep squared
/// distances and cubic volumes in a useful f64 range. Zero-radius atomic balls
/// have zero volume and are ignored; the integration radius must be positive.
pub fn sphere_union_volume(
    region: AtomicSphere,
    atoms: &[AtomicSphere],
    method: ContinuousVolumeMethod,
) -> Result<ContinuousVolumeResult, DescriptorError> {
    for sphere in std::iter::once(&region).chain(atoms) {
        if !sphere.radius.is_finite()
            || sphere.radius < 0.0
            || sphere.radius > 1e100
            || sphere
                .center
                .iter()
                .any(|x| !x.is_finite() || x.abs() > 1e100)
        {
            return Err(error(
                "continuous sphere inputs must be finite, |coordinate| and radius <= 1e100, radius >= 0",
            ));
        }
    }
    if region.radius == 0.0 {
        return Err(error("integration radius must be positive"));
    }
    let atoms: Vec<_> = atoms.iter().copied().filter(|s| s.radius > 0.0).collect();
    match method {
        ContinuousVolumeMethod::MidpointSlices { subdivisions } => {
            if subdivisions == 0 || subdivisions > 8192 {
                return Err(error("midpoint subdivisions must be in 1..=8192"));
            }
            midpoint(region, &atoms, subdivisions, method)
        }
        ContinuousVolumeMethod::Reference {
            absolute_tolerance,
            max_boxes,
        } => {
            if !absolute_tolerance.is_finite() || absolute_tolerance <= 0.0 || max_boxes == 0 {
                return Err(error(
                    "reference tolerance and box budget must be positive and finite",
                ));
            }
            Ok(reference(
                region,
                &atoms,
                absolute_tolerance,
                max_boxes,
                method,
            ))
        }
    }
}

fn error(message: &str) -> DescriptorError {
    DescriptorError(message.into())
}

fn midpoint(
    region: AtomicSphere,
    atoms: &[AtomicSphere],
    n: usize,
    method: ContinuousVolumeMethod,
) -> Result<ContinuousVolumeResult, DescriptorError> {
    let r = region.radius;
    let step = 2.0 * r / n as f64;
    if step == 0.0 {
        return Err(error("midpoint spacing underflows f64"));
    }
    let relative: Vec<_> = atoms
        .iter()
        .map(|a| {
            (
                std::array::from_fn::<_, 3, _>(|i| a.center[i] - region.center[i]),
                a.radius,
            )
        })
        .collect();
    let mut chords = Vec::with_capacity(atoms.len());
    let mut sum = 0.0;
    let mut compensation = 0.0;
    for i in 0..n {
        let x = -r + (i as f64 + 0.5) * step;
        for j in 0..n {
            let y = -r + (j as f64 + 0.5) * step;
            let remaining = r * r - x * x - y * y;
            if remaining <= 0.0 {
                continue;
            }
            let extent = remaining.sqrt();
            chords.clear();
            for &(c, radius) in &relative {
                let dx = x - c[0];
                let dy = y - c[1];
                let squared = radius * radius - dx * dx - dy * dy;
                if squared <= 0.0 {
                    continue;
                }
                let dz = squared.sqrt();
                let lo = (c[2] - dz).max(-extent);
                let hi = (c[2] + dz).min(extent);
                if lo < hi {
                    chords.push((lo, hi));
                }
            }
            chords.sort_unstable_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
            let mut length = 0.0;
            if let Some(&(mut lo, mut hi)) = chords.first() {
                for &(a, b) in &chords[1..] {
                    if a > hi {
                        length += hi - lo;
                        lo = a;
                        hi = b;
                    } else {
                        hi = hi.max(b);
                    }
                }
                length += hi - lo;
            }
            // Compensated summation limits accumulation error; it does not
            // bound the quadrature error at circle/interval-topology boundaries.
            let adjusted = length - compensation;
            let next = sum + adjusted;
            compensation = (next - sum) - adjusted;
            sum = next;
        }
    }
    let volume = (sum * step) * step;
    if !volume.is_finite() {
        return Err(error("continuous midpoint volume overflow"));
    }
    Ok(ContinuousVolumeResult {
        method,
        volume_angstrom3: volume,
        lower_bound_angstrom3: None,
        upper_bound_angstrom3: None,
        absolute_error_bound_angstrom3: None,
        target_met: None,
        work_units: n * n,
        analytical: false,
    })
}

// Adjacent floats without raising the declared Rust 1.85 minimum version.
trait Outward {
    fn up(self) -> Self;
    fn down(self) -> Self;
}
impl Outward for f64 {
    fn up(self) -> Self {
        if self.is_nan() || self == f64::INFINITY {
            return self;
        }
        if self == 0.0 {
            return f64::from_bits(1);
        }
        let bits = self.to_bits();
        f64::from_bits(if self > 0.0 { bits + 1 } else { bits - 1 })
    }
    fn down(self) -> Self {
        -(-self).up()
    }
}

// Every bound operation widens by one adjacent f64 after round-to-nearest.
// Exact zero special cases prevent gratuitous widening of degenerate intervals.
#[derive(Clone, Copy, Debug)]
struct Interval {
    lo: f64,
    hi: f64,
}
impl Interval {
    fn point(x: f64) -> Self {
        Self { lo: x, hi: x }
    }
    fn add(self, b: Self) -> Self {
        Self {
            lo: (self.lo + b.lo).down(),
            hi: (self.hi + b.hi).up(),
        }
    }
    fn sub(self, b: Self) -> Self {
        if self.lo == self.hi && b.lo == b.hi && self.lo == b.lo {
            return Self::point(0.0);
        }
        Self {
            lo: (self.lo - b.hi).down(),
            hi: (self.hi - b.lo).up(),
        }
    }
    fn mul(self, b: Self) -> Self {
        if (self.lo == 0.0 && self.hi == 0.0) || (b.lo == 0.0 && b.hi == 0.0) {
            return Self::point(0.0);
        }
        let products = [
            self.lo * b.lo,
            self.lo * b.hi,
            self.hi * b.lo,
            self.hi * b.hi,
        ];
        if products.iter().any(|p| p.is_nan()) {
            return Self {
                lo: f64::NEG_INFINITY,
                hi: f64::INFINITY,
            };
        }
        Self {
            lo: products
                .into_iter()
                .map(Outward::down)
                .fold(f64::INFINITY, f64::min),
            hi: products
                .into_iter()
                .map(Outward::up)
                .fold(f64::NEG_INFINITY, f64::max),
        }
    }
    fn div_positive(self, b: Self) -> Self {
        debug_assert!(b.lo > 0.0);
        self.mul(Self {
            lo: (1.0 / b.hi).down(),
            hi: (1.0 / b.lo).up(),
        })
    }
    fn square(self) -> Self {
        if self.lo == 0.0 && self.hi == 0.0 {
            return self;
        }
        let min = if self.lo <= 0.0 && self.hi >= 0.0 {
            0.0
        } else {
            self.lo.abs().min(self.hi.abs())
        };
        let max = self.lo.abs().max(self.hi.abs());
        Self {
            lo: (min * min).down().max(0.0),
            hi: (max * max).up(),
        }
    }
    fn sqrt(self) -> Self {
        if self.lo == 0.0 && self.hi == 0.0 {
            return self;
        }
        Self {
            lo: self.lo.max(0.0).sqrt().down().max(0.0),
            hi: self.hi.max(0.0).sqrt().up(),
        }
    }
}

fn ball_volume(radius: f64) -> Interval {
    let r = Interval::point(radius);
    // These neighboring floats enclose π; independently checked using Machin's
    // formula with rational alternating-series bounds in the campaign oracle.
    let pi = Interval {
        lo: std::f64::consts::PI.down(),
        hi: std::f64::consts::PI.up(),
    };
    r.mul(r)
        .mul(r)
        .mul(pi)
        .mul(Interval::point(4.0))
        .div_positive(Interval::point(3.0))
}
fn distance(a: AtomicSphere, b: AtomicSphere) -> Interval {
    let mut squared = Interval::point(0.0);
    for i in 0..3 {
        squared = squared.add(
            Interval::point(a.center[i])
                .sub(Interval::point(b.center[i]))
                .square(),
        );
    }
    squared.sqrt()
}
fn contained(inner: AtomicSphere, outer: AtomicSphere) -> bool {
    if inner.center == outer.center {
        return inner.radius <= outer.radius;
    }
    distance(inner, outer).add(Interval::point(inner.radius)).hi <= outer.radius
}

// For fixed radii, sphere-intersection volume decreases monotonically with d.
// Resolve each endpoint's mathematical branch using interval comparisons.
fn intersection_at_distance(r: f64, s: f64, d: f64) -> Option<Interval> {
    let small = ball_volume(r.min(s));
    let sum = Interval::point(r).add(Interval::point(s));
    let difference = Interval::point(r.max(s)).sub(Interval::point(r.min(s)));
    if d <= difference.lo {
        return Some(small);
    }
    if d >= sum.hi {
        return Some(Interval::point(0.0));
    }
    if d < difference.hi || d > sum.lo || d == 0.0 {
        return None;
    }
    let di = Interval::point(d);
    let height = sum.sub(di);
    let factor = di.add(sum.mul(Interval::point(2.0))).sub(
        difference
            .square()
            .div_positive(di)
            .mul(Interval::point(3.0)),
    );
    let pi = Interval {
        lo: std::f64::consts::PI.down(),
        hi: std::f64::consts::PI.up(),
    };
    let v = height
        .square()
        .mul(factor)
        .mul(pi)
        .div_positive(Interval::point(12.0));
    Some(Interval {
        lo: v.lo.max(0.0),
        hi: v.hi.min(small.hi).max(0.0),
    })
}
fn intersection(a: AtomicSphere, b: AtomicSphere) -> Option<Interval> {
    if a.center == b.center {
        return Some(ball_volume(a.radius.min(b.radius)));
    }
    let d = distance(a, b);
    let lower = intersection_at_distance(a.radius, b.radius, d.hi)?;
    let upper = intersection_at_distance(a.radius, b.radius, d.lo)?;
    Some(Interval {
        lo: lower.lo,
        hi: upper.hi,
    })
}

fn bounded_result(
    bounds: Interval,
    method: ContinuousVolumeMethod,
    work: usize,
    analytical: bool,
    tolerance: f64,
) -> ContinuousVolumeResult {
    let lo = bounds.lo.max(0.0);
    let hi = bounds.hi.max(lo);
    let value = lo + (hi - lo) * 0.5;
    let error = if lo == hi {
        0.0
    } else {
        (value - lo).max(hi - value).up()
    };
    ContinuousVolumeResult {
        method,
        volume_angstrom3: value,
        lower_bound_angstrom3: Some(lo),
        upper_bound_angstrom3: Some(hi),
        absolute_error_bound_angstrom3: Some(error),
        target_met: Some(error <= tolerance),
        work_units: work,
        analytical,
    }
}

#[derive(Clone, Copy)]
struct BoundedSphere {
    center: [Interval; 3],
    radius_squared: Interval,
}
fn relative_sphere(s: AtomicSphere, origin: [f64; 3]) -> BoundedSphere {
    BoundedSphere {
        center: std::array::from_fn(|i| {
            Interval::point(s.center[i]).sub(Interval::point(origin[i]))
        }),
        radius_squared: Interval::point(s.radius).square(),
    }
}
#[derive(Clone, Copy)]
struct BoxCell {
    lo: [f64; 3],
    hi: [f64; 3],
    volume: Interval,
    ordinal: usize,
}
impl PartialEq for BoxCell {
    fn eq(&self, b: &Self) -> bool {
        self.ordinal == b.ordinal
    }
}
impl Eq for BoxCell {}
impl PartialOrd for BoxCell {
    fn partial_cmp(&self, b: &Self) -> Option<Ordering> {
        Some(self.cmp(b))
    }
}
impl Ord for BoxCell {
    fn cmp(&self, b: &Self) -> Ordering {
        self.volume
            .hi
            .total_cmp(&b.volume.hi)
            .then_with(|| b.ordinal.cmp(&self.ordinal))
    }
}
impl BoxCell {
    fn new(lo: [f64; 3], hi: [f64; 3], ordinal: usize) -> Self {
        let mut volume = Interval::point(1.0);
        for i in 0..3 {
            volume = volume.mul(Interval::point(hi[i]).sub(Interval::point(lo[i])));
        }
        volume.lo = volume.lo.max(0.0);
        Self {
            lo,
            hi,
            volume,
            ordinal,
        }
    }
    fn squared_distance(self, s: BoundedSphere) -> Interval {
        let mut result = Interval::point(0.0);
        for i in 0..3 {
            let d = Interval {
                lo: (self.lo[i] - s.center[i].hi).down(),
                hi: (self.hi[i] - s.center[i].lo).up(),
            };
            result = result.add(d.square());
        }
        result.lo = result.lo.max(0.0);
        result
    }
}
enum Classification {
    Outside,
    Inside,
    Unresolved,
}
fn classify(cell: BoxCell, region: BoundedSphere, atoms: &[BoundedSphere]) -> Classification {
    let d = cell.squared_distance(region);
    if d.lo > region.radius_squared.hi {
        return Classification::Outside;
    }
    let in_region = d.hi <= region.radius_squared.lo;
    let mut possible = false;
    for &atom in atoms {
        let distance = cell.squared_distance(atom);
        if in_region && distance.hi <= atom.radius_squared.lo {
            return Classification::Inside;
        }
        possible |= distance.lo <= atom.radius_squared.hi;
    }
    if possible {
        Classification::Unresolved
    } else {
        Classification::Outside
    }
}

fn reference(
    region: AtomicSphere,
    atoms: &[AtomicSphere],
    tolerance: f64,
    max_boxes: usize,
    method: ContinuousVolumeMethod,
) -> ContinuousVolumeResult {
    let exact_case = if atoms.is_empty() {
        Some(Interval::point(0.0))
    } else if atoms.iter().any(|&a| contained(region, a)) {
        Some(ball_volume(region.radius))
    } else if atoms.len() == 1 {
        intersection(region, atoms[0])
    } else if atoms.len() == 2 && atoms.iter().all(|&a| contained(a, region)) {
        intersection(atoms[0], atoms[1]).map(|overlap| {
            ball_volume(atoms[0].radius)
                .add(ball_volume(atoms[1].radius))
                .sub(overlap)
        })
    } else {
        None
    };
    if let Some(bounds) = exact_case {
        return bounded_result(bounds, method, 1, true, tolerance);
    }
    let region_bound = relative_sphere(region, region.center);
    let atom_bounds: Vec<_> = atoms
        .iter()
        .map(|&a| relative_sphere(a, region.center))
        .collect();
    let root = BoxCell::new([-region.radius; 3], [region.radius; 3], 0);
    let mut lower = 0.0;
    let mut upper = root.volume.hi;
    let mut work = 1;
    let mut heap = BinaryHeap::new();
    match classify(root, region_bound, &atom_bounds) {
        Classification::Outside => upper = 0.0,
        Classification::Inside => lower = root.volume.lo,
        Classification::Unresolved => heap.push(root),
    }
    loop {
        let result = bounded_result(
            Interval {
                lo: lower,
                hi: upper,
            },
            method,
            work,
            false,
            tolerance,
        );
        if result.target_met == Some(true) || work > max_boxes.saturating_sub(2) {
            return result;
        }
        let Some(cell) = heap.pop() else {
            return result;
        };
        let axis = (0..3)
            .max_by(|&a, &b| (cell.hi[a] - cell.lo[a]).total_cmp(&(cell.hi[b] - cell.lo[b])))
            .unwrap();
        let middle = cell.lo[axis] + (cell.hi[axis] - cell.lo[axis]) * 0.5;
        if middle == cell.lo[axis] || middle == cell.hi[axis] {
            // Retain its unresolved upper contribution, but do not loop forever
            // on a cell whose endpoints have no representable interior split.
            continue;
        }
        let mut left_hi = cell.hi;
        left_hi[axis] = middle;
        let mut right_lo = cell.lo;
        right_lo[axis] = middle;
        upper = (upper - cell.volume.hi).up();
        for child in [
            BoxCell::new(cell.lo, left_hi, work),
            BoxCell::new(right_lo, cell.hi, work + 1),
        ] {
            work += 1;
            match classify(child, region_bound, &atom_bounds) {
                Classification::Outside => {}
                Classification::Inside => {
                    lower = (lower + child.volume.lo).down().max(0.0);
                    upper = (upper + child.volume.hi).up();
                }
                Classification::Unresolved => {
                    upper = (upper + child.volume.hi).up();
                    heap.push(child);
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn adjacent_float_steps_cover_signed_zero_and_extremes() {
        assert_eq!(0.0_f64.up().to_bits(), 1);
        assert_eq!((-0.0_f64).up().to_bits(), 1);
        assert_eq!(0.0_f64.down().to_bits(), (1_u64 << 63) | 1);
        assert_eq!(f64::MAX.up(), f64::INFINITY);
        assert_eq!(f64::INFINITY.down(), f64::MAX);
        assert_eq!(f64::NEG_INFINITY.up(), -f64::MAX);
        assert_eq!(f64::NEG_INFINITY.down(), f64::NEG_INFINITY);
        assert_eq!(1.0_f64.up().to_bits(), 1.0_f64.to_bits() + 1);
        assert_eq!((-1.0_f64).up().to_bits(), (-1.0_f64).to_bits() - 1);
        assert!(f64::NAN.up().is_nan());
    }

    #[test]
    fn union_order_does_not_change_chord_merging_or_box_results() {
        let region = sphere(0.0, 3.5);
        let atoms = [sphere(-1.0, 1.0), sphere(0.0, 0.75), sphere(1.0, 1.25)];
        let reverse = [atoms[2], atoms[1], atoms[0]];
        for method in [
            ContinuousVolumeMethod::MidpointSlices { subdivisions: 64 },
            ContinuousVolumeMethod::Reference {
                absolute_tolerance: 1e-5,
                max_boxes: 4095,
            },
        ] {
            assert_eq!(
                sphere_union_volume(region, &atoms, method).unwrap(),
                sphere_union_volume(region, &reverse, method).unwrap()
            );
        }
    }

    #[test]
    fn box_distance_bounds_enclose_independent_dyadic_integer_extrema() {
        let mut state = 0x824d_29a1_e723_b18d_u64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (state >> 32) as i64 % 129 - 64
        };
        for index in 0..10000 {
            let mut lo = [0_i64; 3];
            let mut hi = [0_i64; 3];
            let mut center = [0_i64; 3];
            let mut minimum = 0_i128;
            let mut maximum = 0_i128;
            for axis in 0..3 {
                let a = next();
                let b = next();
                lo[axis] = a.min(b);
                hi[axis] = a.max(b);
                center[axis] = next();
                let nearest = if center[axis] < lo[axis] {
                    lo[axis] - center[axis]
                } else if center[axis] > hi[axis] {
                    center[axis] - hi[axis]
                } else {
                    0
                };
                let farthest = (lo[axis] - center[axis])
                    .abs()
                    .max((hi[axis] - center[axis]).abs());
                minimum += i128::from(nearest).pow(2);
                maximum += i128::from(farthest).pow(2);
            }
            // All endpoints are exact eighths; extrema in units of 1/64 are
            // computed solely with integers, independently of interval helpers.
            let cell = BoxCell::new(
                lo.map(|x| x as f64 / 8.0),
                hi.map(|x| x as f64 / 8.0),
                index,
            );
            let sphere = relative_sphere(
                AtomicSphere {
                    center: center.map(|x| x as f64 / 8.0),
                    radius: 1.0,
                },
                [0.0; 3],
            );
            let actual = cell.squared_distance(sphere);
            assert!(actual.lo <= minimum as f64 / 64.0, "{index}");
            assert!(actual.hi >= maximum as f64 / 64.0, "{index}");
            let exact_volume =
                (0..3).map(|i| i128::from(hi[i] - lo[i])).product::<i128>() as f64 / 512.0;
            assert!(cell.volume.lo <= exact_volume && exact_volume <= cell.volume.hi);
        }
    }
    fn sphere(x: f64, r: f64) -> AtomicSphere {
        AtomicSphere {
            center: [x, 0.0, 0.0],
            radius: r,
        }
    }
    fn reference_method() -> ContinuousVolumeMethod {
        ContinuousVolumeMethod::Reference {
            absolute_tolerance: 1e-10,
            max_boxes: 10000,
        }
    }
    fn encloses(result: &ContinuousVolumeResult, value: f64) {
        assert!(
            result.lower_bound_angstrom3.unwrap() <= value
                && value <= result.upper_bound_angstrom3.unwrap(),
            "{result:?}, reference={value}"
        );
    }
    #[test]
    fn analytical_ball_lens_and_union_bounds() {
        let pi = std::f64::consts::PI;
        for (region, atoms, volume) in [
            (sphere(0.0, 2.0), vec![], 0.0),
            (sphere(0.0, 2.0), vec![sphere(0.0, 1.0)], 4.0 * pi / 3.0),
            (sphere(0.0, 1.0), vec![sphere(1.0, 1.0)], 5.0 * pi / 12.0),
            (sphere(0.0, 1.0), vec![sphere(2.0, 1.0)], 0.0),
            (
                sphere(0.0, 3.0),
                vec![sphere(-0.5, 1.0), sphere(0.5, 1.0)],
                9.0 * pi / 4.0,
            ),
        ] {
            let result = sphere_union_volume(region, &atoms, reference_method()).unwrap();
            encloses(&result, volume);
            assert_eq!(result.target_met, Some(true));
        }
    }
    #[test]
    fn adaptive_boxes_report_unresolved_budget_and_tighten() {
        let region = sphere(0.0, 2.0);
        // Three identical balls deliberately take the generic box path.
        let atoms = [sphere(0.0, 1.0); 3];
        let coarse = sphere_union_volume(
            region,
            &atoms,
            ContinuousVolumeMethod::Reference {
                absolute_tolerance: 1e-8,
                max_boxes: 1,
            },
        )
        .unwrap();
        let fine = sphere_union_volume(
            region,
            &atoms,
            ContinuousVolumeMethod::Reference {
                absolute_tolerance: 0.5,
                max_boxes: 20000,
            },
        )
        .unwrap();
        encloses(&coarse, 4.0 * std::f64::consts::PI / 3.0);
        encloses(&fine, 4.0 * std::f64::consts::PI / 3.0);
        assert_eq!(coarse.target_met, Some(false));
        assert!(fine.absolute_error_bound_angstrom3 < coarse.absolute_error_bound_angstrom3);
        assert!(fine.work_units <= 20000);
    }
    #[test]
    fn midpoint_chords_converge_for_centered_sphere() {
        let reference = 4.0 * std::f64::consts::PI / 3.0;
        let calculate = |n| {
            sphere_union_volume(
                sphere(0.0, 2.0),
                &[sphere(0.0, 1.0)],
                ContinuousVolumeMethod::MidpointSlices { subdivisions: n },
            )
            .unwrap()
        };
        let a = calculate(16);
        let b = calculate(128);
        assert!((b.volume_angstrom3 - reference).abs() < (a.volume_angstrom3 - reference).abs());
        assert!(b.absolute_error_bound_angstrom3.is_none());
        assert!(b.target_met.is_none());
    }
    #[test]
    fn continuous_volume_rejects_invalid_inputs() {
        assert!(sphere_union_volume(sphere(0.0, 0.0), &[], reference_method()).is_err());
        assert!(
            sphere_union_volume(
                sphere(0.0, 1.0),
                &[sphere(f64::NAN, 1.0)],
                reference_method()
            )
            .is_err()
        );
        assert!(
            sphere_union_volume(
                sphere(0.0, 1.0),
                &[],
                ContinuousVolumeMethod::MidpointSlices { subdivisions: 0 }
            )
            .is_err()
        );
    }
}
