// Research-only specialization, preserving the original outward-rounded bits.
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
            Self { lo: (self.lo * b.lo).down(), hi: (self.hi * b.hi).up() }
        } else {
            self.mul(b)
        }
    }
    fn chord_div_positive(self, b: Self) -> Self {
        debug_assert!(b.lo > 0.0);
        self.chord_mul(Self { lo: (1.0 / b.hi).down(), hi: (1.0 / b.lo).up() })
    }
}
