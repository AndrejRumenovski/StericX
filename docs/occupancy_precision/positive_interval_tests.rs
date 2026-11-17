#[cfg(test)]
mod positive_interval_tests {
    use super::*;
    fn equal(a: Interval, b: Interval) {
        assert_eq!(a.lo.to_bits(), b.lo.to_bits(), "lower {a:?} {b:?}");
        assert_eq!(a.hi.to_bits(), b.hi.to_bits(), "upper {a:?} {b:?}");
    }
    #[test]
    fn specialized_bounds_match_four_product_reference_across_exponents() {
        let mut state=0xa45e_2689_b213_54df_u64;
        let mut random=|| {state^=state<<13;state^=state>>7;state^=state<<17;f64::from_bits(state & 0x7fef_ffff_ffff_ffff)};
        for _ in 0..100000 {
            let (a,b,c,d)=(random(),random(),random(),random());
            let x=Interval{lo:a.min(b),hi:a.max(b)};let y=Interval{lo:c.min(d),hi:c.max(d)};
            equal(x.chord_mul(y),x.mul(y));
            if y.lo>0.0 {equal(x.chord_div_positive(y),x.div_positive(y));}
        }
        let cases=[Interval::point(0.0),Interval::point(-0.0),Interval::point(f64::from_bits(1)),
            Interval{lo:0.0,hi:f64::MIN_POSITIVE},Interval{lo:-1.0,hi:1.0},Interval{lo:-2.0,hi:-1.0},
            Interval{lo:1.0,hi:f64::MAX},Interval{lo:0.0,hi:f64::INFINITY},Interval{lo:f64::NEG_INFINITY,hi:f64::INFINITY}];
        for x in cases {for y in cases {equal(x.chord_mul(y),x.mul(y));if y.lo>0.0 {equal(x.chord_div_positive(y),x.div_positive(y));}}}
    }
}
