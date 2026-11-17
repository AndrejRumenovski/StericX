// Appended only to a frozen copy of the private certified module for inspection.
pub fn inspect_predicate(q: &serde_json::Value) -> serde_json::Value {
    let origin: [f64; 3] = serde_json::from_value(q["origin"].clone()).unwrap();
    let sphere: AtomicSphere = serde_json::from_value(q["sphere"].clone()).unwrap();
    let lo: [f64; 2] = serde_json::from_value(q["lo"].clone()).unwrap();
    let hi: [f64; 2] = serde_json::from_value(q["hi"].clone()).unwrap();
    let s = relative_sphere(sphere, origin);
    let value = chord_q(s, lo, hi);
    let span = chord_span(0, s, lo, hi).map(|s| serde_json::json!({
        "low":[s.low.lo,s.low.hi],"high":[s.high.lo,s.high.hi],"regular":s.regular}));
    let mid = std::array::from_fn(|i| Interval::point(lo[i])
        .add(Interval::point(hi[i])).mul(Interval::point(0.5)));
    let endpoints: Vec<_> = [-1, 1].into_iter().map(|sign| {
        let v=chord_endpoint_value(ChordEndpoint{sphere:0,sign},&[s],mid);
        [v.lo,v.hi]
    }).collect();
    serde_json::json!({"id":q["id"],"q":[value.lo,value.hi],"span":span,
        "mid":[[mid[0].lo,mid[0].hi],[mid[1].lo,mid[1].hi]],"mid_endpoints":endpoints})
}
