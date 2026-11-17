#[derive(Default)]
struct ChordScratch {
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
    if spans.iter().any(|s| !s.regular || s.lower.is_none() || s.upper.is_none() || s.low.hi > s.high.lo) {
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
            let Some(component) = components.last_mut() else { return false; };
            component.1 = position + 1;
        }
    }
    let mut previous_upper = f64::NEG_INFINITY;
    endpoints.clear();
    for &(start, end) in components.iter() {
        let group = &order[start..end];
        let Some(&lower) = group.iter().min_by(|&&a, &&b| spans[a].low.lo.total_cmp(&spans[b].low.lo)) else {return false;};
        let Some(&upper) = group.iter().max_by(|&&a, &&b| spans[a].high.hi.total_cmp(&spans[b].high.hi)) else {return false;};
        let low = spans[lower];
        let high = spans[upper];
        if low.low.lo < previous_upper { return false; }
        for &i in group {
            if i != lower && low.low.hi > spans[i].low.lo && !same_endpoint(low.lower, spans[i].lower, spheres) { return false; }
            if i != upper && high.high.lo < spans[i].high.hi && !same_endpoint(high.upper, spans[i].upper, spheres) { return false; }
        }
        previous_upper = high.high.hi;
        let (Some(a), Some(b)) = (low.lower, high.upper) else {return false;};
        endpoints.push((a,b));
    }
    true
}
