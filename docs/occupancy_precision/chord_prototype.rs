// Research prototype appended to the frozen interval implementation.
// No midpoint-difference heuristic is used as an error certificate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct ChordEndpoint { sphere: usize, sign: i8 }
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
    lo: [f64; 2], hi: [f64; 2], bounds: Interval,
    candidates: Vec<usize>, ordinal: usize,
}
impl PartialEq for ChordCell { fn eq(&self,b:&Self)->bool {self.ordinal==b.ordinal} }
impl Eq for ChordCell {}
impl PartialOrd for ChordCell { fn partial_cmp(&self,b:&Self)->Option<Ordering>{Some(self.cmp(b))} }
impl Ord for ChordCell {
    fn cmp(&self,b:&Self)->Ordering {
        (self.bounds.hi-self.bounds.lo).total_cmp(&(b.bounds.hi-b.bounds.lo))
            .then_with(||b.ordinal.cmp(&self.ordinal))
    }
}
#[derive(Default)]
struct ChordCounters { evaluations: usize, smooth: usize, range: usize, empty: usize }

fn chord_q(s: BoundedSphere, lo:[f64;2],hi:[f64;2])->Interval {
    let mut q=s.radius_squared;
    for i in 0..2 {q=q.sub(Interval{lo:lo[i],hi:hi[i]}.sub(s.center[i]).square());}
    q
}
fn chord_span(index:usize,s:BoundedSphere,lo:[f64;2],hi:[f64;2])->Option<ChordSpan> {
    let q=chord_q(s,lo,hi);
    if q.hi<=0.0 {return None;}
    let extent=q.sqrt();
    Some(ChordSpan{low:s.center[2].sub(extent),high:s.center[2].add(extent),
        lower:Some(ChordEndpoint{sphere:index,sign:-1}),upper:Some(ChordEndpoint{sphere:index,sign:1}),regular:q.lo>0.0})
}
fn chord_max(a:Interval,b:Interval)->Interval {Interval{lo:a.lo.max(b.lo),hi:a.hi.max(b.hi)}}
fn chord_min(a:Interval,b:Interval)->Interval {Interval{lo:a.lo.min(b.lo),hi:a.hi.min(b.hi)}}
fn chord_length(mut ranges:Vec<(f64,f64)>,upper:bool)->f64 {
    ranges.retain(|&(a,b)|a<b);
    ranges.sort_unstable_by(|a,b|a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
    let Some(&(mut a,mut b))=ranges.first() else {return 0.0;};
    let mut sum=0.0;
    let round=|x:f64|if upper{x.up()}else{x.down()};
    for &(c,d) in &ranges[1..] {
        if c>b {sum=round(sum+round(b-a));a=c;b=d;} else {b=b.max(d);}
    }
    round(sum+round(b-a)).max(0.0)
}
fn same_endpoint(a:Option<ChordEndpoint>,b:Option<ChordEndpoint>,spheres:&[AtomicSphere])->bool {
    match (a,b) {(Some(a),Some(b))=>a.sign==b.sign && spheres[a.sphere]==spheres[b.sphere],_=>false}
}
fn chord_smooth(spans:&[ChordSpan],spheres:&[AtomicSphere])->Option<Vec<(ChordEndpoint,ChordEndpoint)>> {
    if spans.iter().any(|s|!s.regular || s.lower.is_none() || s.upper.is_none() || s.low.hi>s.high.lo) {return None;}
    let mut order:Vec<_>=(0..spans.len()).collect();
    order.sort_unstable_by(|&a,&b|spans[a].low.hi.total_cmp(&spans[b].low.hi).then(a.cmp(&b)));
    let mut components:Vec<Vec<usize>>=Vec::new();
    let mut inner_high=f64::NEG_INFINITY;
    for index in order {
        let s=spans[index];
        if s.low.hi>inner_high {components.push(Vec::new());inner_high=s.high.lo;}
        else {inner_high=inner_high.max(s.high.lo);}
        components.last_mut()?.push(index);
    }
    let mut previous_upper=f64::NEG_INFINITY;
    let mut endpoints=Vec::new();
    for group in components {
        let lower=*group.iter().min_by(|&&a,&&b|spans[a].low.lo.total_cmp(&spans[b].low.lo))?;
        let upper=*group.iter().max_by(|&&a,&&b|spans[a].high.hi.total_cmp(&spans[b].high.hi))?;
        let low=spans[lower];let high=spans[upper];
        if low.low.lo<previous_upper {return None;}
        for &i in &group {
            if i!=lower && low.low.hi>spans[i].low.lo && !same_endpoint(low.lower,spans[i].lower,spheres) {return None;}
            if i!=upper && high.high.lo<spans[i].high.hi && !same_endpoint(high.upper,spans[i].upper,spheres) {return None;}
        }
        previous_upper=high.high.hi;
        endpoints.push((low.lower?,high.upper?));
    }
    Some(endpoints)
}
fn chord_endpoint_value(endpoint:ChordEndpoint,spheres:&[BoundedSphere],mid:[Interval;2])->Interval {
    let s=spheres[endpoint.sphere];
    let mut q=s.radius_squared;
    for (i,m) in mid.into_iter().enumerate() {q=q.sub(m.sub(s.center[i]).square());}
    let extent=q.sqrt();
    if endpoint.sign<0 {s.center[2].sub(extent)} else {s.center[2].add(extent)}
}
fn chord_cell(lo:[f64;2],hi:[f64;2],indices:&[usize],spheres:&[AtomicSphere],
    bounds:&[BoundedSphere],ordinal:usize,counters:&mut ChordCounters)->ChordCell {
    let mut cell=ChordCell{lo,hi,bounds:Interval::point(0.0),candidates:Vec::new(),ordinal};
    counters.evaluations+=1;
    let Some(region)=chord_span(0,bounds[0],lo,hi) else {counters.empty+=1;return cell;};
    let widths=std::array::from_fn::<_,2,_>(|i|Interval::point(hi[i]).sub(Interval::point(lo[i])));
    let area=widths[0].mul(widths[1]);
    let mut spans=Vec::with_capacity(indices.len());
    for &i in indices {
        counters.evaluations+=1;
        let Some(mut s)=chord_span(i,bounds[i],lo,hi) else {continue;};
        let same=spheres[i]==spheres[0];
        s.lower=if same || s.low.lo>=region.low.hi {s.lower}
            else if region.low.lo>=s.low.hi {region.lower} else {None};
        s.upper=if same || s.high.hi<=region.high.lo {s.upper}
            else if region.high.hi<=s.high.lo {region.upper} else {None};
        s.low=chord_max(s.low,region.low);
        s.high=chord_min(s.high,region.high);
        s.regular &= region.regular;
        if s.low.lo>=s.high.hi {continue;}
        cell.candidates.push(i);spans.push(s);
    }
    if spans.is_empty() {counters.empty+=1;return cell;}
    let lower=chord_length(spans.iter().filter(|s|s.regular).map(|s|(s.low.hi,s.high.lo)).collect(),false);
    let upper=chord_length(spans.iter().map(|s|(s.low.lo,s.high.hi)).collect(),true);
    cell.bounds=Interval{lo:lower,hi:upper}.mul(area);
    cell.bounds.lo=cell.bounds.lo.max(0.0);
    let Some(endpoints)=chord_smooth(&spans,spheres) else {counters.range+=1;return cell;};
    let mid=std::array::from_fn(|i|Interval::point(lo[i]).add(Interval::point(hi[i])).mul(Interval::point(0.5)));
    let mut length=Interval::point(0.0);
    let mut curvature=Interval::point(0.0);
    for (a,b) in endpoints {
        length=length.add(chord_endpoint_value(b,bounds,mid).sub(chord_endpoint_value(a,bounds,mid)));
        for endpoint in [a,b] {
            let s=bounds[endpoint.sphere];
            let q=chord_q(s,lo,hi);
            let qmin=Interval::point(q.lo);
            let denominator=qmin.mul(qmin.sqrt());
            if denominator.lo<=0.0 || !denominator.hi.is_finite() {counters.range+=1;return cell;}
            curvature=curvature.add(Interval::point(s.radius_squared.hi).div_positive(denominator));
        }
    }
    let remainder=area.mul(widths[0].square().add(widths[1].square())).mul(curvature)
        .div_positive(Interval::point(24.0));
    let integral=length.mul(area).add(Interval{lo:-remainder.hi,hi:remainder.hi});
    if integral.lo.is_finite() && integral.hi.is_finite() {
        cell.bounds.lo=cell.bounds.lo.max(integral.lo);
        cell.bounds.hi=cell.bounds.hi.min(integral.hi);
        assert!(cell.bounds.lo<=cell.bounds.hi);
        counters.smooth+=1;
    } else {counters.range+=1;}
    cell
}

pub fn chord_probe(input:&serde_json::Value)->serde_json::Value {
    use std::time::Instant;
    let region:AtomicSphere=serde_json::from_value(input["region"].clone()).unwrap();
    let mut atoms:Vec<AtomicSphere>=serde_json::from_value(input["atoms"].clone()).unwrap();
    // Validate through the existing explicit API, without spending a box budget.
    sphere_union_volume(region,&atoms,ContinuousVolumeMethod::Reference{absolute_tolerance:1.0,max_boxes:1}).unwrap();
    let tolerance=input["tolerance"].as_f64().unwrap_or(0.05);
    let budget=input["budget"].as_u64().unwrap_or(262143) as usize;
    assert!(tolerance>0.0 && tolerance.is_finite() && budget>0);
    let start=Instant::now();
    atoms.retain(|a|a.radius>0.0);
    atoms.sort_by(|a,b|a.center[0].total_cmp(&b.center[0]).then(a.center[1].total_cmp(&b.center[1]))
        .then(a.center[2].total_cmp(&b.center[2])).then(a.radius.total_cmp(&b.radius)));
    atoms.dedup();
    let spheres:Vec<_>=std::iter::once(region).chain(atoms).collect();
    let bounds:Vec<_>=spheres.iter().map(|&s|relative_sphere(s,region.center)).collect();
    let mut counters=ChordCounters::default();
    let root=chord_cell([-region.radius;2],[region.radius;2],&(1..spheres.len()).collect::<Vec<_>>(),
        &spheres,&bounds,0,&mut counters);
    let mut lower=root.bounds.lo;let mut upper=root.bounds.hi;
    let mut best_lower=lower;let mut best_upper=upper.min(ball_volume(region.radius).hi);
    let mut heap=BinaryHeap::new();heap.push(root);
    let mut cells=1usize;let mut subdivisions=0usize;let mut peak=1usize;
    let (estimate,error,met,status)=loop {
        best_lower=best_lower.max(lower);best_upper=best_upper.min(upper);
        assert!(best_lower<=best_upper,"independent accumulated bounds must overlap");
        let value=best_lower+(best_upper-best_lower)*0.5;
        let error=if best_lower==best_upper {0.0} else {(value-best_lower).max(best_upper-value).up()};
        if error<=tolerance {break(value,error,true,"tolerance_achieved");}
        if cells>budget.saturating_sub(2) {break(value,error,false,"cell_budget_exhausted");}
        let Some(parent)=heap.pop() else {break(value,error,false,"floating_point_limit");};
        let axis=usize::from(parent.hi[1]-parent.lo[1]>=parent.hi[0]-parent.lo[0]);
        let middle=parent.lo[axis]+(parent.hi[axis]-parent.lo[axis])*0.5;
        if middle<=parent.lo[axis] || middle>=parent.hi[axis] {continue;}
        let mut left=parent.hi;left[axis]=middle;
        let mut right=parent.lo;right[axis]=middle;
        lower=(lower-parent.bounds.lo).down();upper=(upper-parent.bounds.hi).up();
        for (lo,hi) in [(parent.lo,left),(right,parent.hi)] {
            let child=chord_cell(lo,hi,&parent.candidates,&spheres,&bounds,cells,&mut counters);cells+=1;
            lower=(lower+child.bounds.lo).down();upper=(upper+child.bounds.hi).up();
            if child.bounds.hi>child.bounds.lo {heap.push(child);}
        }
        subdivisions+=1;peak=peak.max(heap.len());
    };
    serde_json::json!({"id":input["id"],"volume":estimate,"lower":best_lower,"upper":best_upper,
        "absolute_error_bound":error,"relative_error_bound":if best_lower>0.0 {Some((error/best_lower).up())}else{None},
        "target_met":met,"status":status,"cells":cells,"subdivisions":subdivisions,
        "geometric_evaluations":counters.evaluations,"smooth_cells":counters.smooth,
        "range_cells":counters.range,"empty_cells":counters.empty,"peak_heap_cells":peak,
        "elapsed_ns":start.elapsed().as_nanos()})
}
