// Research-only comparison: exact BVH pruning, row-cell lists, and block parallelism.
struct ProbeNode { low:Vec3, high:Vec3, r2:f32, children:Option<(usize,usize)>, ids:Vec<usize> }
fn probe_tree(atoms:&[AlignedAtom],mut ids:Vec<usize>,nodes:&mut Vec<ProbeNode>)->usize {
    let mut low=Vec3::splat(f32::INFINITY);let mut high=Vec3::splat(f32::NEG_INFINITY);let mut r2=0.0_f32;
    for &i in &ids { low=low.min(atoms[i].position);high=high.max(atoms[i].position);r2=r2.max(atoms[i].radius_squared); }
    let index=nodes.len();nodes.push(ProbeNode{low,high,r2,children:None,ids:Vec::new()});
    if ids.len()<=8 {nodes[index].ids=ids;} else {
        let span=high-low;let axis=if span.x>=span.y && span.x>=span.z {0}else if span.y>=span.z {1}else {2};
        ids.sort_unstable_by(|&a,&b|atoms[a].position[axis].total_cmp(&atoms[b].position[axis]).then(a.cmp(&b)));
        let right=ids.split_off(ids.len()/2);let a=probe_tree(atoms,ids,nodes);let b=probe_tree(atoms,right,nodes);nodes[index].children=Some((a,b));
    } index
}
fn probe_gap(a:f32,b:f32,lo:f32,hi:f32)->f32 {if b<lo {lo-b}else if a>hi {a-hi}else {0.0}}
fn probe_node_mask(index:usize,nodes:&[ProbeNode],grid:&PreparedGrid,block:&GridBlock,atoms:&[AlignedAtom],mask:&mut u64) {
    let full=u64::MAX>>(64-block.len);if *mask==full{return;}
    let n=&nodes[index];let z=&grid.z[block.start..block.start+block.len];
    let dx=probe_gap(block.x,block.x,n.low.x,n.high.x);let dy=probe_gap(block.y,block.y,n.low.y,n.high.y);
    let dz=probe_gap(z[0],z[z.len()-1],n.low.z,n.high.z);
    if (dx*dx+dy*dy)+dz*dz>n.r2{return;}
    if let Some((a,b))=n.children {probe_node_mask(a,nodes,grid,block,atoms,mask);probe_node_mask(b,nodes,grid,block,atoms,mask);} else {
        for &i in &n.ids {*mask|=probe_atom_mask(grid,block,&atoms[i]);if *mask==full {break;}}
    }
}
fn probe_atom_mask(grid:&PreparedGrid,block:&GridBlock,atom:&AlignedAtom)->u64 {
    let dx=block.x-atom.position.x;let dy=block.y-atom.position.y;let xy=dx*dx+dy*dy;
    if xy>atom.radius_squared{return 0;}
    let z=&grid.z[block.start..block.start+block.len];
    if std::is_x86_feature_detected!("avx2") {unsafe {row_mask_avx2(z,atom.position.z,xy,atom.radius_squared)}} else {
        z.iter().enumerate().fold(0,|m,(i,&p)|{let dz=p-atom.position.z;m|((u64::from(xy+dz*dz<=atom.radius_squared))<<i)})
    }
}
pub fn spatial_probe(q:&serde_json::Value)->serde_json::Value {
    use std::{hint::black_box,time::Instant};use rayon::prelude::*;
    let config=BuriedVolumeConfig{sphere_radius:q["radius"].as_f64().unwrap_or(3.5) as f32,density:q["density"].as_f64().unwrap_or(0.01) as f32,..BuriedVolumeConfig::default()};
    let atoms:Vec<_>=q["atoms"].as_array().unwrap().iter().map(|a|AlignedAtom{position:Vec3::from_array(std::array::from_fn(|i|a["position"][i].as_f64().unwrap() as f32)),radius_squared:a["radius_squared"].as_f64().unwrap() as f32}).collect();
    let start=Instant::now();let grid=PreparedGrid::new(integration_grid(config));let grid_ns=start.elapsed().as_nanos();
    let start=Instant::now();let mut nodes=Vec::new();probe_tree(&atoms,(0..atoms.len()).collect(),&mut nodes);let bvh_prepare_ns=start.elapsed().as_nanos();
    let start=Instant::now();let keys:Vec<_>=grid.blocks.iter().map(|b| {
        let axis=|x:f32|((((x+config.sphere_radius)/(2.0*config.sphere_radius))*8.0) as usize).min(7);
        axis(b.x)+8*axis(b.y)
    }).collect();
    let mut extents=vec![(Vec3::splat(f32::INFINITY),Vec3::splat(f32::NEG_INFINITY));64];
    for (b,&k) in grid.blocks.iter().zip(&keys) {let p=Vec3::new(b.x,b.y,0.0);extents[k].0=extents[k].0.min(p);extents[k].1=extents[k].1.max(p);}
    let candidates:Vec<Vec<usize>>=extents.iter().map(|&(lo,hi)|atoms.iter().enumerate().filter_map(|(i,a)|{
        if lo.x>hi.x{return None;}let dx=probe_gap(lo.x,hi.x,a.position.x,a.position.x);let dy=probe_gap(lo.y,hi.y,a.position.y,a.position.y);
        (dx*dx+dy*dy<=a.radius_squared).then_some(i)
    }).collect()).collect();let cells_prepare_ns=start.elapsed().as_nanos();
    let direct=|b:&GridBlock|{let mut m=0;let full=u64::MAX>>(64-b.len);for a in &atoms {m|=probe_atom_mask(&grid,b,a);if m==full{break;}}m};
    let cell=|(i,b):(usize,&GridBlock)|{let mut m=0;let full=u64::MAX>>(64-b.len);for &a in &candidates[keys[i]] {m|=probe_atom_mask(&grid,b,&atoms[a]);if m==full{break;}}m};
    let bvh=|b:&GridBlock|{let mut m=0;probe_node_mask(0,&nodes,&grid,b,&atoms,&mut m);m};
    let expected:Vec<_>=grid.blocks.iter().map(|b|grid.points[b.start..b.start+b.len].iter().enumerate().fold(0,|m,(i,p)|m|(u64::from(atoms.iter().any(|a|p.distance_squared(a.position)<=a.radius_squared))<<i))).collect();
    let pools:Vec<_>=[2,4,6].into_iter().map(|n|rayon::ThreadPoolBuilder::new().num_threads(n).build().unwrap()).collect();
    let modes=["prepared","cell_list","bvh","parallel_2","parallel_4","parallel_6"];
    let evaluate=|mode:usize|->Vec<u64>{match mode {0=>grid.blocks.iter().map(direct).collect(),1=>grid.blocks.iter().enumerate().map(cell).collect(),2=>grid.blocks.iter().map(bvh).collect(),p=>pools[p-3].install(||grid.blocks.par_iter().map(direct).collect())}};
    let mismatches:Vec<_>=(0..modes.len()).map(|m|evaluate(m).iter().zip(&expected).filter(|(a,b)|a!=b).count()).collect();
    let iterations=q["iterations"].as_u64().unwrap_or(20) as usize;let mut samples=Vec::new();
    for rep in 0..5 {for j in 0..modes.len(){let m=if rep%2==0 {j}else {modes.len()-1-j};let start=Instant::now();for _ in 0..iterations {black_box(evaluate(m));}samples.push(serde_json::json!({"mode":modes[m],"rep":rep,"warmup":rep==0,"iterations":iterations,"elapsed_ns":start.elapsed().as_nanos()}));}}
    serde_json::json!({"id":q["id"],"atoms":atoms.len(),"points":grid.len(),"blocks":grid.blocks.len(),"grid_ns":grid_ns,
        "grid_retained_bytes":grid.retained_bytes(),"bvh_prepare_ns":bvh_prepare_ns,"cells_prepare_ns":cells_prepare_ns,
        "bvh_nodes":nodes.len(),"cell_candidate_entries":candidates.iter().map(Vec::len).sum::<usize>(),
        "row_candidate_tests":keys.iter().map(|&k|candidates[k].len()).sum::<usize>(),"full_row_atom_tests":grid.blocks.len()*atoms.len(),
        "modes":modes,"mask_mismatches":mismatches,"samples":samples})
}
